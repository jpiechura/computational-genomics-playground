nextflow.enable.dsl=2

process GWAS_SUM {
tag "gwas_sum_chr${params.chr}_${params.pop}${params.n_samples}"
publishDir "${params.run_outdir}/gwas", mode: 'copy'
label "python311"
input:
    path plink_output


output:
    path "*.gwas.tsv", emit: gwas_tsv



script:
def prefix = "chr${params.chr}.${params.pop}${params.n_samples}"
"""
# Convert PLINK1.9 outputs to a unified TSV:
  #  - For logistic: p19_chr*.assoc.logistic
  #  - For linear:   p19_chr*.assoc.linear
  python3 -m pip install --no-cache-dir -q pandas numpy 
  python3 - <<'PY'
import glob, pandas as pd, numpy as np

trait = "${params.trait_type}"
if trait == "binary":
    files = glob.glob("*.assoc.logistic")
    if not files:
        raise SystemExit("No .assoc.logistic produced")
    dfs = []
    for f in files:
        df = pd.read_csv(f, delim_whitespace=True, comment="#")
        # keep ADD test only
        if "TEST" in df.columns:
            df = df[df["TEST"].eq("ADD")]
        # Standardize
        out = pd.DataFrame({
            "CHR":  df["CHR"],
            "BP":   df["BP"],
            "SNP":  df["SNP"],
            "A1":   df["A1"],
            # PLINK1.9 logistic has OR/BETA depending on flags; when OR present, convert to BETA=ln(OR)
            "BETA": np.log(df["OR"]) if "OR" in df.columns else df.get("BETA", np.nan),
            "STAT":   df.get("STAT", np.nan),
            "P":    df["P"]
        })
        dfs.append(out)
    res = pd.concat(dfs, ignore_index=True)
else:
    files = glob.glob("*.assoc.linear")
    if not files:
        raise SystemExit("No .assoc.linear produced")
    dfs = []
    for f in files:
        df = pd.read_csv(f, delim_whitespace=True, comment="#")
        # keep ADD test only
        if "TEST" in df.columns:
            df = df[df["TEST"].eq("ADD")]
        out = pd.DataFrame({
            "CHR":  df["CHR"],
            "BP":   df["BP"],
            "SNP":  df["SNP"],
            "A1":   df["A1"],
            "BETA": df["BETA"],
            "STAT":   df["STAT"],
            "P":    df["P"]
        })
        dfs.append(out)
    res = pd.concat(dfs, ignore_index=True)

# clean/sort
res = res.dropna(subset=["CHR","BP","SNP","P"])
res["CHR"] = pd.to_numeric(res["CHR"], errors="coerce")
res["BP"]  = pd.to_numeric(res["BP"], errors="coerce")
res["P"]   = pd.to_numeric(res["P"],  errors="coerce")
res = res.sort_values(["CHR","BP"])
res.to_csv("${prefix}.gwas.tsv", sep="\\t", index=False)
PY
  """
}