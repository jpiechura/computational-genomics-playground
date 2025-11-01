nextflow.enable.dsl=2

process QC_HET {
tag "qc_het_${het_stats.baseName}"
publishDir "qc", mode: 'copy'
container "python:3.11"

input:
  path het_stats   // FID IID PHENO
  path pyenv

output:
  path "het_outliers.txt", emit: het_outliers


script:
def prefix = "chr${params.chr}.${params.pop}${params.n_samples}"
"""
# 4) Heterozygosity outliers
./pyenv/bin/python - <<PY
import pandas as pd, numpy as np
df=pd.read_csv("${het_stats}",delim_whitespace=True)
mu,sd=df["F"].mean(),df["F"].std(ddof=1)
m=(df.F<mu-${params.het_sd}*sd)|(df.F>mu+${params.het_sd}*sd)
df.loc[m,["FID","IID"]].to_csv("het_outliers.txt",sep="\\t",header=False,index=False)
print("HET_OUTLIERS", int(m.sum()))
PY

"""
}