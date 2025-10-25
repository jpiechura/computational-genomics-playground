/*
 * main.nf — GWAS demo pipeline: fetch → sample list → subset → plink → gcta simulate
 * DSL2, scripts staged from workflow (no `from file(...)` inside process inputs)
 */

nextflow.enable.dsl = 2

// ---- Defaults ----
def launcher_path = projectDir

params.chr       = params.chr       ?: 22
params.pop       = params.pop       ?: 'EUR'
params.n_samples = params.n_samples ?: 500
params.vcfdir    = "${launcher_path}/data/1kg/vcf"
params.metadir   = "${launcher_path}/data/1kg/meta"
params.plinkdir  = "${launcher_path}/data/1kg/plink"
params.workdir   = "${launcher_path}/work"
params.h2        = params.h2        ?: 0.5
params.n_causal  = params.n_causal  ?: 50
params.seed      = params.seed      ?: 142

// Ensure directories exist
new File(params.vcfdir).mkdirs()
new File(params.metadir).mkdirs()
new File(params.plinkdir).mkdirs()
new File(params.workdir).mkdirs()

workflow.onComplete {
  println "Done. Check ${params.vcfdir}, ${params.metadir}, ${params.plinkdir}, and ${params.workdir}."
}

/* ---------------- 1) Fetch 1000G VCF + panel ---------------- */
process FETCH_1KG {
  tag "chr${params.chr}"

  publishDir "${params.vcfdir}",  mode: 'copy', pattern: "vcf/*"
  publishDir "${params.metadir}", mode: 'copy', pattern: "meta/*"

  input:
    path fetch_script

  output:
    path "vcf/ALL.chr${params.chr}*.vcf.gz",     emit: vcf
    path "vcf/ALL.chr${params.chr}*.vcf.gz.tbi", emit: vcf_tbi
    path "meta/integrated_call_samples_v3.20130502.ALL.panel", emit: panel

  script:
  """
  apt-get update && apt-get install -y --no-install-recommends curl ca-certificates && rm -rf /var/lib/apt/lists/*
  mkdir -p vcf meta
  bash "$fetch_script" ${params.chr} vcf meta
  """
}

/* ---------------- 2) Make random sample list ---------------- */
process MAKE_SAMPLE_LIST {
  tag "${params.pop}_${params.n_samples}"
  container 'debian:bookworm-slim'
  containerOptions '-u 0:0'

  input:
    path panel
    path mklist_script

  publishDir "${params.metadir}", mode: 'copy', pattern: "samples_${params.pop}_${params.n_samples}.txt"

  output:
    path "samples_${params.pop}_${params.n_samples}.txt", emit: samples

  script:
  """
  apt-get update && apt-get install -y --no-install-recommends bash coreutils ca-certificates && rm -rf /var/lib/apt/lists/*
  bash "$mklist_script" "$panel" ${params.pop} ${params.n_samples} "samples_${params.pop}_${params.n_samples}.txt"
  """
}

/* ------------- 3) Subset VCF to those samples -------------- */
process SUBSET_VCF {
  tag "subset_chr${params.chr}_${params.pop}${params.n_samples}"
  container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'

  input:
    path vcf_in
    path samples
    path subset_script

  publishDir "${params.vcfdir}", mode: 'copy', pattern: "chr${params.chr}.${params.pop}${params.n_samples}.vcf.gz*"

  output:
    path "chr${params.chr}.${params.pop}${params.n_samples}.vcf.gz",     emit: vcf_sub
    path "chr${params.chr}.${params.pop}${params.n_samples}.vcf.gz.tbi", emit: vcf_sub_tbi

  script:
  def invcf  = vcf_in.toString()
  def outvcf = "chr${params.chr}.${params.pop}${params.n_samples}.vcf.gz"
  """
  bash "$subset_script" "${invcf}" "${samples}" "${outvcf}"
  """
}

/* --------------- 4) Convert subset VCF → PLINK -------------- */
process PLINKIFY {
  tag "plink_chr${params.chr}_${params.pop}${params.n_samples}"
  container 'quay.io/biocontainers/plink:1.90b6.21--hec16e2b_2'

  input:
    path vcf_file
    path plinkify_script

  publishDir "${params.plinkdir}", mode: 'copy', pattern: "chr${params.chr}.${params.pop}${params.n_samples}.*"

  output:
    path "chr${params.chr}.${params.pop}${params.n_samples}.*", emit: plink_files

  script:
  """
  OUTP=\$PWD/chr${params.chr}.${params.pop}${params.n_samples}
  bash "$plinkify_script" "$vcf_file" "\$OUTP"
  """
}

/* --------- 5) Build GRM + simulate phenotype with GCTA ------- */
process GCTA_SIM {
  tag "gcta_sim_chr${params.chr}_${params.pop}${params.n_samples}"

  input:
    path bed
    path bim
    path fam
    path 'dummy.sh'  // Dummy input to maintain compatibility

  publishDir "${params.workdir}", mode: 'copy', pattern: "sim.*"

  output:
    path "sim.*", emit: sim_outputs

  script:
'''
set -euo pipefail

# Derive PLINK prefix (basename without .bed) from the staged BED path
bed_file="''' + bed + '''"
base="$(basename "$bed_file")"
base="${base%.bed}"

echo "Input files:"
ls -lh

echo "Running GCTA GRM..."
gcta64 --bfile "$base" \
       --make-grm \
       --out sim \
       --thread-num ''' + task.cpus + '''

# Count unique SNP IDs first (BIM col2)
nuniq=$(awk '!seen[$2]++{c++} END{print c+0}' "${base}.bim")
if [ "$nuniq" -lt ''' + params.n_causal + ''' ]; then
  echo "Not enough unique SNP IDs in ${base}.bim: have $nuniq, need ''' + params.n_causal + '''" >&2
  exit 1
fi

# Build a unique, reproducibly shuffled causal list
#  1) extract unique IDs (keep first occurrence)
#  2) assign seeded random key
#  3) sort by key (and by ID as a tiebreaker for determinism)
#  4) take N
awk '!seen[$2]++{print $2}' "${base}.bim" \
  | awk -v seed=''' + params.seed + ''' 'BEGIN{srand(seed)} {printf "%0.12f\t%s\n", rand(), $0}' \
  | sort -k1,1n -k2,2 \
  | cut -f2 \
  | head -n ''' + params.n_causal + ''' > causal.snplist

echo "Selected $(wc -l < causal.snplist) unique causal SNPs."

echo "Simulating quantitative trait..."
gcta64 --bfile "$base" \
       --simu-qt \
       --grm sim \
       --simu-hsq ''' + params.h2 + ''' \
       --simu-causal-loci causal.snplist \
       --simu-seed ''' + params.seed + ''' \
       --out sim
'''
}

/* -------------------- Wiring (workflow) --------------------- */
workflow {
  log.info "Starting pipeline..."

  // Stage 1: Fetch data
  def fetched = FETCH_1KG( file('bin/00_fetch_1000g.sh') )

  // Stage 2: Create sample list
  def picked = MAKE_SAMPLE_LIST(
    fetched.panel,
    file('bin/01_make_sample_list.sh')
  )

  // Stage 3: Subset VCF
  def subset = SUBSET_VCF(
    fetched.vcf,
    picked.samples,
    file('bin/02_subset_vcf.sh')
  )

  // Stage 4: Convert to PLINK format
  def plinked = PLINKIFY(
    subset.vcf_sub,
    file('bin/03_plinkify.sh')
  )

  // Stage 5: Run GCTA simulation
  plinked.plink_files
  .collect()
  .map { files ->
    def base = "chr${params.chr}.${params.pop}${params.n_samples}"
    def bed = files.find { it.name == "${base}.bed" }
    def bim = files.find { it.name == "${base}.bim" }
    def fam = files.find { it.name == "${base}.fam" }
    if (!bed || !bim || !fam) {
      error "Missing PLINK files. Found: BED=${bed}, BIM=${bim}, FAM=${fam}"
    }
    log.info "PLINK files found: BED=${bed}, BIM=${bim}, FAM=${fam}"
    tuple bed, bim, fam
  }
  .set { gcta_inputs }

  GCTA_SIM(
    gcta_inputs.map { it[0] },
    gcta_inputs.map { it[1] },
    gcta_inputs.map { it[2] },
    file('dummy.sh')
  )
}