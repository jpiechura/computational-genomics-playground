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
params.phendir   = "${launcher_path}/pheno"
params.qcdir     = "${launcher_path}/qc"
params.h2        = params.h2        ?: 0.5
params.n_causal  = params.n_causal  ?: 50
params.seed      = params.seed      ?: 142

// Ensure directories exist
new File(params.vcfdir).mkdirs()
new File(params.metadir).mkdirs()
new File(params.plinkdir).mkdirs()
new File(params.workdir).mkdirs()
new File(params.phendir).mkdirs()
new File(params.qcdir).mkdirs()

include { PLINKIFY } from './nf/modules/PLINKIFY.nf'
include { GCTA_SIM } from './nf/modules/GCTA_SIM.nf'
include { QC_PLINK } from './nf/modules/QC_PLINK.nf'
include { QC_HET } from './nf/modules/QC_HET.nf'
include { QC_HETDIF_FILT } from './nf/modules/QC_HETDIF_FILT.nf'
include { MAKE_PY_ENV } from './nf/modules/MAKE_PY_ENV.nf'
include { PCA_PLINK } from './nf/modules/PCA_PLINK.nf'
include { BUILD_COVAR } from './nf/modules/BUILD_COVAR.nf'
include { GWAS_PLINK } from './nf/modules/GWAS_PLINK.nf'


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

  def gcta_results = GCTA_SIM(
    plinked.bed,
    plinked.bim,
    plinked.fam,
  )

  def plink_qc = QC_PLINK(
    plinked.bed,
    plinked.bim,
    plinked.fam,
    gcta_results.pheno_file
  )

  def pyenv = MAKE_PY_ENV()

  def qc_het = QC_HET(
    plink_qc.het_stats,
    pyenv
  )

  def qc_hetfilt = QC_HETDIF_FILT(
    qc_het.het_outliers,
    gcta_results.pheno_file,
    plink_qc.bed,
    plink_qc.bim,
    plink_qc.fam
  )

  def pca_results = PCA_PLINK(
    qc_hetfilt.bed,
    qc_hetfilt.bim,
    qc_hetfilt.fam
  )

Channel.fromPath('bin/append_sex.awk').set { append_sex_script }
  covar_file = BUILD_COVAR(
    pca_results.eigenvec,
    fetched.panel,
    append_sex_script
  )


  def gwas_results = GWAS_PLINK(
    qc_hetfilt.bed,
    qc_hetfilt.bim,
    qc_hetfilt.fam,
    gcta_results.pheno_file,
    covar_file.covar_file
  )
}