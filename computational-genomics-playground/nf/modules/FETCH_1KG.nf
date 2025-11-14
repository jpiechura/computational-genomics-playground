nextflow.enable.dsl=2
/* ---------------- 1) Fetch 1000G VCF + panel ---------------- */
process FETCH_1KG {
  tag "chr${params.chr}"
  label "bcftools"

  publishDir "${params.run_outdir}/1000g",  mode: 'copy', pattern: "*"

  input:
    path fetch_script

  output:
    path "ALL.chr${params.chr}*.vcf.gz",     emit: vcf
    path "ALL.chr${params.chr}*.vcf.gz.tbi", emit: vcf_tbi
    path "integrated_call_samples_v3.20130502.ALL.panel", emit: panel

  script:
  """
  apt-get update && apt-get install -y --no-install-recommends curl ca-certificates && rm -rf /var/lib/apt/lists/*
  bash "$fetch_script" ${params.chr} 
  """
}