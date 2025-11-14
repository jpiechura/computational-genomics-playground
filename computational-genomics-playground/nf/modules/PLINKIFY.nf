/* --------------- 4) Convert subset VCF → PLINK -------------- */
nextflow.enable.dsl=2

process PLINKIFY {
  tag "plink_chr${params.chr}_${params.pop}${params.n_samples}"
  label "plink"

  input:
    path vcf_file
    path plinkify_script

  publishDir "${params.run_outdir}/1000g", mode: 'copy', pattern: "chr${params.chr}.${params.pop}${params.n_samples}.*"

  output:
    path "chr${params.chr}.${params.pop}${params.n_samples}.bim", emit: bim
    path "chr${params.chr}.${params.pop}${params.n_samples}.fam", emit: fam
    path "chr${params.chr}.${params.pop}${params.n_samples}.bed", emit: bed

  script:
  """
  OUTP=\$PWD/chr${params.chr}.${params.pop}${params.n_samples}
  bash "$plinkify_script" "$vcf_file" "\$OUTP"
  """
}