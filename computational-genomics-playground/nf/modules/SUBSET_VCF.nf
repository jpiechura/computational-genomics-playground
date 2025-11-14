nextflow.enable.dsl=2
/* ------------- 3) Subset VCF to those samples -------------- */
process SUBSET_VCF {
  tag "subset_chr${params.chr}_${params.pop}${params.n_samples}"
  label "bcftools"

  input:
    path vcf_in
    path samples
    path subset_script

  publishDir "${params.run_outdir}/1000g", mode: 'copy', pattern: "chr${params.chr}.${params.pop}${params.n_samples}.vcf.gz*"

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
