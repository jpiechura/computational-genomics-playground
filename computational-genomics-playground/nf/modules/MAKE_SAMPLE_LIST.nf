nextflow.enable.dsl=2
/* ---------------- 2) Make random sample list ---------------- */
process MAKE_SAMPLE_LIST {
  tag "${params.pop}_${params.n_samples}"
  label "bcftools"

  input:
    path panel
    path mklist_script

  publishDir "${params.run_outdir}/1000g", mode: 'copy', pattern: "samples_${params.pop}_${params.n_samples}.txt"

  output:
    path "samples_${params.pop}_${params.n_samples}.txt", emit: samples

  script:
  """
  apt-get update && apt-get install -y --no-install-recommends bash coreutils ca-certificates && rm -rf /var/lib/apt/lists/*
  bash "$mklist_script" "$panel" ${params.pop} ${params.n_samples} "samples_${params.pop}_${params.n_samples}.txt"
  """
}