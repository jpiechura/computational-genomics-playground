nextflow.enable.dsl=2
process PLOT_QQ_MANHATTAN {
  tag "gwas_plot_chr${params.chr}_${params.pop}${params.n_samples}"
  container "python:3.11"
  publishDir "gwas", mode: 'copy'
  input:
    path gwas_tsv

  output:
    path "gwas_all.with_q.tsv", emit: gwas_q
    path "lambda.txt",          emit: lambda
    path "top_hits.tsv",        emit: tops
    path "manhattan.png",       emit: manhattan
    path "qq.png",              emit: qq

  

  script:
  """
  set -euo pipefail

  python3 -m pip install --no-cache-dir -q pandas numpy scipy matplotlib

  python3 ${projectDir}/bin/manqq.py "${gwas_tsv}" "${params.p_gws}"
  """
}