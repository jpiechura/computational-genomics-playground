nextflow.enable.dsl=2
process RUN_SUSIE_RSS {
  tag "susie_${locus_id}"
  container 'susier:latest'
  publishDir "finemap", mode: 'copy'

input:
      tuple val(locus_id), path(zfile), path(rfile), path(keepfile)
output:
      tuple val(locus_id),
            path("${locus_id}.pip.tsv"),
            path("${locus_id}.credible_sets.tsv"),
            path("${locus_id}.susie_fit.rds"), emit: susie_results

    script:
    """
    echo "ID: $locus_id"
    echo "Z:  $zfile"
    echo "R:  $rfile"
    echo "K:  $keepfile"
    set -euo pipefail
    susie_rss.R ${zfile} ${rfile} ${keepfile} ${locus_id} 5 ${params.n_samples}
    """
}