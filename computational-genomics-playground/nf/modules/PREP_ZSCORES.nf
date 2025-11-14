process PREP_ZSCORES {
  tag "prepz_${locus_id}"
  label "htslib"
  publishDir "${params.run_outdir}/finemap/${locus_id}", mode: 'copy', overwrite: true

  input:
    tuple val(locus_id), path(zfile), path(prune_in), path(pruned_bim)

  output:
    tuple val(locus_id), path("z.${locus_id}.pruned.tsv"), emit: zscores

  script:
  """
  set -euo pipefail

  tr -d '\\r' < ${prune_in} > prune.nocr.txt

  awk -f ${projectDir}/bin/prep_z.awk \\
      -v BIM="${pruned_bim}" \\
      -v PIN="prune.nocr.txt" \\
      -v OUT="z.${locus_id}.pruned.tsv" \\
      -F'[ \\t]+' \\
      "${zfile}"

  # Rename STAT → Z in header
  awk -F'\\t' 'NR==1{
      for(i=1;i<=NF;i++) if(\$i=="STAT") \$i="Z";
      print;
      next
    }1' OFS='\\t' z.${locus_id}.pruned.tsv > z.${locus_id}.pruned.tmp

  mv z.${locus_id}.pruned.tmp z.${locus_id}.pruned.tsv
  """
}
