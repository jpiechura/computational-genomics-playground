nextflow.enable.dsl=2
process COMPUTE_LD {
  tag "computeld_${locus_id}"
  label "plink"
  publishDir "${params.run_outdir}/finemap/${locus_id}", mode: 'copy'
    input:
      tuple val(locus_id), path(bed), path(bim), path(fam), path(snplist)
    output:
      tuple val(locus_id), path("ld.${locus_id}.ld"), emit : ld_matrix

    script:
    """
    set -euo pipefail

    # LD as correlation (symmetric, SNP x SNP), gzip for R
    plink \\
      --bfile ${bed.baseName} \\
      --extract ${snplist} \\
      --r square \\
      --threads 2 \\
      --out ld.${locus_id}

    """
}