process MAKE_REGION_BFILE {
  tag "makebfile_${locus_id}"
  label "plink"
  publishDir "${params.run_outdir}/finemap/${locus_id}", mode: 'copy', overwrite: true

  input:
    tuple val(locus_id), val(chr), val(start), val(end), path(sumstats)
    path bed
    path bim
    path fam

  output:
    // pruned bfile + prune list (everything Z-prep needs to align)
    tuple val(locus_id),
          path("region_${locus_id}.pruned.bed"),
          path("region_${locus_id}.pruned.bim"),
          path("region_${locus_id}.pruned.fam"),
          path("pruned.${locus_id}.in"),
          emit: pruned_bfile

  script:
  """
  set -euo pipefail

  # 1) SNP ids from sumstats (gz or plain)
  if gzip -t "${sumstats}" >/dev/null 2>&1; then
      zcat "${sumstats}"
  else
      cat "${sumstats}"
  fi | awk 'NR>1{print \$3}' > snps.${locus_id}.txt

  # 2) Build regional bfile filtered by chr/bp + sumstats SNPs + MAF
  plink \\
    --bed ${bed} --bim ${bim} --fam ${fam} \\
    --chr ${chr} --from-bp ${start} --to-bp ${end} \\
    --extract snps.${locus_id}.txt \\
    --maf ${params.min_maf} \\
    --make-bed --out region_${locus_id}

  # 3) Prune within the region bfile
  plink \\
    --bfile region_${locus_id} \\
    --indep-pairwise ${params.indep_kb} ${params.indep_step} ${params.indep_r2} \\
    --out pruned.${locus_id}

  mv pruned.${locus_id}.prune.in pruned.${locus_id}.in

  # 4) Build pruned regional bfile to match the pruned SNP set exactly
  plink \\
    --bfile region_${locus_id} \\
    --extract pruned.${locus_id}.in \\
    --make-bed --out region_${locus_id}.pruned
  """
}
