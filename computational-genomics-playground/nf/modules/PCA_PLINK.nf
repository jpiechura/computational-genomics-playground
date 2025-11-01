nextflow.enable.dsl=2

process PCA_PLINK {
  tag "${bed.baseName}"
  publishDir "pca", mode: 'copy'
  container "${params.plink_container}"

  input:
    path bed
    path bim
    path fam

  output:
    path "*.pca.eigenvec", emit: eigenvec
    path "*.pca.eigenval", emit: eigenval
    path "*.prune.*",      emit: prune

  script:
  def prefix = "chr${params.chr}.${params.pop}${params.n_samples}"
  """
  

  plink --bed $bed --bim $bim --fam $fam --indep-pairwise ${params.ld_kb} ${params.ld_step} ${params.ld_r2} --out ${prefix}.prune
  plink --bed $bed --bim $bim --fam $fam --extract ${prefix}.prune.prune.in --make-bed --out ${prefix}.pruned
  plink --bfile ${prefix}.pruned --pca ${params.n_pcs} header tabs --out ${prefix}.pca
  """
}
