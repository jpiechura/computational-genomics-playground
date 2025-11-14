nextflow.enable.dsl=2

process GWAS_PLINK {
  tag "gwas_chr${params.chr}_${params.pop}${params.n_samples}"
  label "plink"
  input:
    path bed
    path bim
    path fam
    path pheno_file
    path covar_file
  
  publishDir "${params.run_outdir}/gwas", mode: 'copy'

  output:
    path "*.assoc.*", emit: plink_output

  script:
  def prefix = "chr${params.chr}.${params.pop}${params.n_samples}"
  def gwasCmd = (params.trait_type == 'binary')
      ? '--logistic hide-covar'
      : '--linear hide-covar'

  """
  
  plink \\
    --bed $bed --bim $bim --fam $fam \\
    --chr ${params.chr} \\
    --pheno ${pheno_file} \\
    --covar $covar_file \\
    --covar-name ${params.covar_names} \\
    --threads ${task.cpus} \\
    ${gwasCmd} \\
    --allow-no-sex \\
    --out p19_gwas_${prefix}
    """


}