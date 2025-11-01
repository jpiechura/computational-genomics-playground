nextflow.enable.dsl=2

process QC_PLINK {
  tag "qc_plink_${bed.baseName}"
  publishDir "qc", mode: 'copy'
  container "${params.plink_container}"

  input:
    path bed
    path bim
    path fam
    path pheno   // FID IID PHENO
   output:
    path "*.hwe.bed", emit: bed
    path "*.hwe.bim", emit: bim
    path "*.hwe.fam", emit: fam
    path "*.het.het", emit: het_stats

  script:
  def prefix = "chr${params.chr}.${params.pop}${params.n_samples}"
  def hwe_scope = 'controls'  // we'll decide inside script based on data
  """
  set -euo pipefail
  

 
  # 1) Missingness filters
  plink --bed $bed --bim $bim --fam $fam --geno ${params.geno} --mind ${params.mind} --make-bed --out ${prefix}.missfilt

  # 2) MAF
  plink --bfile ${prefix}.missfilt --maf ${params.maf} --make-bed --out ${prefix}.maf

  # 3) HWE (choose controls if phenotype looks binary)
  BIN=\$(awk '{print \$3}' $pheno | awk '(\$1==0||\$1==1||\$1==2){c++} END{print (c>0)?"yes":"no"}')
  if [ "\$BIN" = "yes" ]; then
    plink --bfile ${prefix}.maf --hwe ${params.hwe} midp controls --make-bed --out ${prefix}.hwe
  else
    plink --bfile ${prefix}.maf --hwe ${params.hwe} midp          --make-bed --out ${prefix}.hwe
  fi
  
  plink --bfile ${prefix}.hwe --het --out ${prefix}.het

  """
}