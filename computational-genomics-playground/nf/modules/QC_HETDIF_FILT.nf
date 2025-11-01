nextflow.enable.dsl=2

process QC_HETDIF_FILT {
    tag "qc_het_${bed.baseName}"
    publishDir "qc", mode: 'copy'
    container "${params.plink_container}"
    
    input:
        path het_outliers   
        path pheno
        path bed
        path bim
        path fam
    
    output:
        path "*.qcd.bed", emit: bed
        path "*.qcd.bim", emit: bim
        path "*.qcd.fam", emit: fam
        
    script:
    def prefix = "chr${params.chr}.${params.pop}${params.n_samples}"
    """

  plink --bed $bed --bim $bim --fam $fam --make-bed --out ${prefix}.hetclean --remove ${het_outliers}
  BIN=\$(awk '{print \$3}' $pheno | awk '(\$1==0||\$1==1||\$1==2){c++} END{print (c>0)?"yes":"no"}')
  # 5) Differential missingness (cases/controls)
  if [ "\$BIN" = "yes" ] && [ "${params.diff_missing}" = "true" ]; then
    plink --bfile ${prefix}.hetclean --pheno $pheno --test-missing --out ${prefix}.diffmiss
    awk 'NR>1 && \$NF < ${params.diffmiss_p} {print \$2}' ${prefix}.diffmiss.missing > snps_diffmiss_exclude.txt || true
    EX2=""
    [ -s snps_diffmiss_exclude.txt ] && EX2="--exclude snps_diffmiss_exclude.txt"
    plink --bfile ${prefix}.hetclean \$EX2 --make-bed --out ${prefix}.qcd
  else
    cp ${prefix}.hetclean.bed ${prefix}.qcd.bed
    cp ${prefix}.hetclean.bim ${prefix}.qcd.bim
    cp ${prefix}.hetclean.fam ${prefix}.qcd.fam
  fi
    """

    }
    