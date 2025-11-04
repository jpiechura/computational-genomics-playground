nextflow.enable.dsl=2
process CLUMP_PLINK19 {
  tag "clump_chr${params.chr}_${params.pop}${params.n_samples}"
  container 'quay.io/biocontainers/plink:1.90b6.21--hec16e2b_2'
  publishDir "clump", mode: 'copy'
  input:
    path bed
    path bim
    path fam
    path gwas_all
    

  output:
    path "clump.clumped",  emit: clumped_raw
    path "lead_snps.tsv",  emit: lead_snps

  script:
  """
  # Prepare a minimal file for PLINK clumping (ensures required columns exist)
  awk 'BEGIN{FS=OFS="\\t"} NR==1{for(i=1;i<=NF;i++){h[\$i]=i}; print "CHR","BP","SNP","P"; next} {print \$h["CHR"],\$h["BP"],\$h["SNP"],\$h["P"]}' ${gwas_all} > clump_input.tsv

  plink \\
    --bed $bed --bim $bim --fam $fam \\
    --clump clump_input.tsv \\
    --clump-field P \\
    --clump-snp-field SNP \\
    --clump-p1 ${params.clump_p1} \\
    --clump-r2 ${params.clump_r2} \\
    --clump-kb ${params.clump_kb} \\
    --out clump

  # 3) Parse .clumped to tidy leads (robust to whitespace/indent)
if [ -s clump.clumped ]; then
  awk -f ${projectDir}/bin/parse_clumped.awk clump.clumped > lead_snps.tsv
else
  : > lead_snps.tsv
fi

  """

}