nextflow.enable.dsl=2

process EXTRACT_LOCSTATS {
  tag "extract_${locus_id}"
  container 'quay.io/biocontainers/htslib:1.19--h81da01d_0'
  publishDir "finemap", mode: 'copy'
  input:
    input:
      tuple val(locus_id), val(chr), val(start), val(end)
      path gwas
    
    output:
      tuple val(locus_id), val(chr), val(start), val(end), path("sumstats.${locus_id}.tsv"), emit: loc_sumstats

    script:
    """

set -euo pipefail
IN="${gwas}"

if gzip -t "\$IN" >/dev/null 2>&1; then zcat "\$IN"; else cat "\$IN"; fi \
| awk -v c=${chr} -v s=${start} -v e=${end} 'BEGIN{FS=OFS="\\t"} NR==1 || (\$1==c && \$2>=s && \$2<=e)' \
> sumstats.${locus_id}.tsv


    """ 

}
    
  