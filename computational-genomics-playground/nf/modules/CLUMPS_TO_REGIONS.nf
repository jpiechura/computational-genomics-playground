nextflow.enable.dsl=2
process CLUMPS_TO_REGIONS {
  tag "clump2reg_chr${params.chr}_${params.pop}${params.n_samples}"
  label "python311"
  publishDir "${params.run_outdir}/clump", mode: 'copy'
  input:
    path lead_snps
  output:
    path "loci.tsv", emit: loci

script:
def pad_bp = params.containsKey('pad_kb') ? (params.pad_kb as long) * 1000L : 0L
    """
    set -euo pipefail
    clumps_to_regions.py ${lead_snps} ${pad_bp}
    """
}