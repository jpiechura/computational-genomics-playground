/*
 * main.nf — GWAS demo pipeline: fetch → sample list → subset → plink → gcta simulate
 * DSL2, scripts staged from workflow (no `from file(...)` inside process inputs)
 */

nextflow.enable.dsl = 2

/* create the output directory for the run */
params.run_outdir = "${params.outdir}/${params.runname}"


include { FETCH_1KG } from './nf/modules/FETCH_1KG.nf' /* fixed */
include { MAKE_SAMPLE_LIST } from './nf/modules/MAKE_SAMPLE_LIST.nf' /* fixed */
include { SUBSET_VCF } from './nf/modules/SUBSET_VCF.nf' /* fixed */
include { PLINKIFY } from './nf/modules/PLINKIFY.nf' /* fixed */
include { GCTA_SIM } from './nf/modules/GCTA_SIM.nf' /* fixed */
include { QC_PLINK } from './nf/modules/QC_PLINK.nf'  /* fixed */
include { QC_HET } from './nf/modules/QC_HET.nf' /* fixed */
include { QC_HETDIF_FILT } from './nf/modules/QC_HETDIF_FILT.nf' /* fixed */
include { PCA_PLINK } from './nf/modules/PCA_PLINK.nf' /* fixed */
include { BUILD_COVAR } from './nf/modules/BUILD_COVAR.nf' /* fixed */
include { GWAS_PLINK } from './nf/modules/GWAS_PLINK.nf' /* fixed */
include { GWAS_SUM } from './nf/modules/GWAS_SUM.nf' /* fixed */
include { PLOT_QQ_MANHATTAN } from './nf/modules/PLOT_QQ_MANHATTAN.nf' /* fixed */
include { CLUMP_PLINK19 } from './nf/modules/CLUMP_PLINK19.nf' /* fixed */
include { CLUMPS_TO_REGIONS } from './nf/modules/CLUMPS_TO_REGIONS.nf' /* fixed */
include { EXTRACT_LOCSTATS } from './nf/modules/EXTRACT_LOCSTATS.nf' /* fixed */
include { MAKE_REGION_BFILE } from './nf/modules/MAKE_REGION_BFILE.nf' /* fixed */
include { COMPUTE_LD } from './nf/modules/COMPUTE_LD.nf' /* fixed */
include { PREP_ZSCORES } from './nf/modules/PREP_ZSCORES.nf' /* fixed */
include { RUN_SUSIE_RSS } from './nf/modules/RUN_SUSIE_RSS.nf' /* fixed */


workflow.onComplete {
  println "Done."
}




/* -------------------- Wiring (workflow) --------------------- */
workflow {
  log.info "Starting pipeline..."



  // Stage 1: Fetch data
  def fetched = FETCH_1KG( file('bin/00_fetch_1000g.sh') )

  // Stage 2: Create sample list
  def picked = MAKE_SAMPLE_LIST(
    fetched.panel,
    file('bin/01_make_sample_list.sh')
  )

  // Stage 3: Subset VCF
  def subset = SUBSET_VCF(
    fetched.vcf,
    picked.samples,
    file('bin/02_subset_vcf.sh')
  )

  // Stage 4: Convert to PLINK format
  def plinked = PLINKIFY(
    subset.vcf_sub,
    file('bin/03_plinkify.sh')
  )

  // Stage 5: Run GCTA simulation

  def gcta_results = GCTA_SIM(
    plinked.bed,
    plinked.bim,
    plinked.fam,
  )

  def plink_qc = QC_PLINK(
    plinked.bed,
    plinked.bim,
    plinked.fam,
    gcta_results.pheno_file
  )


  def qc_het = QC_HET(
    plink_qc.het_stats,
  )

  def qc_hetfilt = QC_HETDIF_FILT(
    qc_het.het_outliers,
    gcta_results.pheno_file,
    plink_qc.bed,
    plink_qc.bim,
    plink_qc.fam
  )

  def pca_results = PCA_PLINK(
    qc_hetfilt.bed,
    qc_hetfilt.bim,
    qc_hetfilt.fam
  )

Channel.fromPath('bin/append_sex.awk').set { append_sex_script }
  covar_file = BUILD_COVAR(
    pca_results.eigenvec,
    fetched.panel,
    append_sex_script
  )


  def gwas_results = GWAS_PLINK(
    qc_hetfilt.bed,
    qc_hetfilt.bim,
    qc_hetfilt.fam,
    gcta_results.pheno_file,
    covar_file.covar_file
  )

  def gwas_summary = GWAS_SUM(
    gwas_results.plink_output,
  )

  def gwas_plots = PLOT_QQ_MANHATTAN(
    gwas_summary.gwas_tsv
  )

  def clumped = CLUMP_PLINK19(
    qc_hetfilt.bed,
    qc_hetfilt.bim,
    qc_hetfilt.fam,
    gwas_summary.gwas_tsv
  )

    def regions = CLUMPS_TO_REGIONS(
        clumped.lead_snps
    )

    Channel
        regions.loci
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            if (row.containsKey('start') && row.containsKey('end')) {
                tuple(row.locus_id, row.chr as String, row.start as long, row.end as long)
            } else {
                // lead-SNP row: compute window
                def half = (params.window_kb as long) * 1000L
                tuple(row.locus_id, row.chr as String, (row.bp as long) - half, (row.bp as long) + half)
            }
        }
        .set { loci_ch }

    // 1) subset GWAS to region
    def locstats = EXTRACT_LOCSTATS(loci_ch, gwas_summary.gwas_tsv)

    def region_bfiles = MAKE_REGION_BFILE(
        locstats,
        qc_hetfilt.bed,
        qc_hetfilt.bim,
        qc_hetfilt.fam
    )

    def ld_matrices = COMPUTE_LD(
        region_bfiles
    )
    //make z prep inputs
    // you want val(locus_id), path(zfile), path(prune_in), path(pruned_bim)
    //ld matrices has  tuple val(locus_id),path("ld.${locus_id}.ld"),path("ld.${locus_id}.snplist"),path("region_${locus_id}.pruned.bim"),
    //pruning happens in make_region_bfile
    // Extract (id, sumstats)

    sumstats_ch = locstats.map { id, chr, start, end, sfile -> tuple(id, sfile) }

// (id, prune_in, pruned_bim)
pruned_index_ch = region_bfiles.map { id, bed, pbim, pfam, prune_in ->
  tuple(id, prune_in, pbim)
}

// Join on id → (id, zfile, prune_in, pruned_bim)
prep_inputs = sumstats_ch
  .join(pruned_index_ch)
  .map { id, zfile, prune_in, pruned_bim -> tuple(id, zfile, prune_in, pruned_bim) }



    def z_scores = PREP_ZSCORES(
        prep_inputs
    )
    

  susie_inputs = z_scores
  .map { id, zfile      -> tuple(id, zfile) }   // keep key
  .join(  ld_matrices .map { id, rfile  -> tuple(id, rfile) } )
  .join(  region_bfiles .map { id, x, y, z, keepfile-> tuple(id, keepfile) } )


    def susie = RUN_SUSIE_RSS(
        susie_inputs
    )
        
}