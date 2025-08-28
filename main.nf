/*
============================================
                I l D A
============================================
 Illumina Data Analysis Pipeline
--------------------------------------------
*/

nextflow.enable.dsl=2

// Command line shortcuts, quick entry point
include { clip_reads as clip_R1; clip_reads as clip_R2; merge_fastq_files as merge_R1_files; merge_fastq_files as merge_R2_files; get_baf; get_roh_rg } from './modules/utilities'
include { parallel_gzip as gzip_R1; parallel_gzip as gzip_R2; parallel_gzip as gzip_flash } from './modules/pigz'
include { flash } from './modules/flash'
include { bwa_mem_paired_end; bwa_mem_single_end } from './modules/bwa'
include { sam_to_sorted_bam; index_bam } from './modules/samtools'
include { mark_duplicates; collect_wgs_metrics; collect_insert_size_metrics; collect_alignment_summary_metrics } from './modules/gatk'
include { mosdepth } from './modules/mosdepth'
include { deepvariant } from './modules/deepvariant'
include { bcftools_roh; bcftools_roh_viz } from './modules/bcftools'
include { plot_baf_roh } from './modules/karyoploter'
include { qdnaseq as qdnaseq_5kb; qdnaseq as qdnaseq_10kb; qdnaseq as qdnaseq_50kb; qdnaseq as qdnaseq_100kb } from './modules/qdnaseq'

/*
* Main workflow
*/
workflow {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )
    genomeref_index = Channel.fromPath( params.genomeref_index, checkIfExists: true )
    genomeref_amb = Channel.fromPath( params.genomeref_amb, checkIfExists: true )
    genomeref_ann = Channel.fromPath( params.genomeref_ann, checkIfExists: true )
    genomeref_bwt = Channel.fromPath( params.genomeref_bwt, checkIfExists: true )
    genomeref_pac = Channel.fromPath( params.genomeref_pac, checkIfExists: true )
    genomeref_sa = Channel.fromPath( params.genomeref_sa, checkIfExists: true )

    // Initialise variables based on config file //
    def perform_clipping = "no"
    def perform_merging_files = "no"
    def perform_merging_reads = "no"
    def perform_mapping = "no"
    def perform_snv_calling = "no"
    def perform_roh_calling = "no"
    def perform_cnv_calling = "no"

    // Initialise analysis variables based on config file - to be corrected based on information available in config file
    if (params.clipped_length != "") 
    {
        perform_clipping = "yes"
    }
    if (params.merge_fastq_files == "yes") 
    {
        perform_merging_files = "yes"
    }
    if (params.merge_reads == "yes")
    {
        perform_merging_reads = "yes"
    }
    if (params.mapping == "yes")
    {
        perform_mapping = "yes"
    }
    if (params.snv_calling == "yes")
    {
        perform_snv_calling = "yes"
    }
    if (params.roh_calling == "yes")
    {
        perform_roh_calling = "yes"
    }
    if (params.cnv_calling == "yes")
    {
        perform_cnv_calling = "yes"
    }
   
    // Check input files
    if (params.r1_fastq != "" && params.r2_fastq != "")
    {
        r1_fastq = Channel.fromPath( params.r1_fastq + "*.fastq.gz", followLinks: true, checkIfExists: true ).collect().sort()
        r2_fastq = Channel.fromPath( params.r2_fastq + "*.fastq.gz", followLinks: true, checkIfExists: true ).collect().sort()
    }
    // Get BAM file from config file if provided
    else if (params.bam != "" && params.snv_indel_vcf != "")
    {
        perform_mapping = "no"
        bam = Channel.fromPath( params.bam, followLinks: true, checkIfExists: true ).collect()
        // Index BAM file if necessary
        if (params.bam_index == "")
        {
            index_bam( bam )
            bam_index = index_bam.out.bam_index
        }
        else 
        {
            bam_index = Channel.fromPath( params.bam_index, followLinks: true, checkIfExists: true ).collect()
        }
        mosdepth(bam, bam_index)
        /*collect_wgs_metrics(bam, bam_index, genomeref)
        collect_alignment_summary_metrics(bam, bam_index)
        collect_insert_size_metrics(bam, bam_index)*/
        perform_snv_calling = "no"
        snv_indel_vcf = Channel.fromPath( params.snv_indel_vcf, followLinks: true, checkIfExists: true ).collect()
    }
    else if (params.bam != "")
    {
        perform_mapping = "no"
        bam = Channel.fromPath( params.bam, followLinks: true, checkIfExists: true ).collect()
        // Index BAM file if necessary
        if (params.bam_index == "")
        {
            index_bam( bam )
            bam_index = index_bam.out.bam_index
        }
        else 
        {
            bam_index = Channel.fromPath( params.bam_index, followLinks: true, checkIfExists: true ).collect()
        }
        //mosdepth(bam, bam_index)
        collect_wgs_metrics(bam, bam_index, genomeref)
        collect_alignment_summary_metrics(bam, bam_index)
        collect_insert_size_metrics(bam, bam_index)
    }
    else if (params.snv_indel_vcf != "")
    {
        perform_snv_calling = "no"
        snv_indel_vcf = Channel.fromPath( params.snv_indel_vcf, followLinks: true, checkIfExists: true ).collect()
    }
    else 
    {
        error "No input files provided. Please provide either R1 and R2 FASTQ files, a BAM file or a SNV/indel VCF file."
    }

    // Print assigned variables to stdout
    println "Perform clipping: $perform_clipping"
    println "Perform merging FASTQ files: $perform_merging_files"
    println "Perform merging paired end reads: $perform_merging_reads"
    println "Perform mapping: $perform_mapping"
    println "Perform SNV calling: $perform_snv_calling"
    println "Perform ROH calling: $perform_roh_calling"
    println "Perform CNV calling (hg38 only): $perform_cnv_calling"

    /* Start FASTQ processing */
    // Perform pre-processing
    if (perform_clipping == "yes")
    {
        clip_R1(r1_fastq, params.clipped_length)
        clip_R2(r2_fastq, params.clipped_length)
        r1_fastq = clip_R1.out.clipped_fastq
        r2_fastq = clip_R2.out.clipped_fastq
        gzip_R1(r1_fastq)
        gzip_R2(r2_fastq)
    }
    // Merging automatically done if clipping is done
    if (perform_merging_files == "yes")
    {
        merge_R1_files(r1_fastq,"R1")
        merge_R2_files(r2_fastq,"R2")
        r1_fastq = merge_R1_files.out.merged_fastq
        r2_fastq = merge_R2_files.out.merged_fastq
        gzip_R1(r1_fastq)
        gzip_R2(r2_fastq)
    }
    if (perform_merging_reads == "yes")
    {
        flash(r1_fastq, r2_fastq)
        fastq = flash.out.merged_reads
        gzip_flash
    }
    // Perform mapping
    if (perform_mapping == "yes") 
    {
        if (perform_merging_reads == "yes")
        {
            bwa_mem_single_end(params.sampleid, genomeref, genomeref_amb, genomeref_ann, genomeref_bwt, genomeref_pac, genomeref_sa, fastq)
            sam = bwa_mem_single_end.out.mapped_sam
            sam_log = bwa_mem_single_end.out.runlog
        }
        else 
        {
            bwa_mem_paired_end(params.sampleid, genomeref, genomeref_amb, genomeref_ann, genomeref_bwt, genomeref_pac, genomeref_sa, r1_fastq, r2_fastq)
            sam = bwa_mem_paired_end.out.mapped_sam
            sam_log = bwa_mem_paired_end.out.runlog
        }
        sam_to_sorted_bam(sam, genomeref, sam_log)
        mark_duplicates(params.sampleid, sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index)
        bam = mark_duplicates.out.dedupped_bam
        bam_index = mark_duplicates.out.dedupped_bam_index
        // Perform QC
        /*nanoplot_bam(minimap_align_bamout.out.bam, minimap_align_bamout.out.idx)*/
        mosdepth(bam, bam_index)
        collect_wgs_metrics(bam, bam_index, genomeref)
        collect_alignment_summary_metrics(bam, bam_index)
        if (perform_merging_reads == "no")
        {
            collect_insert_size_metrics(bam, bam_index)
        }
    }
    
    /* Perform SNV and SV calling */
    if (perform_snv_calling == "yes")
    {
        deepvariant(bam, bam_index, genomeref, genomeref_index)
        snv_indel_vcf = deepvariant.out.snv_indel_vcf
    }
    if (perform_roh_calling == "yes")
    {
        bcftools_roh(snv_indel_vcf)
        get_roh_rg(bcftools_roh.out.roh)
        bcftools_roh_viz(snv_indel_vcf,bcftools_roh.out.roh)
        get_baf(snv_indel_vcf)
        plot_baf_roh(get_baf.out.vcf_baf, get_roh_rg.out.roh_rg, mosdepth.out.mosdepth_50kb_regions)
    }

    if (perform_cnv_calling == "yes")
    {
        qdnaseq_5kb(bam, bam_index,5)
        qdnaseq_10kb(bam, bam_index,10)
        qdnaseq_50kb(bam, bam_index,50)
        qdnaseq_100kb(bam, bam_index,100)
    }
}


