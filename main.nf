/*
============================================
                I l D A
============================================
 Illumina Data Analysis Pipeline
--------------------------------------------
*/

nextflow.enable.dsl=2

// Command line shortcuts, quick entry point
include { clip_reads as clip_R1; clip_reads as clip_R2; merge_fastq as merge_R1; merge_fastq as merge_R2 } from './modules/utilities'
include { parallel_gzip as gzip_R1; parallel_gzip as gzip_R2 } from './modules/pigz'
include { bwa_mem } from './modules/bwa'
include { sam_to_sorted_bam; index_bam } from './modules/samtools'
include { mark_duplicates; collect_wgs_metrics; collect_insert_size_metrics; collect_alignment_summary_metrics } from './modules/gatk'
include { deepvariant } from './modules/deepvariant'
include { mosdepth } from './modules/mosdepth'
include { qdnaseq_hg38 as qdnaseq_hg38_5kb; qdnaseq_hg38 as qdnaseq_hg38_10kb; qdnaseq_hg38 as qdnaseq_hg38_50kb; qdnaseq_hg38 as qdnaseq_hg38_100kb } from './modules/qdnaseq'

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
    def perform_merging = "no"
    def perform_mapping = "no"
    def perform_snv_calling = "no"
    def perform_cnv_calling = "no"

    // Initialise analysis variables based on config file - to be corrected based on information available in config file
    if (params.clipped_length != "") 
    {
        perform_clipping = "yes"
    }
    if (params.merge_fastq != "") 
    {
        perform_merging = "yes"
    }
    if (params.mapping == "yes")
    {
        perform_mapping = "yes"
    }
    if (params.snv_calling == "yes")
    {
        perform_snv_calling = "yes"
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
    }

    // Print assigned variables to stdout
    println "Perform mapping: $perform_mapping"
    println "Perform SNV calling: $perform_snv_calling"
    
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
    else if (perform_merging == "yes")
    {
        merge_R1(r1_fastq,"R1")
        merge_R2(r2_fastq,"R2")
        r1_fastq = merge_R1.out.merged_fastq
        r2_fastq = merge_R2.out.merged_fastq
        gzip_R1(r1_fastq)
        gzip_R2(r2_fastq)
    }
    // Perform mapping
    if (perform_mapping == "yes") 
    {
        bwa_mem(params.sampleid, genomeref, genomeref_amb, genomeref_ann, genomeref_bwt, genomeref_pac, genomeref_sa, r1_fastq, r2_fastq)
        sam_to_sorted_bam(bwa_mem.out.mapped_sam, genomeref, bwa_mem.out.runlog)
        mark_duplicates(params.sampleid, sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index)
        bam = mark_duplicates.out.dedupped_bam
        bam_index = mark_duplicates.out.dedupped_bam_index
        // Perform QC
        /*nanoplot_bam(minimap_align_bamout.out.bam, minimap_align_bamout.out.idx)*/
        mosdepth(bam, bam_index)
        collect_wgs_metrics(bam, bam_index, genomeref)
        collect_insert_size_metrics(bam, bam_index)
        collect_alignment_summary_metrics(bam, bam_index)
    }
    
    /* Perform SNV and SV calling */
    if (perform_snv_calling == "yes")
    {
        deepvariant(bam, bam_index, genomeref, genomeref_index)
    }
    if (perform_cnv_calling == "yes")
    {
        qdnaseq_hg38_5kb(bam, bam_index,5)
        qdnaseq_hg38_10kb(bam, bam_index,10)
        qdnaseq_hg38_50kb(bam, bam_index,50)
        qdnaseq_hg38_100kb(bam, bam_index,100)
    }
}


