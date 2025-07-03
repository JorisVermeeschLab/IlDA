
/* 
* Mark duplicates
*/
process mark_duplicates {
    label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'gatk'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    val sample
    path bam
    path bam_index

    output:
    path "${sample}.marked.bam", emit: dedupped_bam
    path "${sample}.marked.bai", emit: dedupped_bam_index
    path "markDuplicates.command.log", emit: runlog

    script:
    """
    gatk MarkDuplicates --INPUT $bam --METRICS_FILE ${sample}_markDuplicates.metrics --CREATE_INDEX true --OUTPUT ${sample}.marked.bam --TMP_DIR ${workDir}
    cp .command.log markDuplicates.command.log
    """
}


/* 
* Collect WGS metrics
*/
process collect_wgs_metrics {
    label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'gatk'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path bam
    path bam_index
    path genomeref

    output:
    path "${params.sampleid}_CollectWgsMetrics.metrics"
    path "collectwgsmetrics.command.log"

    script:
    """
    gatk CollectWgsMetrics --INPUT $bam --OUTPUT ${params.sampleid}_CollectWgsMetrics.metrics --REFERENCE_SEQUENCE $genomeref --TMP_DIR ${workDir}
    cp .command.log collectwgsmetrics.command.log
    """
}


/* 
* Collect insert size metrics
*/
process collect_insert_size_metrics {
    label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'gatk'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path bam
    path bam_index

    output:
    path "${params.sampleid}_CollectInsertSizeMetrics.metrics"
    path "${params.sampleid}_CollectInsertSizeMetrics.png"
    path "collectinsertsizemetrics.command.log"

    script:
    """
    gatk CollectInsertSizeMetrics --INPUT $bam --OUTPUT ${params.sampleid}_CollectInsertSizeMetrics.metrics --Histogram_FILE ${params.sampleid}_CollectInsertSizeMetrics.png --TMP_DIR ${workDir}
    cp .command.log collectinsertsizemetrics.command.log
    """
}



/* 
* Collect alignment summary metrics
*/
process collect_alignment_summary_metrics {
    label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'gatk'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path bam
    path bam_index

    output:
    path "${params.sampleid}_CollectAlignmentSummary.metrics"
    path "collectalignmentsummarymetrics.command.log"

    script:
    """
    gatk CollectAlignmentSummaryMetrics --INPUT $bam --OUTPUT ${params.sampleid}_CollectAlignmentSummary.metrics --TMP_DIR ${workDir}
    cp .command.log collectalignmentsummarymetrics.command.log
    """
}
