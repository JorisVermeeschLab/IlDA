/* 
* Clip reads with awk
*/
process clip_reads {
    label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'

    //publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path fastq
    val clipped_length

    output:
    path "${fastq.getSimpleName()}.clipped${clipped_length}.fastq", emit: clipped_fastq

    script:
    """
    zcat $fastq | awk -v clip=$clipped_length '{new=\$0;flag=NR%4;if(flag==0||flag==2){new=substr(\$0,1,clip)}; print(new)}' > ${fastq.getSimpleName()}.clipped${clipped_length}.fastq
    """
}


/* 
* Merge FASTQ files
*/
process merge_fastq {
    label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'

    //publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path fastqs
    val read
    
    output:
    path "${params.sampleid}.merged.${read}.fastq", emit: merged_fastq

    script:
    """
    zcat $fastqs > ${params.sampleid}.merged.${read}.fastq
    """
}
