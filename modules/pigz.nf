
/* 
* Process (filter + trim) fastq files with fastp
*/
process parallel_gzip {
    //label 'private_node'
    label 'cpu_mid'
    label 'mem_low'
    label 'time_mid'
    label 'pigz'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path fastq

    output:
    path "${fastq}.gz", emit: fastqgz

    script:
    """
    pigz \
        --best \
        --stdout \
        --processes $task.cpus \
        $fastq > ${fastq}.gz
    """
}


process parallel_unzip {
    //label 'private_node'
    label 'cpu_mid'
    label 'mem_low'
    label 'time_low'
    label 'pigz'

    //publishDir path: "${params.outdir}/${params.sampleid}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path fastq

    output:
    path "${params.sampleid}.unzipped.fastq", emit: unzipped_fastq

    script:
    """
    pigz \
        --best \
        --stdout \
        --processes $task.cpus \
        -dc $fastq > ${params.sampleid}.unzipped.fastq
    """
}