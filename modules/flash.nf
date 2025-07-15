
/* 
* Merge overlapping paired-end reads
*/
process flash {
    //label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'flash2'

    //publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path r1_fastq
    path r2_fastq

    output:
    path "${params.sampleid}.merged.fastq", emit: merged_reads
    path "flash.command.log", emit: runlog
    
    script:
    """
    flash2 $r1_fastq $r2_fastq -c > ${params.sampleid}.merged.fastq
    cp .command.log flash.command.log
    """
}
