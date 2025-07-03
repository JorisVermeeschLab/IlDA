
/* 
* Map with BWA-mem
*/
process bwa_mem {
    label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'bwa'

    //publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    val sample
    path genomeref
    path genomeref_amb
    path genomeref_ann
    path genomeref_bwt
    path genomeref_pac
    path genomeref_sa
    path r1_fastq
    path r2_fastq

    output:
    path "${sample}.sam", emit: mapped_sam
    path "bwa.command.log", emit: runlog

    script:
    """
    bwa mem -t $task.cpus -R \'@RG\\tID:RG\\tSM:${sample}\\tPL:Illumina\\tLB:WGS\' $genomeref $r1_fastq $r2_fastq > ${sample}.sam
    cp .command.log bwa.command.log
    """
}
