
/* 
* Convert SAM to sorted BAM using samtools and compute stats
*/
process sam_to_sorted_bam {
    label 'private_node'
    label 'cpu_mid'
    label 'mem_high'
    label 'time_mid'
    label 'samtools'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy' , pattern: "*{stats,log}"
    
    input:
    path mapped_sam
    path genomeref
    path bwalog

    output:
    path "${mapped_sam.getSimpleName()}_sorted.bam", emit: sorted_bam
    path "${mapped_sam.getSimpleName()}_sorted.bam.bai", emit: bam_index
    path "${mapped_sam.getSimpleName()}_sorted.bam.flagstats", emit: flagstats
    path "${mapped_sam.getSimpleName()}_sorted.bam.idxstats", emit: idxstats
    path "${mapped_sam.getSimpleName()}_sorted.bam.stats", emit: stats
    path bwalog, emit: bwalog

    script:
    """
    samtools sort  \
	    -@ $task.cpus \
        --write-index \
	    -m 2G \
        -o ${mapped_sam.getSimpleName()}_sorted.bam##idx##${mapped_sam.getSimpleName()}_sorted.bam.bai \
        --reference $genomeref \
        -T sorttmp_${mapped_sam.getSimpleName()}_sorted \
        $mapped_sam
    samtools flagstat ${mapped_sam.getSimpleName()}_sorted.bam > ${mapped_sam.getSimpleName()}_sorted.bam.flagstats
    samtools idxstats -@ 4 ${mapped_sam.getSimpleName()}_sorted.bam > ${mapped_sam.getSimpleName()}_sorted.bam.idxstats
    samtools stats -@ 4 ${mapped_sam.getSimpleName()}_sorted.bam > ${mapped_sam.getSimpleName()}_sorted.bam.stats
    """
}


/*
* Index BAM file with samtools
*/
process index_bam {
    label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'samtools'

    input:
    path mapped_bam

    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bam_index

    """
    samtools index -b -@ $task.cpus $mapped_bam ${mapped_bam}.bai
    """
}

