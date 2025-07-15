
/* 
* SNV and indel calling on aligned reads using DeepVariant
*/
process deepvariant {
    //label 'private_node'
    label 'deepvariant'
    label 'time_high'
    label 'mem_high'
    //label 'time_low'
    //label 'mem_low'
    label 'cpu_low'
    label ( workflow.profile.contains('qsub') ? null: 'cpu_high' )

    stageInMode 'copy'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path bam
    path bam_index
    path genomeref
    path genomeref_index

    output:
    path "${params.sampleid}_deepvariant.vcf.gz", emit: snv_indel_vcf
    path "${params.sampleid}_deepvariant.gvcf.gz", emit: gvcf
    path "deepvariant.command.log", emit: runlog

    script:
    def localproc = ( workflow.profile.contains('qsub') ? 36: task.cpus )
    """
    export TMPDIR="${workDir}"
    run_deepvariant \
        --model_type ${params.deepvariant_model} \
        --ref $genomeref \
        --reads $bam \
        --output_vcf ${params.sampleid}_deepvariant.vcf.gz \
        --output_gvcf ${params.sampleid}_deepvariant.gvcf.gz \
        --num_shards ${localproc}
    cp .command.log deepvariant.command.log
    """
}
