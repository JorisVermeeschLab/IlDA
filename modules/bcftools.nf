
/* 
* Compute ROH with bcftools
*/
process bcftools_roh {
    //label 'private_node'
    label 'cpu_mid'
    label 'mem_high'
    label 'time_mid'
    label 'bcftools'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'
    
    input:
    path vcf

    output:
    path "${vcf.getSimpleName()}_roh.tsv", emit: roh
    path bcftools_roh_log

    script:
    """
    bcftools roh -G${params.roh_G} --AF-dflt ${params.roh_AF} $vcf > ${vcf.getSimpleName()}_roh.tsv
    cp .command.log bcftools_roh_log
    """
}


/* 
* Visualise ROH with bcftools
*/
process bcftools_roh_viz {
    //label 'private_node'
    label 'cpu_mid'
    label 'mem_high'
    label 'time_mid'
    label 'bcftools'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'
    
    input:
    path vcf
    path roh

    output:
    path "${vcf.getSimpleName()}_roh.html", emit: roh
    path bcftools_rohviz_log

    script:
    """
    /bcftools/misc/roh-viz -i $roh -v $vcf -o ${vcf.getSimpleName()}_roh.html
    cp .command.log bcftools_rohviz_log
    """
}