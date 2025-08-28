
/* 
* Compute coverage with mosdepth
*/
process mosdepth {
    //label 'private_node'
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_mid'
    label 'mosdepth'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'
    
    input:
    path sorted_bam
    path bam_index

    output:
    path "${params.sampleid}.mosdepth.summary.txt"
    path "${params.sampleid}.mosdepth.global.dist.txt"
    path "${params.sampleid}.quantized.bed.gz"
    path "${params.sampleid}.quantized.bed.gz.csi"
    path "${params.sampleid}.dist.html"
    path "${params.sampleid}.5kb.regions.bed.gz"
    path "${params.sampleid}.5kb.regions.bed.gz.csi"
    path "${params.sampleid}.5kb.mosdepth.region.dist.txt"
    path "${params.sampleid}.50kb.regions.bed.gz", emit: mosdepth_50kb_regions
    path "${params.sampleid}.50kb.regions.bed.gz.csi"
    path "${params.sampleid}.50kb.mosdepth.region.dist.txt"
    
    script:
    """
    export MOSDEPTH_Q0=NO_COVERAGE
    export MOSDEPTH_Q1=LOW_COVERAGE
    export MOSDEPTH_Q2=CALLABLE
    export MOSDEPTH_Q3=HIGH_COVERAGE
    mosdepth -n -q 0:1:10:50: ${params.sampleid} $sorted_bam
    mosdepth -n -b 5000 -m ${params.sampleid}.5kb $sorted_bam
    mosdepth -n -b 50000 ${params.sampleid}.50kb $sorted_bam
    python3 /mosdepth-0.3.11/scripts/plot-dist.py -o ${params.sampleid}.dist.html ${params.sampleid}.mosdepth.global.dist.txt
    """
}
