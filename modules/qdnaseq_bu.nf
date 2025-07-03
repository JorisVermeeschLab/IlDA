
/* 
* Copy number calling on aligned reads using QDNAseq (only for hg38)
*/
process qdnaseq_hg38 {
    label 'private_node'
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'qdnaseq'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index

    output:
    path "${params.sampleid}_ReadCounts_5kb.png"
    path "${params.sampleid}_isobarPlot_5kb.png"
    path "${params.sampleid}_noisePlot_5kb.png"
    path "${params.sampleid}_copyNumbersSmooth_5kb.png"
    path "${params.sampleid}_copyNumbersSegmented_5kb.png"
    path "${params.sampleid}_copyNumbersCalled_5kb.png"
    path "${params.sampleid}_CNVcalled_5kb.vcf", emit: cnv_calls

    script:

shell:
    '''
    #!/usr/bin/env Rscript
    library(QDNAseq)
    library(Biobase)
    library("QDNAseq.hg38")

    sorted_bam <- \"!{sorted_bam}\"
    sample_name <- \"!{params.sampleid}\"

    options(future.globals.maxSize= 1048576000)

    bins <- getBinAnnotations(binSize=5, genome="hg38")
    readCounts <- binReadCounts(bins, bamfiles=sorted_bam,chunkSize=100000)

    png(filename=paste0(sample_name,"_ReadCounts_5kb.png"))
    plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
    dev.off()

    readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

    png(filename=paste0(sample_name,"_isobarPlot_5kb.png"))
    isobarPlot(readCountsFiltered)
    dev.off()

    readCountsFiltered <- estimateCorrection(readCountsFiltered)

    png(filename=paste0(sample_name,"_noisePlot_5kb.png"))
    noisePlot(readCountsFiltered)
    dev.off()

    copyNumbers <- correctBins(readCountsFiltered)

    copyNumbersNormalized <- normalizeBins(copyNumbers)

    copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

    png(filename=paste0(sample_name,"_copyNumbersSmooth_5kb.png"))
    plot(copyNumbersSmooth)
    dev.off()

    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")

    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

    png(filename=paste0(sample_name,"_copyNumbersSegmented_5kb.png"))
    plot(copyNumbersSegmented)
    dev.off()

    copyNumbersCalled <- callBins(copyNumbersSegmented)

    png(filename=paste0(sample_name,"_copyNumbersCalled_5kb.png"))
    plot(copyNumbersCalled)
    dev.off()

    exportBins(copyNumbersCalled, paste0(sample_name,"_CNVcalled_5kb.vcf"), format="vcf")
    '''
}



/* 
* Copy number calling on aligned reads using QDNAseq (only for hg38)
*/
process qdnaseq_equCab3 {
    label 'private_node'
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'qdnaseq'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index

    output:
    path "${params.sampleid}_ReadCounts_5kb.png"
    path "${params.sampleid}_isobarPlot_5kb.png"
    path "${params.sampleid}_noisePlot_5kb.png"
    path "${params.sampleid}_copyNumbersSmooth_5kb.png"
    path "${params.sampleid}_copyNumbersSegmented_5kb.png"
    path "${params.sampleid}_copyNumbersCalled_5kb.png"
    path "${params.sampleid}_CNVcalled_5kb.vcf", emit: cnv_calls

    script:

shell:
    '''
    #!/usr/bin/env Rscript
    library(QDNAseq)
    library(Biobase)
    library("QDNAseq.hg38")

    sorted_bam <- \"!{sorted_bam}\"
    sample_name <- \"!{params.sampleid}\"

    options(future.globals.maxSize= 1048576000)

    bins <- getBinAnnotations(binSize=5, genome="hg38")
    readCounts <- binReadCounts(bins, bamfiles=sorted_bam,chunkSize=100000)

    png(filename=paste0(sample_name,"_ReadCounts_5kb.png"))
    plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
    dev.off()

    readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

    png(filename=paste0(sample_name,"_isobarPlot_5kb.png"))
    isobarPlot(readCountsFiltered)
    dev.off()

    readCountsFiltered <- estimateCorrection(readCountsFiltered)

    png(filename=paste0(sample_name,"_noisePlot_5kb.png"))
    noisePlot(readCountsFiltered)
    dev.off()

    copyNumbers <- correctBins(readCountsFiltered)

    copyNumbersNormalized <- normalizeBins(copyNumbers)

    copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

    png(filename=paste0(sample_name,"_copyNumbersSmooth_5kb.png"))
    plot(copyNumbersSmooth)
    dev.off()

    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")

    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

    png(filename=paste0(sample_name,"_copyNumbersSegmented_5kb.png"))
    plot(copyNumbersSegmented)
    dev.off()

    copyNumbersCalled <- callBins(copyNumbersSegmented)

    png(filename=paste0(sample_name,"_copyNumbersCalled_5kb.png"))
    plot(copyNumbersCalled)
    dev.off()

    exportBins(copyNumbersCalled, paste0(sample_name,"_CNVcalled_5kb.vcf"), format="vcf")
    '''
}