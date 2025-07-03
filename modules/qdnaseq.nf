
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
    val qdnaseq_bin_size

    output:
    path "${params.sampleid}_FilteredReadCounts_${qdnaseq_bin_size}kb.bed"
    path "${params.sampleid}_FilteredReadCounts_${qdnaseq_bin_size}kb.igv"
    path "${params.sampleid}_ReadCounts_${qdnaseq_bin_size}kb.png"
    path "${params.sampleid}_isobarPlot_${qdnaseq_bin_size}kb.png"
    path "${params.sampleid}_noisePlot_${qdnaseq_bin_size}kb.png"
    path "${params.sampleid}_copyNumbersSmooth_${qdnaseq_bin_size}kb.png"
    path "${params.sampleid}_copyNumbersSegmented_${qdnaseq_bin_size}kb.png"
    path "${params.sampleid}_copyNumbersCalled_${qdnaseq_bin_size}kb.png"
    path "${params.sampleid}_calledCNVs_${qdnaseq_bin_size}kb.vcf", emit: cnv_calls
    path "${params.sampleid}_calledCNVs_${qdnaseq_bin_size}kb.igv"
    path "${params.sampleid}_copyNumbersSegmented_perChr_${qdnaseq_bin_size}kb.pdf"

    //script:

    shell:
    '''
    #!/usr/bin/env Rscript
    library(QDNAseq)
    library(Biobase)
    library("QDNAseq.hg38")

    sorted_bam <- \"!{sorted_bam}\"
    sample_name <- \"!{params.sampleid}\"
    qdnaseq_bin_size <- as.numeric(\"!{qdnaseq_bin_size}\")

    options(future.globals.maxSize= 1048576000)

    bins <- getBinAnnotations(binSize=qdnaseq_bin_size, genome="hg38")
    readCounts <- binReadCounts(bins, bamfiles=sorted_bam,chunkSize=100000)

    png(filename=paste0(sample_name,"_ReadCounts_",qdnaseq_bin_size, "kb.png"))
    plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
    dev.off()

    readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE) 

    png(filename=paste0(sample_name,"_isobarPlot_",qdnaseq_bin_size, "kb.png"))
    isobarPlot(readCountsFiltered)
    dev.off()

    readCountsFiltered <- estimateCorrection(readCountsFiltered)
    exportBins(readCountsFiltered, file=paste0(sample_name,"_FilteredReadCounts_",qdnaseq_bin_size,"kb.bed"), format="bed")
    exportBins(readCountsFiltered, file=paste0(sample_name,"_FilteredReadCounts_",qdnaseq_bin_size,"kb.igv"), format="igv")

    png(filename=paste0(sample_name,"_noisePlot_",qdnaseq_bin_size, "kb.png"))
    noisePlot(readCountsFiltered)
    dev.off()

    copyNumbers <- correctBins(readCountsFiltered)

    copyNumbersNormalized <- normalizeBins(copyNumbers)

    copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

    png(filename=paste0(sample_name,"_copyNumbersSmooth_",qdnaseq_bin_size, "kb.png"))
    plot(copyNumbersSmooth)
    dev.off()

    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")

    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

    png(filename=paste0(sample_name,"_copyNumbersSegmented_",qdnaseq_bin_size, "kb.png"))
    plot(copyNumbersSegmented)
    dev.off()

    # Plot CNVs per chromosome
    pdf(file=paste0(sample_name,"_copyNumbersSegmented_perChr_",qdnaseq_bin_size, "kb.pdf"))
    f.data <- as(featureData(copyNumbersSegmented), "data.frame")
    for (chromosome in unique(f.data$chromosome)) 
    {
		select <- which(f.data$chromosome == chromosome)
        if (!grepl("X",chromosome) & !grepl("Y",chromosome) & !grepl("Un",chromosome) & !grepl("EBV",chromosome))
        {
            print(paste("Plotting chromosome:", chromosome, sep = " "))
            plot(copyNumbersSegmented[select, 1], ylim=c(-2,2), main=paste0(" Chr",chromosome,"_",qdnaseq_bin_size,"kb"))
        }
	}
	dev.off()

    copyNumbersCalled <- callBins(copyNumbersSegmented)

    png(filename=paste0(sample_name,"_copyNumbersCalled_",qdnaseq_bin_size, "kb.png"))
    plot(copyNumbersCalled)
    dev.off()

    exportBins(copyNumbersCalled, paste0(sample_name,"_calledCNVs_",qdnaseq_bin_size, "kb.vcf"), format="vcf")
    exportBins(copyNumbersCalled, file=paste0(sample_name,"_calledCNVs_",qdnaseq_bin_size,"kb.igv"), format="igv")

    '''
}