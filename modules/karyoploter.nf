
/* 
* Plot BAF with ROH regions
*/
process plot_baf_roh {
    //label 'private_node'
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'karyoploter'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path baf
    path roh

    output:
    path "${params.sampleid}_ROH_L1M.png"
    path "${params.sampleid}_ROH_L1M_BAF_perChr.pdf"

    shell:
    '''
    #!/usr/bin/env Rscript
    library(karyoploteR)

    vaf <- \"!{baf}\"
    roh <- \"!{roh}\"
    sample_name <- \"!{params.sampleid}\"
    genome <- \"!{params.genome}\"
    
    # Settings for plots
    transparency <- 0.85
    trans_sample <- transparent("royalblue",amount=transparency)

    if (genome == "equCab3") {
        horse_cytoband<-read.csv(file="/Horse/equCab3_cytoband.tsv",header=TRUE,sep="\t")
        custom.genome <- toGRanges(horse_cytoband)
        chrs<-c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26","chr27","chr28","chr29","chr30","chr31","chrX")
    } else if (genome == "hg38") {
        custom.genome <- "hg38"
        chrs<-c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
    } else {
        stop("Unsupported genome. Please use 'equCab3' or 'hg38'.")
    }
    
    # Read ROH regions
    rg<-read.csv(file=roh,header=TRUE,sep="\t")
    # Keep only regions longer than 1Mbp
    rg_subset<-rg[as.numeric(rg$Length)>=1000000,]
    rg_subset_bed<-rg_subset[,c("Chr","Start","End")]
    regions<-toGRanges(rg_subset_bed)

    # Plot ROH regions
    png(file=paste0(sample_name,"_ROH_L1M.png"),type="cairo",units="in", width=10, height=10, pointsize=14, res=1000)
    kp<-plotKaryotype(genome = custom.genome, chromosomes = chrs)
    kpPlotRegions(kp, data=regions,col="royalblue",r0=0,r1=0.5)
    dev.off()

    # Read BAF values
    baf_sample<-read.csv(file=vaf,header=FALSE,sep="\t")
    colnames(baf_sample) <- c("chr", "start","end","y") 

    # Make one PDF file with BAF plots per chromosome and ROH regions
    pdf(paste0(sample_name,"_ROH_L1M_BAF_perChr.pdf"), onefile = TRUE)
    for (chr in chrs) { 
        print(chr) 
        baf_sample_chr<-baf_sample[baf_sample$chr==chr,]
        dfBAF_sample <- makeGRangesFromDataFrame(baf_sample_chr, keep.extra.columns = TRUE)
        rg_subset_bed_chr<-rg_subset_bed[rg_subset_bed$Chr==chr,]	
        regions_chr<-toGRanges(rg_subset_bed_chr)
        k <- plotKaryotype(plot.type = 4, genome=custom.genome, chromosomes = chr)
        kpAddBaseNumbers(k, tick.len=4, add.units=TRUE, cex=0.4)
        kpPoints(k, data=dfBAF_sample, cex=0.3, ymin=0, ymax=1, col=trans_sample, pch=1)
        kpPlotRegions(k, data = regions_chr, col=transparent("lightblue1",amount=0.5),r1=1.01)
    }
    dev.off()

    '''
}