
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
    path mosdepth

    output:
    path "${params.sampleid}_DP_ROH_L1M.png"
    path "${params.sampleid}_DP_ROH_L1M_BAF_perChr.pdf"
    path "${params.sampleid}_DP_ROH.png"
    path "${params.sampleid}_DP_ROH_BAF_perChr.pdf"

    shell:
    '''
    #!/usr/bin/env Rscript
    library(karyoploteR)

    vaf <- \"!{baf}\"
    roh <- \"!{roh}\"
    mosdepth <- \"!{mosdepth}\"
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
    # Keep all regions
    rg_bed<-rg[,c("Chr","Start","End")]
    regions<-toGRanges(rg_bed)
    # Keep only regions longer than 1Mbp
    rg_subset<-rg[as.numeric(rg$Length)>=1000000,]
    rg_subset_bed<-rg_subset[,c("Chr","Start","End")]
    regions_subset<-toGRanges(rg_subset_bed)

    # Read BAF values
    baf_sample<-read.csv(file=vaf,header=FALSE,sep="\t")
    colnames(baf_sample) <- c("chr", "start","end","y") 

    # Read depth values by 50 kb bins
    cov<-read.csv(file=mosdepth,header=FALSE,sep="\t")
    colnames(cov) <- c("chr", "start","end","y")
    cov$norm_median<-cov$y/median(cov$y)
    cov$norm_median[cov$norm_median>4] <- 4
    cov$log_median<-log2(cov$norm_median)
    cov$log_median[cov$log_median<(-2)] <- -2
    coverage<-toGRanges(cov)

    # Plot DP & ROH regions (>1MP)
    png(file=paste0(sample_name,"_DP_ROH_L1M.png"),type="cairo",units="in", width=10, height=10, pointsize=14, res=1000)
    pp_2 <- getDefaultPlotParams(plot.type = 2)
    pp_2$data1min <- -2
    pp_2$data1max <- 2
    pp_2$data1height <- 350
    pp_2$data2height <- 50
    kp<-plotKaryotype(plot.type=2, genome = custom.genome, chromosomes = chrs, plot.params = pp_2, cex=0.6)
    kpAxis(kp, cex=0.3)
    kpPoints(kp, coverage, x=coverage$start, y=coverage$log_median, col="royalblue", cex=0.2, ymin=-2, ymax=2)
    kpPlotRegions(kp, data=regions_subset,col="black", data.panel=2)
    dev.off()

    # Plot DP & ROH regions (all))
    png(file=paste0(sample_name,"_DP_ROH.png"),type="cairo",units="in", width=10, height=10, pointsize=14, res=1000)
    pp_2 <- getDefaultPlotParams(plot.type = 2)
    pp_2$data1min <- -2
    pp_2$data1max <- 2
    pp_2$data1height <- 350
    pp_2$data2height <- 50
    kp<-plotKaryotype(plot.type=2, genome = custom.genome, chromosomes = chrs, plot.params = pp_2, cex=0.6)
    kpAxis(kp, cex=0.3)
    kpPoints(kp, coverage, x=coverage$start, y=coverage$log_median, col="royalblue", cex=0.2, ymin=-2, ymax=2)
    kpPlotRegions(kp, data=regions,col="black", data.panel=2)
    dev.off()

    # Make one PDF file with BAF and depth plots per chromosome and ROH regions (>1MP)
    pdf(paste0(sample_name,"_DP_ROH_L1M_BAF_perChr.pdf"), onefile = TRUE)
    for (chr in chrs) { 
        print(chr) 
        baf_sample_chr<-baf_sample[baf_sample$chr==chr,]
        dfBAF_sample <- makeGRangesFromDataFrame(baf_sample_chr, keep.extra.columns = TRUE)
        rg_subset_bed_chr<-rg_subset_bed[rg_subset_bed$Chr==chr,]	
        regions_chr<-toGRanges(rg_subset_bed_chr)
        coverage_chr<-cov[cov$chr==chr,]
        coverages_chr<-toGRanges(coverage_chr)

        pp <- getDefaultPlotParams(plot.type = 1)
        pp$data1height=600
        k <- plotKaryotype(plot.type = 1, genome=custom.genome, chromosomes = chr, plot.params = pp)
        kpAddBaseNumbers(k, tick.len=4, add.units=TRUE, cex=0.4)
        kpDataBackground(k, r0=0, r1=0.3, color="gray97")
        kpAxis(k, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.5)
        kpPoints(k, coverages_chr, x=coverages_chr$start, y=coverages_chr$log_median, col="royalblue", pch=16, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.3)
        kpAxis(k, r0=0.4, r1=1, cex=0.5)
        kpDataBackground(k, r0=0.4, r1=1, color="gray97")
        kpPoints(k, data=dfBAF_sample, cex=0.3, ymin=0, ymax=1, col=trans_sample, pch=1, r0=0.4, r1=1)
        kpPlotRegions(k, data = regions_chr, col=transparent("lightblue1",amount=0.8), border="lightblue1",r0=-0.01, r1=1.01)
    }
    dev.off()

    # Make one PDF file with BAF and depth plots per chromosome and ROH regions (all)
    pdf(paste0(sample_name,"_DP_ROH_BAF_perChr.pdf"), onefile = TRUE)
    for (chr in chrs) { 
        print(chr) 
        baf_sample_chr<-baf_sample[baf_sample$chr==chr,]
        dfBAF_sample <- makeGRangesFromDataFrame(baf_sample_chr, keep.extra.columns = TRUE)
        rg_bed_chr<-rg_bed[rg_bed$Chr==chr,]	
        regions_chr<-toGRanges(rg_bed_chr)
        coverage_chr<-cov[cov$chr==chr,]
        coverages_chr<-toGRanges(coverage_chr)

        pp <- getDefaultPlotParams(plot.type = 1)
        pp$data1height=600
        k <- plotKaryotype(plot.type = 1, genome=custom.genome, chromosomes = chr, plot.params = pp)
        kpAddBaseNumbers(k, tick.len=4, add.units=TRUE, cex=0.4)
        kpDataBackground(k, r0=0, r1=0.3, color="gray97")
        kpAxis(k, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.5)
        kpPoints(k, coverages_chr, x=coverages_chr$start, y=coverages_chr$log_median, col="royalblue", pch=16, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.3)
        kpAxis(k, r0=0.4, r1=1, cex=0.5)
        kpDataBackground(k, r0=0.4, r1=1, color="gray97")
        kpPoints(k, data=dfBAF_sample, cex=0.3, ymin=0, ymax=1, col=trans_sample, pch=1, r0=0.4, r1=1)
        kpPlotRegions(k, data = regions_chr, col=transparent("lightblue1",amount=0.8), border="lightblue1",r0=-0.01, r1=1.01)
    }
    dev.off()

    '''
}