
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
    //path "${params.sampleid}_DP_ROH.png"
    //path "${params.sampleid}_DP_ROH_BAF_perChr.pdf"

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
    #png(file=paste0(sample_name,"_DP_ROH.png"),type="cairo",units="in", width=10, height=10, pointsize=14, res=1000)
    #pp_2 <- getDefaultPlotParams(plot.type = 2)
    #pp_2$data1min <- -2
    #pp_2$data1max <- 2
    #pp_2$data1height <- 350
    #pp_2$data2height <- 50
    #kp<-plotKaryotype(plot.type=2, genome = custom.genome, chromosomes = chrs, plot.params = pp_2, cex=0.6)
    #kpAxis(kp, cex=0.3)
    #kpPoints(kp, coverage, x=coverage$start, y=coverage$log_median, col="royalblue", cex=0.2, ymin=-2, ymax=2)
    #kpPlotRegions(kp, data=regions,col="black", data.panel=2)
    #dev.off()

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
    #pdf(paste0(sample_name,"_DP_ROH_BAF_perChr.pdf"), onefile = TRUE)
    #for (chr in chrs) { 
    #    print(chr) 
    #    baf_sample_chr<-baf_sample[baf_sample$chr==chr,]
    #    dfBAF_sample <- makeGRangesFromDataFrame(baf_sample_chr, keep.extra.columns = TRUE)
    #    rg_bed_chr<-rg_bed[rg_bed$Chr==chr,]	
    #    regions_chr<-toGRanges(rg_bed_chr)
    #    coverage_chr<-cov[cov$chr==chr,]
    #    coverages_chr<-toGRanges(coverage_chr)
    #    pp <- getDefaultPlotParams(plot.type = 1)
    #    pp$data1height=600
    #    k <- plotKaryotype(plot.type = 1, genome=custom.genome, chromosomes = chr, plot.params = pp)
    #    kpAddBaseNumbers(k, tick.len=4, add.units=TRUE, cex=0.4)
    #    kpDataBackground(k, r0=0, r1=0.3, color="gray97")
    #    kpAxis(k, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.5)
    #    kpPoints(k, coverages_chr, x=coverages_chr$start, y=coverages_chr$log_median, col="royalblue", pch=16, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.3)
    #    kpAxis(k, r0=0.4, r1=1, cex=0.5)
    #    kpDataBackground(k, r0=0.4, r1=1, color="gray97")
    #    kpPoints(k, data=dfBAF_sample, cex=0.3, ymin=0, ymax=1, col=trans_sample, pch=1, r0=0.4, r1=1)
    #    kpPlotRegions(k, data = regions_chr, col=transparent("lightblue1",amount=0.8), border="lightblue1",r0=-0.01, r1=1.01)
    #}
    #dev.off()

    '''
}


/* 
* Plot BAF and normalised depth with ROH regions
*/
process plot_baf_roh_gc_normalised {
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
    path "${params.sampleid}_normalisedDP_ROH_L1M.png"
    path "${params.sampleid}_normalisedDP_ROH_L1M_BAF_perChr.pdf"
    //path "${params.sampleid}_DP_ROH.png"
    //path "${params.sampleid}_DP_ROH_BAF_perChr.pdf"

    shell:
    '''
    #!/usr/bin/env Rscript
    library(karyoploteR)

    vaf <- \"!{baf}\"
    roh <- \"!{roh}\"
    mosdepth <- \"!{mosdepth}\"
    sample_name <- \"!{params.sampleid}\"
    genome <- \"!{params.genome}\"
    gc_file <- \"!{params.gc_file}\"

    # Settings for plots
    transparency <- 0.85
    trans_sample <- transparent("royalblue",amount=transparency)

    if (genome == "equCab3") {
        horse_cytoband<-read.csv(file="/Horse/equCab3_cytoband.tsv",header=TRUE,sep="\t")
        custom.genome <- toGRanges(horse_cytoband)
        chrs<-c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26","chr27","chr28","chr29","chr30","chr31","chrX")
        gc<-read.csv(file=gc_file,header=TRUE,sep="\t")
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

    # Normalise depth by GC content
    cov$chr_start <- paste(cov$chr,cov$start,sep="_")
    gc$chr_start <- paste(gc$X.1_usercol,gc$X2_usercol,sep="_")
    cov_gc <- merge(cov, gc, by.x="chr_start", by.y="chr_start")
    cov_gc_chr <- cov_gc[cov_gc$chr=="chr1" | cov_gc$chr=="chr2" | cov_gc$chr=="chr3" | cov_gc$chr=="chr4" | cov_gc$chr=="chr5" | cov_gc$chr=="chr6" | cov_gc$chr=="chr7" | cov_gc$chr=="chr8" | cov_gc$chr=="chr9" | cov_gc$chr=="chr10" | cov_gc$chr=="chr11" | cov_gc$chr=="chr12" | cov_gc$chr=="chr13" | cov_gc$chr=="chr14" | cov_gc$chr=="chr15" | cov_gc$chr=="chr16" | cov_gc$chr=="chr17" | cov_gc$chr=="chr18" | cov_gc$chr=="chr19" | cov_gc$chr=="chr20" | cov_gc$chr=="chr21" | cov_gc$chr=="chr22" | cov_gc$chr=="chr23" | cov_gc$chr=="chr24" | cov_gc$chr=="chr25" | cov_gc$chr=="chr26" | cov_gc$chr=="chr27" | cov_gc$chr=="chr28" | cov_gc$chr=="chr29" | cov_gc$chr=="chr30" | cov_gc$chr=="chr31" | cov_gc$chr=="chrX",]
    # Fit LOESS regression
    loess_model <- loess(norm_median ~ X5_pct_gc, data = cov_gc_chr, span = 0.3)
    # Predict expected coverage from GC content
    cov_gc_chr$expected <- predict(loess_model, cov_gc_chr$X5_pct_gc)
    cov_gc_chr$normalized <- cov_gc_chr$norm_median / cov_gc_chr$expected
    # Rescale to median = 1
    cov_gc_chr$normalized <- cov_gc_chr$normalized / median(cov_gc_chr$normalized, na.rm = TRUE)

    # Add log
    cov_gc_chr$normalized[cov_gc_chr$normalized>4] <- 4
    cov_gc_chr$log_normalized<-log2(cov_gc_chr$normalized)
    cov_gc_chr$log_normalized[cov_gc_chr$log_normalized<(-2)] <- -2

    data_cov<-cov_gc_chr[c("chr","start","end","log_normalized")]
    coverage<-toGRanges(data_cov)

    # Plot DP & ROH regions (>1MP)
    png(file=paste0(sample_name,"_normalisedDP_ROH_L1M.png"),type="cairo",units="in", width=10, height=10, pointsize=14, res=1000)
    pp_2 <- getDefaultPlotParams(plot.type = 2)
    pp_2$data1min <- -2
    pp_2$data1max <- 2
    pp_2$data1height <- 350
    pp_2$data2height <- 50
    kp<-plotKaryotype(plot.type=2, genome = custom.genome, chromosomes = chrs, plot.params = pp_2, cex=0.6)
    kpAxis(kp, cex=0.3)
    kpPoints(kp, coverage, x=coverage$start, y=coverage$log_normalized, col="royalblue", cex=0.2, ymin=-2, ymax=2)
    kpPlotRegions(kp, data=regions_subset,col="black", data.panel=2)
    dev.off()

    # Plot DP & ROH regions (all))
    #png(file=paste0(sample_name,"normalisedDP_ROH.png"),type="cairo",units="in", width=10, height=10, pointsize=14, res=1000)
    #pp_2 <- getDefaultPlotParams(plot.type = 2)
    #pp_2$data1min <- -2
    #pp_2$data1max <- 2
    #pp_2$data1height <- 350
    #pp_2$data2height <- 50
    #kp<-plotKaryotype(plot.type=2, genome = custom.genome, chromosomes = chrs, plot.params = pp_2, cex=0.6)
    #kpAxis(kp, cex=0.3)
    #kpPoints(kp, coverage, x=coverage$start, y=coverage$log_normalized, col="royalblue", cex=0.2, ymin=-2, ymax=2)
    #kpPlotRegions(kp, data=regions,col="black", data.panel=2)
    #dev.off()

    # Make one PDF file with BAF and depth plots per chromosome and ROH regions (>1MP)
    pdf(paste0(sample_name,"_normalisedDP_ROH_L1M_BAF_perChr.pdf"), onefile = TRUE)
    for (chr in chrs) { 
        print(chr) 
        baf_sample_chr<-baf_sample[baf_sample$chr==chr,]
        dfBAF_sample <- makeGRangesFromDataFrame(baf_sample_chr, keep.extra.columns = TRUE)
        rg_subset_bed_chr<-rg_subset_bed[rg_subset_bed$Chr==chr,]	
        regions_chr<-toGRanges(rg_subset_bed_chr)
        coverage_chr<-data_cov[data_cov$chr==chr,]
        print(head(coverage_chr))
        coverages_chr<-toGRanges(coverage_chr)

        pp <- getDefaultPlotParams(plot.type = 1)
        pp$data1height=600
        k <- plotKaryotype(plot.type = 1, genome=custom.genome, chromosomes = chr, plot.params = pp)
        kpAddBaseNumbers(k, tick.len=4, add.units=TRUE, cex=0.4)
        kpDataBackground(k, r0=0, r1=0.3, color="gray97")
        kpAxis(k, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.5)
        kpPoints(k, coverages_chr, x=coverages_chr$start, y=coverages_chr$log_normalized, col="royalblue", pch=16, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.3)
        kpAxis(k, r0=0.4, r1=1, cex=0.5)
        kpDataBackground(k, r0=0.4, r1=1, color="gray97")
        kpPoints(k, data=dfBAF_sample, cex=0.3, ymin=0, ymax=1, col=trans_sample, pch=1, r0=0.4, r1=1)
        kpPlotRegions(k, data = regions_chr, col=transparent("lightblue1",amount=0.8), border="lightblue1",r0=-0.01, r1=1.01)
    }
    dev.off()

    # Make one PDF file with BAF and depth plots per chromosome and ROH regions (all)
    #pdf(paste0(sample_name,"_normalisedDP_ROH_BAF_perChr.pdf"), onefile = TRUE)
    #for (chr in chrs) { 
    #    print(chr) 
    #    baf_sample_chr<-baf_sample[baf_sample$chr==chr,]
    #    dfBAF_sample <- makeGRangesFromDataFrame(baf_sample_chr, keep.extra.columns = TRUE)
    #    rg_bed_chr<-rg_bed[rg_bed$Chr==chr,]	
    #    regions_chr<-toGRanges(rg_bed_chr)
    #    coverage_chr<-data_cov[data_cov$chr==chr,]
    #    coverages_chr<-toGRanges(coverage_chr)
    #    pp <- getDefaultPlotParams(plot.type = 1)
    #    pp$data1height=600
    #    k <- plotKaryotype(plot.type = 1, genome=custom.genome, chromosomes = chr, plot.params = pp)
    #    kpAddBaseNumbers(k, tick.len=4, add.units=TRUE, cex=0.4)
    #    kpDataBackground(k, r0=0, r1=0.3, color="gray97")
    #    kpAxis(k, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.5)
    #    kpPoints(k, coverages_chr, x=coverages_chr$start, y=coverages_chr$log_normalized, col="royalblue", pch=16, r0=0, r1=0.3, ymin=-2, ymax=2, cex=0.3)
    #    kpAxis(k, r0=0.4, r1=1, cex=0.5)
    #    kpDataBackground(k, r0=0.4, r1=1, color="gray97")
    #    kpPoints(k, data=dfBAF_sample, cex=0.3, ymin=0, ymax=1, col=trans_sample, pch=1, r0=0.4, r1=1)
    #    kpPlotRegions(k, data = regions_chr, col=transparent("lightblue1",amount=0.8), border="lightblue1",r0=-0.01, r1=1.01)
    #}
    #dev.off()

    '''
}


/*
* Plot BAF density and violin plots
*/
process plot_density_violin_baf {
    //label 'private_node'
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'karyoploter'

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path baf

    output:
    path "${params.sampleid}_BAF_violin.png"
    path "${params.sampleid}_BAF_violin_lower_099.png"
    path "${params.sampleid}_BAF_violin_lower_095.png"
    path "${params.sampleid}_BAF_density.png"
    path "${params.sampleid}_BAF_density_fill.png"
    path "${params.sampleid}_BAF_density_y2.png"
    path "${params.sampleid}_BAF_density_y2_fill.png"
    
    shell:
    '''
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(RColorBrewer)

    baf_file <- \"!{baf}\"
    sample_name <- \"!{params.sampleid}\"

    n <- 32
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    my_color<-sample(col_vector, n)

    baf<-read.csv(file=baf_file,header=FALSE,sep="\t")
    colnames(baf) <- c("chr","begin","end","BAF") 
    baf_chr <- baf[baf$chr=="chr1" | baf$chr=="chr2" | baf$chr=="chr3" | baf$chr=="chr4" | baf$chr=="chr5" | baf$chr=="chr6" | baf$chr=="chr7" | baf$chr=="chr8" | baf$chr=="chr9" | baf$chr=="chr10" | baf$chr=="chr11" | baf$chr=="chr12" | baf$chr=="chr13" | baf$chr=="chr14" | baf$chr=="chr15" | baf$chr=="chr16" | baf$chr=="chr17" | baf$chr=="chr18" | baf$chr=="chr19" | baf$chr=="chr20" | baf$chr=="chr21" | baf$chr=="chr22" | baf$chr=="chr23" | baf$chr=="chr24" | baf$chr=="chr25" | baf$chr=="chr26" | baf$chr=="chr27" | baf$chr=="chr28" | baf$chr=="chr29" | baf$chr=="chr30" | baf$chr=="chr31" | baf$chr=="chrX",]

    density_name=paste(sample_name, "_BAF_density_fill.png", sep="")
    png(filename=density_name, type="cairo", units="in", width=16, height=8, pointsize=14, res=1000)
    ggplot(baf_chr, aes(x=BAF, color=chr, fill=chr)) + geom_density(alpha=0.5) + scale_color_manual(values=my_color) + scale_fill_manual(values=my_color) + theme_bw()
    dev.off()

    density_name=paste(sample_name, "_BAF_density.png", sep="")
    png(filename=density_name, type="cairo", units="in", width=16, height=8, pointsize=14, res=1000)
    ggplot(baf_chr, aes(x=BAF, color=chr)) + geom_density() + scale_color_manual(values=my_color) + scale_fill_manual(values=my_color) + theme_bw()
    dev.off()

    density_name=paste(sample_name, "_BAF_density_y2_fill.png", sep="")
    png(filename=density_name, type="cairo", units="in", width=16, height=8, pointsize=14, res=1000)
    ggplot(baf_chr, aes(x=BAF, color=chr, fill=chr)) + geom_density(alpha=0.5) + scale_color_manual(values=my_color) + scale_fill_manual(values=my_color) + theme_bw() + ylim(0,2)
    dev.off()

    density_name=paste(sample_name, "_BAF_density_y2.png", sep="")
    png(filename=density_name, type="cairo", units="in", width=16, height=8, pointsize=14, res=1000)
    ggplot(baf_chr, aes(x=BAF, color=chr)) + geom_density() + scale_color_manual(values=my_color) + scale_fill_manual(values=my_color) + theme_bw() + ylim(0,2)
    dev.off()

    violin_name=paste(sample_name, "_BAF_violin.png", sep="")
    png(filename=violin_name, type="cairo", units="in", width=16, height=8, pointsize=14, res=1000)
    ggplot(baf_chr, aes(x=chr, y=BAF, fill=chr)) + geom_violin(alpha=0.5) + theme_bw() + scale_fill_manual(values=my_color)
    dev.off()

    baf_chr_1 <- baf_chr[baf_chr$BAF<0.99,]
    violin_name=paste(sample_name, "_BAF_violin_lower_099.png", sep="")
    png(filename=violin_name, type="cairo", units="in", width=16, height=8, pointsize=14, res=1000)
    ggplot(baf_chr_1, aes(x=chr, y=BAF, fill=chr)) + geom_violin(alpha=0.5) + theme_bw() + scale_fill_manual(values=my_color)
    dev.off()

    baf_chr_2 <- baf_chr[baf_chr$BAF<0.95,]
    violin_name=paste(sample_name, "_BAF_violin_lower_095.png", sep="")
    png(filename=violin_name, type="cairo", units="in", width=16, height=8, pointsize=14, res=1000)
    ggplot(baf_chr_2, aes(x=chr, y=BAF, fill=chr)) + geom_violin(alpha=0.5) + theme_bw() + scale_fill_manual(values=my_color)
    dev.off()

    '''
}