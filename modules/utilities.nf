/* 
* Clip reads with awk
*/
process clip_reads {
    //label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'

    //publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path fastq
    val clipped_length

    output:
    path "${fastq.getSimpleName()}.clipped${clipped_length}.fastq", emit: clipped_fastq

    script:
    """
    zcat $fastq | awk -v clip=$clipped_length '{new=\$0;flag=NR%4;if(flag==0||flag==2){new=substr(\$0,1,clip)}; print(new)}' > ${fastq.getSimpleName()}.clipped${clipped_length}.fastq
    """
}


/* 
* Merge FASTQ files
*/
process merge_fastq_files {
    //label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'

    //publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path fastqs
    val read
    
    output:
    path "${params.sampleid}.merged.${read}.fastq", emit: merged_fastq

    script:
    """
    zcat $fastqs > ${params.sampleid}.merged.${read}.fastq
    """
}


/* 
* Get BAF from VCF file
*/
process get_baf {
    //label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'

    input:
    path vcf
    
    output:
    path "${vcf.getSimpleName()}_GQ15_DP10_QUAL10_PASS_SNVs_AF.bed", emit: vcf_baf

    script:
    """
    zcat $vcf | grep -v "^#" | awk 'BEGIN {FS=OFS="\t"} {m=split(\$9,a,":"); n=split(\$10,b,":"); for (i=1; i<=m; i++) {if (a[i]=="VAF") {af=b[i]}; if (a[i]=="GQ") {gq=b[i]}; if (a[i]=="DP") {dp=b[i]}}; if (gq>=15 && dp >=10 && length(\$4)==1 && length(\$5)==1 && \$6>=10 &&  \$7=="PASS") {print \$1 ":" \$2, af}}' > ${vcf.getSimpleName()}_GQ15_DP10_QUAL10_PASS_SNVs_AF.tsv
    awk 'BEGIN {FS=OFS="\t"} {split(\$1,a,":"); print a[1], a[2]-1, a[2], \$2}' ${vcf.getSimpleName()}_GQ15_DP10_QUAL10_PASS_SNVs_AF.tsv > ${vcf.getSimpleName()}_GQ15_DP10_QUAL10_PASS_SNVs_AF.bed
    """
}


/* 
* Get RG tags from ROH file
*/
process get_roh_rg {
    //label 'private_node'
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'

    input:
    path roh
    
    output:
    path "${roh.getSimpleName()}_RG.tsv", emit: roh_rg

    script:
    """
    echo -e "Type\tSample\tChr\tStart\tEnd\tLength\tMarkers\tQuality" > ${roh.getSimpleName()}_RG.tsv
    grep -v "^#" $roh | grep RG  >> ${roh.getSimpleName()}_RG.tsv
    """
}
