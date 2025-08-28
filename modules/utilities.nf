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
    path "${vcf.getSimpleName()}*_PASS_SNVs_AF.bed", emit: vcf_baf

    shell:
    '''
    mean_DP=`zcat !{vcf} | grep -v "^#" | awk 'BEGIN {FS=OFS="\t"} {m=split(\$9,a,":"); n=split(\$10,b,":"); for (i=1; i<=m; i++) {if (a[i]=="DP") {sum=sum+b[i]}}} END {print sum/NR}'`
    echo "Mean DP: $mean_DP"
    if [[ $(echo "$mean_DP" | bc) < 3 ]]; then
        DP_cutoff=2
    elif [[ $(echo "$mean_DP" | bc) < 6 ]]; then
        DP_cutoff=3
    elif [[ $(echo "$mean_DP" | bc) < 10 ]]; then
        DP_cutoff=5
    else
        DP_cutoff=10
    fi
    echo "DP cutoff: $DP_cutoff"

    mean_GQ=`zcat !{vcf} | grep -v "^#" | awk 'BEGIN {FS=OFS="\t"} {m=split(\$9,a,":"); n=split(\$10,b,":"); for (i=1; i<=m; i++) {if (a[i]=="GQ") {sum=sum+b[i]}}} END {print sum/NR}'`
    echo "Mean GQ: $mean_GQ"
    if [[ $(echo "$mean_GQ" | bc) < 5 ]]; then
        GQ_cutoff=3
    elif [[ $(echo "$mean_GQ" | bc) < 9 ]]; then
        GQ_cutoff=5
    elif [[ $(echo "$mean_GQ" | bc) < 15 ]]; then
        GQ_cutoff=10
    else
        GQ_cutoff=15
    fi
    zcat !{vcf} | grep -v "^#" | awk -v GQ_cutoff=$GQ_cutoff -v DP_cutoff=$DP_cutoff 'BEGIN {FS=OFS="\t"} {m=split(\$9,a,":"); n=split(\$10,b,":"); for (i=1; i<=m; i++) {if (a[i]=="VAF") {af=b[i]}; if (a[i]=="GQ") {gq=b[i]}; if (a[i]=="DP") {dp=b[i]}}; if (gq>=GQ_cutoff && dp>=DP_cutoff && length(\$4)==1 && length(\$5)==1 && \$6>=10 &&  \$7=="PASS") {print \$1,\$2-1, \$2, af}}' > !{vcf.getSimpleName()}_DP${DP_cutoff}_GQ${GQ_cutoff}_QUAL10_PASS_SNVs_AF.bed
    '''
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
