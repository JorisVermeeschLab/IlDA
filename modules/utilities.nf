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
    val read

    output:
    path "${fastq.getSimpleName()}.clipped${clipped_length}.${read}.fastq", emit: clipped_fastq

    script:
    """
    zcat $fastq | awk -v clip=$clipped_length '{new=\$0;flag=NR%4;if(flag==0||flag==2){new=substr(\$0,1,clip)}; print(new)}' > ${fastq.getSimpleName()}.clipped${clipped_length}.${read}.fastq
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

    publishDir path: "${params.outdir}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path vcf
    
    output:
    path "${vcf.getSimpleName()}*_PASS_SNVs_AF.bed", emit: bed_baf
    path "${vcf.getSimpleName()}*_PASS_SNVs_AF.vcf", emit: vcf_baf

    shell:
    '''
    mean_DP=`zcat !{vcf} | grep -v "^#" | awk 'BEGIN {FS=OFS="\t"} {m=split(\$9,a,":"); n=split(\$10,b,":"); for (i=1; i<=m; i++) {if (a[i]=="DP") {sum=sum+b[i]}}} END {print sum/NR}'`
    echo "Mean DP: $mean_DP"
    round_mean_DP=${mean_DP%.*}
    echo $round_mean_DP
    #if [[ $(echo "$mean_DP" | bc) < 3 ]]; then
    if [[ $round_mean_DP -lt 3 ]]; then
        DP_cutoff=2
    elif [[ $round_mean_DP -lt 6 ]]; then
        DP_cutoff=3
    elif [[ $round_mean_DP -lt 10 ]]; then
        DP_cutoff=5
    else
        DP_cutoff=10
    fi
    echo "DP cutoff: $DP_cutoff"

    mean_GQ=`zcat !{vcf} | grep -v "^#" | awk 'BEGIN {FS=OFS="\t"} {m=split(\$9,a,":"); n=split(\$10,b,":"); for (i=1; i<=m; i++) {if (a[i]=="GQ") {sum=sum+b[i]}}} END {print sum/NR}'`
    echo "Mean GQ: $mean_GQ"
    round_mean_GQ=${mean_GQ%.*}
    echo $round_mean_GQ
    if [[ $round_mean_GQ -lt 5 ]]; then
        GQ_cutoff=3
    elif [[ $round_mean_GQ -lt 9 ]]; then
        GQ_cutoff=5
    elif [[ $round_mean_GQ -lt 15 ]]; then
        GQ_cutoff=10
    else
        GQ_cutoff=15
    fi
    echo "GQ cutoff: $GQ_cutoff"

    mean_QUAL=`zcat !{vcf} | grep -v "^#" | awk 'BEGIN {FS=OFS="\t"} {if ($7!="RefCall") {sum=sum+$6}} END {print sum/NR}'`
    echo "Mean QUAL: $mean_QUAL"
    round_mean_QUAL=${mean_QUAL%.*}
    echo $round_mean_QUAL
    if [[ $round_mean_QUAL -lt 5 ]]; then
        QUAL_cutoff=3
    elif [[ $round_mean_QUAL -lt 9 ]]; then
	QUAL_cutoff=5
    elif [[ $round_mean_QUAL -lt 15 ]]; then
        QUAL_cutoff=10
    else
        QUAL_cutoff=15
    fi
    echo "QUAL cutoff: $QUAL_cutoff"

    zcat !{vcf} | grep -v "^#" | awk -v GQ_cutoff=$GQ_cutoff -v DP_cutoff=$DP_cutoff -v QUAL_cutoff=$QUAL_cutoff 'BEGIN {FS=OFS="\t"} {m=split(\$9,a,":"); n=split(\$10,b,":"); for (i=1; i<=m; i++) {if (a[i]=="VAF") {af=b[i]}; if (a[i]=="GQ") {gq=b[i]}; if (a[i]=="DP") {dp=b[i]}}; if (gq>=GQ_cutoff && dp>=DP_cutoff && length(\$4)==1 && length(\$5)==1 && \$6>=QUAL_cutoff &&  \$7=="PASS") {print \$1,\$2-1, \$2, af}}' > !{vcf.getSimpleName()}_DP${DP_cutoff}_GQ${GQ_cutoff}_QUAL${QUAL_cutoff}_PASS_SNVs_AF.bed
    zcat !{vcf} | awk -v GQ_cutoff=$GQ_cutoff -v DP_cutoff=$DP_cutoff -v QUAL_cutoff=$QUAL_cutoff 'BEGIN {FS=OFS="\t"} {if (\$0~/^#/) {print \$0} else {m=split(\$9,a,":"); n=split(\$10,b,":"); for (i=1; i<=m; i++) {if (a[i]=="VAF") {af=b[i]}; if (a[i]=="GQ") {gq=b[i]}; if (a[i]=="DP") {dp=b[i]}}; if (gq>=GQ_cutoff && dp>=DP_cutoff && length(\$4)==1 && length(\$5)==1 && \$6>=QUAL_cutoff &&  \$7=="PASS") {print \$0}}}' > !{vcf.getSimpleName()}_DP${DP_cutoff}_GQ${GQ_cutoff}_QUAL${QUAL_cutoff}_PASS_SNVs_AF.vcf
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
