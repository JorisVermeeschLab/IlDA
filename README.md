# IlDA <br />Illumina Data Analysis pipeline

## Instructions for users

This pipeline uses containers. Upon first use the container given in the config file will be downoaded, eventually converted to singularity and stored in the following directory:
```
/staging/leuven/stg_00019/research/nextflowPipeline/singularity_images/
```

Resources such as reference genome fasta sequence or any other file used in the nextflow pipeline should be stored/found in the following directory:
```
/staging/leuven/stg_00019/research/nextflowPipeline/Resources/
```


To start the pipeline
1. Adapt the  *nextflow.config* file according to the analysis to perform (see below for some details). If the config file is renamed, the full path of the config file should be given when starting the pipeline with *-config* option. <br />Do not forget to modify the parameter *workdir* so that all temporary results are written in your scratch directory, see example below.
```
workDir = '/scratch/leuven/315/vsc31555'
```
2. Adapt the following parameters to suit your environment:
```
module load Java/17.0.6
export APPTAINER_TMPDIR='/scratch/leuven/315/vsc31555/runPipeline/tmp'
export APPTAINERENV_TMPDIR='/scratch/leuven/315/vsc31555/runPipeline/tmp'
export APPTAINER_CACHEDIR='/scratch/leuven/315/vsc31555/runPipeline/tmp'
export SINGULARITY_LOCALCACHEDIR='/scratch/leuven/315/vsc31555/runPipeline'
export TMPDIR='/scratch/leuven/315/vsc31555/runPipeline/tmp'
export SLURM_CLUSTERS='wice'
```
3. Run the following command, with eventually *-bg* option so that it runs in the background:
```
~/nextflow run main.nf -profile singularity,srun_new -config your_nextflow.config
```

In case the pipeline crashes, the config file and/or the code can be modified and the pipeline can be resumed by launching the following command:
```
~/nextflow run main.nf -profile singularity,srun_new -config your_nextflow.config -resume
```

#### Analysis steps (WGS, single sample)
 - Clipping FASTQ files
 - Merging FASTQ files
 - Mapping paired end reads to reference genome with BWA
 - Mark duplicates
 - Evaluate mapping (Picard metrics)
 - Call CNVs with QDNAseq
 - Call SNVs and indels with DeepVariant


## Instructions for developers

Before pushing any change, please make sure to rebase with the following command so that the pushed changes are compatible with the last available version:
```
git pull --rebase
```