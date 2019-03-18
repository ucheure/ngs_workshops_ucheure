#!/usr/bin/env bash
set -o nounset ## QUESTION: Google these settings to find out what they mean.
set -o errexit
set -o pipefail

###############################
## NGS Bash Pipeline Version ##
###############################
VERSION="0.1.0"

###############################
## software requirements     ##
###############################
# conda
# fastqc
# trimmomatic
# bwa
# samtools
# picard
# samblaster
# freebayes
# annovar
# snpeff
# multiqc
# bgzip
# tabix
# parallel
# awscli



###############################
## Some basic compute set up ##
###############################

## PATH to bioconda based ngs software and annovar
## tells the computer where the software is installed
export PATH="$HOME/anaconda3/bin:$PATH"
export PATH="${HOME}/share/software/annovar:$PATH"

## Set number of cpu (cores) to use for various steps in the pipeline
trimmomatic_cpu="12"
bwa_cpu="30"
freebayes_cpu="36"

#########################################
## Set pipeline input/output directory ##
#########################################
## Student data processing to go here
## /home/ubuntu/share/projects
## EG mine would be /home/ubuntu/share/projects/xsjnewhouse
## YOU NEED TO RUN THIS FIRST:
##                            make-project-dir.sh <project dir name>
SHARED_DIRECTORY="${HOME}/share/projects"
MY_PERSONAL_DIRECTORY="xsjnewhouse" ## name of your personal folder
FASTQ_DIR="${SHARED_DIRECTORY}/${MY_PERSONAL_DIRECTORY}/ngs_project/fastq"
ALIGNMENT_DIR="${SHARED_DIRECTORY}/${MY_PERSONAL_DIRECTORY}/ngs_project/alignments"
VCF_DIR="${SHARED_DIRECTORY}/${MY_PERSONAL_DIRECTORY}/ngs_project/variant_calls"
ANNOTATIONS_DIR="${SHARED_DIRECTORY}/${MY_PERSONAL_DIRECTORY}/ngs_project/variant_annotations"
REPORTS_DIR="${SHARED_DIRECTORY}/${MY_PERSONAL_DIRECTORY}/ngs_project/reports"

## Set temp directory used with various tools
temp_directory="${HOME}/share/projects/${MY_PERSONAL_DIRECTORY}/ngs_project/temp"


#######################################
## Symbolinc links to some test data ##
#######################################
## for this workshop
# I have provided some fastq data already
# located on /home/ubuntu/share/ngs_data/fastq/
# to save space we will create a symbolic link to this data 
# e.g ln -s source_file myfile
# google symbolic link
# you only need to do this once

ln -s ${HOME}/share/ngs_data/fastq/workshop_data/WES01_chr22_R1.fastq.gz ${FASTQ_DIR}/WES01_chr22_R1.fastq.gz;
ln -s ${HOME}/share/ngs_data/fastq/workshop_data/WES01_chr22_R2.fastq.gz ${FASTQ_DIR}/WES01_chr22_R2.fastq.gz;

############################
## Some basic data set up ##
############################

## FASTQ FILES
READ1="${FASTQ_DIR}/${1:-WES01_chr22_R1.fastq.gz}" # you can set this here or take it from the command line
READ2="${FASTQ_DIR}/${2:-WES01_chr22_R2.fastq.gz}" # a default value has been set to point to ${FASTQ_DIR}/WES01_chr22_R{1,2}.fastq.gz

FASTQ_BASE_NAME=$(basename ${READ1} .fastq.gz) ## change the extenstion to match your data e.g it might be .fq or .fq.gz

## BASIC READ GROUP INFORMATION 
## https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
RG_SM="SJN" ## Use any unique short name or your initials 
RG_PLATFORM="ILLUMINA"
RG_PLATFORM_UNIT="" ## get this from the fastq files (FLowcell id)
RG_LIBRARY="" ## make something up
RG_DATE=$(date +%Y-%m-%d) ## auto generated
FLOWCELL_LANE="" ## get this from the fastq files 
RG_ID="${RG_SM}.${RG_PLATFORM_UNIT}.${FLOWCELL_LANE}"

## THE REFERENCE GENOME
REF_GENOME_DIR="${HOME}/share/ngs_data/references/Homo_sapiens/GATK/hg19"
REF_GENOME_VERSION="hg19"
REF_GENOME_FASTA="${REF_GENOME_DIR}/ucsc.hg19.fasta"
REF_GENOME_PREFIX=$(basename ${REF_GENOME_FASTA} .fasta)

## ANNOTATION DATABASES
PATH_TO_ANNOVAR_DB="${HOME}/share/software/annovar"

####################
## START PIPELINE ##
####################

## Move to you personal directory
cd ${MY_PERSONAL_DIRECTORY}

## 1.0 fatsqc ---------------------------------------------------=##

fatsqc --threads 2 ${READ1} ${READ2}

## 2.0 trimmomatic -----------------------------------------------##


# set adapter sequence fasta file
ADAPTER_FASTA="${HOME}/anaconda3/share/trimmomatic/adapters"

# set name for trimmomatic paired output
QCD_PE_FASTQ_R1=${FASTQ_BASE_NAME}.paired-trimmed.fq.gz
QCD_PE_FASTQ_R2=${FASTQ_BASE_NAME}.paired-trimmed.fq.gz

# set name for trimmomatic unpaired ouput
QCD_SE_FASTQ_R1=${FASTQ_BASE_NAME}.unpaired-trimmed.fq.gz
QCD_SE_FASTQ_R2=${FASTQ_BASE_NAME}.unpaired-trimmed.fq.gz

# RUN Trimmomatic followed by Fastqc on trimmed files
trimmomatic \
PE \
-phred64 \
-threads ${trimmomatic_cpu} \
-trimlog ${FASTQ_DIR}/${FASTQ_BASE_NAME}_trimmomatic_qc.log \
${FASTQ_DIR}/${READ1} ${FASTQ_DIR}/${READ2} \
${FASTQ_DIR}/${QCD_PE_FASTQ_R1} ${FASTQ_DIR}/${QCD_PE_FASTQ_R2} \
${FASTQ_DIR}/${QCD_SE_FASTQ_R1} ${FASTQ_DIR}/${QCD_SE_FASTQ_R2} \
ILLUMINACLIP:${ADAPTER_FASTA}:2:30:10:5:true \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:50 && \
fastqc \
--threads 4 \
${FASTQ_DIR}/${QCD_PE_FASTQ_R1} \
${FASTQ_DIR}/${QCD_PE_FASTQ_R2} \
${FASTQ_DIR}/${QCD_SE_FASTQ_R1} \
${FASTQ_DIR}/${QCD_SE_FASTQ_R2};

## QUESTION: what does the && mean?

## 3.0 Alingment BWA mem -----------------------------------------------##
# to speed things up we will use lots of pipes
# and take advantage samblaster and sambamba 
# to quickly allow piped processes between software
# and faster indexing of sam and bam files
# align with bwa mem
# mark duplicates with samblaster
# convert to bam and sort with sambamba
# index bam with sambamba

# Name of the BAM FILE
BAM_PREFIX="${RG_SM}.${REF_GENOME_VERSION}.${RG_PLATFORM}.${RG_PLATFORM_UNIT}"

bwa mem \
-M \
-t ${bwa_cpu} \
-R '@RG\tID:${RG_ID}\tSM:${RG_SM}\tPU:${RG_PLATFORM_UNIT}\tPL:${RG_PLATFORM}\tLB:${RG_LIBRARY}\tDT:${RG_DATE}' \
${REF_GENOME_FASTA} \
${FASTQ_DIR}/${QCD_PE_FASTQ_R1} ${FASTQ_DIR}/${QCD_PE_FASTQ_R2} \
| samblaster --addMateTags --excludeDups \
| sambamba view -t ${bwa_cpu} -S -f bam /dev/stdin \
| sambamba sort -t ${bwa_cpu} -m 4GB --tmpdir=${temp_directory} -o ${ALIGNMENT_DIR}/${BAM_PREFIX}.dupemk.bam /dev/stdin && \
sambamba index --threads=${bwa_cpu} ${ALIGNMENT_DIR}/${BAM_PREFIX}.dupemk.bam


## 4.0 Post Alingment Processing -----------------------------------------------##
## You can edit this to include and insert GATK best practices 
## google them and read the paper I provided
## Note: GATK best practices will take a loooong time to run 
## So we choose to run a very basic qc and filtering process on the alignemd bam file

## Basic Filter BAM 
samtools view -bh -q 20 -F 1796 ${ALIGNMENT_DIR}/${BAM_PREFIX}.dupemk.bam > ${ALIGNMENT_DIR}/${BAM_PREFIX}.filtered.bam && \
sambamba index --threads=${bwa_cpu} ${ALIGNMENT_DIR}/${BAM_PREFIX}.filtered.bam;



## 5.0 Variant Calling: Freebayes -----------------------------------------------##

## Custom utils from bcbio 
## Taken from https://github.com/chapmanb/bcbio-nextgen/blob/98c75603907cb22a3e4cd4fd78f7e995b80bddfd/bcbio/variation/vcfutils.py#L76
function fix_ambiguous() {

  awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, "N", $4) } {print}'

}

## Generate regions that are equal in terms of data content, and thus have lower variance
## in runtime.  This will yield better resource utilization.
bamtools coverage -in ${ALIGNMENT_DIR}/${BAM_PREFIX}.filtered.bam \
| coverage_to_regions.py ${REF_GENOME_FASTA} 500 > ${REF_GENOME_FASTA}.500.regions;

## Freebayes Parallel
## usage: freebayes-parallel [regions file] [ncpus] [freebayes arguments]
## lots of piped post processing for downstream analysis in Annovar etc

## Move to dir that contains freebayes-parallel script 
## recommended by software authors
cd /home/ubuntu/anaconda3/bin/

## run freebayes-parallel
freebayes-parallel \
${REF_GENOME_FASTA}.500.regions \
${freebayes_cpu} \
-f ${REF_GENOME_FASTA} \
--standard-filters \
--min-coverage 4 \
${ALIGNMENT_DIR}/${BAM_PREFIX}.filtered.bam \
| vcffilter -f "QUAL > 10 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" -s \
| fix_ambiguous \
| vcfallelicprimitives --keep-geno \
| vcffixup - \
| vt normalize -r ${REF_GENOME_FASTA} -q - 2> /dev/null \
| bgzip -c > ${VCF_DIR}/${BAM_PREFIX}.norm.vcf.gz;

tabix ${VCF_DIR}/${BAM_PREFIX}.norm.vcf.gz;

## move back to your personal directory
cd ${MY_PERSONAL_DIRECTORY}

## 6.0 Variant Annotation -----------------------------------------------##

## 6.1 ANNOVAR

VCF_TO_ANNOTATE="${VCF_DIR}/${BAM_PREFIX}.norm.vcf.gz"

## run annovar
table_annovar.pl \
${VCF_TO_ANNOTATE} \
${PATH_TO_ANNOVAR_DB}/humandb/ \
-buildver hg19 \
-out ${VCF_TO_ANNOTATE}.annovar \
-remove \
-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,clinvar_20180603,ljb26_all \
-operation gx,r,f,f,f,f,f \
-nastring . \
-csvout \
-polish \
-vcfinput

## 6.2 snpeff
snpeff 


## 7.0 QC report generation -----------------------------------------------##
cd ${REPORTS_DIR}

ln -s ${ALIGNMENT_DIR}/${BAM_PREFIX}.filtered.bam ${REPORTS_DIR}/${BAM_PREFIX}.filtered.bam
ln -s ${ALIGNMENT_DIR}/${BAM_PREFIX}.filtered.bam ${REPORTS_DIR}/${BAM_PREFIX}.filtered.bam.bai
ln -s ${ALIGNMENT_DIR}/${BAM_PREFIX}.dupemk.bam ${REPORTS_DIR}/${BAM_PREFIX}.dupemk.bam
ln -s ${ALIGNMENT_DIR}/${BAM_PREFIX}.dupemk.bam ${REPORTS_DIR}/${BAM_PREFIX}.dupemk.bam.bai

for my_bam in ${REPORTS_DIR}/${BAM_PREFIX}.filtered.bam ${REPORTS_DIR}/${BAM_PREFIX}.dupemk.bam;do
	samtools flagstat ${my_bam} && \
	bamtools stats ${my_bam} && \
	picard CollectInsertSizeMetrics && \
	picard CollectMultipleMetrics && \


## MultiQC
multiqc ${SHARED_DIRECTORY}/${MY_PERSONAL_DIRECTORY}/ngs_project/


#############################################################




# NOT DONE
#vep_install -a cf -s homo_sapiens -y GRCh38 -c /output/path/to/GRCh38/vep --CONVERT
#vep_install -a cf -s homo_sapiens -y hg19 -c /output/path/to/GRCh38/vep --CONVERT
#vep_install -a cf -s homo_sapiens -y GRCh37 -c /home/ubuntu/share/software/vep_db/GRCh37 --CONVERT


