#! bin/bash
#ssh into the server and enter password blindly

cd rosalind
cd Amarachukwu
mkdir raw_data
cd raw_data  

#download sample datasets from zendo
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

#download reference genome & unzip reference
mkdir reference
cd reference
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
gunzip hg19.chr5_12_17.fa.gz

cd ..

Removing low quality sequences using Trimmomatic
#dwownlaod trimmomatic & unzip 
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
cp Trimmomatic-039/adapters/TruSeq3-PE.fa $PWD

#create a txt file with nano and input the sequence name
nano trimmed.sh

SLGFSK-N_231335
SLGFSK-T_231336
#save, rename list.txt & exit nano


#Trimming reads

nano trimmed.sh 
mkdir trimmed_reads

for sample in `cat list.txt`
do
       trimmomatic PE -threads 8 ${sample}_r1_chr5_12_17.fastq.gz ${sample}_r2_chr5_12_17.fastq.gz \
               trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
               trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25
       
       fastqc  trimmed_reads/${sample}_r1_paired.fq.gz  trimmed_reads/${sample}_r2_paired.fq.gz \
                 -o trimmed_reads/Fastqc_results
done 

multiqc  trimmed_reads/Fastqc_results  -o trimmed_reads/Fastqc_results
#save and exit nano
bash trimmed.sh

#Mapping or Aligning the reads

nano mapping.sh

#Reads Mapping
mkdir Mapping

#Index reference
bwa index reference/hg19.chr5_12_17.fa

#Perform alignment with bwa mem
bwa mem -R '@RG\tID:231335\tSM:Normal' ref/hg19.chr5_12_17.fa \
	trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' ref/hg19.chr5_12_17.fa \
	trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam	


#Conversion SAM file to BAM file sorting and indexing the sorted BAM file 

nano sort.sh

for sample in `cat list.txt`
do
        Convert SAM to BAM and sort it 
        samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -n -@ 32 > Mapping/${sample}.sorted.bam
        
        Index BAM file
        samtools index Mapping/${sample}.sorted.bam
done
#save and exit nano
bash sort.sh

#Mapped reads filtering

#nano filtering.sh

for sample in `cat list.txt`
do
	#Filter BAM files
        samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > Mapping/${sample}.filtered1.bam
done
#save and exit nano
bash filtering.sh


#Duplicate removal

#nano markdup.sh
for sample in `cat list.txt`
do
	samtools collate ${sample}.filtered1.bam ${sample}.namecollate
        samtools fixmate -m ${sample}.namecollate.bam ${sample}.fixmate.bam
        samtools sort -@ 32 -o ${sample}.positionsort.bam ${sample}.fixmate.bam
        samtools markdup -@32 -r ${sample}.positionsort.bam ${sample}.clean.bam
done

#save and exit nano
bash markdup.sh

#Left Align BAM

nano leftalign.sh
for sample in `cat list.txt`
do      
        cat ${sample}.clean.bam  | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > ${sample}.leftAlign.bam    
#save and exit nano
bash leftalign.sh

#Recalibrate read mapping qualities 

nano recalibrate.sh
for sample in `cat list.txt`
do
        samtools calmd -@ 32 -b ${sample}.leftAlign.bam hg19.chr5_12_17.fa > ${sample}.recalibrate.bam
done 
#save and exit nano
bash recalibrate.sh

#Refilter read mapping qualities

nano refilter.sh
for sample in `cat list.txt`
do
        bamtools filter -in ${sample}.recalibrate.bam -mapQuality "<=254" > ${sample}.refilter.bam
done
#save and exit nano
bash refilter.sh

cd ..

#Variant calling and classification

wget  https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar	   

#convert variant to pileup
nano variant.sh

mkdir Variants

for sample in `cat list.txt`
do
        samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                > Variants/${sample}.pileup
done

#Variant calling
	varscan somatic Variants/SLGFSK-N_231335.pileup \
        Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
        --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 
       

#merge vcf using bcftools
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf

#save and exit nano
bash variant.sh

#Download snpEff database
snpEff download hg19

#variant annotation
 snpEff hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vc
 
