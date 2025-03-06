#### Variants calling from transcriptome data

```bash
##get the index file for reference genome
#/lustre/nobackup/WUR/ABGC/ni010/Software/STAR/STAR-2.7.11a/bin/Linux_x86_64/STAR --runThreadN 32 --runMode genomeGenerate --genomeSAindexNbases 13 --genomeDir /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/index --sjdbGTFfile /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.109.chr.gtf --genomeFastaFiles /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.dna.toplevel.fa --sjdbOverhang 100
##mapping to reference genome
for i in $(cat id); do
/lustre/nobackup/WUR/ABGC/ni010/Software/STAR/STAR-2.7.11a/bin/Linux_x86_64/STAR --runThreadN 32 --genomeDir /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/index --readFilesCommand gunzip -c --readFilesIn /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/raw_data/${i}_1.fq.gz /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/raw_data/${i}_2.fq.gz --sjdbGTFfile /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.109.chr.gtf --outSAMtype BAM Unsorted --outFileNamePrefix /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/1_mapping/${i}/${i} --limitSjdbInsertNsj 2000000 --outSAMattrRGline ID:3YW3 SM:ov PL:ILLUMINA
done


##data cleanup
for name in $(cat id); do
java -Xmx32g -jar /home/jwyuan/software/picard-tools-1.119/SortSam.jar \
INPUT=/home/jwyuan/nax/chapter2/eQTL/out/${name}Aligned.out.bam \
OUTPUT=/home/jwyuan/nax/chapter2/eQTL/sort_bam/${name}.sort.bam \
SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/home/jwyuan/nax/tmp \
&& java -Xmx32g -jar /home/jwyuan/software/picard-tools-1.119/MarkDuplicates.jar \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
INPUT=/home/jwyuan/nax/chapter2/eQTL/sort_bam/${name}.sort.bam \
OUTPUT=/home/jwyuan/nax/chapter2/eQTL/mark_dup/${name}.sort.dedup.bam \
METRICS_FILE=/home/jwyuan/nax/chapter2/eQTL/${name}.sort.dedup.metrics \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/home/jwyuan/nax/tmp \
&& java -Xmx32g -jar /home/jwyuan/software/picard-tools-1.119/ReorderSam.jar \
I=/home/jwyuan/nax/chapter2/eQTL/mark_dup/${name}.sort.dedup.bam \
O=/home/jwyuan/nax/chapter2/eQTL/mark_dup_reorder/${name}.sort.dedup.reorder.bam \
REFERENCE=/home/jwyuan/nax/ref/GRCg7w_dna_toplevel.fa \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/home/jwyuan/nax/tmp \
&& rm -f /home/jwyuan/nax/chapter2/eQTL/mark_dup/${name}.sort.dedup.bam \
&& rm -f /home/jwyuan/nax/chapter2/eQTL/${name}.sort.dedup.metrics \
&&
/home/jwyuan/software/samtools-1.12/samtools index \
/home/jwyuan/nax/chapter2/eQTL/mark_dup_reorder/${name}.sort.dedup.reorder.bam \
/home/jwyuan/software/gatk-4.2.2.0/gatk SplitNCigarReads \
-R /home/jwyuan/nax/ref/GRCg7w_dna_toplevel.fa \
-I /home/jwyuan/nax/chapter2/eQTL/mark_dup_reorder/${name}.sort.dedup.reorder.bam \
-O /home/jwyuan/nax/chapter2/eQTL/split/${name}_split.bam

##Variants calling
/home/jwyuan/software/gatk-4.2.2.0/gatk --java-options "-Xmx32G" HaplotypeCaller \
-R /home/jwyuan/nax/ref/GRCg7w_dna_toplevel.fa \
-I /home/jwyuan/nax/chapter2/eQTL/recal_bam/${name}.recal.bam \
-ERC GVCF \
-O /home/jwyuan/nax/chapter2/eQTL/raw_vcf/${name}.g.vcf \
--standard-min-confidence-threshold-for-calling 20 \
--create-output-variant-index FALSE \

&& gzip /home/jwyuan/nax/chapter2/eQTL/raw_vcf/${name}.g.vcf
done

##combine gVCFs
gatk CombineGVCFs -R /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.dna.toplevel.fa --variant *.g.vcf -O /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/6_raw_vcf/rna.g.vcf

#variants calling
gatk GenotypeGVCFs -R /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.dna.toplevel.fa -V rna.g.vcf -O rna_raw.vcf

#extract SNPs
gatk SelectVariants -R /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.dna.toplevel.fa -V rna_raw.vcf --select-type SNP -O rna_snp.vcf

#filter SNPs
gatk VariantFiltration -R /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.dna.toplevel.fa -V rna_snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'SNP_filter' -O rna_filter_snp.vcf

#extract filtered SNPs
gatk SelectVariants  -R /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.dna.toplevel.fa -V rna_filter_snp.vcf --exclude-filtered  -O rna_filtered_snp.vcf

##transfer to PLINK format
module load plink/2.0 beagle
plink --vcf rna_filtered_snp.vcf --allow-extra-chr --chr-set 41 --recode --transpose --keep-allele-order --out rna_filtered_snp
mv rna_filtered_snp.tped rna_filtered_snp_raw.tped
sed -n '/*/!p' rna_filtered_snp_raw.tped > rna_filtered_snp.tped
plink --tfile rna_filtered_snp --allow-extra-chr --chr-set 41 --make-bed --out rna_filtered_snp
plink2 --bfile rna_filtered_snp --allow-extra-chr --chr-set 41 --keep-allele-order --set-all-var-ids @:# --make-bed --out rna_snpid_raw

##filter SNPs
awk '$1 ~ /^J|^M/ {print $2}' rna_snpid.bim > snpid_nochr
awk '$5=="." {print $2}' rna_snpid.bim > snpid_0 
cat snpid_0 snpid_nochr > snpid_remove
sort snpid_remove | unique > snpid_remove_1
plink --bfile rna_snpid_raw --allow-extra-chr --chr-set 41 --keep-allele-order --exclude snpid_remove_1 --make-bed --out rna_snpid

##impute with beagle
plink2 --bfile /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/6_raw_vcf/rna_snpid_ref --keep-allele-order --allow-extra-chr --chr-set 41 --recode vcf-iid --chr $i --out chr$i
plink2 --bfile /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref --keep-allele-order --allow-extra-chr --chr-set 41 --recode vcf-iid --chr $i --out ref_chr$i
sed 's/[/]/|/g' ref_chr$i.vcf > ref_chr${i}_phased.vcf
java -Xmx50g -Xss5m -jar beagle.jar gt=chr$i.vcf ref=ref_chr${i}_phased.vcf out=chr${i}_impu nthreads=8

ls */*.vcf.gz  | xargs vcf-concat > rna_imput.vcf

##merge the PLINK files from the transcriptome and WGS
module load plink
plink --vcf rna_imput.vcf --allow-extra-chr --chr-set 41 --make-bed --keep-allele-order --out rna_imput
plink --bfile rna_imput --allow-extra-chr --chr-set 41 --keep id_ovary_plink --keep-allele-order --make-bed --out rna_ovary_imput
plink --bfile /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref --allow-extra-chr --chr-set 41 --keep id_ovary_noimput --make-bed --keep-allele-order --out rna_ovary_noimput
plink --bfile rna_ovary_imput --bmerge rna_ovary_noimput --chr-set 41 --keep-allele-order --allow-extra-chr --nonfounders --make-bed --out rna_ovary
plink --bfile rna_ovary --chr-set 41 --allow-extra-chr --nonfounders --keep-allele-order --maf 0.05 --make-bed --out rna_ovary_maf
plink --bfile rna_ovary_maf --chr-set 41 --allow-extra-chr --keep-allele-order --make-bed --update-ids change_ids --out rna_ovary_maf_changeid
plink --bfile rna_ovary_maf_changeid --chr-set 41 --allow-extra-chr --keep-allele-order --indiv-sort n --make-bed --out rna_ovary_maf_sort
plink --bfile rna_ovary_maf_sort --chr-set 41 --allow-extra-chr --keep-allele-order --recode A-transpose --out rna_ovary_maf_sort
##headline
head -n 1 rna_ovary_maf_sort.traw | sed 's/\t/ /g' | sed 's/_/ /g' > headline
tr ' ' '\n' < headline > headline_1
sed -n '1~2p' headline_1 | sed '1,3d' | sed '1i snpid' > headline_2
tr '\n' ' ' < headline_2 > headline_3
sed -e '$a\' headline_3 > headline_4
###eQTL-mapping input data, genotype, SNP.txt
sed 's/\t/ /g' rna_ovary_maf_sort.traw | cut --complement -d' ' -f1,3,4,5,6 | sed '1d' > SNP_raw.txt
cat headline_4 SNP_raw.txt > SNP.txt
```

#### TPM for each gene when mapped to the chicken reference genome (GRCg7W)

```bash
trait=(5YW1 5YW2 5YW3 5YW4 5YW5 5YW6)
i=0
while (($i<=5))
do
mkdir ${trait[$i]}
hisat2 -p 16 --dta --rna-strandness RF \
-x /lustre/nobackup/WUR/ABGC/ni010/seq/ref/index/GRCg7w \
-1 /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/raw_data/${trait[$i]}_1.fq.gz \
-2 /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/raw_data/${trait[$i]}_2.fq.gz \
-S /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.sam \
--summary-file /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.summary \
--met-file /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.metrics \
#&& rm /home/jwyuan/jyuan/output/1fe/${trait[$i]}.sam

#####sequence summary
awk '{print $2,$3}' /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.metrics > /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/QC/${trait[$i]}_seq.txt
#awk -v OFS='\t' 'NR==1{print "file", $0};FNR==2{print FILENAME, $0}' *.txt | sed -i 's/.txt/ /g' out.txt

######mapping summary
tail -2 /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.summary > /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/QC/${trait[$i]}_mapping.txt
#awk -v OFS='\t' 'NR==1{print "file", $0};FNR==2{print FILENAME, $0}' *.txt > out &&  sed -i 's/.txt/ /g' out && awk '{print $1,$2}' out > mapping.txt && rm -rf out
module load samtools
samtools sort -@ 16 \
-o /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.bam \
/lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.sam

#========gene counts================
featureCounts \
-a /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.109.chr.gtf \
-T 6 --countReadPairs -p -g gene_id -t exon \
-o /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.counts.txt \
/lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.bam

#==========gene fpkm===========#
mkdir /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/gene/${trait[$i]}
stringtie -p 32 -e -B \
-A //lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/gene/${trait[$i]}/${trait[$i]}.tab \
-G /lustre/nobackup/WUR/ABGC/ni010/seq/ref/GRCg7w.109.chr.gtf \
-o /lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.gtf \
/lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/8_fpkm/${trait[$i]}/${trait[$i]}.bam

let "i++"
done
```
#### Estimate PEERs effects

```R
library("peer")
library("preprocessCore")
TPM_tmp0<-read.table("seq/6_eQTL/9_eQTL/data/GE.txt",header=T)
row.names(TPM_tmp0)<-TPM_tmp0$gene_id
TPM_tmp0<-TPM_tmp0[,-1]
expr_matrix00<-TPM_tmp0
expr_matrix00_qn<-normalize.quantiles(as.matrix(TPM_tmp0))
covs = read.table("seq/EP_revise/peer/Cov.txt",header=T,row.names = 1)
covs_1<-apply(t(covs),2,as.double)

INT<-function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
expr_matrix00_qn_ind<-t(apply(expr_matrix00_qn,MARGIN=1,FUN=INT))
colnames(expr_matrix00_qn_ind)<-colnames(expr_matrix00)
rownames(expr_matrix00_qn_ind)<-rownames(expr_matrix00)
expr_peer<-expr_matrix00_qn_ind
expr_matrix00_qn_ind<-as.data.frame(expr_matrix00_qn_ind)
expr_matrix00_qn_ind$id<-row.names(expr_matrix00)
model=PEER()
PEER_setPhenoMean(model,as.matrix(t(expr_peer)))
PEER_setCovariates(model, as.matrix(covs_1))
dim(PEER_getPhenoMean(model))
PEER_setNk(model,30)
PEER_setNmax_iterations(model,1000)
PEER_getNk(model)
PEER_update(model)
factors=PEER_getX(model)
factors
Alpha = PEER_getAlpha(model)
pdf(paste0(30,".r_demo_covs.pdf"),width=8,height=8)
plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance", main="")
dev.off()
write.table(factors,"120.peer.covariance.txt",sep="\t",row.names=T,quote =FALSE)
write.table(expr_matrix00_qn_ind,"120.expression.qn_ind.txt",sep="\t",row.names=T,quote =FALSE)
````

#### eQTL-mapping analysis

```bash
module load R
Rscript eQTL_mapping.R
````
