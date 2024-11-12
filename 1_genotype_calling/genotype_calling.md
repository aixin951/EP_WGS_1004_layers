#### Using fastqc to get the quality of raw sequencing data

```bash
for i in $(cat sample.list); do
fastqc /lustre/nobackup/WUR/ABGC/ni010/seq/0_raw_data/${i}/${i}_1.fq.gz -o /lustre/nobackup/WUR/ABGC/ni010/seq/1_fastqc/ -t 8
done
```

#### Trim to clean reads

```bash
for i in $(cat sample.list); do
trimmomatic PE -phred33 -threads 10 /lustre/nobackup/WUR/ABGC/ni010/seq/0_raw_data/${i}/${i}_1.fq.gz /lustre/nobackup/WUR/ABGC/ni010/seq/0_raw_data/${i}/${i}_2.fq.gz /lustre/nobackup/WUR/ABGC/ni010/seq/2_clean_data/${i}/${i}_1.clean.fq.gz /lustre/nobackup/WUR/ABGC/ni010/seq/2_clean_data/unpaired/${i}/${i}_1.unpaired.fq.gz /lustre/nobackup/WUR/ABGC/ni010/seq/2_clean_data/${i}/${i}_2.clean.fq.gz /lustre/nobackup/WUR/ABGC/ni010/seq/2_clean_data/unpaired/${i}/${i}_2.unpaired.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
done
```

#### Using the population mapping pipeline mapping to the reference genome

We used the GRCg.7w as the reference genome, and here is the link to the pipeline: https://wiki.anunna.wur.nl/index.php/Population_mapping_pipeline

Edit the config.yaml file:

```bash
ASSEMBLY: /lustre/nobackup/WUR/ABGC/ni010/seq/GRCg7w.dna.toplevel.fa.gz
OUTDIR: /lustre/nobackup/WUR/ABGC/ni010/seq/3_mapping
PATHS_WITH_FILES: /lustre/nobackup/WUR/ABGC/ni010/seq/2_clean_data/sub/
```

#### Sort the bam files, add index and get the quality of mapping files

```bash
##sort
for i in $(cat sample.list); do
module load samtools
samtools sort -m 2G -@ 7 -O bam ${i}.bam > ${i}.sorted.bam
samtools index -@ 16 ${i}.sorted.bam
done

##mapping quality
for i in $(cat sample.list); do
qualimap bamqc -bam /lustre/nobackup/WUR/ABGC/ni010/seq/3_mapping/sort_bam/${i}.sorted.bam --java-mem-size=16G -nt 8 -outformat PDF -outdir /lustre/nobackup/WUR/ABGC/ni010/seq/3_mapping/mapping_stats/${i}
done
```

#### Call variants through population variants calling pipeline

Call variants for each chromosome using freebayes

Here is the link to the pipeline: https://wiki.anunna.wur.nl/index.php/Population_variant_calling_pipeline

```bash
module load freebayes bcftools vcflib python/2.7.15 samtools
/lustre/nobackup/WUR/ABGC/ni010/seq/3_mapping/processed_reads/freebayes-parallel.sh <(/lustre/nobackup/WUR/ABGC/ni010/seq/3_mapping/processed_reads/fasta_generate_regions.py /lustre/nobackup/WUR/ABGC/ni010/seq/3_mapping/GRCg7w_2.dna.toplevel.fa.fai 100000) 16 \
-f /lustre/nobackup/WUR/ABGC/ni010/seq/3_mapping/GRCg7w.dna.toplevel.fa \
--use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 \
-L bam_list.txt | vcffilter -f 'QUAL > 20' | bgzip -c > /lustre/nobackup/WUR/ABGC/ni010/seq/4_SNP_calling/chr2.vcf.gz
tabix -p vcf /lustre/nobackup/WUR/ABGC/ni010/seq/4_SNP_calling/chr2.vcf.gz
```

#### Combine VCF files, transfer to PLINK format, and filter

attach the calc_grm_WW.inp 

```
12528961
free_WW.bed PED.txt
genotypes
1
vanraden
grm
G ASReml
print_giv=asc
print_geno=no
8
```



```bash
##combine VCF files
module load vcftools plink
ls *.vcf.gz | xargs vcf-concat > free.vcf

module load plink R plink/2.0
plink --vcf free.vcf --chr-set 41 --allow-extra-chr --make-bed --out free
##extract snps only
plink2 --bfile free --chr-set 41 --allow-extra-chr --snps-only --make-bed --out free1_1
##rename the variants
plink2 --bfile free1_1 --chr-set 41 --allow-extra-chr --set-all-var-ids @:# --make-bed --out free1_2
##remove the missing individual and snps 
plink --bfile free1_2 --chr-set 41 --allow-extra-chr --geno 0.05 --mind 0.05 --make-bed --out free_2
##remove the maf
plink --bfile free_2 --chr-set 41 --allow-extra-chr --maf 0.005 --make-bed --out free_3

##only for pure lines individuals a test for Hardy Weinberg equilibrium P < 1e–4
plink --bfile free_3 --chr-set 41 --allow-extra-chr --keep WW_list ---hardy --out WW
plink --bfile free_3 --chr-set 41 --allow-extra-chr --keep YY_list ---hardy --out YY
awk '$9 < 1e-4' YY.hwe | awk '{print $2}' > YY_list
awk '$9 < 1e-4' WW.hwe | awk '{print $2}' > WW_list
sed -n '/Z/!p' WW_list | sed -n '/W/!p' > WW_list_remWZ
sed -n '/Z/!p' YY_list | sed -n '/W/!p' > YY_list_remWZ
cat WW_list_remWZ YY_list_remWZ | sort | uniq > WW_YY_list_remWZ_union
plink --bfile free_3 --allow-extra-chr --chr-set 41 --nonfounders --exclude WW_YY_list_remWZ_union --make-bed --out free_final

##check for inconsistencies between pedigree and genomic information within WW, YY, and crossbreds
module load SHARED/calc_grm/main
module load intel/compiler/64/2017
calc_grm --gvsa --graph

##remove conflicting animals 
plink --bfile free_final --allow-extra-chr --chr-set 41 --remove remove_an --make-bed --out free_final_reAn

##remove SNPs which is homozygous for both purelines
plink --bfile free_final_reAn --allow-extra-chr --chr-set 41 --keep YY_id --nonfounders --freq --out YY
plink --bfile free_final_reAn --allow-extra-chr --chr-set 41 --keep WW_id --nonfounders --freq --out WW
awk '$5==0 {print $2}' WW.frq > WW_MAF_0_id
awk '$5==0 {print $2}' YY.frq > YY_MAF_0_id
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' WW_MAF_0_id YY_MAF_0_id > hom_WW_YY
plink --bfile free_final_reAn --allow-extra-chr --chr-set 41 --exclude hom_WW_YY --make-bed --out free_final_reAn_reHom
```

#### Impute with beagle

Beagle was used to impute the missing genotype for each chromosome

```bash
plink --bfile /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/beagle/free_final_reAn_reHom --allow-extra-chr --chr-set 41 --recode vcf-iid --chr $i --out chr$i
module load beagle
java -Xmx50g -Xss5m -jar beagle.jar gt=chr$i.vcf out=chr${i}_impu.vcf nthreads=8

##combine imputed VCF files
module load vcftools plink
ls /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/beagle/*/*.vcf.gz  | xargs vcf-concat > free_final_reAn_reHom_imput.vcf
plink --vcf free_final_reAn_reHom_imput.vcf --chr-set 41 --allow-extra-chr --make-bed --out free_final_reAn_reHom_imput

##set the A1 to alternative allele
plink2 --bfile initial_new_allorder --allow-extra-chr --chr-set 41 --snps-only --make-bed --set-all-var-ids @:# --out initial_new_allorder_snps
awk '{print $2,$5}' initial_new_allorder_snps.bim > ALT_allele
plink --bfile c --allow-extra-chr -chr-set 41 --a1-allele ALT_allele --make-bed --out free_final_reAn_reHom_imput_ref
```
