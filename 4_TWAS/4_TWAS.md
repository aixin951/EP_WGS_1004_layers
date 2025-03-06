#### Correct TPM files with breeds, ages and PEERs

```R
##Rscript correct age and breed
# Read the data from the file
data <- read.table("GE_TPM0_3_transpose.txt", header = TRUE)
# Extract relevant columns
IID <- data$geneid
breed <- data$Breed
age <- data$Age
PEER1 <- data$P1
PEER2 <- data$P2
PEER3 <- data$P3
PEER4 <- data$P4
PEER5 <- data$P5
PEER6 <- data$P6
PEER7 <- data$P7
PEER8 <- data$P8
PEER9 <- data$P9
PEER10 <- data$P10
gene_expr <- data[, -(1:13)]  # Exclude columns IID, Breed, Age
# List to store model results and residuals
residuals_list <- list()
# Iterate through each gene expression level
for (gene_index in 1:ncol(gene_expr)) {
  gene_name <- colnames(gene_expr)[gene_index]
  # Fit a linear model with age and breed as fixed effects
  model <- lm(gene_expr[, gene_index] ~ age + breed + PEER1 + PEER2 + PEER3 + PEER4 + PEER5 + PEER6 + PEER7 + PEER8 + PEER9 + PEER10)
  # Extract residuals
  residuals <- residuals(model)

  # Store residuals in the list
  residuals_list[[gene_name]] <- residuals
}
# Create a data frame with IID, Breed, Age, and residuals for each gene
#residuals_df <- data.frame(IID = IID, Breed = breed, Age = age, P1 = PEER1, P2 = PEER2, P3 = PEER3, P4 = PEER4, P5 = PEER5, P6 = PEER6, P7 = PEER7, P8 = PEER8, P9 = PEER9, P10 = PEER10)
residuals_df <- data.frame(IID = IID)

for (gene_index in 1:length(residuals_list)) {
  gene_name <- names(residuals_list)[gene_index]
  residuals_df[[gene_name]] <- residuals_list[[gene_name]]
}
residuals_df_1 <- t(residuals_df)
# Save the residuals data frame to a file (e.g., "residuals_data.csv")
write.table(residuals_df_1, file = "residuals_data_TPM0.txt", col.names = FALSE,quote=FALSE)
```

#### TWAS model training

```bash
module load plink
main_dir=seq/7_TWAS/2_attp_tpmCorrectAgeandBreed
sub_job=$i
#step 1. model training
#ml GCC OpenMPI R
Rscript \${main_dir}/src/predixcan_r.r \\
       --model_training \\
        --main_dir \${main_dir} \\
        --plink_file_name \${main_dir}/bfile/free_final_reAn_reHom_imput_ref \\
        --expression_file_name \${main_dir}/expression/residuals_data_TPM0.txt \\
        --subjob_id \${sub_job} \\
        --n_genes_for_each_subjob 1000 \\
        --annotation_file_name \${main_dir}/anno/GRCg7w_110_twas.gft
# If a '--parallel' flag is added, max(n-1,1) core(s) will be used for parallel model training, where n is the number of available cores.

main_dir=seq/7_TWAS/2_attp_tpmCorrectAgeandBreed
model_name=attp2_correct_TPM0
Rscript ${main_dir}/src/predixcan_r.r \
         --generate_db_and_cov \
         --main_dir ${main_dir} \
         --plink_file_name ${main_dir}/bfile/free_final_reAn_reHom_imput_ref \
         --expression_file_name ${main_dir}/expression/residuals_data_TPM0.txt \
         --annotation_file_name ${main_dir}/anno/GRCg7w_110_twas.gft \
         --output_file_name ${model_name}
```

#### TWAS analysis

```bash
##Model A
module load R
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
j=0
while(($j<=12))
do
cp seq/5_GWAS/4_association/wombat/chicken7w_dbSNP_ref.dat ${trait[$j]}
cp comb_rsid.R ${trait[$j]}
cd ${trait[$j]}
awk '{print $1,$4,$5}' add_q_new | sed 's/"//g' | sed '1i snpid beta se' > add_twas_changebeta
Rscript comb_rsid.R
awk '$4!="NA"' com_add_changebeta.dat | sed 's/"//g' | sed '1d' | awk '{print $2,$5,$6,$3,$4}' | sed '1i rsid beta se eff_allele ref_allele' | sed 's/:/./g' > JTI_${trait[$j]}_input_add_changebeta.dat

Rscript seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/src/predixcan_r.r \
         --asso_test \
         --db_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.db \
         --cov_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.cov \
         --gwas_path JTI_${trait[$j]}_input_add.dat \
         --gwas_variant_col rsid \
         --gwas_beta_col beta \
         --gwas_se_col se \
         --gwas_eff_allele_col eff_allele \
         --gwas_ref_allele_col ref_allele \
         --asso_out_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/output_twas_correct_changebeta_rm900/JTI_ovary_${trait[$j]}.txt \
         --parallel
cd ..
let "j++"
done

##comb_rsid.R
library(data.table)
library(dplyr)
pvalue_add <- fread("add_twas_changebeta", header=T)
rsid <- fread("chicken7w_dbSNP_ref.dat",header=T)
com_add <- inner_join(rsid,pvalue_add,by="snpid")
write.table(com_add,"com_add_changebeta.dat",row.names=FALSE)


##Model AD
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
j=0
while(($j<=12))
do
cp seq/5_GWAS/4_association/wombat/chicken7w_dbSNP_ref.dat ${trait[$j]}
cp comb_rsid.R ${trait[$j]}
cd ${trait[$j]}
awk '{print $1,$4,$5}' add_q_new | sed 's/"//g' | sed '1i snpid beta se' > add_twas_changebeta
awk '{print $1,$4,$5}' dom_q_new | sed 's/"//g' | sed '1i snpid beta se' > dom_twas_changebeta
Rscript comb_rsid.R
awk '$4!="NA"' com_add_changebeta.dat | sed 's/"//g' | sed '1d' | awk '{print $2,$5,$6,$3,$4}' | sed '1i rsid beta se eff_allele ref_allele' | sed 's/:/./g' > JTI_${trait[$j]}_input_add.dat
awk '$4!="NA"' com_dom_changebeta.dat | sed 's/"//g' | sed '1d' | awk '{print $2,$5,$6,$3,$4}' | sed '1i rsid beta se eff_allele ref_allele' | sed 's/:/./g' > JTI_${trait[$j]}_input_dom.dat

Rscript seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/predixcan_r.r \
         --asso_test \
         --db_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.db \
         --cov_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.cov \
         --gwas_path JTI_${trait[$j]}_input_add.dat \
         --gwas_variant_col rsid \
         --gwas_beta_col beta \
         --gwas_se_col se \
         --gwas_eff_allele_col eff_allele \
         --gwas_ref_allele_col ref_allele \
         --asso_out_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/output_twas_correct_changebeta_rm900/EW/JTI_ovary_${trait[$j]}_add.txt \
         #--parallel

Rscript seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/predixcan_r.r \
         --asso_test \
         --db_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.db \
         --cov_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.cov \
         --gwas_path JTI_${trait[$j]}_input_dom.dat \
         --gwas_variant_col rsid \
         --gwas_beta_col beta \
         --gwas_se_col se \
         --gwas_eff_allele_col eff_allele \
         --gwas_ref_allele_col ref_allele \
         --asso_out_path seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/output_twas_correct_changebeta_rm900/EW/JTI_ovary_${trait[$j]}_dom.txt \
         #--parallel

cd ..
let "j++"
done

```
