##Model A
module load R
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
j=0
while(($j<=12))
do
cp /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/chicken7w_dbSNP_ref.dat ${trait[$j]}
cp comb_rsid.R ${trait[$j]}
cd ${trait[$j]}
awk '{print $1,$4,$5}' add_q_new | sed 's/"//g' | sed '1i snpid beta se' > add_twas_changebeta
Rscript comb_rsid.R
awk '$4!="NA"' com_add_changebeta.dat | sed 's/"//g' | sed '1d' | awk '{print $2,$5,$6,$3,$4}' | sed '1i rsid beta se eff_allele ref_allele' | sed 's/:/./g' > JTI_${trait[$j]}_input_add_changebeta.dat

Rscript /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/src/predixcan_r.r \
         --asso_test \
         --db_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.db \
         --cov_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.cov \
         --gwas_path JTI_${trait[$j]}_input_add.dat \
         --gwas_variant_col rsid \
         --gwas_beta_col beta \
         --gwas_se_col se \
         --gwas_eff_allele_col eff_allele \
         --gwas_ref_allele_col ref_allele \
         --asso_out_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/output_twas_correct_changebeta_rm900/JTI_ovary_${trait[$j]}.txt \
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
cp /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/chicken7w_dbSNP_ref.dat ${trait[$j]}
cp comb_rsid.R ${trait[$j]}
cd ${trait[$j]}
awk '{print $1,$4,$5}' add_q_new | sed 's/"//g' | sed '1i snpid beta se' > add_twas_changebeta
awk '{print $1,$4,$5}' dom_q_new | sed 's/"//g' | sed '1i snpid beta se' > dom_twas_changebeta
Rscript comb_rsid.R
awk '$4!="NA"' com_add_changebeta.dat | sed 's/"//g' | sed '1d' | awk '{print $2,$5,$6,$3,$4}' | sed '1i rsid beta se eff_allele ref_allele' | sed 's/:/./g' > JTI_${trait[$j]}_input_add.dat
awk '$4!="NA"' com_dom_changebeta.dat | sed 's/"//g' | sed '1d' | awk '{print $2,$5,$6,$3,$4}' | sed '1i rsid beta se eff_allele ref_allele' | sed 's/:/./g' > JTI_${trait[$j]}_input_dom.dat

Rscript /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/predixcan_r.r \
         --asso_test \
         --db_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.db \
         --cov_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.cov \
         --gwas_path JTI_${trait[$j]}_input_add.dat \
         --gwas_variant_col rsid \
         --gwas_beta_col beta \
         --gwas_se_col se \
         --gwas_eff_allele_col eff_allele \
         --gwas_ref_allele_col ref_allele \
         --asso_out_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/output_twas_correct_changebeta_rm900/EW/JTI_ovary_${trait[$j]}_add.txt \
         #--parallel

Rscript /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/predixcan_r.r \
         --asso_test \
         --db_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.db \
         --cov_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/attp2_correct_TPM0.cov \
         --gwas_path JTI_${trait[$j]}_input_dom.dat \
         --gwas_variant_col rsid \
         --gwas_beta_col beta \
         --gwas_se_col se \
         --gwas_eff_allele_col eff_allele \
         --gwas_ref_allele_col ref_allele \
         --asso_out_path /lustre/nobackup/WUR/ABGC/ni010/seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/output_twas_correct_changebeta_rm900/EW/JTI_ovary_${trait[$j]}_dom.txt \
         #--parallel

cd ..
let "j++"
done
````
