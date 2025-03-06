#### Extract TWAS significant genes from eQTL-mapping results

```bash
for i in $(cat seq/7_TWAS/2_attp_tpmCorrectAgeandBreed/output_twas_correct_changebeta_rm900/2_attp_tpmCorrectAgeandBreed_changebeta_re900_all_sigGene); do
        grep -w "${i}" seq/6_eQTL/9_eQTL_new/cis-eqtl-all > seq/9_MR/output_eQTL_correct/${i}_snp.dat
done

##For additive SNP effects in Model A
module load R plink
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
j=0
while(($j<=12))
do
  mkdir seq/9_MR/output/${trait[$j]}
  awk '$5<=0.05 {print $2}' seq/9_MR/output_twas_correct/JTI_ovary_${trait[$j]}.txt > seq/9_MR/output/${trait[$j]}/sign_gene
  cp seq/9_MR/ld_effect_allele_rmZW_ref_re900 seq/9_MR/output/${trait[$j]}
  cp seq/9_MR/*.r seq/9_MR/output/${trait[$j]}
  cp seq/9_MR/gwas_add/MRJTI_gwas_${trait[$j]}_changebeta.dat seq/9_MR/output/${trait[$j]}/gwas_input
  cd seq/9_MR/output/${trait[$j]}
  for i in $(cat sign_gene); do
    awk '{print $0, ($3 / $4)}' seq/9_MR/output_eQTL_correct/${i}_snp.dat | awk '{print $1,$3,$7,$5}' | sed '1i snpid eqtl_beta eqtl_se eqtl_p' > eqtl_raw.dat
    Rscript clump.r
    awk '{print $1,$5,$6,$7,$8,$2,$3,$4}' clum_dat.dat | sed '1d' | sed '1i SNP CHR bp A1 A2 beta se P' > clum.dat
    plink --bfile seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref --allow-extra-chr --chr-set 41 --clump clum.dat --clump-p1 1 --clump-r2 0.2 --out clumped_${i}
    awk '{print $3}' clumped_${i}.clumped > clumped_${i}_snpid
    sed '1d' clumped_${i}_snpid | sed '1i snpid' > clumped_snpid
    Rscript combin.r
    sed 's/"//g' mrjti_input_raw | sed '1d' | sed '1i rsid effect_allele ldscore eqtl_beta eqtl_se eqtl_p gwas_beta gwas_se gwas_p' > mrjti_input_${i}.dat
    Rscript MR-JTI.r \
    --df_path mrjti_input_${i}.dat \
    --n_genes 1 \
    --out_path mrjti_result_${i}.csv
    #rm -f clumped_${i}.clumped clumped_${i}_snpid
  done
let "j++"
done


#For additive SNP effects in Model AD
module load R plink
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
j=0
while(($j<=12))
do
  mkdir seq/9_MR/output_add/${trait[$j]}
  awk '$5<=0.05 {print $2}' seq/9_MR/output_twas_correct/JTI_ovary_${trait[$j]}_add.txt > seq/9_MR/output_add/${trait[$j]}/sign_gene
  cp seq/9_MR/ld_effect_allele_rmZW_ref_re900 seq/9_MR/output_add/${trait[$j]}
  cp seq/9_MR/*.r seq/9_MR/output_add/${trait[$j]}
  cp seq/9_MR/gwas_dom/MRJTI_gwas_${trait[$j]}_add_changebeta.dat seq/9_MR/output_add/${trait[$j]}/gwas_input
  cd seq/9_MR/output_add/${trait[$j]}
  for i in $(cat sign_gene); do
    awk '{print $0, ($3 / $4)}' seq/9_MR/output_eQTL_correct/${i}_snp.dat | awk '{print $1,$3,$7,$5}' | sed '1i snpid eqtl_beta eqtl_se eqtl_p' > eqtl_raw.dat
    Rscript clump.r
    awk '{print $1,$5,$6,$7,$8,$2,$3,$4}' clum_dat.dat | sed '1d' | sed '1i SNP CHR bp A1 A2 beta se P' > clum.dat
    plink --bfile seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref --allow-extra-chr --chr-set 41 --clump clum.dat --clump-p1 1 --clump-r2 0.2 --out clumped_${i}
    awk '{print $3}' clumped_${i}.clumped > clumped_${i}_snpid
    sed '1d' clumped_${i}_snpid | sed '1i snpid' > clumped_snpid
    Rscript combin.r
    sed 's/"//g' mrjti_input_raw | sed '1d' | sed '1i rsid effect_allele ldscore eqtl_beta eqtl_se eqtl_p gwas_beta gwas_se gwas_p' > mrjti_input_${i}.dat
    Rscript MR-JTI.r \
    --df_path mrjti_input_${i}.dat \
    --n_genes 1 \
    --out_path mrjti_result_${i}.csv
    #rm -f clumped_${i}.clumped clumped_${i}_snpid
  done
let "j++"
done

##For dominance SNP effects in Model AD
module load R plink
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
j=0
while(($j<=12))
do
  mkdir seq/9_MR/output_dom/${trait[$j]}
  awk '$5<=0.05 {print $2}' seq/9_MR/output_twas_correct/JTI_ovary_${trait[$j]}_dom.txt > seq/9_MR/output_dom/${trait[$j]}/sign_gene
  cp seq/9_MR/ld_effect_allele_rmZW_ref_re900 seq/9_MR/output_dom/${trait[$j]}
  cp seq/9_MR/*.r seq/9_MR/output_dom/${trait[$j]}
  cp seq/9_MR/gwas_dom/MRJTI_gwas_${trait[$j]}_dom_changebeta.dat seq/9_MR/output_dom/${trait[$j]}/gwas_input
  cd seq/9_MR/output_dom/${trait[$j]}
  for i in $(cat sign_gene); do
    awk '{print $0, ($3 / $4)}' seq/9_MR/output_eQTL_correct/${i}_snp.dat | awk '{print $1,$3,$7,$5}' | sed '1i snpid eqtl_beta eqtl_se eqtl_p' > eqtl_raw.dat
    Rscript clump.r
    awk '{print $1,$5,$6,$7,$8,$2,$3,$4}' clum_dat.dat | sed '1d' | sed '1i SNP CHR bp A1 A2 beta se P' > clum.dat
    plink --bfile seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref --allow-extra-chr --chr-set 41 --clump clum.dat --clump-p1 1 --clump-r2 0.2 --out clumped_${i}
    awk '{print $3}' clumped_${i}.clumped > clumped_${i}_snpid
    sed '1d' clumped_${i}_snpid | sed '1i snpid' > clumped_snpid
    Rscript combin.r
    sed 's/"//g' mrjti_input_raw | sed '1d' | sed '1i rsid effect_allele ldscore eqtl_beta eqtl_se eqtl_p gwas_beta gwas_se gwas_p' > mrjti_input_${i}.dat
    Rscript MR-JTI.r \
    --df_path mrjti_input_${i}.dat \
    --n_genes 1 \
    --out_path mrjti_result_${i}.csv
    #rm -f clumped_${i}.clumped clumped_${i}_snpid
  done
let "j++"
done
```
