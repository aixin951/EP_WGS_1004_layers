#### Make SNPCounts file

generate additive and additive-dominance SNPCounts files

```bash
##additive SNPCounts
plink --bfile free_final_reAn_reHom_imput_ref --allow-extra-chr --chr-set 41 --keep-allele-order --recode A --out free_final_reAn_reHom_imput_ref_SNPCounts_add
sort -k 1n free_final_reAn_reHom_imput_ref_SNPCounts_add.raw > free_final_reAn_reHom_imput_ref_SNPCounts_add_sort.raw
cut --complement -d' ' -f1,2,3,4,5,6 free_final_reAn_reHom_imput_ref_SNPCounts_add_sort.raw | sed '1d' > free_final_reAn_reHom_imput_ref_SNPCounts_add_sort_1.raw
module load R
Rscript trans.R
sed 's/ //g' free_final_reAn_reHom_imput_ref_SNPCounts_add_sort_1.dat > SNPCounts_add_new.dat
##trans.R
#!/usr/bin/env Rscript
library(data.table)
data <- fread("free_final_reAn_reHom_imput_ref_SNPCounts_add_sort_1.raw")
data1 <- t(data)
write.table(data1,"free_final_reAn_reHom_imput_ref_SNPCounts_add_sort_1.dat",row.names=FALSE,col.names=FALSE)

##additive dominance SNPCounts file
plink --bfile free_final_reAn_reHom_imput_ref --allow-extra-chr --chr-set 41 --recode AD --keep-allele-order --out free_final_reAn_reHom_imput_ref_SNPCounts_add_dom
sort -k 1n free_final_reAn_reHom_imput_ref_SNPCounts_add_dom.raw | cut --complement -d' ' -f1,2,3,4,5,6 | sed '1d' > free_final_reAn_reHom_imput_ref_SNPCounts_add_dom_sort.raw
module load R
Rscript trans_dom.R
sed 's/ //g' free_final_reAn_reHom_imput_ref_SNPCounts_add_dom_1.dat > SNPCounts_add_dom_new.dat
##trans_dom.R
#!/usr/bin/env Rscript
library(data.table)
data <- fread("free_final_reAn_reHom_imput_ref_SNPCounts_add_dom_sort.raw")
data1 <- t(data)
write.table(data1,"free_final_reAn_reHom_imput_ref_SNPCounts_add_dom_1.dat",row.names=FALSE,col.names=FALSE)
```

#### Carry out GWAS with additive model and additive-dominance model

Additive model (Model A)

```bash
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
i=6
j=0
m=1
while(($i<=18))&&(($j<=12))&&(($m<=13))

do
mkdir ${trait[$j]}
cp /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref.* ${trait[$j]}
cp calc_grm.inp ${trait[$j]}
cp trans.R ${trait[$j]}
cp wombat_1.par ${trait[$j]}
cp wombat_2.par ${trait[$j]}
cp id_933_sort ${trait[$j]}

awk -v col=$i '$col!=-99 {print $1,$2,$3,$4,$5,$col}' wombat_raw_933.dat | sed '1i animal PE dom breed cage trait' > ${trait[$j]}/wombat_trait.dat

cd ${trait[$j]}
awk '{print $1}' wombat_trait.dat | sed '1d' | uniq > trait_id
awk '$0=NR" "$0' trait_id > animal.codes
paste trait_id trait_id | sed 's/\t/ /g' > ${trait[$j]}_id_plink
module load plink
plink --bfile free_final_reAn_reHom_imput_ref --allow-extra-chr --chr-set 41 --keep-allele-order --keep ${trait[$j]}_id_plink --make-bed --out ${trait[$j]}

##prepare calc_grm input
awk '$1 ~ /^Z|^W/ {print $2}' ${trait[$j]}.bim > WZ_SNP_id
plink --bfile ${trait[$j]} --allow-extra-chr --chr-set 41 --exclude WZ_SNP_id --keep-allele-order --recode A --out ${trait[$j]}_reWZ
sort -k 1n ${trait[$j]}_reWZ.raw | cut --complement -d' ' -f1,3,4,5,6 | sed '1d' > reWZ_calc.raw
rm ${trait[$j]}_reWZ.raw

module unload gcc
module load SHARED/calc_grm/main
module load intel/compiler/64/2017
calc_grm

awk '{print $1,$2,$3}' G.grm > G_3.grm
awk '{print $2}' reWZ_calc.raw > snp_t
paste  trait_id snp_t | sed 's/\t/ /g' | sed '1i animal snp'> data_snp
rm snp_t

module load R
Rscript trans.R

awk '{print $2,$1,$3}' G_asreml.giv > tmp
cat deter tmp > animal.gin
rm tmp

awk '$3=="NA" {print FNR}' com_id > abs_id
tr '\n' ' ' <  abs_id > abs_id_t
sed -e '$a\' abs_id_t | sed 's/\s*$//g' | sed 's/ /,/g' > abs_id_t_2
remo=$(cat abs_id_t_2)
cut -b $remo --complement /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/SNPCounts_add_new.dat > SNPCounts.dat

var=$(cat var)
sed -i "24c $var" wombat_1.par
sed -i "26c $var" wombat_1.par
/lustre/nobackup/WUR/ABGC/ni010/Software/wombat/wombat --dense --threads8 wombat_1.par

animal=$(grep "2  COVS" SumEstimates.out | awk -F ' ' '{print $6}')
re=$(grep "1  COVS" SumEstimates.out | awk -F ' ' '{print $6}')
sed -i "26c $animal" wombat_2.par
sed -i "28c $re" wombat_2.par
/lustre/nobackup/WUR/ABGC/ni010/Software/wombat/wombat --snap --dense --threads8 wombat_2.par
let "i++"
let "j++"
let "m++"
done
```

Additive-dominance model (Model AD)

```bash
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
i=6
j=0
m=1
while(($i<=18))&&(($j<=12))&&(($m<=13))
do
mkdir ${trait[$j]}
cp /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref.* ${trait[$j]}
cp /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/calc_grm.inp ${trait[$j]}
cp /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/id_933_sort ${trait[$j]}
cp trans.R ${trait[$j]}
cp wombat_1.par ${trait[$j]}
cp wombat_2.par ${trait[$j]}

awk -v col=$i '$col!=-99 {print $1,$2,$3,$4,$5,$col}' wombat_raw_933.dat | sed '1i animal PE dom breed cage trait' > ${trait[$j]}/wombat_trait.dat

cd ${trait[$j]}
awk '{print $1}' wombat_trait.dat | sed '1d' | uniq > trait_id
awk '$0=NR" "$0' trait_id > animal.codes
cp animal.codes dom.codes
paste trait_id trait_id | sed 's/\t/ /g' > ${trait[$j]}_id_plink
module load plink
plink --bfile free_final_reAn_reHom_imput_ref --allow-extra-chr --chr-set 41 --keep-allele-order --keep ${trait[$j]}_id_plink --make-bed --out ${trait[$j]}_raw
awk '$1 ~/^30/' ${trait[$j]}_id_plink > ${trait[$j]}_id_plink_WW
awk '$1 ~/^50/' ${trait[$j]}_id_plink > ${trait[$j]}_id_plink_YY
plink --bfile ${trait[$j]}_raw --allow-extra-chr --chr-set 41 --keep ${trait[$j]}_id_plink_WW --keep-allele-order --hardy --out freq_WW
plink --bfile ${trait[$j]}_raw --allow-extra-chr --chr-set 41 --keep ${trait[$j]}_id_plink_YY --keep-allele-order --hardy --out freq_YY
awk '{print $2,$6}' freq_WW.hwe | sed '1d' | sed 's/\// /g' | awk '$2<=5 || $3<=5 || $4<=5 {print $1}' > freq_WW_maf
awk '{print $2,$6}' freq_YY.hwe | sed '1d' | sed 's/\// /g' | awk '$2<=5 || $3<=5 || $4<=5 {print $1}' > freq_YY_maf
comm -12 <(sort freq_WW_maf) <(sort freq_YY_maf) > freq_YY_WW_id
plink --bfile ${trait[$j]}_raw --allow-extra-chr --chr-set 41 --exclude freq_YY_WW_id --keep-allele-order --make-bed --out ${trait[$j]}

##prepare cacl_grm input
awk '$1 ~ /^Z|^W/ {print $2}' ${trait[$j]}.bim > WZ_SNP_id
plink --bfile ${trait[$j]} --allow-extra-chr --chr-set 41 --exclude WZ_SNP_id --recode A --keep-allele-order --out ${trait[$j]}_reWZ
sort -k 1n ${trait[$j]}_reWZ.raw | cut --complement -d' ' -f1,3,4,5,6 | sed '1d' > reWZ_calc.raw

awk '$1 !~/^Z|^W/' freq_YY_WW_id > freq_YY_WW_id_ZW
wc -l freq_YY_WW_id_ZW | awk '{print $1}' > rm_number
rm_number_1=$(cat rm_number)
expr 11892646 - ${rm_number_1} > snp_number
snp_number_1=$(cat snp_number)
sed -i "1c $snp_number_1" calc_grm.inp

##G-matrix
module unload gcc
module load SHARED/calc_grm/main
module load intel/compiler/64/2017
calc_grm
awk '{print $1,$2,$3}' G.grm > G_3.grm

mv G.grm G_add.gom
mv G_asreml.giv G_asreml_add.giv
rm calculated_all_freq.dat
rm high_grm_coefs.log
rm ID_vs_row_number_G.txt
rm genomic_inbr_coef.dat
sed -i "7c D ASReml" calc_grm.inp

##D-matrix
calc_grm
awk '{print $1,$2,$3}' D.grm > D_3.grm

awk '{print $2,$3}' reWZ_calc.raw > snp_t
paste  trait_id snp_t | sed 's/\t/ /g' | sed '1i animal snp znp'> data_snp
rm snp_t

module load R
Rscript trans.R

awk '{print $2,$1,$3}' G_asreml.giv > tmp
cat deter_G tmp > animal.gin
rm tmp
awk '{print $2,$1,$3}' D_asreml.giv > tmp
cat deter_D tmp > dom.gin
rm tmp

awk 'FNR==NR{a[$0]++};FNR!=NR{if($2 in a){print}}' freq_YY_WW_id_ZW /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/snpid_add_dom_transpose_rmWZ_new | awk '{print $1}' > line_number
awk 'NR==FNR{l[$0];next;} !(FNR in l)' line_number /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/SNPCounts_add_dom_new_rmWZ.dat > SNPCounts_rmSNP.dat
awk '$3=="NA" {print FNR}' com_id > abs_id
tr '\n' ' ' <  abs_id > abs_id_t
sed -e '$a\' abs_id_t | sed 's/\s*$//g' | sed 's/ /,/g' > abs_id_t_2
remo=$(cat abs_id_t_2)
cut -b $remo --complement SNPCounts_rmSNP.dat > SNPCounts.dat

var=$(cat var)
var_dom=$(cat var_dom)
sed -i "26c $var" wombat_1.par
sed -i "28c $var_dom" wombat_1.par
sed -i "30c $var" wombat_1.par

echo -ne '\n' | /lustre/nobackup/WUR/ABGC/ni010/Software/wombat/wombat --dense --threads8 wombat_1.par

animal=$(grep "2  COVS" SumEstimates.out | awk -F ' ' '{print $6}')
dom=$(grep "3  COVS" SumEstimates.out | awk -F ' ' '{print $6}')
re=$(grep "1  COVS" SumEstimates.out | awk -F ' ' '{print $6}')
sed -i "29c $animal" wombat_2.par
sed -i "31c $dom" wombat_2.par
sed -i "33c $re" wombat_2.par

echo -ne '\n' | /lustre/nobackup/WUR/ABGC/ni010/Software/wombat/wombat --snap --dense --threads8 wombat_2.par" > ${trait[$j]}.sh

let "i++"
let "j++"
let "m++"
done
```

#### plot the mahattan plot

Additive model (Model A)
```bash
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
j=0
while(($j<=12))
do

cp *.R ${trait[$j]}
cd ${trait[$j]}
awk '$1 !~ /^Z|^W/ {print $2,$1,$4}' /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref.bim | sed '1i rs chr ps' > pos_raw
grep -v -wFf freq_YY_WW_id_ZW pos_raw > pos
awk '{print $1 }' pos | sed '1d' > pos_id
awk '{print $1,$2,$3}' SNPSolutions.dat | sed -n '1d' | sed '1i Beta SE t_value' > add_t
paste pos add_t | sed 's/\t/ /g' > add
snp_number_1=$((2 * $(cat snp_number)))
sed -i "7c add_p\$p_value = 2*pt(q=abs(add_p\$t_value),df=$snp_number_1,lower.tail=FALSE)" q_m_new_0.01.R
Rscript q_m_new_0.01.R

cd ..
let "j++"
done
```

Additive-dominance model (Model AD)
```bash
trait=(CEN200 CEN300 CEN400 CEN500 CEN600 CEN700 EN300 EN400 EN500 EN600 EN700 EN300_500 EN500_700)
j=0
while(($j<=12))
do

cp *.R ${trait[$j]}
cd ${trait[$j]}
awk '$1 !~ /^Z|^W/ {print $2,$1,$4}' /lustre/nobackup/WUR/ABGC/ni010/seq/5_GWAS/4_association/wombat/free_final_reAn_reHom_imput_ref.bim | sed '1i rs chr ps' > pos_raw
grep -v -wFf freq_YY_WW_id_ZW pos_raw > pos
awk '{print $1 }' pos | sed '1d' > pos_id
awk '{print $1,$2,$3}' SNPSolutions.dat | sed -n '1~2p' | sed '1i Beta SE t_value' > add_t
awk '{print $1,$2,$3}' SNPSolutions.dat | sed -n '2~2p' | sed '1i Beta SE t_value' > dom_t
paste pos add_t | sed 's/\t/ /g' > add
paste pos dom_t | sed 's/\t/ /g' > dom
snp_number_1=$((2 * $(cat snp_number)))
sed -i "7c add_p\$p_value = 2*pt(q=abs(add_p\$t_value),df=$snp_number_1,lower.tail=FALSE)" q_m_new_0.01.R
sed -i "31c dom_p\$p_value = 2*pt(q=abs(dom_p\$t_value),df=$snp_number_1,lower.tail=FALSE)" q_m_new_0.01.R
Rscript q_m_new_0.01.R

cd ..
let "j++"
done
```
