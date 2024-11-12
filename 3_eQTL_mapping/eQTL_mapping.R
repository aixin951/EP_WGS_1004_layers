````R
#!/usr/bin/env Rscrip/lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/9_eQTL/eqtl_mappingt
library(MatrixEQTL) 

base.dir = "/lustre/nobackup/WUR/ABGC/ni010/seq/6_eQTL/9_eQTL_new"

useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="")
snps_location_file_name = paste(base.dir, "/data/snpsloc.txt", sep="")

expression_file_name = paste(base.dir, "/data/GE.txt", sep="")
gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="")

covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="") 


output_file_name_cis = "cis-eqtl-all" 
output_file_name_tra = "trans-eqtl"


pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 0

errorCovariance = numeric()
#errorCovariance = read.table("Sample_Data/errorCovariance.txt");

#Distance for local gene-SNP pairs
cisDist = 1e6

snps = SlicedData$new()
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1       
snps$fileSkipColumns = 1    
snps$fileSliceSize = 2000;     
snps$LoadFile(SNP_file_name); 


gene = SlicedData$new()
gene$fileDelimiter = " "
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)


cvrt = SlicedData$new()
cvrt$fileDelimiter = " "
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}


snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name,  header = TRUE, stringsAsFactors = FALSE);


me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name  = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
````
