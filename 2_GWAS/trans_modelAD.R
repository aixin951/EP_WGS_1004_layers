library(dplyr)
id_934 <- read.table("id_933_sort",header=F)
id_trait <- read.table("trait_id",header=F)
id_934$value <- 1
id_trait$value <- 1
com <- left_join(id_934,id_trait,by="V1")
write.table(com,"com_id",row.names=FALSE,col.names=FALSE)

#pca <- read.table("pca.txt",heaer=T)
snp <- read.table("data_snp",header=T)
wombat_raw <- read.table("wombat_trait.dat",header=T)
wombat <- left_join(wombat_raw,snp,by="animal")
#wombat <- left_join(wombat_1,pca,by="IID")
write.table(wombat,"wombat.dat",row.names=FALSE,col.names=FALSE)

var <- var(wombat[,6])
var_1 <- var*0.03
var_2 <- var*0.002
write.table(var_1,"var",row.names=FALSE,col.names=FALSE)
write.table(var_2,"var_dom",row.names=FALSE,col.names=FALSE)


library(tidyr)
makeSymm <- function(m) {
 m[upper.tri(m)] <- t(m)[upper.tri(m)]
 return(m)
 }
G_m <- read.table("G_3.grm",header=F)
Gm <- spread(G_m,V2,V3)
Gm_1 <- Gm[,-1]
full_G <- makeSymm(as.matrix(Gm_1))
deter_G <- log(det(full_G),10)
write.table(deter_G,"deter_G",row.names=FALSE,col.names=FALSE)
D_m <- read.table("D_3.grm",header=F)
Dm <- spread(D_m,V2,V3)
Dm_1 <- Dm[,-1]
full_D <- makeSymm(as.matrix(Dm_1))
deter_D <- log(det(full_D),10)
write.table(deter_D,"deter_D",row.names=FALSE,col.names=FALSE)
