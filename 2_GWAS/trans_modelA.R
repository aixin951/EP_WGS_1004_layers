library(dplyr)
id_934 <- read.table("id_933_sort",header=F)
id_trait <- read.table("trait_id",header=F)
id_934$value <- 1
id_trait$value <- 1
com <- left_join(id_934,id_trait,by="V1")
write.table(com,"com_id",row.names=FALSE,col.names=FALSE)

snp <- read.table("data_snp",header=T)
wombat_raw <- read.table("wombat_trait.dat",header=T)
wombat <- left_join(wombat_raw,snp,by="animal")
write.table(wombat,"wombat.dat",row.names=FALSE,col.names=FALSE)

var <- var(wombat[,6])
var_1 <- var*0.03
write.table(var_1,"var",row.names=FALSE,col.names=FALSE)

library(tidyr)
makeSymm <- function(m) {
 m[upper.tri(m)] <- t(m)[upper.tri(m)]
 return(m)
 }
G_m <- read.table("G_3.grm",header=F)
Gm <- spread(G_m,V2,V3)
Gm_1 <- Gm[,-1]
full_G <- makeSymm(as.matrix(Gm_1))
deter <- log(det(full_G),10)
write.table(deter,"deter",row.names=FALSE,col.names=FALSE)
