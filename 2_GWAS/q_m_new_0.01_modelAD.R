#!/usr/bin/env Rscript
library(data.table)
library(qvalue)
library(CMplot)
add_p <- fread("add")
add_p$p_value = 2*pt(q=abs(add_p$t_value),df=11892645,lower.tail=FALSE)
q_value <- qvalue(add_p$p_value)
add_q <- cbind(add_p,q_value$qvalues)
min_value_add <- min(add_q$V2)
add_q$abs_0.01 <- abs(add_q$V2-0.01)
if(min_value_add<0.015){
        print('sign_0.01_add')
thre_0.01_add <- mean(add_q[add_q$abs_0.01 == min(add_q$abs_0.01),]$p_value)
} else {print('no_sign_0.01_add')}
write.table(add_q,"add_q_new_new",row.names=FALSE,col.names=FALSE)
if(min_value_add<0.015){
write.table(thre_0.01_add,"thre_0.01_add_new",row.names=FALSE,col.names=FALSE)
} else {print('no_sign_0.01_add_thre')}

CMplot(add_p[,-c(4,5,6)],plot.type="q",conf.int=TRUE,box=FALSE,file="jpg",file.name="add",dpi=300,file.output=TRUE,verbose=TRUE,width=5,height=5)
if(min_value_add<0.015)
{
       CMplot(add_p[,-c(4,5,6)],plot.type="m",LOG10=TRUE,col=c("#B1B0DD","#B3E2B4"),threshold=thre_0.01_add,threshold.lty=2,threshold.lwd=1.5, signal.col="red",threshold.col="grey",file="tif",file.name="add",dpi=500,file.output=TRUE,verbose=TRUE,width=18,height=8)
} else {
     CMplot(add_p[,-c(4,5,6)],plot.type="m",LOG10=TRUE,col=c("#B1B0DD","#B3E2B4"),file="tif",file.name="add",dpi=500,file.output=TRUE,verbose=TRUE,width=18,height=8)}

dom_p <- fread("dom")
dom_p$p_value = 2*pt(q=abs(dom_p$t_value),df=11892645,lower.tail=FALSE)
q_value <- qvalue(dom_p$p_value)
dom_q <- cbind(dom_p,q_value$qvalues)
min_value_dom <- min(dom_q$V2)
dom_q$abs_0.01 <- abs(dom_q$V2-0.01)
if(min_value_dom<0.015){
        print('sign_0.01_dom')
thre_0.01_dom <- mean(dom_q[dom_q$abs_0.01 == min(dom_q$abs_0.01),]$p_value)
} else {print('no_sign_0.01_dom')}
write.table(dom_q,"dom_q_new_new",row.names=FALSE,col.names=FALSE)
if(min_value_dom<0.015){
write.table(thre_0.01_dom,"thre_0.01_dom_new",row.names=FALSE,col.names=FALSE)
} else {print('no_sign_0.01_dom_thre')}

CMplot(dom_p[,-c(4,5,6)],plot.type="q",conf.int=TRUE,box=FALSE,file="jpg",file.name="dom",dpi=100,file.output=TRUE,verbose=TRUE,width=5,height=5)
if(min_value_dom<0.015)
{
      CMplot(dom_p[,-c(4,5,6)],plot.type="m",LOG10=TRUE,col=c("#B1B0DD","#B3E2B4"),threshold=thre_0.01_dom,threshold.lty=2,threshold.lwd=1.5, signal.col="red",threshold.col="grey",file="tif",file.name="dom",dpi=500,file.output=TRUE,verbose=TRUE,width=18,height=8)
} else {
      CMplot(dom_p[,-c(4,5,6)],plot.type="m",LOG10=TRUE,col=c("#B1B0DD","#B3E2B4"),file="tif",file.name="dom",dpi=500,file.output=TRUE,verbose=TRUE,width=18,height=8)}
