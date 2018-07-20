#Statistical analysis related to the article : 
#--- Molecular and immune classification of glioblastoma define patient prognosis by Klopfenstein et al.



#---  First step is to load the different datasets and correct the batch effect between different affymetrix platforms



#Loading the glioblastoma cell lines samples 

setwd("~/Recherche/Glioblastome/Glio/Glio/Donnees_exportees_analysees")
cell_line=read.table("cell_line_processed.txt",sep="\t",header=TRUE,row.names=1)

#Loading of the TCGA glioblastoma samples 
setwd("~/Recherche/Glioblastome/Glio/Glio/Donnees_exportees_analysees/TCGA")
TCGA=read.table("TCGA_gbm.txt",sep="\t",header=TRUE,row.names=1)


#Loading of the clinical data from TCGA samples (discovery cohort)
info_clin=read.table("TCGA_clinical.txt",sep="\t",header=TRUE,row.names=1)

colnames(TCGA)=gsub("\\.","-",colnames(TCGA))
TCGA=TCGA[,which(colnames(TCGA)%in%info_clin[,1])]
info_clin=info_clin[which(info_clin[,1]%in%colnames(TCGA)),]

colnames(TCGA)==info_clin[,1]


#we get rid of the patients who have G-CIMP subtypes 

sample=which(info_clin$EXPRESSION_SUBTYPE=="G-CIMP")


TCGA=TCGA[,-sample]


#Loading of the Rembrandt glioblastoma samples (validation cohort)
setwd("~/Recherche/Glioblastome/Glio/Glio/Donnees_exportees_analysees/Rembrandt")
Rembrandt=read.table("Rembrandt_genes.txt",sep="\t",header=TRUE,row.names=1)




#keeping commun set of genes between the three datasets : cell line, TCGA, Rembrandt
cell_line=cell_line[rownames(cell_line)%in%rownames(TCGA),]
TCGA=TCGA[rownames(TCGA)%in%rownames(cell_line),]

Rembrandt=Rembrandt[rownames(Rembrandt)%in%rownames(TCGA),]





#--- Batch effect correction between Affymetrix platforms
data=cbind(TCGA,Rembrandt,cell_line)

library(sva)
batch=c(rep("h133a",ncol(TCGA)),rep("h133p",(ncol(Rembrandt)+ncol(cell_line))))
my_data=ComBat(data, batch, mod=NULL, par.prior = TRUE,
               prior.plots = FALSE)
library(ade4)
my_pca=dudi.pca(t(my_data))
s.class(my_pca$li,as.factor(c(rep("tumor",(ncol(Rembrandt)+ncol(TCGA))),rep("cell_line",ncol(cell_line)))))


TCGA_corrected=my_data[,1:ncol(TCGA)]
Rembrandt_corrected=my_data[,(ncol(TCGA)+1):(ncol(TCGA)+ncol(Rembrandt))]
cell_line_corrected=my_data[,(ncol(TCGA)+ncol(Rembrandt)+1):790]



ematrix=cell_line_corrected


#Clustering of the 129 samples of glioblastoma cell lines 
#corresponding to the figure 1A of the article 


library(FactoMineR)
library(factoextra)
library(ade4)


res.pca=PCA(t(ematrix),scale.unit = TRUE,ncp=57)

res.hcpc=HCPC(res.pca,graph=FALSE,nb.clust=-1)

fviz_dend(res.hcpc, 
          cex = 0.7,                     # Taille du text
          palette = "jco",               # Palette de couleur ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Rectangle autour des groupes
          rect_border = "jco",           # Couleur du rectangle
          labels_track_height = 0.8      # Augment l'espace pour le texte
)


fviz_cluster(res.hcpc,
             repel = TRUE,            # Evite le chevauchement des textes
             show.clust.cent = TRUE, # Montre le centre des clusters
             palette = "jco",         # Palette de couleurs, voir ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)


tiff("plot.tiff",width=3000,height=1856,res=600)
plot(res.hcpc, choice = "3D.map",ind.names=FALSE,cex.axis=3.5,cex.symbols=3.5,xlab="",ylab="")
dev.off()
table(res.hcpc$data.clust$clust)





#--- Estimation of each cell line cluster in each patient 
#--- TCGA cohort 

#Loading the estimations

setwd("~/Recherche/Glioblastome/Glio/Glio/Donnees_exportees_analysees/TCGA")

res=read.table("cell_line_estimated_TCGA.txt",sep="\t")



#checking the association between the clinical data and transcriptome data

GBM=TCGA_corrected
info_clin=read.table("TCGA_clinical.txt",sep="\t",header=TRUE,row.names=1)

colnames(GBM)=gsub("\\.","-",colnames(GBM))
GBM=GBM[,which(colnames(GBM)%in%info_clin[,1])]
info_clin=info_clin[which(info_clin[,1]%in%colnames(GBM)),]

colnames(GBM)==info_clin[,1]



#clustering the patients based on the cell line estimations
#--- Figure 1B


library(cluster)

d=dist(scale(res[,1:3]),method="euclidean")


H.model=hclust(d,method="ward.D")

fviz_nbclust(scale(res[,1:3]),diss=NULL, hcut, method = "gap_stat",k.max=6,hc_method = "ward.D",nboot = 200)

sub_grp <- cutree(H.model, k = 3)

par(mfrow=c(1,1))
plot(H.model, cex = 0.6)
rect.hclust(H.model, k =3,border = 2:9)




fviz_cluster(list(data = res[,1:3], cluster = sub_grp))



table(sub_grp)
cell_line_classification=sub_grp
# 
# sub_grp
# 1   2   3 
# 198  81 202 


#Building the model to apply on validation cohort
#--- supervised model (random Forest based on discovery cohort)
library(randomForest)
model_lignee=randomForest(res[,1:3],as.factor(sub_grp))



#Loading Rembrandt datasets and related clinical data
#Figure 1C

setwd("~/Recherche/Glioblastome/Glio/Glio/Donnees_exportees_analysees/Rembrandt")
res_test=read.table("cell_line_estimated.txt",sep="\t",header=TRUE,row.names=1)
colnames(res_test)=colnames(res)

molecular_test=read.table("molecular_Rembrandt.txt",sep="\t")


info_clin_Rembrandt=read.table("info_clin_gbm.txt",sep="\t",header=TRUE,row.names=1)

#Predicting the group for each patient of the validation group
test=predict(model_lignee,res_test[,1:3])



library(survival)
library(survminer)
dss_t=info_clin$OS_MONTHS
dss_s=info_clin$OS_STATUS
dss_s=ifelse(dss_s=="DECEASED",1,0)

dss_t=as.numeric(as.character(dss_t))
sample=which(dss_t>24)

dss_t[sample]=24
dss_s[sample]=0




#Prognostic role of Classical classification on patients on TCGA cohort (OS)
#--- Figure 1D



top=data.frame(dss_t=dss_t,dss_s=dss_s,molecular=as.character(info_clin$EXPRESSION_SUBTYPE))

fit <- survfit(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")



#Prognostic role of cell line classification on patients on TCGA cohort (OS)
#--- Figure 1F

top=data.frame(dss_t=dss_t,dss_s=dss_s,molecular=sub_grp)

fit <- survfit(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")




dss_t=info_clin_Rembrandt$OVERALL_SURVIVAL_MONTHS
dss_s=info_clin_Rembrandt$EVENT_OS
dss_s=ifelse(dss_s=="EVENT",1,0)

dss_t=as.numeric(as.character(dss_t))
sample=which(dss_t>24)

dss_t[sample]=24
dss_s[sample]=0



#Prognostic role of Classical classification on patients on Rembrandt cohort (OS)
#--- Figure 1E

top=data.frame(dss_t=dss_t,dss_s=dss_s,molecular=molecular_test[,1])

fit <- survfit(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")


#Prognostic role of cell line classification on patients on Rembrandt cohort (OS)
#--- Figure 1G

top=data.frame(dss_t=dss_t,dss_s=dss_s,molecular=test)

fit <- survfit(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")



#DFS analysis only available for TCGA cohort



dss_t=info_clin$DFS_MONTHS
dss_s=info_clin$DFS_STATUS

dss_s=ifelse(dss_s=="DiseaseFree",0,1)

fit <- survfit(Surv(dss_t,dss_s)~1,data=info_clin)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

sample=which(dss_t>12)

dss_s[sample]=0
dss_t[sample]=12

#Prognostic role of classical classification on patients on TCGA cohort (DFS)
#--- Figure S2A
top=data.frame(dss_t=dss_t,dss_s=dss_s,molecular=info_clin$EXPRESSION_SUBTYPE)

fit <- survfit(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")


#Prognostic role of cell line classification on patients on TCGA cohort (DFS)
#--- Figure S2B

top=data.frame(dss_t=dss_t,dss_s=dss_s,molecular=sub_grp)

fit <- survfit(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")





#Testing classical and cell-lines classification on patients from both cohorts pulled together


classical_classif=c(as.character(info_clin$EXPRESSION_SUBTYPE),as.character(molecular_test[,1]))
OS_time=c(info_clin$OS_MONTHS,info_clin_Rembrandt$OVERALL_SURVIVAL_MONTHS)
OS_event=c(as.character(info_clin$OS_STATUS),as.character(info_clin_Rembrandt$EVENT_OS))
OS_event=ifelse((OS_event=="DECEASED"|OS_event=="EVENT"),1,0)

sample=which(OS_time>24)
OS_time[sample]=24
OS_event[sample]=0


top=data.frame(dss_t=OS_time,dss_s=OS_event,molecular=classical_classif)

fit <- survfit(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

classical_cox=coxph(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
summary(classical_cox)




new_classif=c(cell_line_classification,test)


top=data.frame(dss_t=OS_time,dss_s=OS_event,molecular=new_classif)

fit <- survfit(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

new_cox=coxph(Surv(dss_t,dss_s)~as.factor(molecular),data=top)
summary(new_cox)





classical_lp=predict(classical_cox,type="lp")
new_lp=predict(new_cox,type="lp")



temp=coxph(Surv(OS_time[-t],OS_event[-t])~classical_lp)
tem2=coxph(Surv(OS_time[-t],OS_event[-t])~new_lp)
AIC(temp)
AIC(tem2)













#--- Metagene analysis


setwd("~/tools/metagene")

#Loadings the list of genes related to each immune cells

Bcells=read.table("Bcells.txt",sep="\t",row.names=1,header=TRUE)
CD4=read.table("CD4.txt",sep="\t",row.names=1,header=TRUE)
CD8=read.table("CD8.txt",sep="\t",row.names=1,header=TRUE)
dendritic=read.table("dendritic.txt",sep="\t",row.names=1,header=TRUE)
macro=read.table("macro.txt",sep="\t",row.names=1,header=TRUE)
monocytes=read.table("monocytes.txt",sep="\t",row.names=1,header=TRUE)
NKcells=read.table("NKcells.txt",sep="\t",row.names=1,header=TRUE)
plasmacytoids=read.table("plasmacytoids.txt",sep="\t",row.names=1,header=TRUE)
Tfh=read.table("Tfh.txt",sep="\t",row.names=1,header=TRUE)
Tgd=read.table("Tgd.txt",sep="\t",row.names=1,header=TRUE)
Th1=read.table("Th1.txt",sep="\t",row.names=1,header=TRUE)
Th17=read.table("Th17.txt",sep="\t",row.names=1,header=TRUE)
Th2=read.table("Th2.txt",sep="\t",row.names=1,header=TRUE)
Treg=read.table("Treg.txt",sep="\t",row.names=1,header=TRUE)





#Computing the metagene value

Bcells_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(Bcells),])
CD4_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(CD4),])
CD8_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(CD8),])
dendritic_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(dendritic),])
macro_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(macro),])
mono_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(monocytes),])
NK_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(NKcells),])
plasmacytoids_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(plasmacytoids),])
Tfh_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(Tfh),])
Tgd_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(Tgd),])
Th1_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(Th1),])
Th17_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(Th17),])
Th2_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(Th2),])
Treg_TCGA=colMeans(2^TCGA_corrected[rownames(TCGA_corrected)%in%rownames(Treg),])



#Molecular subtypes for each TCGA patient

Mesenchymal=which(info_clin$EXPRESSION_SUBTYPE=="Mesenchymal")
Proneural=which(info_clin$EXPRESSION_SUBTYPE=="Proneural")
Neural=which(info_clin$EXPRESSION_SUBTYPE=="Neural")
Classical=which(info_clin$EXPRESSION_SUBTYPE=="Classical")



metagene=data.frame(Bcells=Bcells_TCGA,CD4=CD4_TCGA,CD8=CD8_TCGA,dendritic=dendritic_TCGA,macro=macro_TCGA,monocytes=mono_TCGA,NKcells=NK_TCGA,plasmacytoids=plasmacytoids_TCGA,Tfh=Tfh_TCGA,Tgd=Tgd_TCGA,Th1=Th1_TCGA,Th17=Th17_TCGA,Th2=Th2_TCGA,Treg=Treg_TCGA)


#Average value by molecular subtypes 

colMeans(metagene[Mesenchymal,])
colMeans(metagene[Proneural,])
colMeans(metagene[Neural,])
colMeans(metagene[Classical,])



#Average value for each subtype

colMeans(metagene[cell_line_classification==1,])
colMeans(metagene[cell_line_classification==2,])
colMeans(metagene[cell_line_classification==3,])



#¨Preparing Overall Survival data



dss_t=info_clin$OS_MONTHS
dss_s=info_clin$OS_STATUS

dss_s=ifelse(dss_s=="DECEASED",1,0)

fit <- survfit(Surv(dss_t,dss_s)~1,data=info_clin)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

sample=which(dss_t>24)

dss_s[sample]=0
dss_t[sample]=24




#OS on molecular classification 

#--- We delete patients with os time = 0

delete=which(dss_t==0)




normalisation=function(X)
{
  (X-mean(X))/sd(X)
}


continous_survival=function(OS.time,OS.events,biomarker)
{
  
  my_cox=coxph(Surv(OS.time,OS.events)~biomarker)
  hazard_ratio=exp(coef(my_cox))
  conf_int=exp(confint(my_cox))
  p.value=summary(my_cox)$coefficients[5]
  
  results=c(hazard_ratio,conf_int,p.value)
  names(results)=c("HR","lower .95","upper .95","p.value")
  
  return(results)
}
#We compute univariate cox model for each metagene and for each molecular or cell line subtype 

metagene_normalized=apply(metagene,2,normalisation)
res=apply(metagene_normalized[-delete,],2,OS.time=dss_t[-delete],OS.events=dss_s[-delete],continous_survival)
write.table(res,"results.txt",sep="\t")



res=apply(metagene_normalized[Mesenchymal,],2,OS.time=dss_t[Mesenchymal],OS.events=dss_s[Mesenchymal],continous_survival)

write.table(res,"results.txt",sep="\t")




res=apply(metagene_normalized[Neural,],2,OS.time=dss_t[Neural],OS.events=dss_s[Neural],continous_survival)

write.table(res,"results_Neural.txt",sep="\t")



res=apply(metagene_normalized[Proneural,],2,OS.time=dss_t[Proneural],OS.events=dss_s[Proneural],continous_survival)

write.table(res,"results_Proneural.txt",sep="\t")


res=apply(metagene_normalized[Classical,],2,OS.time=dss_t[Classical],OS.events=dss_s[Classical],continous_survival)

write.table(res,"results_Classical.txt",sep="\t")


res=apply(metagene_normalized[which(cell_line_classification==1),],2,OS.time=dss_t[which(cell_line_classification==1)],OS.events=dss_s[which(cell_line_classification==1)],continous_survival)

write.table(res,"results_1.txt",sep="\t")


res=apply(metagene_normalized[which(cell_line_classification==2),],2,OS.time=dss_t[which(cell_line_classification==2)],OS.events=dss_s[which(cell_line_classification==2)],continous_survival)

write.table(res,"results_2.txt",sep="\t")


res=apply(metagene_normalized[which(cell_line_classification==3),],2,OS.time=dss_t[which(cell_line_classification==3)],OS.events=dss_s[which(cell_line_classification==3)],continous_survival)

write.table(res,"results_3.txt",sep="\t")


#Clustering on metagene value

d=dist(scale(metagene),method="euclidean")


H.model=hclust(d,method="ward.D")

fviz_nbclust(scale(metagene), hcut, method = "gap_stat",
             hc_method = "ward.D2",k.max=6,nboot=200)


#4 clusters found
sub_grp <- cutree(H.model, k = 4)

par(mfrow=c(1,1))
plot(H.model, cex = 0.6)
rect.hclust(H.model, k=4,border = 2:9)


results=data.frame(metagene,molecular=info_clin$EXPRESSION_SUBTYPE,group=sub_grp,cell_line=as.character(cell_line_classification))
write.table(results,"results.txt",sep="\t")

model_metagene=randomForest(metagene,as.factor(sub_grp))



#Disease Free survival 



dss_t=info_clin$DFS_MONTHS
dss_s=info_clin$DFS_STATUS

dss_s=ifelse(dss_s=="DiseaseFree",0,1)

fit <- survfit(Surv(dss_t,dss_s)~1,data=info_clin)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

sample=which(dss_t>12)

dss_s[sample]=0
dss_t[sample]=12



#DFS on molecular classification 

#--- On enleve les patients pour lesquels la DFS est égale 0

delete=which(dss_t==0)


res=apply(metagene_normalized[-delete,],2,OS.time=dss_t[-delete],OS.events=dss_s[-delete],continous_survival)
write.table(res,"results.txt",sep="\t")



res=apply(metagene_normalized[Mesenchymal,],2,OS.time=dss_t[Mesenchymal],OS.events=dss_s[Mesenchymal],continous_survival)

write.table(res,"results.txt",sep="\t")




res=apply(metagene_normalized[Neural,],2,OS.time=dss_t[Neural],OS.events=dss_s[Neural],continous_survival)

write.table(res,"results_Neural.txt",sep="\t")



res=apply(metagene_normalized[Proneural,],2,OS.time=dss_t[Proneural],OS.events=dss_s[Proneural],continous_survival)

write.table(res,"results_Proneural.txt",sep="\t")


res=apply(metagene_normalized[Classical,],2,OS.time=dss_t[Classical],OS.events=dss_s[Classical],continous_survival)

write.table(res,"results_Classical.txt",sep="\t")


res=apply(metagene_normalized[which(cell_line_classification==1),],2,OS.time=dss_t[which(cell_line_classification==1)],OS.events=dss_s[which(cell_line_classification==1)],continous_survival)

write.table(res,"results_1.txt",sep="\t")


res=apply(metagene_normalized[which(cell_line_classification==2),],2,OS.time=dss_t[which(cell_line_classification==2)],OS.events=dss_s[which(cell_line_classification==2)],continous_survival)

write.table(res,"results_2.txt",sep="\t")


res=apply(metagene_normalized[which(cell_line_classification==3),],2,OS.time=dss_t[which(cell_line_classification==3)],OS.events=dss_s[which(cell_line_classification==3)],continous_survival)

write.table(res,"results_3.txt",sep="\t")




#--- Metagene analysis on Rembrandt cohort


Bcells_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(Bcells),])
CD4_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(CD4),])
CD8_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(CD8),])
dendritic_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(dendritic),])
macro_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(macro),])
mono_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(monocytes),])
NK_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(NKcells),])
plasmacytoids_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(plasmacytoids),])
Tfh_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(Tfh),])
Tgd_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(Tgd),])
Th1_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(Th1),])
Th17_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(Th17),])
Th2_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(Th2),])
Treg_Rembrandt=colMeans(2^Rembrandt_corrected[rownames(Rembrandt_corrected)%in%rownames(Treg),])





Mesenchymal=which(molecular_test=="Mesenchymal")
Proneural=which(molecular_test=="Proneural")
Neural=which(molecular_test=="Neural")
Classical=which(molecular_test=="Classical")


metagene=data.frame(Bcells=Bcells_Rembrandt,CD4=CD4_Rembrandt,CD8=CD8_Rembrandt,dendritic=dendritic_Rembrandt,macro=macro_Rembrandt,monocytes=mono_Rembrandt,NKcells=NK_Rembrandt,plasmacytoids=plasmacytoids_Rembrandt,Tfh=Tfh_Rembrandt,Tgd=Tgd_Rembrandt,Th1=Th1_Rembrandt,Th17=Th17_Rembrandt,Th2=Th2_Rembrandt,Treg=Treg_Rembrandt)


colMeans(metagene[Mesenchymal,])
colMeans(metagene[Proneural,])
colMeans(metagene[Neural,])
colMeans(metagene[Classical,])


cell_line_Rembrandt=test



colMeans(metagene[cell_line_Rembrandt==1,])
colMeans(metagene[cell_line_Rembrandt==2,])
colMeans(metagene[cell_line_Rembrandt==3,])




dss_t=info_clin_Rembrandt$OVERALL_SURVIVAL_MONTHS
dss_s=info_clin_Rembrandt$EVENT_OS
dss_s=ifelse(dss_s=="EVENT",1,0)

dss_t=as.numeric(as.character(dss_t))
sample=which(dss_t>24)

dss_t[sample]=24
dss_s[sample]=0







metagene_normalized=apply(metagene,2,normalisation)
res=apply(metagene_normalized,2,OS.time=dss_t,OS.events=dss_s,continous_survival)
write.table(res,"results.txt",sep="\t")



res=apply(metagene_normalized[Mesenchymal,],2,OS.time=dss_t[Mesenchymal],OS.events=dss_s[Mesenchymal],continous_survival)

write.table(res,"results.txt",sep="\t")




res=apply(metagene_normalized[Neural,],2,OS.time=dss_t[Neural],OS.events=dss_s[Neural],continous_survival)

write.table(res,"results.txt",sep="\t")



res=apply(metagene_normalized[Proneural,],2,OS.time=dss_t[Proneural],OS.events=dss_s[Proneural],continous_survival)

write.table(res,"results.txt",sep="\t")


res=apply(metagene_normalized[Classical,],2,OS.time=dss_t[Classical],OS.events=dss_s[Classical],continous_survival)

write.table(res,"results.txt",sep="\t")






res=apply(metagene_normalized[which(cell_line_Rembrandt==1),],2,OS.time=dss_t[which(cell_line_Rembrandt==1)],OS.events=dss_s[which(cell_line_Rembrandt==1)],continous_survival)

write.table(res,"results.txt",sep="\t")


res=apply(metagene_normalized[which(cell_line_Rembrandt==2),],2,OS.time=dss_t[which(cell_line_Rembrandt==2)],OS.events=dss_s[which(cell_line_Rembrandt==2)],continous_survival)

write.table(res,"results.txt",sep="\t")


res=apply(metagene_normalized[which(cell_line_Rembrandt==3),],2,OS.time=dss_t[which(cell_line_Rembrandt==3)],OS.events=dss_s[which(cell_line_Rembrandt==3)],continous_survival)

write.table(res,"results.txt",sep="\t")


class_immuno=predict(model_metagene,metagene)





setwd("~/Recherche/Glioblastome/Glio/Glio/Donnees_exportees_analysees/TCGA")
res_abs=read.table("TCGA_res_abs.txt",sep="\t",header=TRUE,row.names = 1)


Mesenchymal=which(info_clin$EXPRESSION_SUBTYPE=="Mesenchymal")
Classical=which(info_clin$EXPRESSION_SUBTYPE=="Classical")
Proneural=which(info_clin$EXPRESSION_SUBTYPE=="Proneural")
Neural=which(info_clin$EXPRESSION_SUBTYPE=="Neural")




colMeans(res_abs[Mesenchymal,])
colMeans(res_abs[Classical,])
colMeans(res_abs[Proneural,])
colMeans(res_abs[Neural,])



colMeans(res_abs[cell_line_classification==1,])
colMeans(res_abs[cell_line_classification==2,])
colMeans(res_abs[cell_line_classification==3,])



d=dist(scale(res_abs),method="euclidean")


H.model=hclust(d,method="ward.D")

fviz_nbclust(scale(res_abs), hcut, method = "gap_stat",
             hc_method = "ward.D",k.max=6,nboot=200)

sub_grp <- cutree(H.model, k = 6)

par(mfrow=c(1,1))
plot(H.model, cex = 0.6)
rect.hclust(H.model, k=6,border = 2:9)




fviz_cluster(list(data = scale(res_abs), cluster = sub_grp))



table(sub_grp)

#we get rid of the patient that is unique in group 6 it seems to be an outlier for this clustering step

del=which(sub_grp==6)
sub_grp=sub_grp[-del]
res_abs=res_abs[-del,]

#we used this table to produce figure 3 (A-C-E) of the article

results=data.frame(res_abs,molecular=info_clin$EXPRESSION_SUBTYPE[-del],group=sub_grp,cell_line=cell_line_classification[-del])







#Survival analysis on Overall Survival 


dss_t=info_clin$OS_MONTHS
dss_s=info_clin$OS_STATUS

dss_s=ifelse(dss_s=="DECEASED",1,0)

fit <- survfit(Surv(dss_t,dss_s)~1,data=info_clin)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

sample=which(dss_t>24)

dss_s[sample]=0
dss_t[sample]=24

delete=which(dss_t==0)

res_abs_normalized=apply(res_abs,2,normalisation)
res=apply(res_abs_normalized[-delete,],2,continous_survival,OS.time=dss_t[-delete],OS.events=dss_s[-delete])
write.table(res,"results.txt",sep="\t")



res=apply(res_abs_normalized[Mesenchymal,],2,continous_survival,OS.time=dss_t[Mesenchymal],OS.events=dss_s[Mesenchymal])
write.table(res,"results.txt",sep="\t")




res=apply(res_abs_normalized[Classical,],2,continous_survival,OS.time=dss_t[Classical],OS.events=dss_s[Classical])
write.table(res,"results.txt",sep="\t")


res=apply(res_abs_normalized[Neural,],2,continous_survival,OS.time=dss_t[Neural],OS.events=dss_s[Neural])
write.table(res,"results.txt",sep="\t")


res=apply(res_abs_normalized[Proneural,],2,continous_survival,OS.time=dss_t[Proneural],OS.events=dss_s[Proneural])
write.table(res,"results.txt",sep="\t")



FAC=which(cell_line_classification==1)
Glioma=which(cell_line_classification==2)
Metabolic=which(cell_line_classification==3)


res=apply(res_abs_normalized[FAC,],2,continous_survival,OS.time=dss_t[FAC],OS.events=dss_s[FAC])
write.table(res,"results.txt",sep="\t")




res=apply(res_abs_normalized[Glioma,],2,continous_survival,OS.time=dss_t[Glioma],OS.events=dss_s[Glioma])
write.table(res,"results.txt",sep="\t")



res=apply(res_abs_normalized[Metabolic,],2,continous_survival,OS.time=dss_t[Metabolic],OS.events=dss_s[Metabolic])
write.table(res,"results.txt",sep="\t")


#Survival analysis on Disease Free Survival 


dss_t=info_clin$DFS_MONTHS
dss_s=info_clin$DFS_STATUS

dss_s=ifelse(dss_s=="DiseaseFree",0,1)

fit <- survfit(Surv(dss_t,dss_s)~1,data=info_clin)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

sample=which(dss_t>12)

dss_s[sample]=0
dss_t[sample]=12



delete=which(dss_t==0)


res=apply(res_abs_normalized[-delete,],2,continous_survival,OS.time=dss_t[-delete],OS.events=dss_s[-delete])
write.table(res,"results.txt",sep="\t")



res=apply(res_abs_normalized[Mesenchymal,],2,continous_survival,OS.time=dss_t[Mesenchymal],OS.events=dss_s[Mesenchymal])
write.table(res,"results.txt",sep="\t")




res=apply(res_abs_normalized[Classical,],2,continous_survival,OS.time=dss_t[Classical],OS.events=dss_s[Classical])
write.table(res,"results.txt",sep="\t")


res=apply(res_abs_normalized[Neural,],2,continous_survival,OS.time=dss_t[Neural],OS.events=dss_s[Neural])
write.table(res,"results.txt",sep="\t")


res=apply(res_abs_normalized[Proneural,],2,continous_survival,OS.time=dss_t[Proneural],OS.events=dss_s[Proneural])
write.table(res,"results.txt",sep="\t")



FAC=which(cell_line_classification==1)
Glioma=which(cell_line_classification==2)
Metabolic=which(cell_line_classification==3)


res=apply(res_abs_normalized[FAC,],2,continous_survival,OS.time=dss_t[FAC],OS.events=dss_s[FAC])
write.table(res,"results.txt",sep="\t")




res=apply(res_abs_normalized[Glioma,],2,continous_survival,OS.time=dss_t[Glioma],OS.events=dss_s[Glioma])
write.table(res,"results.txt",sep="\t")



res=apply(res_abs_normalized[Metabolic,],2,continous_survival,OS.time=dss_t[Metabolic],OS.events=dss_s[Metabolic])
write.table(res,"results.txt",sep="\t")



#Same analysis on Rembrandt cohort



res_lvl0=immune_cell_express(sig_lvl0,2^Rembrandt_corrected,QN=FALSE)
res_lvl1=immune_cell_express(sig_lvl1,2^Rembrandt_corrected,QN=FALSE)
res_lvl2=immune_cell_express(sig_lvl2,2^Rembrandt_corrected,QN=FALSE)
res_lvl3=immune_cell_express(sig_lvl3,2^Rembrandt_corrected,QN=FALSE)
res_lvl4=immune_cell_express(sig_lvl4,2^Rembrandt_corrected,QN=FALSE)
res_myelo=immune_cell_express(sig_myelo,2^Rembrandt_corrected,QN=FALSE)
res_CD4=immune_cell_express(sig_CD4,2^Rembrandt_corrected,QN=FALSE)
res_macro=immune_cell_express(sig_macro,2^Rembrandt_corrected,QN=FALSE)


#calcul des quantités absolues 

res1_abs=sweep(res_lvl1[,1:4],1,res_lvl0[,2],"*")

res2_abs=sweep(res_lvl2[,1:4],1,res_lvl1[,2],'*')
res2_abs=sweep(res2_abs,1,res_lvl0[,2],'*')

res3_abs=sweep(res_lvl3[,1:3],1,res_lvl2[,4],'*')
res3_abs=sweep(res3_abs,1,res_lvl1[,2],'*')
res3_abs=sweep(res3_abs,1,res_lvl0[,2],'*')


res4_abs=sweep(res_lvl4[,1:4],1,res_lvl3[,2],'*')
res4_abs=sweep(res4_abs,1,res_lvl2[,4],'*')
res4_abs=sweep(res4_abs,1,res_lvl1[,2],'*')
res4_abs=sweep(res4_abs,1,res_lvl0[,2],'*')

res_myelo_abs=sweep(res_myelo[,1:5],1,res_lvl1[,3],'*')
res_myelo_abs=sweep(res_myelo_abs,1,res_lvl0[,2],'*')



res_macro_abs=sweep(res_macro[,1:2],1,res_myelo[,3],'*')
res_macro_abs=sweep(res_macro_abs,1,res_lvl1[,3],'*')
res_macro_abs=sweep(res_macro_abs,1,res_lvl0[,2],'*')


res_CD4_abs=sweep(res_CD4[,1:5],1,res_lvl3[,1],'*')
res_CD4_abs=sweep(res_CD4_abs,1,res_lvl2[,4],'*')
res_CD4_abs=sweep(res_CD4_abs,1,res_lvl1[,2],'*')
res_CD4_abs=sweep(res_CD4_abs,1,res_lvl0[,2],'*')


res=data.frame(cancer=res_lvl1[,1],lympho=res_lvl1[,2],myelo=res_lvl1[,3],stroma=res_lvl1[,4],Bcells=res_lvl2[,1],NK=res_lvl2[,2],PlasmaCells=res_lvl2[,3],Tcells=res_lvl2[,4],CD8=res_lvl3[,1],CD4=res_lvl3[,2],Tgd=res_lvl3[,3],CD8_naive=res_lvl4[,4],CM_CD8=res_lvl4[,1],EM_CD8=res_lvl4[,2],EMRA_CD8=res_lvl4[,3],Th1=res_CD4[,1],Th2=res_CD4[,2],Th17=res_CD4[,3],Treg=res_CD4[,4],Tfh=res_CD4[,5],dendritic=res_myelo[,1],granulocytes=res_myelo[,2],macro=res_myelo[,3],monocytes=res_myelo[,4],plasmacytoid=res_myelo[,5],M1=res_macro[,1],M2=res_macro[,2])
res_abs_Rembrandt=data.frame(cancer_abs=res1_abs[,1],lympho_abs=res1_abs[,2],myelo_abs=res1_abs[,3],stroma_abs=res1_abs[,4],Bcells_abs=res2_abs[,1],NK_abs=res2_abs[,2],PlasmaCells_abs=res2_abs[,3],Tcells_abs=res2_abs[,4],CD8_abs=res3_abs[,1],CD4_abs=res3_abs[,2],Tgd_abs=res3_abs[,3],CD8_naive_abs=res4_abs[,4],CM_CD8_abs=res4_abs[,1],EM_CD8_abs=res4_abs[,2],EMRA_CD8_abs=res4_abs[,3],Th1_abs=res_CD4_abs[,1],Th2_abs=res_CD4_abs[,2],Th17_abs=res_CD4_abs[,3],Treg_abs=res_CD4_abs[,4],Tfh_abs=res_CD4_abs[,5],dendritic_abs=res_myelo_abs[,1],granulocytes_abs=res_myelo_abs[,2],macro_abs=res_myelo_abs[,3],monocytes_abs=res_myelo_abs[,4],plasmacytoid_abs=res_myelo_abs[,5],M1_abs=res_macro_abs[,1],M2_abs=res_macro_abs[,2])





Mesenchymal=which(molecular_test[,1]=="Mesenchymal")
Classical=which(molecular_test[,1]=="Classical")
Proneural=which(molecular_test[,1]=="Proneural")
Neural=which(molecular_test[,1]=="Neural")




colMeans(res_abs[Mesenchymal,])
colMeans(res_abs[Classical,])
colMeans(res_abs[Proneural,])
colMeans(res_abs[Neural,])



colMeans(res_abs[cell_line_Rembrandt==1,])
colMeans(res_abs[cell_line_Rembrandt==2,])
colMeans(res_abs[cell_line_Rembrandt==3,])






#Survival analysis on Overall Survival 


dss_t=info_clin_Rembrandt$OVERALL_SURVIVAL_MONTHS
dss_s=info_clin_Rembrandt$EVENT_OS

dss_s=ifelse(dss_s=="EVENT",1,0)

fit <- survfit(Surv(dss_t,dss_s)~1,data=info_clin)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

sample=which(dss_t>24)

dss_s[sample]=0
dss_t[sample]=24



res_abs_normalized=apply(res_abs,2,normalisation)
res=apply(res_abs_normalized,2,continous_survival,OS.time=dss_t,OS.events=dss_s)
write.table(res,"results.txt",sep="\t")



res=apply(res_abs_normalized[Mesenchymal,],2,continous_survival,OS.time=dss_t[Mesenchymal],OS.events=dss_s[Mesenchymal])
write.table(res,"results.txt",sep="\t")




res=apply(res_abs_normalized[Classical,],2,continous_survival,OS.time=dss_t[Classical],OS.events=dss_s[Classical])
write.table(res,"results.txt",sep="\t")


res=apply(res_abs_normalized[Neural,],2,continous_survival,OS.time=dss_t[Neural],OS.events=dss_s[Neural])
write.table(res,"results.txt",sep="\t")


res=apply(res_abs_normalized[Proneural,],2,continous_survival,OS.time=dss_t[Proneural],OS.events=dss_s[Proneural])
write.table(res,"results.txt",sep="\t")



FAC=which(cell_line_Rembrandt==1)
Glioma=which(cell_line_Rembrandt==2)
Metabolic=which(cell_line_Rembrandt==3)


res=apply(res_abs_normalized[FAC,],2,continous_survival,OS.time=dss_t[FAC],OS.events=dss_s[FAC])
write.table(res,"results.txt",sep="\t")




res=apply(res_abs_normalized[Glioma,],2,continous_survival,OS.time=dss_t[Glioma],OS.events=dss_s[Glioma])
write.table(res,"results.txt",sep="\t")



res=apply(res_abs_normalized[Metabolic,],2,continous_survival,OS.time=dss_t[Metabolic],OS.events=dss_s[Metabolic])
write.table(res,"results.txt",sep="\t")




#--- Multivariate model with LASSO shrinkage 

setwd("~/Recherche/Glioblastome/Glio/Glio/Donnees_exportees_analysees/TCGA")

res_abs=read.table("TCGA_res_abs.txt",sep="\t",header=TRUE,row.names=1)


dss_t=info_clin$OS_MONTHS
dss_s=info_clin$OS_STATUS

dss_s=ifelse(dss_s=="DECEASED",1,0)

fit <- survfit(Surv(dss_t,dss_s)~1,data=info_clin)
ggsurvplot(fit, risk.table = TRUE, surv.median.line = "hv")

sample=which(dss_t>24)

dss_s[sample]=0
dss_t[sample]=24





delete=which(dss_t==0)

res_abs_normalized=apply(res_abs,2,normalisation)
molecular=as.factor(as.character(info_clin$EXPRESSION_SUBTYPE))
cell_line=as.factor(cell_line_classification)
age=as.numeric(info_clin$AGE)
sexe=as.factor(info_clin$GENDER)


#we dont take in the matrix the variable that can be obtained as linear combination of others/ they don't bring any information
top=data.frame(dss_t=dss_t,dss_s=dss_s,molecular=molecular,cell_line=cell_line,age=age,sexe=sexe,res_abs_normalized[,-c(1,4,2,3,8,9,10,15,23)])
sample= unique(which(is.na(top),arr.ind=TRUE)[,1])
top=top[-sample,]


form=Surv(dss_t,dss_s)~.+cell_line:Bcells_abs+cell_line:NK_abs+cell_line:PlasmaCells_abs+cell_line:Tgd_abs+cell_line:CD8_naive_abs+cell_line:CM_CD8_abs+cell_line:EM_CD8_abs+cell_line:Th1_abs+cell_line:Th2_abs+cell_line:Th17_abs+cell_line:Treg_abs+cell_line:Tfh_abs+cell_line:dendritic_abs+cell_line:granulocytes_abs+cell_line:monocytes_abs+cell_line:plasmacytoid_abs+cell_line:M1_abs+cell_line:M2_abs
Mat=model.matrix(form,data=top[,-c(3,5,6)])

#LASSO model 
library(glmnet)
penal_factor=rep(1,ncol(Mat))
penal_factor[c(2,3)]=c(0,0)
coxnet=cv.glmnet(Mat,Surv(top$dss_t,top$dss_s),family="cox",maxit=10000,nfolds=479,penalty.factor=penal_factor,standardize=FALSE)
plot(coxnet)
coxnet$lambda


Coefficients <- coef(coxnet, s =coxnet$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients  <- Coefficients[Active.Index]



val_pred=predict(coxnet,newx=Mat,s= coxnet$lambda.min,type="link")




fit <- survfit(Surv(dss_t, dss_s)~as.factor(ifelse(val_pred>median(val_pred),1,0)), data = top)
ggsurvplot(fit, risk.table = TRUE, break.time.by = 1, surv.median.line = "hv")

summary(coxph(Surv(dss_t, dss_s)~as.factor(ifelse(val_pred>median(val_pred),1,0)),data=top))



write.table(cbind(top$dss_t,top$dss_s,as.factor(ifelse(val_pred>median(val_pred),1,0))),"results.txt",sep="\t")



###################################################################################################################
##---------------------------------Clinical variables on OS


#age 

my_cox=coxph(Surv(dss_t,dss_s)~age,data=top)
summary(my_cox)


#sex 

my_cox=coxph(Surv(dss_t,dss_s)~as.factor(sexe),data=top)
summary(my_cox)


#IDH

IDH=as.character(info_clin$IDH1_STATUS)
IDH=IDH[-sample]

temp=which(IDH=="R132G"|IDH=="R132H")


IDH[temp]="Mutated"

table(IDH)


my_cox=coxph(Surv(dss_t,dss_s)~as.factor(IDH),data=top)
summary(my_cox)


#MGMT

MGMT=info_clin$MGMT_STATUS
MGMT=MGMT[-sample]

my_cox=coxph(Surv(dss_t,dss_s)~as.factor(MGMT),data=top)
summary(my_cox)


#Treatment

treatment=info_clin$Treatment
treatment=treatment[-sample]

my_cox=coxph(Surv(dss_t,dss_s)~as.factor(treatment),data=top)
summary(my_cox)



#composite biomarker


my_cox=coxph(Surv(dss_t,dss_s)~val_pred,data=top)
summary(my_cox)







###################################################################################################################
##---------------------------------Multivariate model 

library(MASS)


my_cox=coxph(Surv(dss_t,dss_s)~age+sexe+as.factor(MGMT)+as.factor(treatment),data=top)
summary(my_cox)


AIC(my_cox)




my_cox2=coxph(Surv(dss_t,dss_s)~age+sexe+as.factor(MGMT)+as.factor(treatment)+val_pred,data=top)
summary(my_cox2)

AIC(my_cox)



anova(my_cox,my_cox2)



# 
# Analysis of Deviance Table
# Cox model: response is  Surv(dss_t, dss_s)
# Model 1: ~ age + sexe + as.factor(MGMT) + as.factor(treatment)
# Model 2: ~ age + sexe + as.factor(MGMT) + as.factor(treatment) + val_pred
# loglik  Chisq Df P(>|Chi|)    
# 1 -897.33                        
# 2 -890.88 12.908  1 0.0003271 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1













model_final=coxph(Surv(dss_t,dss_s)~val_pred+age,data=top)
summary(model_final)


val_pred2=predict(model_final,type="lp")


grouped=rep(0,length(val_pred2))
quantiles=quantile(val_pred2,probs=c(0,0.25,0.5,0.75,1))

grouped[val_pred2<quantiles[2]]="Low"
grouped[val_pred2>=quantiles[2]&val_pred2<quantiles[4]]="Medium"
grouped[val_pred2>=quantiles[4]]="High"



table(grouped)


fit <- survfit(Surv(dss_t, dss_s)~as.factor(grouped), data = top)
ggsurvplot(fit, risk.table = TRUE, break.time.by = 1, surv.median.line = "hv")


summary(coxph(Surv(dss_t[-sample],dss_s[-sample])~as.factor(grouped)))


write.table(cbind(top$dss_t,top$dss_s,grouped),"results.txt",sep="\t")










#-- Validation coxboost


setwd("~/Recherche/Glioblastome/Glio/Glio/Donnees_exportees_analysees/Rembrandt")

info_clin_Rembrandt=read.table("info_clin_gbm.txt",sep="\t",header=TRUE,row.names=1)
res_abs_Rembrandt=read.table("res_abs_rembrandt.txt",sep="\t",header=TRUE,row.names=1)
cell_line_Rembrandt=read.table("cell_line_classification.txt",sep="\t",header=TRUE,row.names=1)
molecular_Rembrandt=read.table("molecular_Rembrandt.txt",sep="\t",header=TRUE,row.names=1)

temp_Rembrandt=apply(res_abs_Rembrandt,2,normalisation)

dss_t=info_clin_Rembrandt$OVERALL_SURVIVAL_MONTHS
dss_s=info_clin_Rembrandt$EVENT_OS
dss_s=ifelse(dss_s=="EVENT",1,0)
sample=which(dss_t>24)

dss_t[sample]=24
dss_s[sample]=0


age=info_clin_Rembrandt$AGE_RANGE
age=as.character(age)
age=substr(age,1,2)
age=as.numeric(age)
sexe=info_clin_Rembrandt$GENDER
sexe=as.character(sexe)
sexe[sexe==""]=NA



top=data.frame(dss_t=dss_t,dss_s=dss_s,molecular=as.factor(molecular_Rembrandt[,1]),cell_line=as.factor(cell_line_Rembrandt[,1]),age=age,sexe=sexe,temp_Rembrandt[,-c(1,4,2,3,8,9,10,15,23)])
sample=unique(which(is.na(top),arr.ind=TRUE)[,1])
top=top[-sample,]


form=Surv(dss_t,dss_s)~.+cell_line:Bcells_abs+cell_line:NK_abs+cell_line:PlasmaCells_abs+cell_line:Tgd_abs+cell_line:CD8_naive_abs+cell_line:CM_CD8_abs+cell_line:EM_CD8_abs+cell_line:Th1_abs+cell_line:Th2_abs+cell_line:Th17_abs+cell_line:Treg_abs+cell_line:Tfh_abs+cell_line:dendritic_abs+cell_line:granulocytes_abs+cell_line:monocytes_abs+cell_line:plasmacytoid_abs+cell_line:M1_abs+cell_line:M2_abs
# form=Surv(dss_t,dss_s)~.
Mat=model.matrix(form,data=top[,-c(3,5,6)])








coxnet_predicted=predict(coxnet,newx=Mat,type="link",s=coxnet$lambda.min)


fit <- survfit(Surv(dss_t, dss_s)~as.factor(ifelse(coxnet_predicted>median(val_pred),1,0)), data = top)
ggsurvplot(fit, risk.table = TRUE, break.time.by = 1, surv.median.line = "hv")

summary(coxph(Surv(dss_t, dss_s)~as.factor(ifelse(coxnet_predicted>median(val_pred),1,0)),data=top))



write.table(cbind(top$dss_t,top$dss_s,as.factor(ifelse(coxnet_predicted>median(val_pred),1,0))),"results.txt",sep="\t")


t=data.frame(val_pred=coxnet_predicted,age=top$age)
colnames(t)=c("val_pred","age")
coxboost_predicted2=predict(model_final,newdata=t,type="lp")


grouped=rep(0,length(coxnet_predicted))

grouped[coxboost_predicted2<quantiles[2]]="Low"
grouped[coxboost_predicted2>=quantiles[2]&coxboost_predicted2<quantiles[4]]="Medium"
grouped[coxboost_predicted2>=quantiles[4]]="High"

write.table(cbind(top$dss_t,top$dss_s,grouped),"results.txt",sep="\t")
