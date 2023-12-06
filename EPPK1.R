setwd("C:\\Users\\qliul\\Desktop\\TCPA\\总")
library(dplyr)
library(limma)
a=read.csv('TCPA_TCGA.csv',header = T)
table(!duplicated(a$id))
rownames(a)=a$id
table(a$PR)
fivenum(a$PR)
d=ifelse(a$PR<=-0.053778,'normal','tumor')
table(d)
dat=a
dat=dat[,c(-1)]
dat=t(dat)
dat=as.data.frame(dat)
list=d
list=factor(d,levels = c('normal','tumor'),ordered = F)
head(list)
list <- model.matrix(~0+factor(list))  
colnames(list) <- c("normal", "tumor")
df.fit <-lmFit(dat,list)  
df.matrix <- makeContrasts(tumor - normal , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
head(tempOutput)
nrDEG = na.omit(tempOutput) 
diffsig <- nrDEG  
write.csv(diffsig, "all.limmaOut-4.csv")
b <- read.table("all.limmaOut_4.txt",sep = "\t",header = T,row.names = 1)
a <- read.table("retu.txt",sep = "\t",header = T,row.names = 1)
library(dplyr)
deg.genes <- b %>% 
  filter(P.Value < 0.001)
test=a[rownames(deg.genes),]
group <- data.frame(a=c(rep('PR_',75),rep('PR+',225)))
rownames(group) <- colnames(test)
library(pheatmap) 
outFile <- "my_heatmap1.pdf"      
pdf(file=outFile,width=6,height=5.5)
pheatmap(test,
         annotation=group,
         cluster_cols = F,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames = T,
         show_colnames = F,
         scale="row",  
         #border_color ="NA",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=6)
dev.off()
td=read.table("testdata.txt",header=T) 
library(survival) 
pFilter=0.05 
outResult=data.frame() 
sigGenes=c("surstat","surtime") 
for(i in colnames(td[,3:ncol(td)])){ 
  tdcox <- coxph(Surv(surtime, surstat) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox) 
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] 
  if(pvalue<pFilter){ 
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outResult,file = "univariateCOX1.xls",sep = "/t",row.names = F,quote = F)
write.table(outResult,file="UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)
UniCoxSurSigGeneExp=td[,sigGenes] 
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)
write.table(UniCoxSurSigGeneExp,file="UniCoxSurSigGeneExp.txt",sep="\t",row.names=F,quote=F)
tducs <- read.table("UniCoxSurvival.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(tducs)#提取基因名
hr <- sprintf("%.3f",tducs$"HR")
hrLow  <- sprintf("%.3f",tducs$"L95CI")
hrHigh <- sprintf("%.3f",tducs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))
n <- nrow(tducs)
nRow <- n+1 
ylim <- c(1,nRow) 
layout(matrix(c(1,2),nc=2),width=c(2,2)) 
xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8 
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.4,n:1,pValue,adj=1,cex=text.cex);text(1.4,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh))) 
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="black",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
rm(list = ls())
df <- read.csv("huoshantu.csv",row.names = 1) 
df$threshold = factor(ifelse(df$P_Value  < 0.05 & abs(df$fd) >= 0.5849625, ifelse(df$fd >= 0.5849625 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df) 
ggplot(df,aes(x=fd,y= -log10(P_Value),fill = threshold))+
  geom_point(colour = "black", shape = 21, stroke = 0.5)+
  scale_fill_manual(values=c("#d90424","#0066CC","#bdbdbd"))+
  geom_text_repel(
    data = df[df$P_Value<0.05&abs(df$fd)>0.5849625,],
    aes(label = gene),
    size = 4.5,
    color = "black",
    segment.color = "black", show.legend = FALSE )+
  ylab('-log10 (Pvalue)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-0.5849625,0.5849625),lty=2,col="black",lwd=0.5) +#添加横线|logFoldChange|>0.25
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +
  theme_classic(  
    base_line_size = 1 
  )+
  guides(fill=guide_legend(override.aes = list(size=5)))+ 
  theme(axis.title.x = element_text(size = 10, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 10,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.text = element_text(color="black", 
                                   size = 7, 
                                   face = "bold")
  )


ggsave("volcano_2.0.pdf", height = 4, width = 5)
setwd("C:\\Users\\qliul\\Desktop\\TCPA\\总\\KM")
a=read.table("Age≤60.txt",sep="\t",check.names=F,header=T) #读取文件
library(survival)     
library(ggplot2)
library(ggpubr)
library(survminer)
median(a$EPPK1)
group <- ifelse(a$EPPK1>0.282415,"high","low")
a <- cbind(a,group)
fit <- survfit(Surv(futime, fustat)~group, data=a)  
ggsurvplot(fit,pval=TRUE)
ggsurvplot(
  fit,
  data = a,
  size = 1,                
  conf.int = F,          
  pval = TRUE,              
  risk.table = F,       
  risk.table.col = "strata",
  surv.median.line = "hv",
  xlab="Time(days)",
  legend.labs =
    c("High experssion","Low expression"),    
  legend.title = "Age≤60",
  risk.table.height = 0.25, 
  risk.table.size = 2,
  ggtheme = theme_bw()      
)
b=read.table("stageIV.txt",sep="\t",check.names=F,header=T)
group1 <- ifelse(b$EPPK1>0.282415,"high","low")
fit <- survfit(Surv(futime, fustat)~group1, data=b)  
ggsurvplot(fit,pval=TRUE)
ggsurvplot(
  fit,
  data = b,
  size = 1,                
  conf.int = F,          
  pval = TRUE,           
  risk.table = F,        
  risk.table.col = "strata",
  surv.median.line = "hv",
  xlab="Time(days)",
  legend.labs =
    c("High experssion","Low expression"),    
  legend.title = "StageIV",
  risk.table.height = 0.25, 
  risk.table.size = 2,
  ggtheme = theme_bw()      
)
setwd("C:\\Users\\qliul\\Desktop\\TCPA\\总")
aa=read.table("COX.txt",header=T,sep="\t",row.names=1) 
library(survival)
library(plyr)
str(aa)
aa$surstat<-factor(aa$surstat)
summary(aa$surstat)
y <- Surv(time=aa$surtime,event=aa$surstat==1)
Uni_cox_model<- function(x){
  FML <- as.formula(paste0 ("y~",x))
  cox<- coxph(FML,data=aa)
  cox1<-summary(cox)
  HR <- round(cox1$coefficients[,2],2)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <-round(cox1$conf.int[,3],2)
  CI95 <-round(cox1$conf.int[,4],2)
  Uni_cox_model<- data.frame('Characteristics' = x,
                             'HR' = HR,
                             'CI5' = CI5,
                             'CI95' = CI95,
                             'p' = PValue)
  return(Uni_cox_model)}  
names(aa)
variable.names<- colnames(aa)[c(3:7)] #例：这里选择了3-10号变量
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox <- ldply(Uni_cox,data.frame)
Uni_cox$CI<-paste(Uni_cox$CI5,'-',Uni_cox$CI95)
Uni_cox<-Uni_cox[,-3:-4]
write.table(Uni_cox,file="Uni_cox.txt",sep="\t",row.names=F,quote=F)
View(Uni_cox)
Uni_cox$Characteristics[Uni_cox$p<0.05]
mul_cox_model<- as.formula(paste0 ("y~",
                                   paste0(Uni_cox$Characteristics[Uni_cox$p<0.05],
                                          collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=aa)
cox4<-summary(mul_cox) 
summary(mul_cox)
mul_HR<- round(cox4$coefficients[,2],2) 
mul_PValue<- round(cox4$coefficients[,5],4) 
mul_CI1<-round(cox4$conf.int[,3],2)
mul_CI2<-round(cox4$conf.int[,4],2)
mul_CI<-paste(mul_CI1,'-',mul_CI2)
mul_cox1<- data.frame("HR"=mul_HR,"CI"=mul_CI, "P"=mul_PValue)
write.table(mul_cox1,file="mul_cox.txt",sep="\t",row.names=T,quote=F)
aa=read.table("COX+PR_.txt",header=T,sep="\t",row.names=1) 
library(survival)
library(plyr)
str(aa)
aa$surstat<-factor(aa$surstat)
summary(aa$surstat)
y <- Surv(time=aa$surtime,event=aa$surstat==1)
Uni_cox_model<- function(x){
  FML <- as.formula(paste0 ("y~",x))
  cox<- coxph(FML,data=aa)
  cox1<-summary(cox)
  HR <- round(cox1$coefficients[,2],2)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <-round(cox1$conf.int[,3],2)
  CI95 <-round(cox1$conf.int[,4],2)
  Uni_cox_model<- data.frame('Characteristics' = x,
                             'HR' = HR,
                             'CI5' = CI5,
                             'CI95' = CI95,
                             'p' = PValue)
  return(Uni_cox_model)}  
names(aa)
variable.names<- colnames(aa)[c(3:7)] 
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox <- ldply(Uni_cox,data.frame)
Uni_cox$CI<-paste(Uni_cox$CI5,'-',Uni_cox$CI95)
Uni_cox<-Uni_cox[,-3:-4]
write.table(Uni_cox,file="Uni_coxPR_.txt",sep="\t",row.names=F,quote=F)
View(Uni_cox)
setwd("C:\\Users\\qliul\\Desktop\\TCPA\\总")
library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
a=read.table("GO前20.txt",sep="\t",header=T,check.names=F)
ggplot(a, aes(GeneRatio, Description)) +
  geom_point(aes(y=reorder(Description,GeneRatio),color=p.adjust, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
rt=read.table("KEGG前20.txt",sep="\t",header=T,check.names=F)
library(ggplot2)
ggplot(rt, aes(GeneRatio, Description)) +
  geom_point(aes(y=reorder(Description,GeneRatio),color=p.adjust, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
a=read.table("DO前10.txt",sep="\t",header=T,check.names=F)
ggplot(a, aes(GeneRatio, Description)) +
  geom_point(aes(y=reorder(Description,GeneRatio),color=p.adjust, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(vioplot)
library(ggExtra)
gene="EPPK1"
expFile="mRNAmatrix.txt"    
immFile="CIBERSORT-Results.txt"      
pFilter=0.05    
immune=read.table("CIBERSORT-Results.txt",header=T,sep="\t",check.names=F,row.names=1)
data=read.table("EPPK1表达量PR_.txt",header=T,sep="\t",check.names=F,row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])
sameSample=intersect(row.names(immune),row.names(data))
rt=cbind(immune[sameSample,,drop=F],data[sameSample,,drop=F])
data=rt[,-(ncol(rt)-1)]
data=melt(data,id.vars=c("gene"))
colnames(data)=c("gene","Immune","Expression")
group=levels(factor(data$gene))
data$gene=factor(data$gene,levels=c("low","high"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data,x="Immune",y="Expression",fill="gene",
                  xlab="",
                  ylab="Fraction",
                  legend.title=gene,
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=gene),symnum.args=list(cutpoints=c(0,0.001,0.01,0.05,1),symbols=c("***","**","*","")),label="p.signif")
print(boxplot)
write.table(rt,file="EPPK1+免疫细胞表达量1.txt",sep="\t",quote=F,col.names=T)
a=read.table("T.txt",sep="\t",check.names=F,header=T) #读取文件
library(survival)     
library(ggplot2)
library(ggpubr)
library(survminer)
median(a$T)
group <- ifelse(a$T>0.02599932,"high","low")
a <- cbind(a,group)
fit <- survfit(Surv(futime, fustat)~group, data=a)  
ggsurvplot(fit,pval=TRUE)
ggsurvplot(
  fit,
  data = a,
  size = 1,               
  conf.int = F,          
  pval = TRUE,             
  risk.table = F,       
  risk.table.col = "strata",
  surv.median.line = "hv",
  xlab="Time(days)",
  legend.labs =
    c("High experssion","Low expression"),    
  legend.title = "γδT cells",
  risk.table.height = 0.25, 
  risk.table.size = 2,
  ggtheme = theme_bw()      
)
res.cut <- surv_cutpoint(a,time = "futime",event = "fustat",
                         variables = c("Monocytes"))
summary(res.cut)
plot(res.cut,"Monocytes",palette ="npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(futime,fustat)~Monocytes,data = res.cat)
ggsurvplot(
  fit,
  data = a,
  size = 1,                
  conf.int = F,         
  pval = TRUE,             
  risk.table = F,       
  risk.table.col = "strata",
  surv.median.line = "hv",
  xlab="Time(days)",
  legend.labs =
    c("High experssion","Low expression"),    
  legend.title = "Monocytes",
  risk.table.height = 0.25, 
  risk.table.size = 2,
  ggtheme = theme_bw()     
)
library(ggpubr)             
riskFile="PR.txt"     
tciaFile="TCIA-ClinicalData.tsv"     
tcia=read.csv(tciaFile, header=T, sep="\t", check.names=F, row.names=1)
tcia=tcia[,c("ips_ctla4_neg_pd1_neg","ips_ctla4_neg_pd1_pos","ips_ctla4_pos_pd1_neg","ips_ctla4_pos_pd1_pos")]
tcia<-tcia[complete.cases(tcia),]
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(row.names(risk),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
risk=risk[group==0, ,drop=F]
row.names(risk)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(risk))
risk=cbind(id=row.names(risk), risk)
sameSample=intersect(row.names(tcia), row.names(risk))
tcia=tcia[sameSample, , drop=F]
risk=risk[sameSample, c("PR","EPPK1"), drop=F]
data=cbind(tcia, risk)
fivenum(data$PR)
median(data$EPPK1)
data$EPPK1=ifelse(data$EPPK1>2.659, "high", "low")
data$PR=ifelse(data$PR>9.55690, "PR+", "PR_")
data$EPPK1=ifelse(data$EPPK1=="high", "high", "low")
group=levels(factor(data$EPPK1))
data$EPPK1=factor(data$EPPK1, levels=c("low", "high"))
group=levels(factor(data$EPPK1))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
for(i in colnames(data)[1:(ncol(data)-1)]){
  rt=data[,c(i, "EPPK1")]
  colnames(rt)=c("IPS", "EPPK1")
  gg1=ggviolin(rt, x="EPPK1", y="IPS", fill = "EPPK1", 
               xlab="", ylab=i,
               legend.title="EPPK1",
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  print(gg1)
}
library(limma)
library(plyr)
library(ggplot2)
library(ggpubr)
tideFile="免疫治疗.csv"          
riskFile="C+EPPK1+PR+.txt"      
tide=read.csv(tideFile, header=T, sep=",", check.names=F, row.names=1)
tide$Responder=ifelse(tide$Responder=="TRUE", "responder", "Non-responder")
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0, ,drop=F]
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=cbind(id=row.names(tide), tide)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$EPPK1=factor(risk$EPPK1, levels=c("low", "high"))
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "EPPK1", drop=F]
data=cbind(tide, risk)
p=ggbarplot(tide, x = "id", y = "TIDE", fill = "Responder",
            color = "white",
            palette = c("#6666FF", "#FF6600"),
            sort.by.groups = FALSE,
            xlab = FALSE,
            ylab = "TIDE value",
            width=1,
            legend.title = ""
)+
  theme( axis.ticks.x = element_blank(), axis.text.x = element_blank())
print(p)
dev.off()
rt1=data[,c("Responder", "EPPK1")]
df=as.data.frame(table(rt1))
df=ddply(df, .(EPPK1), transform, percent = Freq/sum(Freq) * 100)
df=ddply(df, .(EPPK1), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
p=ggplot(df, aes(x = factor(EPPK1), y = percent, fill = Responder)) +
  geom_bar(position = position_stack(), stat = "identity", width = .3) +
  scale_fill_manual(values=c("#6666FF","#FF6600"))+
  xlab("EPPK1")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Responder"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_bw()
print(p)
dev.off()

