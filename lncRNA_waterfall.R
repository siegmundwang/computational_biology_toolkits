library(NMF)
raw <- data.frame(read.table("isoforms.fpkm_table_exp.known.txt", sep="\t", row.names=1, header=T))
annotation_col <- data.frame(read.table("Annotation_all_4.txt", sep="\t", row.names=1, header=T))
dim(raw);

## filter1
ldata <- log2(raw+1)
std <- apply(ldata, 1, sd)
std_cf <- quantile(std,seq(0.1,1,0.01))[62] # 71%
ind <- which(std>=std_cf)
data <- ldata[ind,]
dim(data)
write.table(data,file="isoform_lncrna_exp_sd0.71_log.txt", sep="\t")

## consense est
# estim.r <- nmf(data, 2:8, nrun = 60, seed = 123456, .opt="vp8")
# pdf(width=15,height=11,file="AML_isoform_lncrna_exp_NMF_sd0.71_log-est2-8.pdf");
# par(mfcol=c(1,2))
# plot(estim.r)
# consensusmap(estim.r,annCol = annotation_col)
# dev.off()

res <- nmf(data,6,"brunet",nrun=60, seed= 123456, .opt="vp8") # indis
pdf(width=10,height=8,file="AML_lncRNA_exp_NMF_0.71.pdf");
consensusmap(res, Rowv = runif(ncol(res)),annCol = annotation_col)
si <- silhouette(res, what = 'features')
basismap(res,Rowv = si)
coefmap(res)
plot(silhouette(res, what='consensus'))
dev.off()

## predict
#predict(res, prob=TRUE)
#predict(res, 'rows', prob=TRUE)


## extractFeatures
fea_e <- extractFeatures(res,0.5)
fea_cls1 <- fea_e[[1]]
fea_cls2 <- fea_e[[2]]
fea_cls3 <- fea_e[[3]]
fea_cls4 <- fea_e[[4]]
fea_cls5 <- fea_e[[5]]
fea_cls6 <- fea_e[[6]]
list_fea <- c(fea_cls1,fea_cls2,fea_cls3,fea_cls4,fea_cls5,fea_cls6)
fea_edata <- data[list_fea,]
dim(fea_edata)
# res_feae <- nmf(fea_edata, 6, "brunet", nrun=60, seed= 123456, .opt="vp6")
# pdf(width=10,height=8,file="AML_lncRNA_exp_NMF_0.71_feae0.5.pdf");
# consensusmap(res_feae, Rowv = runif(ncol(res_feae)),annCol = annotation_col)
# si <- silhouette(res_feae, what = 'features')
# basismap(res_feae,Rowv = si)
# coefmap(res_feae)
# plot(silhouette(res_feae, what='consensus'))
# dev.off()


## features choose
con <- silhouette(res,what = 'consensus')
fea <- silhouette(res,what = 'features')

fea_cdata <- data[which(fea[,3]>0.55),] # sil>0.55
#estim.r1 <- nmf(fea_cdata, 2:8, nrun = 60, seed = 123456, .opt="vp8")

# feac.random <- randomize(fea_cdata)
# estim.r1.random <- nmf(feac.random, 2:8, nrun = 60, seed = 123456, .opt="vp6")


# pdf(width=15,height=11,file="AML_isoform_lncrna_exp_NMF_sd0.71_feac0.55_log-est2-8.pdf");
# par(mfcol=c(1,2))
# plot(estim.r1, estim.r1.random)
# consensusmap(estim.r1,annCol = annotation_col)
# dev.off()

res_feac <- nmf(fea_cdata, 6, "brunet", nrun=60, seed=123456, .opt="vp8") # ultimate
pdf(width=10,height=8,file="AML_lncRNA_exp_NMF_0.71_feac0.55.pdf");
consensusmap(res_feac, Rowv = runif(ncol(res_feac)),annCol = annotation_col)
si_ft <- silhouette(res_feac, what = 'features')
si_spl <- silhouette(res_feac, what = 'consensus')
basismap(res_feac,Rowv = si)
coefmap(res_feac)
plot(si_spl)
dev.off()


## Fertig alles
order_cluster <- as.data.frame(si_spl[,]) # ordered cluster crucial for following any
##### 这一步不甚精通，silhouette本质为带有属性的list，具体面向对象的细节还要在下边掌握

## gradually vectorization
options(stringsAsFactors = FALSE)
setwd("D:/R/project/lncRNA_cluster") # if in my personal computer

sur_raw <- read.table("clinical_patient_laml.txt", stringsAsFactors = F, header = T, sep = "\t")
sur_raw_for_anl <- data.frame(sur_raw$bcr_patient_barcode, sur_raw$time, sur_raw$dailt)
names(sur_raw_for_anl) <- c("ID", "time", "status")
sur_raw_for_anl$ID <-gsub("-AB-", "_", sur_raw_for_anl$ID)
sur_raw_for_anl <- na.omit(sur_raw_for_anl)
rownames(sur_raw_for_anl) <- sur_raw_for_anl$ID
# sur_raw_for_anl <- sur_raw_for_anl[names(raw_mRNA_unique)[-1],] 
dim(sur_raw_for_anl)
# sample_cluster_ordered$ID <- rownames(sample_cluster_ordered)
order_cluster$ID <- rownames(order_cluster)
sur_cluster <- merge(order_cluster, sur_raw_for_anl)
sur_cluster <- sur_cluster[,c(1,2,5,6)] # 完成了数据的整理
dim(sur_cluster)
sur_cluster$ID <- as.character(sur_cluster$ID) # sucks to autotransfer to factor!

require(GGally)
require(survival)
sur.fit <- survfit(Surv(time, as.numeric(status))~cluster, data = sur_cluster)
#autoplot(sur.fit, conf.int = FALSE, surv.linetype = 'solid', lwd = 5) # 最基本的绘图
#autoplot(sur.fit, surv.linetype = 'dashed', conf.int = FALSE, censor.shape = '*', censor.size = 5, facets = TRUE, ncol = 2)
ggplt <- ggsurv(sur.fit, CI = "def", plot.cens = TRUE, surv.col = "gg.def",
       cens.col = "gg.def", lty.est = 1, lty.ci = 2, size.est = 2,
       size.ci = size.est, cens.size = 2, cens.shape = 3, back.white = FALSE,
       xlab = "Time", ylab = "Survival", main = "") + geom_text(aes(x = 2500, y = 1), label = "p= 0.0186",size = 5)
ggplt
#plot(sur.fit, lty = 1:6, lwd = 1:6, ylab = "S(t)", xlab = "t", main = "Survival Functions")
#legend(x = "topright", legend = c("cluster1", "cluster2","cluster3" ,"cluster4", "cluster5", "cluster6"), lty = 1:6)
survdiff(Surv(time, as.numeric(status))~cluster, data = sur_cluster)
# 

sur_raw_new <- read.table('TCGA_AML_clinic_new.txt', sep = '\t', header = T)
sur_raw_new$ID <-gsub("-AB-", "_", sur_raw_new$Index)
sur_raw_new <- sur_raw_new[,-1]
sur_cluster_new <- merge(order_cluster, sur_raw_new)
sur_cluster_new <- sur_cluster_new[,-c(3,4)]
sur.fit.new <- survfit(Surv(time, as.numeric(status))~cluster, data = sur_cluster)
survdiff(Surv(time, as.numeric(status))~cluster, data = sur_cluster)
ggplt <- ggsurv(sur.fit.new, CI = "def", plot.cens = TRUE, surv.col = "gg.def",
                cens.col = "gg.def", lty.est = 1, lty.ci = 2, size.est = 2,
                size.ci = size.est, cens.size = 2, cens.shape = 3, back.white = FALSE,
                xlab = "Time", ylab = "Survival", main = "") + geom_text(aes(x = 2500, y = 1), label = "p= 0.0186",size = 5)
ggplt

#实打实的cox
require(mfp)
require(hdnom)
require(survival)
suppressMessages(require("doParallel"))
registerDoParallel(detectCores())

x <- as.matrix(sur_cluster_final1_any[, -c(1,2)])
y <- with(sur_cluster_final1_any, Surv(time, status))
cluster1.aenetfit <- hdcox.aenet(x, y, nfolds = 10, rule = "lambda.1se", seed = c(5, 7), parallel = T) # 目前的问题是特征太多

fit1 <- cluster1.aenetfit$aenet_model
alpha1 <- cluster1.aenetfit$aenet_best_alpha
lambda1 <- cluster1.aenetfit$aenet_best_lambda
adapen1 <- cluster1.aenetfit$pen_factor

suppressMessages(require("rms"))
x.df <- as.data.frame(x)
dd <- datadist(x.df)
options(datadist = "dd")
time1 <- sur_cluster_final1_any$time
status1 <- sur_cluster_final1_any$status
nom1 <- hdnom.nomogram(fit1, model.type = "aenet", x, time1, status1, x.df, lambda = lambda, pred.at = 365 * 2,funlabel = "2-Year Overall Survival Probability")
plot(nom1)

val.int1 <- hdnom.validate(x, time1, status1, model.type = "aenet", alpha = alpha, lambda = lambda, pen.factor = adapen, method = "bootstrap", boot.times = 10, tauc.type = "UNO", tauc.time = seq(1, 5, 0.5) * 365, seed = 42, trace = FALSE)
val.int1
summary(val.int1)
plot(val.int1)

#cal.int1 <- hdnom.calibrate(x, time1, status1, model.type = "aenet", alpha = alpha, lambda = lambda, pen.factor = adapen, method = "bootstrap", boot.times = 10, pred.at = 365 * 5, ngroup = 3, seed = 42, trace = FALSE) #group的意义不良定???
#cal.int1
#plot(cal.int1)

cmp.val1 <- hdnom.compare.validate(x, time1, status1, model.type = c("lasso", "alasso"), method = "cv", nfolds = 5, tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365, seed = 42, trace = FALSE)
cmp.val1
plot(cmp.val1)
plot(cmp.val, interval = TRUE)

cox.cluster2 <- mfp(Surv(time, as.numeric(status)) ~ ., family = cox, data = sur_cluster_final2_any, select = 0.05, verbose = T) # 目前的问题是特征太多
summary(cox.cluster2)
plot(survfit(cox.cluster2))

imp_cluster3 <- read.table("cluster3.for.network.txt", header = T, sep = "\t")
imp3 <- as.data.frame(table(imp_cluster3$lncRNA))
imp3 <- imp3[order(imp3$Freq, decreasing = T),]
imp_fl3 <- paste(imp3[imp3$Freq > 10,1], collapse = "+") # 过滤机制
cox.cluster3 <- mfp(Surv(time, as.numeric(status)) ~ ., family = cox, data = sur_cluster_final3_any, select = 0.05, verbose = T) # 无模型？
summary(cox.cluster3)

imp_cluster4 <- read.table("cluster4.for.network.txt", header = T, sep = "\t")
imp4 <- as.data.frame(table(imp_cluster4$lncRNA))
imp4 <- imp4[order(imp4$Freq, decreasing = T),]
imp_fl4 <- paste(imp4[imp4$Freq > 10,1], collapse = "+")
cox.cluster4 <- mfp(Surv(time, as.numeric(status)) ~ ., family = cox, data = sur_cluster_final4_any, select = 0.05, verbose = T) # 无模型？
summary(cox.cluster4)

imp_cluster5 <- read.table("cluster5.for.network.txt", header = T, sep = "\t")
imp5 <- as.data.frame(table(imp_cluster5$lncRNA))
imp5 <- imp5[order(imp5$Freq, decreasing = T),]
imp_fl5 <- paste(imp5[imp5$Freq > 10,1], collapse = "+")
cox.cluster5 <- mfp(Surv(time, as.numeric(status)) ~ ., family = cox, data = sur_cluster_final5_any, select = 0.05, verbose = T)
summary(cox.cluster5) # 同样是不限制

imp_cluster6 <- read.table("cluster6.for.network.txt", header = T, sep = "\t")
imp6 <- as.data.frame(table(imp_cluster6$lncRNA))
imp6 <- imp6[order(imp6$Freq, decreasing = T),]
imp_fl6 <- paste(imp6[imp6$Freq > 1,1], collapse = "+")
cox.cluster6 <- mfp(Surv(time, as.numeric(status)) ~ ., family = cox, data = sur_cluster_final6_any, select = 0.05, verbose = T)
summary(cox.cluster6)# 显著???

# 下面绘制waterfall
require(GenVisR)
setwd("D:\\R\\project\\lncRNA_cluster")
amlgmi.raw <- read.table("genvis.tsv", header = T, stringsAsFactors = F, sep = '\t')
amlgmi <- amlgmi.raw[, c(2, 10, 17)]
names(amlgmi)[1] <- "sample"
amlgmi[amlgmi$trv_type=='-', 3] <- NA
pdf("waterfall.pdf")
waterfall(amlgmi, fileType = "MGI")
dev.off()
