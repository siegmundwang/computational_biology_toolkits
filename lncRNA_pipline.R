require(readr) # for read files faster
require(dplyr) # for data manipulation
require(NMF) # for cluster 
require(biomaRt) # for ID conversion
require(stringr) # for string manipulation
require(survival) # for survial analysis
require(missForest) # for survival data imputation
require(survMisc) # for cutpoint selection
require(GGally) # for beautiful survival curves
require(hdnom) # for high dimensional cox analysis
require(doParallel) # for parallel computing
require(ggplot2) # the basis of graphics
require(GenVisR) # for waterfall map 
setwd("D:/R/project/lncRNA_cluster") # working directory
options(stringsAsFactors = FALSE) # hideous option

# read the raw lncRNA data, raw mRNA data, expression data for each cluster and
# relationships between mRNA cluster and lncRNA cluster

FastRead <- function(filename) {
  raw <- read_tsv(filename)
  raw <- as.data.frame(raw)
  row.names(raw) <- raw$sample_id
  raw <- raw[, -1]
  raw
}

raw.lncRNA <- FastRead("isoforms.fpkm_table_lncRNA.txt")
raw.mRNA <- FastRead("genes.fpkm_table_mrna.txt")
all.files <- list.files(path = getwd(), full.names = F, pattern = "significant.txt")
exp.cluster.value <- lapply(all.files, function(i) read.table(i, header = T, row.names = 1))
sample.fab.cluster <- FastRead("Annotation_all_5.txt")

# filter

Filter <- function(dt) {
  filtered.dt <- dt[apply(dt > 0.1, 1, sum)/length(dt[1,]) > 0.25, ] 
  dim(filtered.dt) ## 5774
  ldata <- log2(filtered.dt + 1)
  std <- apply(filtered.dt, 1, sd)
  std_cf <- quantile(std,seq(0.1,1,0.1))[8] # top 20%
  ind <- which(std >= std_cf)
  ldata[ind,]
}

filtered.lnc <- Filter(raw.lncRNA)
write.table(filtered.lnc, file = "AML_isoform_k_lncrna_0.10.25sd0.85_log.txt", sep = "\t")
sample.fab <- sample.fab.cluster[, 1:2]

# model fitting
## hardcore nmf
fit.result <- nmf(filtered.lnc, 6, "brunet", nrun = 60, .opt = "vp8")  
pdf(width = 10, height = 8, file = "AML_isoform_k_lncRNA_NMF0.10.250.85.pdf");
consensusmap(fit.result, annCol = sample.fab)
basismap(fit.result)
coefmap(fit.result)
dev.off()

## features choose
ft <- silhouette(fit.result, what = 'features')
data_ft <- filtered.lnc[which(ft[, 3] > 0.1),]
fit.final <- nmf(data_ft, 6, "brunet", nrun = 200, .opt = "vp8")

## plot
si.spl <- silhouette(fit.final, what = "consensus")
pdf(width = 10, height = 8, file = "AML_isoform_k_lncrna_0.10.25sd0.85_fea0.1.pdf");
consensusmap(fit.final, Rowv = runif(ncol(fit.final)),annCol = sample.fab)
si.ft <- silhouette(fit.final, what = 'features')
basismap(fit.final, Rowv = si.ft)
coefmap(fit.final)
plot(si.spl)
dev.off()

## sample order 
sample.order <- rownames(si.spl)

# calculate correlation of mRNA　and lncRNA expression
for (i in 1:6) {
  num_1 <- paste("lnc.cluster.", i, sep = '')
  num_2 <- paste("m.cluster.", i, sep = '')
  judge <- paste("lncRNA_Cluster", i, sep = '')
  assign(num_1, t(na.omit(raw.lncRNA[row.names(exp.cluster.value[[i]]), rownames(sample.fab.cluster[sample.fab.cluster$lncRNA_Cluster == judge, ])])))
  assign(num_2, t(na.omit(raw.mRNA[, rownames(sample.fab.cluster[sample.fab.cluster$lncRNA_Cluster == judge, ])])))
}
rm(i) 
lnc.cluster.3 <- lnc.cluster.3[, colSums(lnc.cluster.3) != 0] # some are all zeros

cor.cluster.1 <- t(na.omit(t(cor(lnc.cluster.1, m.cluster.1))))
cor.cluster.2 <- t(na.omit(t(cor(lnc.cluster.2, m.cluster.2))))
cor.cluster.3 <- t(na.omit(t(cor(lnc.cluster.3, m.cluster.3))))
cor.cluster.4 <- t(na.omit(t(cor(lnc.cluster.4, m.cluster.4))))
cor.cluster.5 <- t(na.omit(t(cor(lnc.cluster.5, m.cluster.5))))
cor.cluster.6 <- t(na.omit(t(cor(lnc.cluster.6, m.cluster.6))))

# get input matrix, output format: lncRNA, mRNA, correlation
input.cluster.1 <- data.frame()
for (i in seq(nrow(cor.cluster.1))) {
  for (j in seq(ncol(cor.cluster.1))) {
    if (abs(cor.cluster.1[i,j]) >= 0.85) {
      tmp.vector <- c(rownames(cor.cluster.1)[i],
                      colnames(cor.cluster.1)[j], cor.cluster.1[i,j])
      input.cluster.1 <- rbind(input.cluster.1, tmp.vector)
    }
  }
}
colnames(input.cluster.1) <- c('lncRNA', 'mRNA', 'correlation')
write.table(input.cluster.1,'cluster1.for.network.txt', row.names = F, sep = '\t') 

input.cluster.2 <- data.frame()
for (i in seq(nrow(cor.cluster.2))) {
  for (j in seq(ncol(cor.cluster.2))) {
    if (abs(cor.cluster.2[i,j]) >= 0.85) {
      tmp.vector <- c(rownames(cor.cluster.2)[i], 
                      colnames(cor.cluster.2)[j], cor.cluster.2[i,j])
      input.cluster.2 <- rbind(input.cluster.2, tmp.vector)
    }
  }
}
colnames(input.cluster.2) <- c('lncRNA', 'mRNA', 'correlation')
write.table(input.cluster.2,'cluster2.for.network.txt', row.names = F, sep = '\t')

input.cluster.3 <- data.frame()
for (i in seq(nrow(cor.cluster.3))) {
  for (j in seq(ncol(cor.cluster.3))) {
    if (abs(cor.cluster.3[i,j]) >= 0.85) {
      tmp.vector <- c(rownames(cor.cluster.3)[i], 
                      colnames(cor.cluster.3)[j], cor.cluster.3[i,j])
      input.cluster.3 <- rbind(input.cluster.3, tmp.vector)
    }
  }
}
colnames(input.cluster.3) <- c('lncRNA', 'mRNA', 'correlation')
write.table(input.cluster.3,'cluster3.for.network.txt', row.names = F, sep = '\t')

input.cluster.4 <- data.frame()
for (i in seq(nrow(cor.cluster.4))) {
  for (j in seq(ncol(cor.cluster.4))) {
    if (abs(cor.cluster.4[i,j]) >= 0.85) {
      tmp.vector <- c(rownames(cor.cluster.4)[i], 
                      colnames(cor.cluster.4)[j], cor.cluster.4[i,j])
      input.cluster.4 <- rbind(input.cluster.4, tmp.vector)
    }
  }
}
colnames(input.cluster.4) <- c('lncRNA', 'mRNA', 'correlation')
write.table(input.cluster.4,'cluster4.for.network.txt', row.names = F, sep = '\t')

input.cluster.5 <- data.frame()
for (i in seq(nrow(cor.cluster.5))) {
  for (j in seq(ncol(cor.cluster.5))) {
    if (abs(cor.cluster.5[i,j]) >= 0.85) {
      tmp.vector <- c(rownames(cor.cluster.5)[i], 
                      colnames(cor.cluster.5)[j], cor.cluster.5[i,j])
      input.cluster.5 <- rbind(input.cluster.5, tmp.vector)
    }
  }
}
colnames(input.cluster.5) <- c('lncRNA', 'mRNA', 'correlation')
write.table(input.cluster.5,'cluster5.for.network.txt', row.names = F, sep = '\t')

input.cluster.6 <- data.frame()
for (i in seq(nrow(cor.cluster.6))) {
  for (j in seq(ncol(cor.cluster.6))) {
    if (abs(cor.cluster.6[i,j]) >= 0.85) {
      tmp.vector <- c(rownames(cor.cluster.6)[i], 
                      colnames(cor.cluster.6)[j], cor.cluster.6[i,j])
      input.cluster.6 <- rbind(input.cluster.6, tmp.vector)
    }
  }
}
colnames(input.cluster.6) <- c('lncRNA', 'mRNA', 'correlation')
write.table(input.cluster.6,'cluster6.for.network.txt', row.names = F, sep = '\t')

# convert XLOC of coexression table to Hugo gene name
coexpression.files <- list.files(path = getwd(), full.names = F, pattern = "network.txt")
coexpression.list <- lapply(coexpression.files, function(i) read.table(i, header = T)) # elements are 'input.cluster.n'
gene.name <- read.table("genes.attr_table", header = T)
gene.name <- gene.name[,c(4,5)] # XLOC and ENST, but one XLOC many ENST
coxnm.cluster.1 <- inner_join(input.cluster.1, gene.name, by = c("mRNA" = "gene_id")) # function in dplyr
coxnm.cluster.2 <- inner_join(input.cluster.2, gene.name, by = c("mRNA" = "gene_id"))
coxnm.cluster.3 <- inner_join(input.cluster.3, gene.name, by = c("mRNA" = "gene_id"))
coxnm.cluster.4 <- inner_join(input.cluster.4, gene.name, by = c("mRNA" = "gene_id"))
coxnm.cluster.5 <- inner_join(input.cluster.5, gene.name, by = c("mRNA" = "gene_id"))
coxnm.cluster.6 <- inner_join(input.cluster.6, gene.name, by = c("mRNA" = "gene_id"))

cluster1.gene <- coxnm.cluster.1$gene_short_name
cluster2.gene <- coxnm.cluster.2$gene_short_name
cluster3.gene <- coxnm.cluster.3$gene_short_name
cluster4.gene <- coxnm.cluster.4$gene_short_name
cluster5.gene <- coxnm.cluster.5$gene_short_name
cluster6.gene <- coxnm.cluster.6$gene_short_name

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
GetNameAndEntrez <- function(genelist) {
  ans <- getBM(attributes=c('external_gene_name', 'entrezgene'), 
               filters = "ensembl_transcript_id", values = genelist, mart = ensembl)
  ans <- unique(ans$external_gene_name)
  ans <- ans[!is.na(ans)]
}
cluster1.gene.hugo <- GetNameAndEntrez(cluster1.gene)
cluster2.gene.hugo <- GetNameAndEntrez(cluster2.gene)
cluster3.gene.hugo <- GetNameAndEntrez(cluster3.gene)
cluster4.gene.hugo <- GetNameAndEntrez(cluster4.gene)
cluster5.gene.hugo <- GetNameAndEntrez(cluster5.gene)
cluster6.gene.hugo <- GetNameAndEntrez(cluster6.gene)

# prepare for GSEA
rownames(gene.name) <- gene.name$gene_id # get gene XLOC
gene.name.filtered <- gene.name[rownames(raw.mRNA),] # get gene XLOC in raw.mRNA
raw.mRNA.gsea <- data.frame("NAME" = gene.name.filtered[rownames(raw.mRNA),2], raw.mRNA) # add ENTS to raw.mRNA
gene.list.simple <- sapply(gene.name.filtered$gene_short_name, function(x) str_sub(x, 1, 15)) # get the total result of ENTS, considering that one XLOC many ENTS
{
  
  gene.symbol.all <- list()
  gene.symbol.all <- getBM(attributes = c('hgnc_symbol', 'ensembl_transcript_id'),
                           filters = "ensembl_transcript_id", 
                           values = gene.list.simple, mart = ensembl)
  gene.symbol.all <- gene.symbol.all[gene.symbol.all$hgnc_symbol != "", ]
  
}

raw.mRNA.gsea$NAME <- sapply(raw.mRNA.gsea$NAME, function(x) str_sub(x, 1, 15)) # only leave one ENTS
raw.mRNA.gsea <- merge(raw.mRNA.gsea, gene.symbol.all, by.x = "NAME", by.y = "ensembl_transcript_id" )
raw.mRNA.gsea$NAME <- raw.mRNA.gsea$hgnc_symbol
raw.mRNA.gsea  <- raw.mRNA.gsea[, -173] # all above convert raw.mRNA name to hgnc_symbol


raw.mRNA.unique <- aggregate(raw.mRNA.gsea, by = list(raw.mRNA.gsea$NAME), FUN = mean) 
raw.mRNA.unique$NAME <- raw.mRNA.unique$Group.1
raw.mRNA.unique <- raw.mRNA.unique[, -1] # all above convert raw.mRNA to unique rows
sur.raw <- raw.mRNA.unique[, c("NAME", sample.order)] # ordered raw mRNA data
write.table(raw.mRNA.gesa.final, file = "for_enrichment_ordered.txt", sep = '\t', row.names = F) # well formated mRNA expression tsv for GSEA
cluster.number <- si.spl[, 1] - 1 # to confirm the cls format for GSEA
sink("gsea_2.cls", append = T)
cat(cluster.number)
sink()

# now for survial analysis
## data manipulation
sur.raw <- as.data.frame(read_tsv("clinical_patient_laml.txt") %>%
  select(ID = bcr_patient_barcode, time, status = dailt))
sur.raw$ID <- gsub("-AB-", "_", sur.raw$ID) # trim the sample name
rownames(sur.raw) <- sur.raw$ID
sur.raw <- sur.raw[sample.order,] # make the sample as same as that in mRNA data
dim(sur.raw)
sur.add.cluster <- cbind(sur.raw, si.spl[]) %>% # same rownames so just using cbind
  select(ID, time, status, cluster)
dim(sur.add.cluster) # find some NAs
imputed.sur.cluster <- missForest(sur.add.cluster[,-1])[[1]]
imputed.sur.cluster$ID <- rownames(imputed.sur.cluster)
imputed.sur.cluster$status <- round(imputed.sur.cluster$status) # some status are not integer

## actually KM survival anylasis
sur.fit <- survfit(Surv(time, status) ~ cluster, 
                   data = imputed.sur.cluster)
survdiff(Surv(time, status) ~ cluster, data = imputed.sur.cluster)
fit.curve <- ggsurv(sur.fit, CI = "def", plot.cens = TRUE, surv.col = "gg.def",
                cens.col = "gg.def", lty.est = 1, lty.ci = 2, size.est = 2,
                size.ci = size.est, cens.size = 2, cens.shape = 3, 
                back.white = FALSE, xlab = "Time", ylab = "Survival",
                main = "Survival Curves") + 
  geom_text(aes(x = 2500, y = 1),label = "p = 0.000207", size = 5)
fit.curve
ggsave(plot = fit.curve,
       filename = "survival_curve.pdf",
       height = 8,
       width = 8)

## actual cox survial analysis based on lncRNA expression profile
novel.rna.data <- read_tsv("isoforms.fpkm_table_exp_k+n.ordered.txt")
survial.final.data <- read_tsv("Survival_final.txt")
all.lncRNA <- Reduce(union, list(colnames(lnc.cluster.1), 
                                 colnames(lnc.cluster.2), 
                                 colnames(lnc.cluster.3), 
                                 colnames(lnc.cluster.4), 
                                 colnames(lnc.cluster.5), 
                                 colnames(lnc.cluster.6)))
sur.quant <- data.frame(t(raw.lncRNA[all.lncRNA,]))
sur.quant$ID <- rownames(sur.quant) # all the lncRNA expression for all samples
sur.quant <- merge(sur.quant, imputed.sur.cluster, by = "ID") %>%
  filter(time > 0)# get the expression and time
## penalized Cox model
x <- as.matrix(sur.quant %>% select(-ID, -time, -status, -cluster)) 
y <- with(sur.quant, Surv(time, status))
model.sur.mnet <- hdcox.mnet(x, y, nfolds = 10, parallel = T) # 
lnc.mnet <- as.data.frame(model.sur.mnet$mnet_model$beta)
lnc.mnet$ID <- rownames(lnc.mnet)
lnc.mnet <- lnc.mnet %>% filter(`2.4745` != 0)
key.lnc <- lnc.mnet$ID
sur.quant$lnc.score <- as.vector(as.matrix(sur.quant[, key.lnc]) %*% 
  as.matrix(lnc.mnet$`2.4745`)) # must be rigid matrix  

## lnc.score KM survival anylasis
sur.lnc.score <- sur.quant %>% 
  select(ID, time, status, lnc.score) %>% arrange(lnc.score)
write_tsv(sur.lnc.score, path = "sur_lnc_score.txt")
cut.point <- cutp(survfit(Surv(time, status) ~ lnc.score,
                          data = sur.lnc.score))$lnc.score$lnc.score[1]
sur.lnc.score <- bind_rows((sur.lnc.score %>% 
                             filter(lnc.score > cut.point) %>% mutate(risk = 1)), 
                           (sur.lnc.score %>%  filter(lnc.score <= cut.point) 
                            %>% mutate(risk = 0))) 
sur.fit.lnc <- survfit(Surv(time, status) ~ risk, data = sur.lnc.score)
survdiff(Surv(time, status) ~ risk, data = sur.lnc.score)
plot(sur.fit.lnc)
lnc.curve <- ggsurv(sur.fit.lnc, CI = "def", plot.cens = TRUE, surv.col = "gg.def",
                    cens.col = "gg.def", lty.est = 1, lty.ci = 2, size.est = 2,
                    size.ci = size.est, cens.size = 2, cens.shape = 3, 
                    back.white = FALSE, xlab = "Time", ylab = "Survival",
                    main = "Survival Curves") + 
  geom_text(aes(x = 2500, y = 1),label = "p = 2.63e-10", size = 5)
lnc.curve
ggsave(plot = lnc.curve,
       filename = "lnc_survival_curve.pdf",
       height = 8,
       width = 8)

# draw waterfall map
raw.mutaion <- read.delim("genvis.tsv") %>% select(sample = TCGA_id, 2, 10 ,17)
raw.mutaion[raw.mutaion$trv_type == '-', 3] <- NA
pdf("waterfall.pdf", width = 10, height = 10)
waterfall(raw.mutaion, fileType = "MGI", mainRecurCutoff = 0.06)
dev.off()
