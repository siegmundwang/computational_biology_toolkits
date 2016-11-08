# 载入所需的各种包, 默认各类数据在当前目录 
require(readr) # 更快地读取数据
require(dplyr) # 数据操作
require(survival) # 生存分析
require(survMisc) # 选择生存分析中的cutpoint
require(GGally) # 更好的生存曲线
require(hdnom) # 高维数据生存分析
require(doParallel) # 并行运算
require(ggplot2) # 绘图基础

# 进行数据预处理
novel.rna.data <- read.delim("isoforms.fpkm_table_exp_k+n.ordered.txt", 
                            row.names = 1, 
                            stringsAsFactors = F) # 读取所有lncRNA表达数据, header默认为True, 选取第一列为行名
novel.rna.data <- as.data.frame(t(novel.rna.data)) # 转置从而方便同生存数据的合并
novel.rna.data$sample_id <- rownames(novel.rna.data)
survival.raw.data <- read_tsv("Survival_final.txt")[c("sample_id", 
                                                      "OS_months", "Status")]
survival.data <- merge(novel.rna.data, survival.raw.data, by = "sample_id") # 将表达谱数据和生存数据合并

# 进行生存分析
x <- survival.data %>% select(-sample_id, -OS_months, -Status) # 选中表达数据作为生存分析的变量
y <- with(survival.data, Surv(OS_months, Status)) # 建立生存对象
model.sur.mnet <- hdcox.mnet(x, y, nfolds = 10, parallel = T) # 利用mnet方法进行cox生存分析, 方法可以自行改变, 运算时间较长
## 有多种方法可进行高维数据生存分析, 
## "lasso", "alasso", "flasso", "enet", "aenet", "mcp", "mnet", "scad", 
## "snet"等, 经过多次试验, 发现`mnet`方法最好. 
lnc.mnet <- as.data.frame(model.sur.mnet$mnet_model$beta) # 获得具体的线性模型
lnc.mnet$ID <- rownames(lnc.mnet)
lnc.mnet <- lnc.mnet %>% filter(`0.6313` != 0) # 0.6313为截距项, 同时是系数列的列名, 筛得拟合后的特征
key.lnc <- lnc.mnet$ID # 得到关键的lncRNA名
lnc.score <- novel.rna.data[lnc.mnet$ID] # 结合原有的表达数据, 取其子集
lnc.score$score <- as.vector(as.matrix(lnc.score) %*% lnc.mnet$`0.6313`) # 获得每一个样本的lncRNA分数, 用以危险分层, 矩阵相乘后返回矩阵, 故需显式转换为向量
lnc.score$sample_id <- rownames(lnc.score)
survival.lnc.score <- merge(lnc.score, survival.raw.data, by = "sample_id")
survival.lnc.score <- survival.lnc.score[order(survival.lnc.score$score), ] # 以lncRNA分数进行排序

# 分割lncRNA分数, 危险分层
cut.point <- cutp(survfit(Surv(OS_months, Status) ~ score, data = survival.lnc.score))$score$score[1] # 原有cutp返回复杂的数据结构, 此命令一行获取cut点
survival.lnc.score$risk <- 0 # 定义风险变量 
survival.lnc.score[survival.lnc.score$score > cut.point, "risk"] <- 1 # 大于cut.point为1
sur.fit.lnc <- survfit(Surv(OS_months, Status) ~ risk, data = survival.lnc.score)
survdiff(Surv(OS_months, Status) ~ risk, data = survival.lnc.score) # 结果比较理想, p值可视为0
plot(sur.fit.lnc) # 绘制生存曲线的草图
lnc.score.curve <- ggsurv(sur.fit.lnc, CI = "def", plot.cens = TRUE, surv.col = "gg.def",
                    cens.col = "gg.def", lty.est = 1, lty.ci = 2, size.est = 2,
                    size.ci = size.est, cens.size = 2, cens.shape = 3, 
                    back.white = FALSE, xlab = "Time", ylab = "Survival",
                    main = "Survival Curves") + 
  geom_text(aes(x = 120, y = 1),label = "p ~ 0", size = 5) # 绘制生存曲线
lnc.score.curve
ggsave(plot = lnc.score.curve,
       filename = "lnc_survival_curve.pdf",
       height = 8,
       width = 8) # 保存生存曲线
