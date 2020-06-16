######this script is mainly used for run k-means based on SFM signautre#######
############_K_means_GSE39582_________#############################
setwd("path")
tar=read.csv("SFM.signautre.csv",head=T,row.names = 1)
head(tar)
dim(tar)
#tar=na.omit(tar)
###########
###########expression matrix can be obtained from curatedCRCData R package
library(curatedCRCData)
library(affy)
data(GSE39582_eset)
#write.csv(exprs(GSE39582_eset), file="GSE39582_eset_exprs.csv")
exp=exprs(GSE39582_eset)
clin=pData(GSE39582_eset)
sel_exp=exp[rownames(exp) %in% rownames(tar),]
dim(sel_exp)
##############################
library(factoextra)
sel_exp=t(sel_exp)

fviz_nbclust(sel_exp, kmeans, method = "gap_stat")
res.hk <-hkmeans(sel_exp, 6)
fviz_dend(res.hk, cex = 0.6, palette = "lancet", 
          rect = TRUE, rect_border = "lancet", rect_fill = TRUE)
group=as.data.frame(res.hk$cluster)
head(group)
colnames(group)="cluster_group"
sub_clin=clin[,c(31,32)]
sub_clin$DFS_status=ifelse(sub_clin$dfs_status=="deceased_or_recurrence", 1, 0)
colnames(sub_clin)[2]="DFS"
head(sub_clin)
com_clin=cbind(group,sub_clin)
head(com_clin)
#################################################################################################
#######################_____heatmap___cluster____################################################
sel_exp=exp[rownames(tar),]
sel_exp[1:5,1:5]
sel_exp=(sel_exp - rowMeans(sel_exp))/apply(sel_exp,1,sd)
getwd()

data=com_clin
newdata=data[order(data$cluster_group),]
SEL_EXP=sel_exp[,rownames(newdata)]
#write.csv(SEL_EXP,"script_out_GSE39582/for_heatmap/sel_exp_for_heatmap_GSE39582.csv")
#write.csv(newdata,"script_out_GSE39582/for_heatmap/cluster_group_GSE39582.csv")
library(ComplexHeatmap)
library(circlize)
mycol <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
annot_df <- data.frame(group=newdata$cluster_group)
# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
#col = list(group = c("3" = "forestgreen", "2" = "darkred","1"="blue4"))

col = list(group = c("3" = "forestgreen", "2" = "darkred","1"="blue4",
                     "4"="purple","5"="tomato2","6"="yellow"))

# Create the heatmap annotation
ha <- HeatmapAnnotation(annot_df, col = col)
# Combine the heatmap and the annotation
Heatmap(SEL_EXP, name = "Group", col = mycol, top_annotation = ha,cluster_columns = FALSE ,
        show_row_names = T,show_column_names = T) 
###########################################________survival analysis
library(survival)
library(survminer)
#data$group <- ifelse(data$score >= 0.296168007, "high", "low")
data=com_clin
head(data)
fit2 <- survfit(Surv(data$DFS, data$DFS_status) ~ cluster_group, data = data)
ggsurvplot(fit2, data = data,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","  logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "group",               # Change legend titles
           legend.labs = c("1","2","3","4","5","6"),  # Change legend labels
           palette = "lancet",                    # Use JCO journal color palette
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE,             # Hide tables y axis text
           ylab=" Disease free survival probability",
           xlab="Time(days)"
)
ggsurvplot(fit2, data = data,
           #title =paste("HR","=", hr," (", hr_left, " - ", hr_right, ")","  logrank","P = ", pvalue, sep=""),
           pval = TRUE, pval.method = TRUE, # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "group",               # Change legend titles
           legend.labs = c("1","2","3","4","5","6"),  # Change legend labels
           palette = c("blue4","darkred","forestgreen","purple","tomato2","yellow"),# Use JCO journal color palette
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE,             # Hide tables y axis text
           ylab=" Disease free survival probability",
           xlab="Time(days)"
)
table(data$cluster_group)
#write.csv(data,"TCGA/script_out_GSE39582/cluster_group_GSE39582.csv")