library(DiffBind)
Control_K27ame3_cluster <- dba(sampleSheet = "Control_K27me3_1-6_7-19.csv",dir=system.file("extra", package="DiffBind"))
Control_K27ame3_cluster_count <- dba.count(Control_K27ame3_cluster)

Control_K27ame3_cluster_contrast <- dba.contrast(Control_K27ame3_cluster_count,categories = DBA_CONDITION, minMembers = 2)
Control_K27ame3_cluster_results <- dba.analyze(Control_K27ame3_cluster_contrast)
Control_K27ame3_cluster_results.DB <- dba.report(Control_K27ame3_cluster_results)
write.table(Control_K27ame3_cluster_results.DB, file = "Control_K27ame3_cluster_Diff.csv", sep = ",")

Control_K27ame3_cluster_plot <- plot(Control_K27ame3_cluster_count)
write.csv(Control_K27ame3_cluster_plot,file="E:/R files/Control_K27ame3_cluster_plot.csv")