library(umap)
Control_K27me3_all_peaks<-read.csv("E:/R files/Control_K27me3_all_peaks.csv",row.names = 1)
Control_K27me3_all_peaks<-as.matrix(Control_K27me3_all_peaks)
Control_K27me3_all_peaks<-t(Control_K27me3_all_peaks)
Control_K27me3_all_peaks_umap<-umap(Control_K27me3_all_peaks)
#UMAP dimension reduction to 50
# Control_K27me3_all_peaks_umap_50<-umap(Control_K27me3_all_peaks,n_components=50)
Control_K27me3_all_peaks_umap_matrix<-as.matrix(Control_K27me3_all_peaks_umap[["layout"]])
Control_K27me3_all_peaks_umap_50_project<-umap(Control_K27me3_all_peaks_umap_50_matrix)

#k-means optimal cluster number (Elbow method)
wss_Control_K27me3_all_peaks<-sapply(1:15, function(k){kmeans(Control_K27me3_all_peaks,k,nstart = 10,iter.max = 100)$tot.withinss})
plot(wss_Control_K27me3_all_peaks,type="l")
#k-means clustering
Control_K27me3_all_peaks_k2<-kmeans(Control_K27me3_all_peaks,2,nstart = 10,iter.max = 100)
write.csv(Control_K27me3_all_peaks_k2["cluster"],file="E:/R files/Control_K27me3_all_peaks_k2.csv")
#Louvain clustering
library(igraph)
Control_K27me3_all_peaks_graph<-graph_from_data_frame(Control_K27me3_all_peaks,directed = FALSE)
Control_K27me3_all_peaks_simple<-simplify(Control_K27me3_all_peaks_graph)
Control_K27me3_all_peaks_Louvain<-cluster_louvain(Control_K27me3_all_peaks_graph)
sizes(Control_K27me3_all_peaks_Louvain)
unique(membership(Control_K27me3_all_peaks_Louvain))
as.factor(membership(Control_K27me3_all_peaks_Louvain))