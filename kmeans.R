K27ac_EB_enhancer<-read.csv("E:/R files/K27ac_EB_enhancer.csv",row.names = 1)
K27ac_EB_enhancer<-as.matrix(K27ac_EB_enhancer)
wss_K27ac_EB_enhancer<-sapply(1:15, function(k){kmeans(K27ac_EB_enhancer,k,nstart = 10,iter.max = 100)$tot.withinss})
plot(wss_K27ac_EB_enhancer,type="l")

K27ac_EB_enhancer_cluster_k5<-kmeans(K27ac_EB_enhancer,5,nstart = 10,iter.max = 100)
write.csv(K27ac_EB_enhancer_cluster_k5["cluster"],file="E:/R files/K27ac_EB_enhancer_cluster_k5.csv")

# sort the original matrix based on the cluster number assignment from "K27ac_EB_enhancer_cluster_k5.csv" and save the re-ordered matrix to "K27ac_EB_enhancer_k5.csv"

K27ac_EB_enhancer_k5<-read.csv("E:/R files/K27ac_EB_enhancer_k5.csv",row.names = 1)
K27ac_EB_enhancer_k5<-as.matrix(K27ac_EB_enhancer_k5)

break_K27ac_EB_enhancer<-seq(min(K27ac_EB_enhancer),max(K27ac_EB_enhancer),by=0.01)
palette_K27ac_EB_enhancer<- c(colorRampPalette(c("white", "darkblue"))(220),colorRampPalette(c("darkblue", "darkblue"))(length(break_K27ac_EB_enhancer)-220))
heatmap(K27ac_EB_enhancer_k5,Rowv = NA,Colv = NA,scale = "none",col=palette_K27ac_EB_enhancer,break(break_K27ac_EB_enhancer),labRow = NA,labCol = NA)
image(break_K27ac_EB_enhancer,1,as.matrix(1:length(break_K27ac_EB_enhancer)),col=palette_K27ac_EB_enhancer,ylab="",yaxt="n",xlab="")
