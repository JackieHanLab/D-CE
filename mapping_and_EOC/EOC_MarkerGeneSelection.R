#Gene expression matrix, TPM or FPKM.
TPM <- read.delim("./TPM.txt", stringsAsFactors=FALSE)

#Reconstructed coordinates
coords <- read.delim("./Reconstructed_coordinates.txt", header=FALSE)
coords <- as.matrix(coords)

#The genes with most significant p-value is the best marker gene.The process can be a little slow, so if you have 10,000 genes, it takes about an hour.
gradient_cor <- matrix(nrow = nrow(TPM),ncol = 3)
gradient_cor <- as.data.frame(gradient_cor)
colnames(gradient_cor) <- c("gene","MarkerGene_RCC","pvalue")
gradient_cor$gene <- rownames(TPM)

for (m in 1:nrow(gradient_cor)) {
  genei <- gradient_cor$gene[m]
  Color1 <- TPM[rownames(TPM) ==genei,]
  coord<-Color1
  coord<-as.matrix(t(coord))
  
  corx<-matrix(ncol = 30,nrow = 30)
  corp<-matrix(ncol = 30,nrow = 30)
  for(i in 1:30){
    for(j in 1:30){
      alaph<-(pi/15)*i
      bet<-(pi/15)*j
      transmat<-matrix(nrow = 3,ncol = 3)
      transmat[1,1]<-(cos(alaph)*cos(bet))
      transmat[1,2]<-sin(alaph)
      transmat[1,3]<-(cos(alaph)*sin(bet))
      transmat[2,1]<-((-cos(bet))*sin(alaph))
      transmat[2,2]<-cos(alaph)
      transmat[2,3]<-((-sin(bet))*sin(alaph))
      transmat[3,1]<-(sin(bet)*(-1))
      transmat[3,2]<-0
      transmat[3,3]<-cos(bet)
      coordstrans<-coords%*%transmat
      Cor <- cor.test(coordstrans[,1],coord[,1],method = "spearman")
      corp[i,j]<-Cor$p.value
      corx[i,j]<-Cor$estimate
    }
  }
  
  gradient_cor$MarkerGene_RCC[m] <- corx[which(corp==min(corp))]
  gradient_cor$pvalue[m] <- min(corp)
  print(m)
}