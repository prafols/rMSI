#Mask selection RAMAN
#Gerard Baquer Gómez 16/04/2021

#1. Load Datasets
folder<-"/home/gbaquer/msidata/3. Coregistration/Zebra_Fish"
RAMAN_tags<-c("CTRL","NS","S","MS","HS")
RAMAN_imzml<-c("/imzml/Ctrl_50x_514nm_2s_100__map_7um_center1500_1_CRR_BC_Sm_Copy.imzML",
               "/imzml/NS_50x_514nm_2s_6um_100__map_3_CRR_BC_Sm_Copy.imzML",
               "/imzml/S_50x_514nm_2s_100__map_5um_center1500_1_CRR_BC_Sm_Copy.imzML",
               "/imzml/MS_50x_514nm_2s_6um_100__map_1_CRR_BC_Sm_Copy.imzML",
               "/imzml/HS_514nm_50x_100__2s_5um_center1500_map1_CRR_BC_Sm_Copy.imzML")
RAMAN_dataset<-lapply(paste(folder,RAMAN_imzml,sep=""),rMSI::LoadMsiData)
RAMAN_dataset_normalized<-lapply(RAMAN_dataset,function(x)rMSI::NormalizeTIC(rMSI::NormalizeRMS(x)))
RAMAN_datacube<-lapply(RAMAN_dataset,function(x)rMSIproc::loadDataCube(x,cubeSel = 1))
RAMAN_datacube_normalized <- lapply(seq_along(RAMAN_dataset),function(i)RAMAN_datacube[[i]]/RAMAN_dataset_normalized[[i]]$normalizations$TIC)
RAMAN_pos <- lapply(RAMAN_dataset,function(x)x$pos)

#2. Load Masks
RAMAN_mask_png<-c("/masks/CTRL_mask1.png",
               "/masks/NS_mask1.png",
               "/masks/S_mask1.png",
               "/masks/MS_mask1.png",
               "/masks/HS_mask1.png")

RAMAN_mask<-lapply(paste(folder,RAMAN_mask_png,sep=""),png::readPNG)

lapply(RAMAN_mask,mmand::display)

#3. Check dimensions

a<-sapply(RAMAN_dataset,function(x)x$size)
b<-sapply(RAMAN_mask,function(x)dim(x)[c(2,1)])
a/b #dimensions don't match

#4. PCA with stochastic resonance (kmeans)
RAMAN_datacube_clusered_ions<-lapply(RAMAN_datacube_normalized,function(x)apply(x,2,function(y)kmeans(y,5)$clus))
RAMAN_datacube_pca<-lapply(RAMAN_datacube_clusered_ions,prcomp)
RAMAN_datacube_pca_n<-lapply(RAMAN_datacube_normalized,prcomp)

#It doesn't improve the normal pca

#5. PCA with stochastic resonance (thresholding)

RAMAN_datacube_autoscaled<-lapply(RAMAN_datacube,function(x)apply(x,2,function(y)(y-min(y))/(max(y)-min(y))))
RAMAN_datacube_SR<-lapply(RAMAN_datacube_autoscaled,function(x)x>0.1)

RAMAN_TIC_log<-lapply(RAMAN_datacube,function(x)log(apply(x,1,sum)))
for(i in seq_along(RAMAN_dataset))
  RAMAN_TIC_log[[i]][!is.finite(RAMAN_TIC_log[[i]])]<-0

RAMAN_mask_kmeans<-lapply(RAMAN_TIC_log,function(x)kmeans(x,2)$clus)

thresholds<-c(0.7,0.68,0.56,0.54,0.56)#Hand selected thresholds for each dataset
RAMAN_mask_threshold<-lapply(seq_along(RAMAN_dataset),function(i)RAMAN_TIC_log[[i]]>thresholds[i]*max(RAMAN_TIC_log[[i]]))


rMSI::PlotValues(RAMAN_pos[[4]],RAMAN_mask_threshold[[4]])

#6. Convolution
convolve<-function(pos,val){
  new_val=val
  pos_tag<-apply(pos,1,function(x)paste(x,collapse=" "))
  for(i in 1:nrow(pos)){
    n<-expand.grid(pos[i,"x"]+c(-1,0,1),pos[i,"y"]+c(-1,0,1))
    n_i <- which(pos_tag %in% apply(n,1,function(x)paste(x,collapse=" ")))
    new_val[i]<-sum(val[n_i])>length(n_i)/2
  }
  return(new_val)
}

RAMAN_mask_convolved<-lapply(seq_along(RAMAN_dataset),function(i)convolve(RAMAN_pos[[i]],RAMAN_mask_threshold[[i]]))

rMSI::PlotValues(RAMAN_pos[[5]],RAMAN_mask_convolved[[5]])


#7. Filter Datasets

RAMAN_datacube_filtered<-lapply(seq_along(RAMAN_dataset),function(i)RAMAN_datacube[[i]][RAMAN_mask_convolved[[i]],])
RAMAN_pos_filtered<-lapply(seq_along(RAMAN_dataset),function(i)RAMAN_pos[[i]][RAMAN_mask_convolved[[i]],])

rMSI::PlotValues(RAMAN_pos_filtered[[5]],rep(1,nrow(RAMAN_pos_filtered[[5]])))
rMSI::PlotValues(RAMAN_pos_filtered[[5]],RAMAN_datacube_filtered[[5]][,1])
