#rMSIbenford
library(benford.analysis) #Dependencia opcional
pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/2020_AgPaper_Dataset_1.zip")

#EXploration of the validity of Benfort Law and its applications to MSI data

#TOF study

#Helper functions
#1.Load all datasets
files <- list.files("/home/gbaquer/msidata/2. Benford Law/TOF Datasets",pattern = "\\.zip$",full.names = T)
for(f in files)
  generate_benfort_report(rMSIproc::LoadPeakMatrix(f),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_10_GeneralResults/",images = F,magnitude="intensity",digits=2)
#2. Generate Benfort Report [Graphs, MAD image, MAD spectrum]
generate_benfort_report <- function(pks,path,split=T,images=F,extra_name="",magnitude="area",digits=3,norm="",percentage=0.2)
{
  #Create subfolder
  n<-unlist(lapply(strsplit(pks$names,".",fixed = T),function(x)x[1]))
  #Remove 0s
  #pks$intensity[which(pks$intensity==0)]=10e-6
  #Compute benford on mean spectra
  for(i in 1:length(pks$numPixels)){
    print(i)
    tiff(paste(path,n[i],extra_name,"_1_benford_mean_spectra.tiff",sep=""))
    # bf<-benford.analysis::benford(apply((rMSIcleanup:::get_one_peakMatrix(pks,i))$area,2,mean),1,discrete = F)
    bf<-benford.analysis::benford((rMSIcleanup:::get_one_peakMatrix(pks,i))[[magnitude]],discrete = F,number.of.digits = digits)
    bf$info$data.name=n[i]
    plot(bf)
    dev.off()
  }

  if(images){
    #Compute benford on kmeans
    #Compute benford pixel-wise
    tiff(paste(path,n[1],extra_name,"_2_benford_pixel_MAD.tiff",sep=""),width=1500,height=1500)
    pixel_bf<-pixel_benford(pks,magnitude,digits,norm,percentage)
    MAD<-unlist(lapply(pixel_bf,function(x)x$MAD))
    rMSIproc::plotValuesImage(pks,MAD,ncol=3) 
    dev.off()
    
    tiff(paste(path,n[1],extra_name,"_3_benford_pixel_distortion.tiff",sep=""),width=1500,height=1500)
    dist<-unlist(lapply(pixel_bf,function(x)x$distortion.factor))
    dist[which(is.na(dist))]=max(dist,na.rm = T)
    rMSIproc::plotValuesImage(pks,dist,ncol=3) 
    dev.off()
    #Compute benford ion-wise
  }
}
#3. Add noise
add_noise<-function(pks,size=0.3,type="gaussian",level=0.1){
  noisy_pks<-pks
  #Generate area
  min<-apply(pks$pos,2,min)
  max<-apply(pks$pos,2,max)
  range<-(max-min)
  mid<-min+range/2
  d<-size*range/2
  area<-tmp<-apply(pks$pos,1,function(p)(p[1]<(mid[1]+d))&&(p[1]>(mid[1]-d))&&(p[2]<(mid[2]+d))&&(p[2]>(mid[2]-d)))
  
  #Add noise to area
  mean_TIC<-mean(apply(pks$intensity,1,sum))/length(pks$mass)
  if(type=="gaussian"){
    noise=rnorm(length(area),sd=level*mean_TIC)
  }
  if(type=="poisson"){
    noise=rpois(length(area),lambda = level*mean_TIC)
  }
  noisy_pks$intensity[area,]<-noisy_pks$intensity[area,]+noise
  noisy_pks$intensity[which(noisy_pks$intensity<0)]<-0
  noisy_pks$area[area,]<-noisy_pks$area[area,]+noise
  noisy_pks$area[which(noisy_pks$area<0)]<-0
  return(noisy_pks)
}


non_empty_pixels<-function(pks,magnitude="intensity",percentage=0.1)
{
  TICs<-apply(pks[[magnitude]],1,sum)
  return((TICs>(median(TICs)*percentage)))
  
}
ion_benford<-function(pks,magnitude="area",digits=2,norm="",percentage=0){
  return(apply(normalize(pks,magnitude,norm)[non_empty_pixels(pks,magnitude,percentage),],2,function(x)benford.analysis::benford(x,number.of.digits = digits,discrete=F)))
}
pixel_benford<-function(pks,magnitude="area",digits=3,norm="",percentage=0){
  return(apply(normalize(pks,magnitude,norm)[non_empty_pixels(pks,magnitude,percentage),],1,function(x)benford.analysis::benford(x,number.of.digits = digits,discrete=F)))
}
normalize<-function(pks,magnitude="area",norm="")
{
  if(norm=="")
    return(pks[[magnitude]])
  else
    return(pks[[magnitude]]/pks$normalizations[[norm]])
}
###########
#Bad datasets
library(benford.analysis)

files <- list.files("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Bad Datasets",pattern = "\\.zip$",full.names = T)
pks_list<-lapply(files,rMSIproc::LoadPeakMatrix)
lapply(pks_list,function(x)generate_benfort_report(x,"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201020_11_BadResults/",images = F,magnitude="area",digits=3,percentage = 0.5))

a<-lapply(pks_list,function(x)lapply(c(0,0.2,0.5,0.8,1),function(y)benford(x$area[non_empty_pixels(x,percentage = y)],discrete = F,number.of.digits = 3)))

###########

#Absolute metric
#Dependene Problems of MAD
#1. Number of datapoints
#2. Number of digits chosen

#Explore the number of datapoints

digits<-1:3
pixels<-c(1:9%o%10^(1:5),10^6)

res<-lapply(pks_list,function(x)lapply(digits,function(d)lapply(pixels,function(p)if(p<length(x$area)) benford(x$area[sample(length(x$area),p)],number.of.digits = d,discrete=F))))
#res<-lapply(pks_list,function(x)lapply(digits,function(d)lapply(pixels,function(p)if(p<nrow(x$pos)) benford(x$area[sample(nrow(x$pos),p),],number.of.digits = d,discrete=F))))
MAD<-unlist(sapply(res,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$MAD))))
ideal_res<-lapply(digits,function(d)lapply(pixels,function(p)benford(BenfordTests::rbenf(p),number.of.digits = d,discrete=F)))
ideal_MAD<-unlist(sapply(ideal_res,function(x)sapply(x,function(y)y$MAD)))

p<-length(pixels)
d<-length(digits)
x<-length(pks_list)

df<-data.frame(MAD=c(ideal_MAD,MAD),dataset=c(rep("Ideal Benford",p*d),rep(sapply(pks_list,function(y)paste(y$names,collapse = "\n")),each=p*d)),pixels=rep(pixels,(x+1)*d),digits=rep(rep(digits,each=p),x+1))
df$normalized_MAD <- rep(NA,nrow(df))
df$normalized_MAD[which(!is.na(df$MAD))] <- abs(log10(df$MAD[!is.na(df$MAD)]))
   
ggplot(subset(df,digits==3),aes(x=pixels,y=MAD,group=dataset,color=dataset))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()

#Explore the empty pixel percentage
digits<-1:3
percentage<-c(0,1:9%o%10^c(-2,-1),1,1.5,2)

res<-lapply(pks_list,function(x)lapply(digits,function(d)lapply(percentage,function(p)benford(x$area[non_empty_pixels(x,"area",p),],number.of.digits = d,discrete=F)$MAD)))

#res<-lapply(pks_list,function(x)lapply(digits,function(d)lapply(pixels,function(p)if(p<nrow(x$pos)) benford(x$area[sample(nrow(x$pos),p),],number.of.digits = d,discrete=F))))
MAD<-unlist(sapply(res,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$MAD))))
ideal_res<-lapply(digits,function(d)lapply(pixels,function(p)benford(BenfordTests::rbenf(p),number.of.digits = d,discrete=F)))
ideal_MAD<-unlist(sapply(ideal_res,function(x)sapply(x,function(y)y$MAD)))
ideal_Chisq<-unlist(sapply(ideal_res,function(x)sapply(x,function(y)y$stats$chisq$parameter)))
ideal_pval<-unlist(sapply(ideal_res,function(x)sapply(x,function(y)y$stats$chisq$p.value)))

p<-length(percentage)
d<-length(digits)
x<-length(pks_list)

tmp_df<-data.frame(MAD=unlist(res),dataset=rep(sapply(pks_list,function(y)paste(y$names,collapse = "\n")),each=p*d),percentage=rep(percentage,(x)*d),digits=rep(rep(digits,each=p),x))



ggplot(subset(tmp_df,digits==1),aes(x=percentage,y=MAD,group=dataset,color=dataset))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()

#Global results

files <- list.files(c("/home/gbaquer/msidata/2. Benford Law/TOF Datasets","/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Bad Datasets"),pattern = "\\.zip$",full.names = T)
pks_list<-lapply(files,function(x)rMSIcleanup:::get_one_peakMatrix(rMSIproc::LoadPeakMatrix(x),1))

digits<-1:3
pixels<-c(1:9%o%10^(1:6),10^7)

res<-lapply(pks_list,function(x)lapply(digits,function(d)lapply(pixels,function(p)if(p<length(x$area)) benford(x$area[sample(length(x$area),p)],number.of.digits = d,discrete=F)[c("MAD","stats")])))
#res<-lapply(pks_list,function(x)lapply(digits,function(d)lapply(pixels,function(p)if(p<nrow(x$pos)) benford(x$area[sample(nrow(x$pos),p),],number.of.digits = d,discrete=F))))
MAD<-unlist(sapply(res,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$MAD))))
Chisq<-unlist(sapply(res,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$stats$chisq$parameter))))
pval<-unlist(sapply(res,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$stats$chisq$p.value))))
p<-length(pixels)
d<-length(digits)
x<-length(pks_list)

df_all<-data.frame(MAD=c(ideal_MAD,MAD),dataset=c(rep("Ideal Benford",p*d),rep(sapply(pks_list,function(y)paste(y$names,collapse = "\n")),each=p*d)),pixels=rep(pixels,(x+1)*d),digits=rep(rep(digits,each=p),x+1))
df_all$type<-rep("OK",nrow(df_all))
df_all$type[grepl("DHB",df_all$dataset)]<-"DHB_LUMC"
df_all$type[grepl("AuBSi|garlic",df_all$dataset)]<-"Alex"
ggplot(subset(df_all,digits==3),aes(x=pixels,y=MAD,group=dataset,color=dataset))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()

#After filtering
areas_list_filtered<-lapply(pks_list,function(x)x$area[non_empty_pixels(x,"area",0.1),])
res_filtered<-lapply(areas_list_filtered,function(x)lapply(digits,function(d)lapply(pixels,function(p)if(p<length(x)) benford(x[sample(length(x),p)],number.of.digits = d,discrete=F)[c("MAD","stats")])))
MAD<-unlist(sapply(res_filtered,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$MAD))))
MAD<-unlist(sapply(res_filtered,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$MAD))))
Chisq<-unlist(sapply(res_filtered,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$stats$chisq$parameter))))
pval<-unlist(sapply(res_filtered,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$stats$chisq$p.value))))

p<-length(pixels)
d<-length(digits)
x<-length(pks_list)

df_all_filtered<-data.frame(MAD=c(MAD),Chisq=c(Chisq),pval=c(pval),dataset=rep(sapply(pks_list,function(y)paste(y$names,collapse = "\n")),each=p*d),pixels=rep(pixels,(x)*d),digits=rep(rep(digits,each=p),x))
df_all_filtered$type<-rep("OK",nrow(df_all))
df_all_filtered$type[grepl("DHB",df_all$dataset)]<-"DHB_LUMC"
df_all_filtered$type[grepl("AuBSi|garlic",df_all$dataset)]<-"Alex"
ggplot(subset(df_all_filtered,digits==3),aes(x=pixels,y=MAD,group=dataset,color=type))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()

#After filterning normalization
digits<-1:3
pixels<-c(1:9%o%10^(1:6),10^7)
normalizations<-c("","TIC","RMS","AcqTic")

pks_list_norm <- pks_list[sapply(pks_list,function(x)"AcqTic" %in% names(x$normalizations))]
for(i in 1:length(pks_list_norm))
{
  seleceted_pixels<-non_empty_pixels(pks_list_norm[[i]],"area",0.1)
  pks_list_norm[[i]]$area<-pks_list_norm[[i]]$area[seleceted_pixels,]
  pks_list_norm[[i]]$intensity<-pks_list_norm[[i]]$intensity[seleceted_pixels,]
  pks_list_norm[[i]]$pos<-pks_list_norm[[i]]$pos[seleceted_pixels,]
  new_norm<-list()
  for(n in names(pks_list_norm[[i]]$normalizations)){
    new_norm[[n]]<-pks_list_norm[[i]]$normalizations[[n]][seleceted_pixels]
  }
  pks_list_norm[[i]]$normalizations<-new_norm
}

res_filtered_norm<-lapply(normalizations,function(n)lapply(pks_list_norm,function(x)lapply(digits,function(d)lapply(pixels,function(p)if(p<length(x$area)) benford(normalize(x,"area",n)[sample(length(x$area),p)],number.of.digits = d,discrete=F)[c("MAD","stats")]))))
MAD<-c(unlist(sapply(res_filtered_norm,function(x)sapply(x,function(y)sapply(y,function(z)sapply(z,function(k)if(is.null(k)) NA else k$MAD))))))
p<-length(pixels)
d<-length(digits)
x<-length(pks_list_norm)
n<-length(normalizations)

df_all_filtered_norm<-data.frame(MAD=MAD,norm=rep(normalizations,each=x*d*p),dataset=rep(rep(sapply(pks_list_norm,function(y)paste(y$names,collapse = "\n")),each=p*d),n),pixels=rep(rep(pixels,(x)*d),n),digits=rep(rep(rep(digits,each=p),x),n))
df_all_filtered_norm$type<-rep("OK",nrow(df_all_filtered_norm))
df_all_filtered_norm$type[grepl("DHB",df_all_filtered_norm$dataset)]<-"DHB_LUMC"
df_all_filtered_norm$type[grepl("AuBSi|garlic",df_all_filtered_norm$dataset)]<-"Alex"
df_all_filtered_norm$name=paste("[d]",df_all_filtered_norm$dataset,"[n]:",df_all_filtered_norm$norm)

for(i in seq_along(pks_list_norm))
{
  tiff(paste("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201021_12_Normalization/",pks_list_norm[[i]]$names[1],".tiff",sep=""))
  plot(ggplot(subset(df_all_filtered_norm,digits==3&dataset==levels(dataset)[i]),aes(x=pixels,y=MAD,group=name,color=norm))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10())
  dev.off()
}

#Plot REmoval percentage 2%
for(p in pks_list){
  tiff(paste("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/",p$names[1],".tiff",sep=""))
  rMSIproc::plotPeakImage(p,column=1)
  dev.off()
}
for(p in pks_list){
  tiff(paste("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/",p$names[1],"_MED20percent.tiff",sep=""))
  rMSIproc::plotValuesImage(p,non_empty_pixels(p,magnitude="area",percentage=0.2)+1)
  dev.off()
}
#Noise Experiment
pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/20200713_Au_Arterias.zip")
pks_13<-rMSIcleanup:::get_one_peakMatrix(pks,13)
pks_15<-rMSIcleanup:::get_one_peakMatrix(pks,15)

generate_benfort_report(add_noise(pks_13,level=0.1),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_6_Noise_All/",extra_name="_noise_0_1",images=F)
generate_benfort_report(add_noise(pks_13,level=1),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_6_Noise_All/",extra_name="_noise_1",images=F)
generate_benfort_report(add_noise(pks_13,level=5),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_6_Noise_All/",extra_name="_noise_5",images=F)


generate_benfort_report(add_noise(pks_15,level=0.1),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_noise_0_1",images=T)
generate_benfort_report(add_noise(pks_15,level=1),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_noise_1",images=T)
generate_benfort_report(add_noise(pks_15,level=5),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_noise_5",images=T)


generate_benfort_report(add_noise(pks_13,type="poisson",level=0.1),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_poisson_noise_0_1",images=T)
generate_benfort_report(add_noise(pks_13,type="poisson",level=1),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_poisson_noise_1",images=T)
generate_benfort_report(add_noise(pks_13,type="poisson",level=5),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_poisson_noise_5",images=T)


generate_benfort_report(add_noise(pks_15,type="poisson",level=0.1),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_poisson_noise_0_1",images=T)
generate_benfort_report(add_noise(pks_15,type="poisson",level=1),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_poisson_noise_1",images=T)
generate_benfort_report(add_noise(pks_15,type="poisson",level=5),"/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_Noise/",extra_name="_poisson_noise_5",images=T)


#Bin size experiments

img<-rMSI::LoadMsiData("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_56140149 002-proc.tar")

for(b in c(1,5,10,50)){
  tmp<-rMSIproc::ProcessImage(img,EnablePeakPicking = T,EnableAlignment = F,EnableSmoothing = F,EnableSpectraNormalization = F,EnableCalibration = F
                              ,BinTolerance = b,BinToleranceUsingPPM = F,SNR=5,UseBinning = T)
  rMSIproc::StorePeakMatrix(paste("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_BinTol_",b,".zip",sep=""),tmp$peakMat)
}

generate_benfort_report(rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_BinTol_1.zip"),path = "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_3_BinTol/",images = T,extra_name = "_BinTol_1")
generate_benfort_report(rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_BinTol_5.zip"),path = "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_3_BinTol/",images = T,extra_name = "_BinTol_5")
generate_benfort_report(rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_BinTol_10.zip"),path = "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_3_BinTol/",images = T,extra_name = "_BinTol_10")
generate_benfort_report(rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_BinTol_50.zip"),path = "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_3_BinTol/",images = T,extra_name = "_BinTol_50")

for(s in c(2,5,10,50)){
  tmp<-rMSIproc::ProcessImage(img,EnablePeakPicking = T,EnableAlignment = F,EnableSmoothing = F,EnableSpectraNormalization = F,EnableCalibration = F
                              ,BinTolerance = 5,BinToleranceUsingPPM = F,SNR=s,UseBinning = T)
  rMSIproc::StorePeakMatrix(paste("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_SNR_",s,".zip",sep=""),tmp$peakMat)
}

generate_benfort_report(rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_SNR_2.zip"),path = "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_4_SNR/",images = T,extra_name = "_SNR_2")
generate_benfort_report(rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_SNR_5.zip"),path = "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_4_SNR/",images = T,extra_name = "_SNR_5")
generate_benfort_report(rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_SNR_10.zip"),path = "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_4_SNR/",images = T,extra_name = "_SNR_10")
generate_benfort_report(rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_SNR_50.zip"),path = "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_4_SNR/",images = T,extra_name = "_SNR_50")

sweep_preprocessing<-function(img,path,name=get_name(img$name),BinTol=6,SNR=5){
  for(bt in BinTol){
    for(snr in SNR){
      tmp<-rMSIproc::ProcessImage(img,EnablePeakPicking = T,EnableAlignment = F,EnableSmoothing = F,EnableSpectraNormalization = F,EnableCalibration = F
                                  ,BinTolerance = bt,BinToleranceUsingPPM = F,SNR=snr,UseBinning = T)
      rMSIproc::StorePeakMatrix(paste(path,name,"_BinTolerance",bt,"_SNR",snr,".zip",sep=""),tmp$peakMat)
    }
  }
}
get_name<-function(str){
  strsplit(str,".",fixed=T)[[1]][1]
}

path="/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/"
sweep_preprocessing(rMSI::LoadMsiData("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END28_56150310_002-proc.tar"),path=path,SNR=c(2,5,10,50))
sweep_preprocessing(rMSI::LoadMsiData("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200717_END46_56160849_002-proc.tar"),path=path,SNR=c(2,5,10,50))

sweep_preprocessing(rMSI::LoadMsiData("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200714_END9_56140350_002-proc.tar"),path=path,SNR=c(2,5,10,50))
sweep_preprocessing(rMSI::LoadMsiData("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END25_56150211_002-proc.tar"),path=path,SNR=c(2,5,10,50))
sweep_preprocessing(rMSI::LoadMsiData("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END31_56160011_002-proc.tar"),path=path,SNR=c(2,5,10,50))

#Graphs
df <- data.frame(stringsAsFactors = F, files = c("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END28_56150310_002-proc_BinTolerance6_SNR2.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END28_56150310_002-proc_BinTolerance6_SNR5.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END28_56150310_002-proc_BinTolerance6_SNR10.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END28_56150310_002-proc_BinTolerance6_SNR50.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200717_END46_56160849_002-proc_BinTolerance6_SNR2.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200717_END46_56160849_002-proc_BinTolerance6_SNR5.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200717_END46_56160849_002-proc_BinTolerance6_SNR10.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200717_END46_56160849_002-proc_BinTolerance6_SNR50.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_56140149 002-proc_BinTolerance6_SNR2.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_56140149 002-proc_BinTolerance6_SNR5.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_56140149 002-proc_BinTolerance6_SNR10.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200713_END2_56140149 002-proc_BinTolerance6_SNR50.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200714_END9_56140350_002-proc_BinTolerance6_SNR2.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200714_END9_56140350_002-proc_BinTolerance6_SNR5.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200714_END9_56140350_002-proc_BinTolerance6_SNR10.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200714_END9_56140350_002-proc_BinTolerance6_SNR50.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END31_56160011_002-proc_BinTolerance6_SNR2.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END31_56160011_002-proc_BinTolerance6_SNR5.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END31_56160011_002-proc_BinTolerance6_SNR10.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END31_56160011_002-proc_BinTolerance6_SNR50.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END25_56150211_002-proc_BinTolerance6_SNR2.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END25_56150211_002-proc_BinTolerance6_SNR5.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END25_56150211_002-proc_BinTolerance6_SNR10.zip",
                                                 "/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Preprocessing Experiment/20200715_END25_56150211_002-proc_BinTolerance6_SNR50.zip"))

df$dataset<-rep(c("0_END28","0_END46","0_END2","1_END9","1_END31","1_END25"),each=4)
df$variable<-rep(c(2,5,10,50),6)
pks_list <- lapply(df$files,rMSIproc::LoadPeakMatrix)

mean_spectra <- lapply(pks_list,function(x)(x$area)*1000)
benford_results <- lapply(mean_spectra,function(x)benford.analysis::benford(x,number.of.digits = 3,discrete = F))

df$MADs <- unlist(lapply(benford_results,function(x)x$MAD))

lapply(benford_results,plot)
as_pks_list<-function(pks){
  return(lapply(1:length(pks$numPixels),function(i)rMSIcleanup:::get_one_peakMatrix(pks,i)))
}
generate_MAD_plot_ordered <- function(benford_results,sort=T,dataset=1:length(benford_results),variable=rep(1,length(benford_results))){
  df<-data.frame(dataset=dataset,variable=variable)
  df$MAD <- unlist(lapply(benford_results,function(x)x$MAD))
  ggplot(data=df, aes(x=dataset, y=MAD, group=variable), mapping = aes(x = reorder(dataset, MAD), MAD)) +
    geom_line()+
    geom_point()
}
generate_MAD_plot <- function(benford_results,sort=T,dataset=1:length(benford_results),variable=rep(1,length(benford_results))){
  df<-data.frame(dataset=dataset,variable=as.factor(variable))
  df$MAD <- unlist(lapply(benford_results,function(x)x$MAD))
  ggplot(data=df, aes(x=dataset, y=MAD, group=variable,col=variable)) +
    geom_line()+
    geom_point()
}
compute_mean_benford <- function(pks_list,magnitude="intensity"){
  mean_spectra <- lapply(pks_list,function(x)(x[[magnitude]])*1000)
  benford_results <- lapply(mean_spectra,function(x)benford.analysis::benford(x,number.of.digits = 1,discrete = F))
  return(benford_results)
}

generate_MAD_violin<-function(benford_results,sort=T,dataset=1:length(benford_results),variable=rep(1,length(benford_results))){
  df<-data.frame(dataset=dataset,variable=as.factor(variable))
  df$MAD <- unlist(lapply(benford_results,function(x)x$MAD))
  ggplot(df, aes(x=variable, y=MAD)) + 
    geom_violin()
}

pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/20200713_Au_Arterias.zip")
#Store MAD ARterias
tiff("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_1_General_Arterias/20200713_END1_56140138002_4_benford_MAD_ALL.tiff",width=1500,height=500)
tmp<-compute_mean_benford(as_pks_list(pks))
generate_MAD_plot_ordered(tmp,dataset = unlist(lapply(strsplit(pks$names,"_"),function(x)x[2])))
dev.off()

#Store MAD SNR
tiff("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_4_SNR/SNR_sweep_area.tiff")
generate_MAD_plot(compute_mean_benford(pks_list,"area"),dataset=df$dataset,variable=df$variable)
dev.off()

tiff("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201014_4_SNR/SNR_sweep_intensity.tiff")
generate_MAD_plot(compute_mean_benford(pks_list,"intensity"),dataset=df$dataset,variable=df$variable)
dev.off()

#Changing Laser Power
bd<-"/home/gbaquer/msidata/1. In-source Fragmentation/1.3. AuBSi Wells/20200701_AuBSi_wells_std/20200701_AuBSi_wells_std/"
pks<-rMSIproc::LoadPeakMatrix(paste(bd,"20200701_AuBSi_wells_std-peaks.zip",sep=""))
pks_list<-lapply(1:3,function(i)rMSIcleanup:::get_one_peakMatrix(pks,i))
xml_files<-paste(bd,c("20200701_AuBSi_wells_std_minEnergy_Regions.xml","20200701_AuBSi_wells_std_medEnergy_Regions.xml","20200701_AuBSi_wells_std_maxEnergy_Regions.xml"),sep="")
ROIs_list <-lapply(seq_along(xml_files),function(i)rMSI::ReadBrukerRoiXML(pks_list[[i]],xml_files[i]))
ROIs<-list()


for(ROI in ROIs_list)
  ROIs<-append(ROIs,ROI)

metadata<-data.frame(name=unlist(lapply(ROIs,function(x)x$name)))

#Correct a mistake 
metadata$name[which(metadata$name=="310")]="Pool_5C"

#Process the rest of the metadata
metadata$type<-as.factor(unlist(lapply(metadata$name,function(x)strsplit(as.character(x),"_")[[1]][1])))
metadata$concentration<-as.numeric(gsub("[^0-9]","", metadata$name))
metadata$replicate<-as.factor(gsub("[!0-9]","", unlist(lapply(metadata$name,function(x)strsplit(as.character(x),"_")[[1]][2]))))
metadata$power<-rep(c("MIN","MED","MAX"),each=41)


#Compute benford
digits<-1:3
pixels<-c(1:9%o%10^(1:3),10^4,unique(sapply(ROIs,function(x)length(x$id)*length(pks$mass))))

areas_list <- lapply(ROIs,function(x)(pks$area[x$id,]))

res<-lapply(areas_list,function(x)lapply(digits,function(d)lapply(pixels,function(p)if(p<length(x)) benford(x[sample(length(x),p)],number.of.digits = d,discrete=F)[c("MAD","stats")])))
MAD<-unlist(sapply(res,function(x)sapply(x,function(y)sapply(y,function(z)if(is.null(z)) NA else z$MAD))))

p<-length(pixels)
d<-length(digits)
x<-length(ROIs)

df<-data.frame(MAD=c(MAD),name=rep(metadata$name,each=p*d),type=rep(metadata$type,each=p*d),replicate=rep(metadata$replicate,each=p*d),power=rep(metadata$power,each=p*d),concentration=rep(metadata$concentration,each=p*d),pixels=rep(pixels,(x)*d),digits=rep(rep(digits,each=p),x))
df$dataset<-paste(df$name,df$power,"_")

ggplot(subset(df,digits==1&power=="MED"&replicate=="B"),aes(x=pixels,y=MAD,group=dataset,color=name))+geom_line()+geom_point()+scale_x_log10()+scale_y_log10()


#New Ideas

#Compare the conformity of all normalizations [DONE]
b<-lapply(cbind(1,pks_13$normalizations),function(x)benford.analysis::benford(pks_13$intensity/x,discrete=F))


#Another application
#How to determine the TIC threshold to remove a empty pixels [DOUBT]
plot(benford.analysis::benford(pks_15$area[-which(apply(pks_15$intensity,1,sum)<0.03*max(apply(pks_15$intensity,1,sum))),]))
b<-lapply(c(0.03,0.05,0.1,0.2,0.3,0.5,0.8),function(x)benford.analysis::benford(pks_15$area[-which(apply(pks_15$intensity,1,sum)<x*max(apply(pks_15$intensity,1,sum))),]))

#Another application
#Sort MAD of ions. What differences are there between the most conforming and the least conforming? [DONE]



#EXoplore effects of sample size in p-value of the chi-squared
a<-lapply(c(1:9 %o% 10^(1:4)),function(x)BenfordTests::chisq.benftest(pks_13$intensity[sample(x)],digits=3))
plot(c(1:9 %o% 10^(1:4)),unlist(lapply(a,function(x)x$p.value)))

BenfordTests::chisq.benftest(pks_13$intensity[sample(length(pks_13$intensity),1000)],digits=2)

a<-lapply(c(1:9 %o% 10^(1:4)),function(x)benford.analysis::benford(pks_13$intensity[sample(length(pks_13$intensity),x)],number.of.digits =3))

#Splitting between good and bad

#Compute ion_MAD
#Assess if gaussian

#If not gaussian split with kmeans
tmp<-kmeans(ion_MAD,2)
plot(density(ion_MAD,bw=0.0001))
lines(density(ion_MAD[tmp$cluster==1],bw=0.0001))
lines(density(ion_MAD[tmp$cluster==2],bw=0.0001))


#check correlation

pixel_benford<-apply(pks_13$area,1,function(x)benford.analysis::benford(x,discrete=F,number.of.digits = 1))
ion_benford<-apply(pks_13$area,2,function(x)benford.analysis::benford(x,discrete=F))


pixel_conformity<-unlist(lapply(pixel_benford,function(x)x$MAD.conformity))

n_pixels<-5
grouping<-rep(1:ceiling(nrow(pks_13$area)/n_pixels),each=n_pixels)[1:nrow(pks_13$area)]

grouping_pixel_benford<-lapply(1:max(grouping),function(x)benford.analysis::benford(pks_13$area[which(grouping==x),],discrete=F,number.of.digits = 1))
grouping_pixel_conformity<-unlist(lapply(grouping_pixel_benford,function(x)x$MAD.conformity))



########
#Ionwise study

pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/20200713_Au_Arterias.zip")
pks_list<-lapply(seq_along(pks$names),function(i)rMSIcleanup:::get_one_peakMatrix(pks,i))

ion_results<-lapply(pks_list,function(x)ion_benford(x,magnitude = "area",digits=2,percentage = 0.1))
ion_TIC<-sapply(pks_list,function(x)apply(x$area,2,sum))
ion_MAD<-sapply(ion_results,function(x)sapply(x,function(y)y$MAD))

m<-length(pks$mass)
d<-length(pks$names)
df<-data.frame(MAD=c(ion_MAD),TIC=c(ion_TIC),dataset=rep(pks$names,each=m),mz=rep(pks$mass,d),i=rep(seq_along(pks$mass),d),stringsAsFactors = T)
df$name=sapply(strsplit(as.character(df$dataset),"_"),function(x)x[2])
df$log2TIC<-log2(df$TIC)

ggplot(subset(df,dataset==levels(dataset)[10]),aes(x=TIC,y=MAD))+geom_point()+scale_x_log10()
for(i in 1:6){
  tiff(paste("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/20201021_13_IonWise/Group_",i,".tiff",sep=""))
  p1<-ggplot(subset(df,dataset==levels(dataset)[i*3+1:3]),aes(x=MAD,group=name,color=name))+geom_density()
  p2<-ggplot(subset(df,dataset==levels(dataset)[i*3+1:3]),aes(x=TIC,group=name,color=name))+geom_density()+scale_x_log10()
  plot(grid.arrange(p1,p2))
  dev.off()
}

worst_ions<-sort(table(apply(ion_MAD,2,function(x)order(x,decreasing =T)[1:10])),decreasing = T)
best_ions<-sort(table(apply(ion_MAD,2,function(x)order(x,decreasing =F)[1:10])),decreasing = T)

df$type<-rep("NORMAL",nrow(df))
df$type[df$i%in%as.numeric(names(worst_ions))]<-"BAD"
df$type[df$i%in%as.numeric(names(best_ions))]<-"GOOD"

ggplot(subset(df,i%in%as.numeric(c(names(best_ions)[1:3],names(worst_ions)[1:3]))),aes(x=name,y=MAD,group=i,color=type))+geom_point()+geom_line()
ggplot(subset(df,i%in%as.numeric(c(names(best_ions)[1:3],names(worst_ions)[1:3]))),aes(x=name,y=TIC,group=i,color=type))+geom_point()+geom_line()
ggplot(df[sample(nrow(df),nrow(df)),],aes(x=TIC,y=MAD,color=type))+geom_point()+scale_x_log10()
ggplot(subset(df[sample(nrow(df),nrow(df)),],dataset==levels(dataset)[3]),aes(x=TIC,y=MAD,color=type))+geom_point()+scale_x_log10()
  
i=3     
tiff(paste("/home/gbaquer/msidata/2. Benford Law/TOF Datasets/Results/BEST_",i,"_in_",best_ions[i],"_datasets.tiff",sep=""),width=1000,height=1500)
rMSIproc::plotPeakImage(pks,matrix="area",ncol=4,column=as.numeric(names(best_ions[i])),labels=unique(df$name))
dev.off()

ggplot(subset(df,i==as.numeric(names(best_ions[5]))),aes(x=TIC,y=MAD,group=i))+geom_point()

#PCA
bfd_1<-(sapply(ion_results[[7]],function(x)x$bfd$data.dist))
pca<-prcomp(bfd_1)
df<-data.frame(PC=pca$x)
ggplot(df,aes(x=PC.PC1,y=PC.PC2))+geom_point()+
  xlab(paste("PC1",round(100*summary(pca)$importance[2,1],2),"%"))+
  ylab((paste("PC2",round(100*summary(pca)$importance[2,2],2),"%")))
