# LIPID FRAGMENTATION IDENTIFICATION WORKFLOW
# Gerard Baquer Gómez
# 29/09/20

# 0. LITERATURE DATA
RefMet<-read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/refmet.csv")
load("/home/gbaquer/msidata/sp_MS2ID_20191002_125914.RData")
pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.2. Standard Study /Munster_Au_Res_140k_10um.zip")

files_neg<-list.files(path = "/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/imzml",pattern=".zip",full.names = T)
pks_neg<-lapply(files_neg,rMSIproc::LoadPeakMatrix)

results_neg<-lapply(pks_neg,function(x)identify_fragments(x,20,"neg",F,F))
T1<-do.call(rbind,lapply(results_neg,function(x)cbind(x$T1,data.frame(dataset=x$parameters$names[1]))))
T2<-do.call(rbind,lapply(results_neg,function(x)cbind(x$T2,data.frame(dataset=x$parameters$names[1]))))
T1$complete_id<-paste(T1$id,T1$dataset,sep="_")
T2$complete_parental_id<-paste(T2$parental_id,T2$dataset,sep="_")
T2$complete_fragment_id<-paste(T2$fragment_id,T2$dataset,sep="_")
T2_filters<-filter_fragments(T2,mode="neg")

MS2ID_lipids<-do.call(rbind,lapply(1:nrow(lipids),function(i)generate_lipid_list(lipids$name[i],lipids$synonyms[i])))
#T1$dataset<-rep(sapply(pks_neg,function(x)x$names[1]),each=sapply(results_neg,function(x)nrow(x$T1)))
#T2$dataset<-rep(sapply(pks_neg,function(x)x$names[1]),each=sapply(results_neg,function(x)nrow(x$T2)))

lipids<-data.frame(name=c("DG","TG","PA","PE","PC","PG","PI","PS","CL","Cer","CerP","SM","HexCer","SFT","GM3","FRAG"),stringsAsFactors = F)
lipids$standard <- c("DG 18:1","TG 18:1","PA 18:1", "PE 16:0/18:1 & PE P-18:1","PC 16:0/18:1 & PC P-18:1","PG 18:0","PI 16:0/18:1","PS 16:0","CL 16:1","Cer d18:/18:0","CerP d18:1/16:0","SM d18:1/12:0","GlcCer d18:1/18:0","SFT d18:1/12:0","Ganglioside GM3")
lipids$synonyms <-c("Dag|Diacylglycerol","Tag|Triacylglycerol","Phosphatidic acid|Phosphatidate",
                    "Phophatidylethanolamine|GPEtn","Phosphatidylcholine|GPCho","Phosphatidylglycerol|GPG",
                    "Phosphatidylinositol|Pino","Phosphatidylserine|Pser","Cardiolipin","DHCer|Ceramide",
                    "Ceramide phosphate","Sphingomyelin","","","","")
lipid_dict<-as.list(lipids$name)
names(lipid_dict)<-lipids$name
lipid_dict$PEtOH <- "PE"
lipid_dict$DGTS <- "DG"
lipid_dict$SHexCer <- "HexCer"
lipid_dict$exCer<-"HexCer"

adducts_pos<-data.frame(name=c("M+Na-COOH","M-H2O+H","M+H","M-H2O+Na","M-H20+K","M+Na","M+K","M-H+2Na","M-H+Na+K","M-2H+3Na","M-H+2K","M-2H+2Na+K","M-2H+Na+2K","M-2H+3K","M+"),stringsAsFactors = F)
adducts_pos$mass<-c(-21.0006,-14.9876,1.0073,6.9943,22.9682,22.9892,38.9632,44.9712,60.9451,66.9531,76.9190,82.9271,98.9010,114.8749,0)
adducts_pos$mode<-rep("pos",nrow(adducts_pos))
adducts_pos$DG<-c(0,1,0,1,1,3,1,0,0,0,0,0,0,0,1)
adducts_pos$TG<-c(0,0,0,0,0,3,1,0,0,0,0,0,0,0,1)
adducts_pos$PA<-c(0,0,1,0,0,1,1,3,1,1,1,1,1,1,1)
adducts_pos$PE<-c(0,0,1,0,0,2,1,3,1,0,1,0,0,0,1)
adducts_pos$PC<-c(0,0,2,0,0,3,1,0,0,0,0,0,0,0,1)
adducts_pos$PG<-c(0,0,0,0,0,1,1,3,1,0,1,0,0,0,1)
adducts_pos$PI<-c(0,0,0,0,0,1,1,3,1,0,1,0,0,0,1)
adducts_pos$PS<-c(0,0,0,0,0,1,1,1,1,3,1,1,1,1,1)
adducts_pos$CL<-c(0,0,0,0,0,1,1,2,1,3,1,1,1,1,1)
adducts_pos$Cer<-c(0,2,1,1,1,3,1,0,0,0,0,0,0,0,1)
adducts_pos$CerP<-c(0,0,0,0,0,1,1,3,1,0,1,0,0,0,1)
adducts_pos$SM<-c(0,0,2,0,0,3,1,0,0,0,0,0,0,0,1)
adducts_pos$HexCer<-c(0,0,0,0,0,3,2,0,0,0,0,0,0,0,1)
adducts_pos$SFT<-c(0,0,0,0,0,0,0,3,0,0,1,0,0,0,1)
adducts_pos$GM3<-c(1,0,0,0,0,1,1,3,1,0,1,0,0,0,1)
adducts_pos$FRAG<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)

adducts_neg<-data.frame(name=c("M-2(H2O)","M-oh","M-Ch3","M-H","M+Na-2h","M+K-2h","M+Dan-H","M-"),stringsAsFactors = F)
adducts_neg$mass<-c(-36.0212,-17.0027,-15.0229,-1.0073,20.9747,36.9486,157.0771,0)
adducts_neg$mode<-rep("neg",nrow(adducts_neg))
adducts_neg$DG<-c(0,0,0,0,0,0,0,1)
adducts_neg$TG<-c(0,0,0,0,0,0,0,1)
adducts_neg$PA<-c(0,0,0,3,0,0,1,1)
adducts_neg$PE<-c(0,0,0,3,0,0,0,1)
adducts_neg$PC<-c(0,0,3,0,0,0,0,1)
adducts_neg$PG<-c(0,0,0,3,0,0,0,1)
adducts_neg$PI<-c(0,0,0,3,0,0,0,1)
adducts_neg$PS<-c(0,0,0,3,1,1,0,1)
adducts_neg$CL<-c(0,0,0,1,3,1,0,1)
adducts_neg$Cer<-c(0,0,0,3,0,0,0,1)
adducts_neg$CerP<-c(0,0,0,3,0,0,1,1)
adducts_neg$SM<-c(0,0,3,0,0,0,0,1)
adducts_neg$HexCer<-c(0,0,0,3,0,0,0,1)
adducts_neg$SFT<-c(1,1,0,3,0,0,0,1)
adducts_neg$GM3<-c(0,0,0,3,0,0,0,1)
adducts_neg$FRAG<-c(0,0,0,0,0,0,0,1)
  
adducts<-rbind(adducts_pos,adducts_neg)
rownames(adducts)<-adducts$name

adducts_fragment<-adducts
adducts_fragment[which(adducts$mode=="pos"),c('DG','Cer')]<-0
adducts_fragment[which(adducts$mode=="pos"&adducts$name=="M-H20+H"),c('DG','Cer')]<-1

fragmentation<-data.frame(matrix("",ncol=nrow(lipids),nrow=nrow(lipids)),stringsAsFactors = F)
colnames(fragmentation)<-lipids$name
rownames(fragmentation)<-lipids$name

fragmentation_pos<-fragmentation
fragmentation_neg<-fragmentation
rm(fragmentation)

fragmentation_pos['PC','PA']<-"-N(CH3)3-Cho"
fragmentation_pos['PE','PA']<-"-EtnA"
fragmentation_pos['PI','PA']<-"-Inositol"
fragmentation_pos['PS','PA']<-"-Serine"

fragmentation_pos['PC','DG']<-"-P_Cho"
fragmentation_pos['PE','DG']<-"-P_EtnA"
fragmentation_pos['PI','DG']<-"-P_Inositol"

fragmentation_pos['SM','CerP']<-"-N(CH3)3-Cho"
fragmentation_pos['SM','Cer']<-"-P_Cho"

fragmentation_pos['GM3','FRAG']<-"-Sialic_Acid"

fragmentation_pos['CL','FRAG']<-"-DG+H2O | -DG | -RCOOH"

fragmentation_neg['PC','PA']<-"-CH2N(CH3)3 | -N(CH3)3 | -Cho"
fragmentation_neg['PS','PA']<-"-Serine"

fragmentation_neg['PI','FRAG']<-"-C4H10O5"

fragmentation_neg['SM','CerP']<-"-CH2N(CH3)3 | -N(CH3)3 | -Cho"

fragmentation_neg['GM3','Cer']<-"-GM3HG"
fragmentation_neg['GM3','HexCer']<-"-Sialic_Acid-Hex"
fragmentation_neg['GM3','FRAG']<-"-Sialic_Acid"
fragmentation_neg['GM3','Cer']<-"-Sialic_Acid-Lac"

fragmentation_neg['HexCer','Cer']<-"-Hex"

fragmentation_neg['SFT','FRAG']<-"-OH | 2H2O"

fragmentation_neg['CL','FRAG']<-"-RCOOH-RCOOH | -RCOOH-RCOH | -RCOOH | -RCOH | -DG+H2O | -DG"


losses<-NULL

losses['P_Cho']<-183.0660444
losses['P_EtnA']<-141.0190942
losses['P_Inositol']<-259.0218935
losses['N(CH3)3']<-59.0734993
losses['Cho']<-85.08914936
losses['EtnA']<-43.04219917
losses['Inositol']<-162.0528234
losses['Serine']<-87.03257699
losses['Sialic_Acid']<-291.0954165
#losses['CH2N(CH3)3']<-73.08914936
losses['CH2N(CH3)3']<-71.0734993 #Modified after meeting with Lucía. Itshould be this CN(CH3)3
losses['Hex']<-162.0528234
losses['Lac']<-324.1056468
losses['OH']<-17.00273965
losses['2H2O']<-36.02112937
losses['C4H10O5']<-138.0528234
losses['GM3HG']<-615.2010634

# 1. HELPER FUNCTIONS
parse_fragmentation<- function(s)
{
  s<- gsub(" ", "", s, fixed = T)
  n<-unlist(strsplit(s,"|",fixed=T))
  v<-sapply(strsplit(n,"-",fixed = T),function(x)(-1)*sum(losses[x],na.rm = T))
  names(v)<-n
  return(v)
}

find_c_delta <- function(n)
{
  pattern<-"([0-9]+:[a-z]?[0-9]+)"
  matches<-gsub("[a-z]","",unlist(regmatches(n, gregexpr(pattern, n))))
  if(length(matches)==0)
    result=c(NA,NA)
  else
    result<-apply(matrix(as.numeric(unlist(strsplit(matches,":"))),ncol=length(matches)),1,sum)
  return(result)
}
find_lipid <- function(n,s)
{
  pattern<-paste("(",s,")+",sep="") #pattern<-paste("[",s,"]\\w+",sep="")
  return(unlist(regmatches(n, gregexpr(pattern, n)))[1])
}
find_lipid_c_delta <- function(n,s=paste(lipids$name,collapse ="|")){
  return(paste(find_lipid(n,s),paste(find_c_delta(n),collapse = ":")))
}
find_adduct <- function(n)
{
    pattern<-"([+|-][^0-9][^]]+)"
    return(gsub(" ","",paste("M",unlist(regmatches(n, gregexpr(pattern, n)))[1],sep="")))
}
get_tol<-function(m1,m2){
  return(abs((m1-m2)/m2)*10^6)
}
withintol<-function(m1,m2,tol){
  return(get_tol(m1,m2)<tol)
}

searchMS2ID<-function(mz,m_adduct,tol=5,substring=""){
  matches<-subset(subset(df_metametabolits,withintol(monoisotopic_molecular_weight,mz-m_adduct,tol)),grepl(substring,nommetabolit))
  unique_result=NULL
  if(nrow(matches)>0){
    c_delta<-matrix(unlist(lapply(matches$nommetabolit,find_c_delta)),nrow=2)
    
    result<-data.frame(metabolite_id=matches$idmetabolit,mz=matches$monoisotopic_molecular_weight,stringsAsFactors = F)
    result$original_lipid<-unlist(lapply(matches$nommetabolit,function(n)find_lipid(n,substring)))
    result$lipid<-unlist(lapply(result$original_lipid,function(x)if(is.null(unlist(lipid_dict[x])))"FRAG" else unlist(lipid_dict[x])))
    result$c<-c_delta[1,]
    result$delta<-c_delta[2,]
    result$ppm_error<-get_tol(result$mz,(mz-m_adduct))
    
    unique_result<-result[!duplicated(result[,c('lipid','c','delta')]),]
  }
  return(unique_result)
}


#Import export functions
load_gt<-function(path="/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_neg.csv"){
  table<-read.csv(path,sep=";")
    #lipid
  table$lipid<-sapply(table$assignement,function(s)find_lipid(s,paste(names(adducts)[-(1:3)],collapse="|")))
  #C_delta
  c_delta<-sapply(table$assignement,find_c_delta)
  table$c<-c_delta[1,]
  table$delta<-c_delta[2,]
  #adduct
  table$adduct<-sapply(table$assignement,find_adduct)
  return(table)
}

# 2. IDENTIFICATION FUNCTIONS
identify_fragments<-function(pks,tol=5, m="pos",filter_monoisotopic=T,perform_T2=T){
  #Get Monoisotopic
  if(filter_monoisotopic){
    monoisotopic_i<-(rMSIproc:::isotopeAnnotation(pks,tolerance = tol))$monoisotopicPeaks
  }else{
    monoisotopic_i<-1:length(pks$mass)
  }
    
  a<-subset(adducts,mode==m)
  if(m=="pos"){
    fragmentation<-fragmentation_pos
  }else{
    fragmentation<-fragmentation_neg
  }
    
  #Perform all DB searches [It can be optimized]
 
  T1_db_matches<- data.frame(NULL,stringsAsFactors = F)
  for(i in monoisotopic_i){
    for(j in 1:nrow(a)){
      #before: (l in colnames(a)[which(a[j,-(1:3)]!=0)+3])
      for(l in colnames(a)[-(1:3)]){
        fs=c(0,parse_fragmentation(fragmentation[l,]))
        for(k in seq_along(fs))
        {
          r<-searchMS2ID(pks$mass[i],a$mass[j]+as.numeric(fs[k]),tol,substring=l)
          if(!is.null(r)){
            r$i<-i
            r$adduct<-a$name[j]
            r$experimental_mz<-pks$mass[i]
            r$fragmentation<-names(fs)[k]
          }
          T1_db_matches<-rbind(T1_db_matches,r)
        }
      }
    }
  }
  
  # T1_db_matches<- data.frame(NULL,stringsAsFactors = F)
  # for(i in monoisotopic_i){
  #   for(j in 1:nrow(a)){
  #     r<-searchMS2ID(pks$mass[i],a$mass[j],tol,substring=paste(colnames(a)[which(a[j,-(1:3)]!=0)+3],collapse = "|"))
  #     if(!is.null(r)){
  #       r$i<-i
  #       r$adduct<-a$name[j]
  #       r$experimental_mz<-pks$mass[i]
  #     }
  #     T1_db_matches<-rbind(T1_db_matches,r)
  #   }
  # }
  
  T1_db_matches<-T1_db_matches[order(-T1_db_matches$mz),]
  T1_db_matches$id<-1:nrow(T1_db_matches)
  
  #Update adduct likelyhood and times found
  T1_db_matches$adduct_likelyhood<-unlist(lapply(1:nrow(T1_db_matches),function(i){a[T1_db_matches[i,]$adduct,T1_db_matches[i,]$lipid]}))
  T1_db_matches$times_found<-0
  times_found<-aggregate(list(numdup=rep(1,nrow(T1_db_matches))), T1_db_matches[c('lipid','c','delta')], length)
  for(i in 1:nrow(times_found)){
    tf<-times_found[i,]
    T1_db_matches[rownames(subset(T1_db_matches,lipid==tf$lipid&c==tf$c&delta==tf$delta)),]$times_found<-tf$numdup
  }

  #Perform all correlation checks
  T2_possible_fragments<- data.frame(NULL,stringsAsFactors = F)
  if(perform_T2)
  {
    c<-cor(pks$intensity[monoisotopic_i,monoisotopic_i])
    
    
    for(l in lipids$name){
      parentals<-subset(T1_db_matches,lipid==l)
      fragment_lipids<-names(fragmentation)[which(fragmentation[l,]!="")]
      fragments<-subset(T1_db_matches,lipid%in%fragment_lipids)
      if(nrow(parentals)!=0&nrow(fragments)!=0){
        for(i in 1:nrow(parentals)){
          parental=parentals[i,]
          for(j in 1:nrow(fragments)){
            fragment=fragments[j,]
            if(parental$mz>fragment$mz){
              f<-data.frame(parental_id=parental$id,parental_experimental_mz=parental$experimental_mz,parental_lipid=parental$lipid,parental_c=parental$c,parental_delta=parental$delta,
                            fragment_id=fragment$id,fragment_experimental_mz=fragment$experimental_mz,fragment_lipid=fragment$lipid,fragment_c=fragment$c,fragment_delta=fragment$delta,
                            cor=c[match(parental$i,monoisotopic_i),match(fragment$i,monoisotopic_i)],diff=parental$mz-fragment$mz, combined_adduct_likelyhood=NA,
                            combined_ppm_error=parental$ppm_error+fragment$ppm_error,adduct_match=(parental$adduct==fragment$adduct),stringsAsFactors = F)
              T2_possible_fragments<-rbind(T2_possible_fragments,f)
            }
          }
        }
      }
      
    }
    #Update putative fragments and putative parentals
    
    T1_db_matches$putative_fragments<-unlist(lapply(T1_db_matches$id,function(x)sum(T2_possible_fragments$parental_id==x)))
    T1_db_matches$putative_parents<-unlist(lapply(T1_db_matches$id,function(x)sum(T2_possible_fragments$fragment_id==x)))
    
    
  }
  
  
  
  

  return(list(T1=T1_db_matches,T2=T2_possible_fragments, parameters=list(names=pks$names,tol=tol,mode=mode,filter_monoisotopic=filter_monoisotopic,monoisotopic_i=monoisotopic_i)))
}

#It probably should be included in the previous function

filter_fragments<-function(T2,tol=5,mode="pos")
{
  if(mode=="pos")
    fragmentation<-fragmentation_pos
  else
    fragmentation<-fragmentation_neg
  
  f<-NULL
  for(i in 1:nrow(T2)){
    if(any(withintol(T2[i,]$diff,abs(parse_fragmentation(fragmentation[T2[i,]$parental_lipid,T2[i,]$fragment_lipid])),tol)))
      f<-c(f,i)
  }
    
  return(T2[f,])
}

#Generate T3
generate_lipid_list<-function(lipid,synonyms="")
{
  if(synonyms!="")
    matches<-subset(df_metametabolits,grepl(lipid,nommetabolit)|(grepl(tolower(synonyms),tolower(nommetabolit))))
  else
    matches<-subset(df_metametabolits,grepl(lipid,nommetabolit))
  
  l<-data.frame(t(sapply(matches$nommetabolit,find_c_delta)))
  if(ncol(l)==2){
    names(l)<-c("c","delta")
    l$lipid<-lipid
    l$mz<-matches$monoisotopic_molecular_weight
  }
  else{
    l=NULL
  }
  return(unique(l))
  
}
generate_T3<-function(T1,T2)
{
  
  #Match paper mz
  MS2ID_lipids<-do.call(rbind,lapply(lipids$name,generate_lipid_list))
  
  tols<-outer(T1$experimental_mz,table$mz,get_tol)
  closest<-apply(tols,1,which.min)
  T1$paper_mz<-sapply(seq_along(closest),function(i) if(tols[i,closest[i]]>5) NA else table$mz[closest[i]])
  
  tols<-outer(T2_filters$fragment_experimental_mz....fragment.experimental_mz,table$mz,get_tol)
  closest<-apply(tols,1,which.min)
  T2_filters$paper_fragment_mz<-sapply(seq_along(closest),function(i) if(tols[i,closest[i]]>5) NA else table$mz[closest[i]])
  
  table$matched_mz<-table$mz%in%T1$paper_mz
  
  table$matched_id<-paste(table$mz,table$lipid,table$c,table$delta,tolower(table$adduct))%in%paste(T1$paper_mz,T1$lipid,T1$c,T1$delta,tolower(T1$adduct))
  
  table$matched_lipid<-paste(table$mz,table$lipid,table$c,table$delta)%in%paste(T1$paper_mz,T1$lipid,T1$c,T1$delta)
  
  table$lipid_in_MS2ID_old<-table$lipid_in_MS2ID
 
  table$lipid_in_MS2ID<-paste(table$lipid,table$c,table$delta) %in% paste(MS2ID_lipids$lipid,MS2ID_lipids$c,MS2ID_lipids$delta)
  
  subset(table,!matched_id&matched_mz)
  
  table$adduct_considered<-tolower(table$adduct)%in%tolower(subset(adducts,mode=="neg")$name)
  
  #Solving the missed matched that had an adduct considered (This is mostly a tolerance issue)
  s1<-subset(table,!matched_id&adduct_considered&matched_mz)
  subset(T1,paper_mz==s1$mz[2])
  
  table$matched_fragment_mz<-table$mz%in%T2_filters$paper_fragment_mz
  table$matched_fragment_id<-paste(table$mz,table$lipid,table$c,table$delta)%in%paste(T2_filters$paper_fragment_mz,T2_filters$parental_lipid,T2_filters$parental_c,T2_filters$parental_delta)
  
  write.csv(table,"/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results.csv")
  
  #Total Cases: 227
  nrow(table)
  
  #Different YELLOW cases 37 (16.3%)
  #YELLOW Case 1: Masses not accounted for 11 (4.8%)
  nrow(subset(table,!matched_mz))
  #YELLOW Case 2: Lipid not present in DB 26 (11.5%)
  nrow(subset(table,!lipid_in_MS2ID_old))
  
  #Case 2.1: No matches 2
  nrow(subset(table,is.na(lipid)))
  #Case 2.2: Found synonym 6
  nrow(subset(table,!lipid_in_MS2ID_old& lipid_in_MS2ID))
  #Case 2.3: Truly not found 18
  nrow(subset(table,!lipid_in_MS2ID&!is.na(lipid)))
  tmp<-subset(table,lipid_in_MS2ID&!is.na(lipid))
  tmp2<-(paste(tmp$lipid," ",tmp$c,":",tmp$delta,sep=""))
  length(tmp2)
  
  
  #Different GREEN cases:  
  
  
  #Different RED cases
  #RED Case 1: Adduct or fragment not accounted for
  #RED Case 2: Tolerance
  
  return(table)
}


generate_T1<-function(mz,tol=5,m="neg")
{
  ix<-1:length(mz)

  a<-subset(adducts,mode==m)
  if(m=="pos"){
    fragmentation<-fragmentation_pos
  }else{
    fragmentation<-fragmentation_neg
  }
  
  target_mz<-outer(outer(mz,adducts$mass,function(x,y)x+y),losses,function(x,y)x-y)
  results<-lapply(MS2ID_lipids$mz,function(x)withintol(x,target_mz,tol))
  return(results)
}
#Perform all DB searches [It can be optimized]


#IT'S A WHOLE NEW WORLD

RefMet<-read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/refmet.csv")

#Trim RefMet
lipids_pattern<-paste(lipids$name,collapse ="|")
RefMet_lipids<-subset(RefMet,grepl(lipids_pattern,refmet_name))
RefMet_lipids[c("c","delta")]<-t(sapply(RefMet_lipids$refmet_name,find_c_delta))
RefMet_lipids<-subset(RefMet_lipids,!is.na(c))
RefMet_lipids$lipid<-sapply(RefMet_lipids$refmet_name,function(x)find_lipid(x,lipids_pattern))
RefMet_lipids<-RefMet_lipids[rownames(unique(RefMet_lipids[c("lipid","c","delta","exactmass")])),]
RefMet_lipids<-RefMet_lipids[order(RefMet_lipids$exactmass),]

#Temporary fix
RefMet_lipids_min_oxigen<-RefMet_lipids
RefMet_lipids_min_oxigen$exactmass<-RefMet_lipids_min_oxigen$exactmass-15.99491462
RefMet_lipids<-rbind(RefMet_lipids,RefMet_lipids_min_oxigen)

sc<-lipids$name
sc[which(!lipids$name%in%unique(RefMet$sub_class))]<-c("DAG","TAG","Cer-1-P", "", "Gangliosides", "")

RefMet_lipids<-subset(RefMet_lipids,sub_class%in%sc)
#Generate a and f
a<-adducts$mass
names(a)<-adducts$name
f<-unlist(lapply(unique(unlist(c(fragmentation_neg,fragmentation_pos))),parse_fragmentation))
f<-f[!(f==0)]
f<-c(f,0)
f["-NH2"]<-16.01872407
f["-NH3"]<-17.02654911 #MODIFICATION LUCIA added NH3 instead of NH2

offset<-outer(a,f,function(x,y)x+y)

#Generate target

table<-load_gt("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_neg.csv")
mz<-unique(table$mz)

target<-outer(mz,offset,function(x,y)x-y)

#Generate matches

tol=20
db_mz<-RefMet_lipids$exactmass
target_min<-target-target*tol*10^-6
target_max<-target+target*tol*10^-6

matches<-lapply(db_mz,function(x)which((target_min<x)&(target_max>x),arr.ind = T))

#Generate final data.table
results<-do.call(rbind,lapply(seq_along(matches),function(i,x=matches[[i]])if(nrow(x)==0) NULL else data.frame(RefMet_lipids[i,c("sub_class","lipid","c","delta","exactmass")],mz_i=x[,1],adduct_i=x[,2],fragment_i=x[,3])))
results$experimental_mz<-mz[results$mz_i]
results$adduct<-names(a)[results$adduct_i]
results$adduct_mode<-adducts$mode[results$adduct_i]
results$adduct_score<-sapply(1:nrow(results),function(i)adducts[results$adduct_i[i],results$lipid[i]])
results$fragmentation<-names(f)[results$fragment_i]
results$fragmentation_possible<-sapply(1:nrow(results),function(i)frag_possible(results$adduct_mode[i],results$lipid[i],results$fragmentation[i]))
results$ppm_error=get_tol(results$exactmass,results$experimental_mz-a[results$adduct_i]-f[results$fragment_i])


r<-subset(results,adduct_mode=="neg")

frag_possible<-function(m,l,frag){
  if(m=="neg")
    f=fragmentation_neg
  else
    f=fragmentation_pos
  s<- paste(f[l,],collapse="|")
  s<- gsub(" ", "", s, fixed = T)
  n<-unique(unlist(strsplit(s,"|",fixed=T)))
  return(frag %in% n)
}

#Validate results

gt <- table

results$matched_lipid<-paste(results$lipid,results$c,results$delta,results$experimental_mz)%in%paste(gt$lipid,gt$c,gt$delta,gt$mz)
paste(results$lipid,results$c,results$delta,results$experimental_mz)%in%paste(gt$lipid,gt$c,gt$delta,gt$mz)
gt$matched_lipid<-paste(gt$lipid,gt$c,gt$delta,gt$mz)%in%paste(r$lipid,r$c,r$delta,r$experimental_mz)
gt$in_RefMet<-paste(gt$lipid,gt$c,gt$delta)%in%paste(RefMet_lipids$lipid,RefMet_lipids$c,RefMet_lipids$delta)
gt$matched_lipid_MS2ID<-table$matched_lipid


#Cross Validation with METASPACE
results_pos<-read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results_pos.csv")
results_neg<-read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results_neg.csv")

metaspace_pos_LM <- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_pos_LM.csv",sep=";")
metaspace_pos_HMDB <- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_pos_HMDB.csv",sep=";")
metaspace_pos_ChEBI<- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_pos_ChEBI.csv",sep=";")
metaspace_neg_LM <- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_neg_LM.csv",sep=";")
metaspace_neg_HMDB <- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_neg_HMDB.csv",sep=";")
metaspace_neg_ChEBI<- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_neg_ChEBI.csv",sep=";")

metaspace_pos_LM$moleculeNames[1:10]


get_lipids <- function(s){
  return(unique(sapply(strsplit(as.character(s),",")[[1]],find_lipid_c_delta)))
}

metaspace_pos_LM_lipids<-lapply(metaspace_pos_LM$moleculeNames,get_lipids)
metaspace_neg_LM_lipids<-lapply(metaspace_neg_LM$moleculeNames,get_lipids)

metaspace_pos_LM$paper_mz<-sapply(metaspace_pos_LM$mz,function(x)results_pos$mz[which.min(abs(x-results_pos$mz))])
metaspace_neg_LM$paper_mz<-sapply(metaspace_neg_LM$mz,function(x)results_neg$mz[which.min(abs(x-results_neg$mz))])

results_pos$matched_METASPACE<-paste(results_pos$lipid,paste(results_pos$c,results_pos$delta,sep=":"),results_pos$mz)%in%unlist(lapply(metaspace_pos_LM_lipids,function(x)paste(x,metaspace_pos_LM$paper_mz)))
results_neg$matched_METASPACE<-paste(results_neg$lipid,paste(results_neg$c,results_neg$delta,sep=":"),results_neg$mz)%in%unlist(lapply(metaspace_neg_LM_lipids,function(x)paste(x,metaspace_neg_LM$paper_mz)))


write.csv(results_pos,"/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results_pos_metaspace.csv")
write.csv(results_neg,"/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results_neg_metaspace.csv")

# 3. DISPLAY FUNCTIONS

display_fragmentation_network<-function(T1,T2){
  library('igraph')
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
  
  library(igraph)
  
  # Create data
  data <- matrix(sample(0:1, 400, replace=TRUE, prob=c(0.8,0.2)), nrow=20)
  data <- (c>0.9)
  losses<-round(c(59.0734993,85.08914936 ),2)
  T2_visualisation<-subset(T2_possible_fragments,cor>0.7&round(diff,2)%in%losses)
  
  edgelist<-matrix("",nrow(T2_visualisation),2)
  edgelist[,1]<-unlist(lapply(1:nrow(T2_visualisation),function(i)paste(T2_visualisation[i,c('parental_lipid','parental_c','parental_delta')],collapse=":")))
  edgelist[,2]<-unlist(lapply(1:nrow(T2_visualisation),function(i)paste(T2_visualisation[i,c('fragment_lipid','fragment_c','fragment_delta')],collapse=":")))

  T1_db_matches[unique(match(T2_visualisation$parental_id,T1_db_matches$id)),c("i","adduct")]
  T1_db_matches[unique(match(T2_visualisation$fragment_id,T1_db_matches$id)),c("i","adduct")]
  
  
  network<-graph_from_edgelist(edgelist)
  plot(network, layout=layout.auto, main="fruchterman.reingold")
  
  # random graph
  net = rgraph(10, mode = "graph", tprob = 0.5)
  net = network(net, directed = FALSE)
  
  # vertex names
  network.vertex.names(net) = letters[1:10]
  ggnet2(net)
}

bd<-"/home/gbaquer/msidata/1. In-source Fragmentation/1.3. AuBSi Wells/20200701_AuBSi_wells_std/20200701_AuBSi_wells_std/"
pks<-rMSIproc::LoadPeakMatrix(paste(bd,"20200701_AuBSi_wells_std-peaks.zip",sep=""))
pks_list<-lapply(1:3,function(i)rMSIcleanup:::get_one_peakMatrix(pks,i))


#Edit position

pks_neg[[1]]$pos[1:1800,1]<-rep(1:90,each=20)
pks_neg[[1]]$pos[1:1800,2]<-rep(1:20,each=90)
pks_neg[[1]]$pos[1801,]<-c(91,1)
pks_neg[[1]]$pos[1802,]<-c(91,2)
