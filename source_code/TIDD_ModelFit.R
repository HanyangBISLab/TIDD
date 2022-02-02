#   TIDD: peptide rescoring s/w
#
#   Written by H. Li <hllee@hanyang.ac.kr>
#
#   Copyright (C) 2022 BIS Labs, Hanyang Univ. Korea
#
#   TIDD is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   TIDD is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with TIDD.  If not, see <http://www.gnu.org/licenses/>.



library(ROCR)
library(readr)
library(e1071)

###############################################################

#write input/output file 

#inputfile<-"/home/hllee/Documents/TIDD/modif/phosphp/hela/Hela_phospho_comet/Hela_format_Feature.tsv"
#outputfile<-"/home/hllee/Documents/TIDD/modif/phosphp/hela/Hela_phospho_comet/Hela_format_Feature_MF.tsv"
inputfile<-"/home/hllee/Documents/TIDD/modif/phosphp/AML/AML_phospho_comet/AML_phospho_format_Feature.tsv"
outputfile<-"/home/hllee/Documents/TIDD/modif/phosphp/AML/AML_phospho_comet/AML_phospho_format_Feature_MF_impute.tsv"
MaxNum_iteration<-10
minPLEN<-6
trainsize<-10000
k_cross<-5
out_ycol<-"TorD"
init_score<-"Xcorr"
features=c("TorD","Xcorr","AbsDeltaMass","Charge","PeptideLength","Triptic","MissedCleavage","PrecursorM","DeltaMass","CalMW","TIC","SumYmatchInt","SumBmatchInt","FracYmatchInt","FracBmatchInt","MaxIntALL","MaxYionInt","MaxBionInt","SeqCoverYion","SeqCoverBion","ConsecutiveYion","ConsecutiveBion","MassErrMean","MassErrSD","NumofAnnoPeaks") 
#features=c("TorD","CalXcorr","NumofAnnoPeaks","ConsecutiveYion","SeqCoverYion")
###############################################################

# read file

  data<-read.delim(sprintf("%s",inputfile))
 
  naid<-which(is.na(data$MassErrMean)==TRUE)
  data[naid,]$MassErrMean<-mean(na.omit(data$MassErrMean))

  naid<-which(is.na(data$MassErrSD)==TRUE)
  data[naid,]$MassErrSD<-mean(na.omit(data$MassErrSD))

  data<-na.omit(data)
  plen<-which(data[,"PeptideLength"]>=minPLEN)
  data<-data[plen,]
  names(data)
  
  data[,init_score]<-as.numeric(data[,init_score])  
  
  dec_id<-which(data[,out_ycol]=="D")
  
  scores<-unique(round(data[dec_id, init_score],3))
  
  #data$lnExpect<- (-1)*data$lnExpect
  
  #tmpdata<-cbind(data$Xcorr,data$TorD)
  
  IDs<-rep(0,as.numeric(MaxNum_iteration)+1)
  IDs_name<-rep(0,as.numeric(MaxNum_iteration)+1)
  
  #naid<-which(is.na(data$MassErrMean)==TRUE)
  #data[naid,]$MassErrMean<-mean(na.omit(data$MassErrMean))
  
  
###############################################################
  
# init putative target identification
  
  Results_FDR<- -1000
  score_threshold<-1000
  CIDS<-0
  
  for(k in 1:length(scores))
  {
    TP<-which(data[,init_score]>=scores[k] & data[,out_ycol]=="T")
    FP<-which(data[,init_score]>=scores[k] & data[,out_ycol]=="D")
    
    FDR=length(FP)/length(TP)
    
    if((0.01-FDR)>=0 & abs(0.01-FDR)<=abs(0.01-Results_FDR))
    {
      
      if(CIDS<length(TP))
      {
        Results_FDR<-FDR
        score_threshold<-scores[k]
        CIDS<-length(TP)
      }
    }
  }
  
  Results_FDR
  score_threshold
  print({paste("Init ID:",CIDS)})
 
###############################################################
  
  # change categorical features to "factor"
  
  IDs[1]<-CIDS 
  IDs_name[1]<-"Init"
  
  tar_candidate_id<-which(data[,init_score]>=score_threshold & data[,out_ycol]=="T")
  
  ndata<-data[,features]
  ndata[,out_ycol]<-as.factor(as.character(ndata[,out_ycol]))
  
  if("Charge" %in% features)
  {
    id<-which(ndata$Charge>=4)
    
    ndata[id,]$Charge<-4
    
    ndata$Charge<-as.factor(ndata$Charge)
  }
  
  if("Triptic" %in% features)
  {
    ndata$Triptic<-as.factor(ndata$Triptic)
  }
  
  colnames(ndata)[which(names(ndata)==out_ycol)]<-"TorD"

###############################################################
  
  # fitting 
  
  for(iteration in 1: as.numeric(MaxNum_iteration))
  {
   
    minsize<-min(length(dec_id),length(tar_candidate_id),trainsize)
    
    tar_c<-sample(tar_candidate_id,minsize,replace=FALSE)
    dec_c<-sample(dec_id,minsize,replace=FALSE)
    
    CV<-as.numeric(k_cross)
    cv.error<-rep(0,CV)
    auc.error<-rep(0,CV)
    
    x_tar<-ndata[tar_c,]
    x_dec<-ndata[dec_c,]
    
    set.seed(11)
    t_fold<-cut(seq(1,nrow(x_tar)),breaks=CV,labels=FALSE)
    d_fold<-cut(seq(1,nrow(x_dec)),breaks=CV,labels=FALSE)
    
    # cross validation
    for(i in 1:CV)
    {
      
      t_testid<-which(t_fold==i,arr.ind=TRUE)
      d_testid<-which(d_fold==i,arr.ind=TRUE)
      
      Tdata<-rbind(x_tar[t_testid,],x_dec[d_testid,])
      TrainData<-rbind(x_tar[-t_testid,],x_dec[-d_testid,])
      
      set.seed(1)
      svm.fit=svm(formula = TrainData$TorD~ ., data = TrainData, probability=TRUE, kernel="linear", type="C-classification")
      
      pred.prob<-predict(svm.fit,Tdata,probability=TRUE)
      pred_1<-prediction(attr(pred.prob,"probabilities")[,1],Tdata$TorD)
      perf<-performance(pred_1,"auc")
      auc.error[i]<-as.numeric(as.character(unlist(perf@y.values[[1]])))
    }
    
    #          print(auc.error)
    SiM<-100
    index<-1
    meanAUC<-mean(auc.error)
    for(k in 1:CV)
    {
      if(abs(auc.error[k]-meanAUC)<SiM)
      {
        SiM<-abs(auc.error[k]-meanAUC)
        index<-k
      }
    }
    #TEST
    t_testid<-which(t_fold==index,arr.ind=TRUE)
    d_testid<-which(d_fold==index,arr.ind=TRUE)
    
    TrainData<-rbind(x_tar[-t_testid,],x_dec[-d_testid,])
    svm.fit=svm(formula = TrainData$TorD ~ ., data = TrainData, type = "C",probability=TRUE, kernel="linear")
    pred.prob<-predict(svm.fit,ndata,probability=TRUE)
    
    pred_1<-prediction(attr(pred.prob,"probabilities")[,1],ndata$TorD)
    newTestData<-cbind(data,attr(pred.prob,"probabilities")[,1])   
    colnames(newTestData)[ncol(newTestData)]<-c(paste0 ("SVM_Prob_",iteration))
    
    write.table(newTestData,sprintf("%s",outputfile),sep="\t")
    scorename<-colnames(newTestData)[ncol(newTestData)]
    
    scores<-unique(round(newTestData[dec_id,scorename],3))
    length(scores)
    
    Results_FDR<- -1000
    score_threshold<-1000
    CIDS<-0
    
    for(k in 1:length(scores))
    {
      
      TP<-which(newTestData[,scorename]>=scores[k] & newTestData[,out_ycol]=="T")
      FP<-which(newTestData[,scorename]>=scores[k] & newTestData[,out_ycol]=="D")
      
      if(length(TP)>0)
      {
        FDR=length(FP)/length(TP)
        
        if(0.01-FDR>=0 & abs(0.01-FDR)<=abs(0.01-Results_FDR))
        {
          
          if(CIDS<-length(which(newTestData[,scorename]>=scores[k] & newTestData[,out_ycol]=="T")))
          {
            Results_FDR<-FDR
            score_threshold<-scores[k]
            CIDS<-length(which(newTestData[,scorename]>=scores[k] & newTestData[,out_ycol]=="T"))
          }
        }
      }
      
    }
    
    Results_FDR
    score_threshold
    
    tar_candidate_id<-which(newTestData[,scorename]>=score_threshold & newTestData[,out_ycol]=="T")
    print(paste(paste("iteration ",iteration),": ", length(tar_candidate_id)))
    IDs[(iteration+1)]<-length(tar_candidate_id)
    IDs_name[(iteration+1)]<-paste0("Iter",iteration)
    
    write.table(newTestData[tar_candidate_id,],sprintf("%s_ID.tsv",outputfile),sep="\t")
    
  }  

  