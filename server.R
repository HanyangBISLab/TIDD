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



library(shiny)
library(e1071)
library(ROCR)

# Shiny Server Side -------
options(shiny.maxRequestSize = 10240*1024^2)

server <- function(input, output, session) {
  
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste("D", Sys.Date(), "_",input$downloadpsm, sep="")
    },
    content = function(file) {
      tmp<-read.delim(paste0('data/PSM/',input$downloadpsm))
      
      write.table(tmp, file,sep="\t")
    })
  output$downloadmgfData <- downloadHandler(
    filename = function() { 
      paste("D", Sys.Date(), "_",input$downloadmgf, sep="")
    },
    content = function(file) {
      
      tmp<-read.delim(paste0('data/MGF/test/',input$downloadmgf))
      
      write.table(tmp, file,sep="\t")
    })
  
  observe({
    
    if (is.null(input$upload)) 
      return()
    else
    {
      path="data/"
      if(input$filetype=="1")
        path<-paste0(path,"PSM/")
      else
        path<-paste0(path,"MGF/Local/")
      
      file.copy(input$upload$datapath, paste0(path,input$upload$name), recursive=TRUE)
      
    }
  })
  
  
  observeEvent(input$loadButton,{
    
    psms_part<-read.delim(paste0('data/PSM/',input$psmfile))
    psms<-renderTable(psms_part)
    
    print(typeof(psms_part))
    print(dim(psms))
    
    output$psm_filename<-renderText(input$psmfile)
    output$psmlist_preview<-renderDataTable(psms_part,options = list(scrollX = TRUE,
                                                                     lengthChange=FALSE,
                                                                     info=FALSE,
                                                                     searching=FALSE))
    tmp=rep(0,length(names(psms_part)))
    
    for(i in 1:length(names(psms_part)))
    {
      tmp[i]<-paste0(i,".",names(psms_part)[i]) 
    }
    
    updateRadioButtons(session=session,inputId="psmf",choices = tmp)
    
    observeEvent(input$featureButton,{
      output$tidd0<-renderText({"Extracing features: finished [result file: *feature_*.tsv and *_Xcorr.tsv"})
      
      output$tidd1<-renderText({
        
        progress<-shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message="Extracting features now, it takes a few minites...",value=0)
        
        
        specdir<-paste0('data/MGF/',input$msfile,"/")
        psm<-paste0('data/PSM/',input$psmfile)
        
        cmd<-sprintf("java -jar JAVA/TIDD.jar -i %s -dir %s -id %s -dec %s -ft %s", psm,specdir,input$required_feature,input$dec_pre,input$frag_tol) 
        
        
        system(cmd,intern=TRUE)
        
        psm_feature<-paste0('data/PSM/',strsplit(input$psmfile,".tsv")[1],'_Feature.tsv')
        info<-strsplit(input$required_feature,",")[[1]]
        
        mgffileCol<-as.numeric(info[1])-1
        chargeCol<-as.numeric(info[3])-1
        peptideCol<-as.numeric(info[5])-1
        
        psm_num_Col<-length(psms_part)-1
        
        titleCol<-psm_num_Col+1
        precursorMZ_Col<-psm_num_Col+3
        
        #Xcorr Usage: java -jar xcorr.jar spectrumFile resultFile binWidth binOffset titleCol peptideCol precursorMZCol chargeCol mgffilenameCol header");
        cmd2<-sprintf("java -jar JAVA/Xcorr_v2.jar %s %s 0.02 0.0 %s %s %s %s %s true",specdir,psm_feature,titleCol,peptideCol,precursorMZ_Col,chargeCol,mgffileCol)
        
        system(cmd2,intern=TRUE)
        
        a<-"DONE!"
        b<-"1._feature.tsv: contain annotated peaks - related features"
        c<-"2._xcorrColAdded.tsv: add calculated xcorr features"
        
        print(paste(a,b,c,sep="\n"))
      })
    })
    
    observeEvent(input$psmf,{
      
      index<-""
      b<-""
      if(!is.null(psms_part))
      {
        index<-as.character(input$psmf)
        index<-as.integer(as.numeric(strsplit(index,"\\.")[[1]][1]))
        
        
        output$summary_feature<-renderPrint({
          
          
          print(dim(psms_part))
          print(names(psms_part)[index])
          print(summary(psms_part[c(1:nrow(psms_part)),index]))
          
        })
        
        if(is.numeric(psms_part[1,index]))
        {
          output$feature_plot<-renderPlot({
            plot(density(psms_part[c(1:nrow(psms_part)),index]), main=names(psms_part)[index])
          })  
        }
        
      }
      
    })
  })
  
  feature_part<-observeEvent(input$load2Button,{
    
    tmp<-read.delim(paste0('data/PSM/',input$psm_feature_file),nrows=1)
    
    output$feature_filename<-renderText(paste0("Input file:",input$psm_feature_file))
    
    lst<-names(tmp)
    
    updateCheckboxGroupInput(session,"check1",choices=lst, selected=c("DeltaXcorr","AbsDeltaMass","Charge","PeptideLength","Triptic","MissedCleavage","PrecursorM","DeltaMass","CalMW","TIC","SumYmatchInt","SumBmatchInt","FracYmatchInt","FracBmatchInt","MaxIntALL","MaxYionInt","MaxBionInt","SeqCoverYion","SeqCoverBion","ConsecutiveYion","ConsecutiveBion","MassErrMean","MassErrSD","NumofAnnoPeaks","TorD","CalXcorr") )
    
    updateSelectInput(session=session, inputId = "init_score",choices=c(" ", input$check1), selected=c(" "))
    updateSelectInput(session=session, inputId = "y", choices=c(" ", input$check1), selected=c(" "))
    
  })
  
  observe({
    x<-input$check1
    if(is.null(x))
      x<-character(0)
    
    updateSelectInput(session=session, inputId = "init_score",choices=x,)
    updateSelectInput(session=session, inputId = "y",choices=x)
  })
  
  observeEvent(input$modelbutton,{
    
    output$confirm<-renderText({
      
      IDs<-rep(0,as.numeric(input$iteration)+1)
      IDs_name<-rep(0,as.numeric(input$iteration)+1)
      
      if(is.null(input$check1))
      {  return("Select features and parameters!")}
      else
      {
        tmp<-read.delim(paste0('data/PSM/',input$psm_feature_file))
        
        output$static<-renderText(
          {
            paste("# instances: ", nrow(tmp))
          })
        
        data<-na.omit(tmp)
        dec_id<-which(data[,input$y]=="D")
        
        scores<-unique(round(data[dec_id,input$init_score],3))
        Results_FDR<- -1000
        score_threshold<-1000
        CIDS<-0
        
        
        #  print(length(scores))
        
        for(k in 1:length(scores))
        {
          TP<-which(data[,input$init_score]>=scores[k] & data[,input$y]=="T")
          FP<-which(data[,input$init_score]>=scores[k] & data[,input$y]=="D")
          
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
        
        output$initID<-renderText({paste("Init ID:",CIDS)})
        IDs[1]<-CIDS 
        IDs_name[1]<-"Init"
        tar_candidate_id<-which(data[,input$init_score]>=scores[k] & data[,input$y]=="T")
        
        ndata<-data[,input$check1]
        ndata[,input$y]<-as.factor(as.character(ndata[,input$y]))
        
        if("Charge" %in% input$check1)
        {
          id<-which(ndata$Charge>=4)
          
          ndata[id,]$Charge<-4
          
          ndata$Charge<-as.factor(ndata$Charge)
        }
        if("Triptic" %in% input$check1)
        {
          ndata$Triptic<-as.factor(ndata$Triptic)
        }
        
        #      print(summary(ndata))
        
        #      print(typeof(ndata[,1]))
        #      print(typeof(ndata[,2]))
        
        colnames(ndata)[which(names(ndata)==input$y)]<-"TorD"
        
        fitting<-reactive({
          
          for(iteration in 1: as.numeric(input$iteration))
          {
            
            minTrain_size<-min(tar_candidate_id,dec_id, as.numeric(input$tain_size))
            tar_c<-sample(tar_candidate_id,minTrain_size,replace=FALSE)
            dec_c<-sample(dec_id,minTrain_size,replace=FALSE)
            
            CV<-as.numeric(input$k_cross)
            cv.error<-rep(0,CV)
            auc.error<-rep(0,CV)
            
            x_tar<-ndata[tar_c,]
            x_dec<-ndata[dec_c,]
            
            set.seed(11)
            t_fold<-cut(seq(1,nrow(x_tar)),breaks=CV,labels=FALSE)
            d_fold<-cut(seq(1,nrow(x_dec)),breaks=CV,labels=FALSE)
            
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
            #########TEST
            t_testid<-which(t_fold==index,arr.ind=TRUE)
            d_testid<-which(d_fold==index,arr.ind=TRUE)
            
            TrainData<-rbind(x_tar[-t_testid,],x_dec[-d_testid,])
            svm.fit=svm(formula = TrainData$TorD ~ ., data = TrainData, type = "C",probability=TRUE, kernel="linear")
            pred.prob<-predict(svm.fit,ndata,probability=TRUE)
            
            pred_1<-prediction(attr(pred.prob,"probabilities")[,1],ndata$TorD)
            newTestData<-cbind(ndata,attr(pred.prob,"probabilities")[,1])   
            colnames(newTestData)[ncol(newTestData)]<-c(paste0 ("SVM_Prob_",iteration))
            
            write.table(newTestData,sprintf("data/PSM/%s_ModelFitting_Results.tsv",input$psm_feature_file),sep="\t")
            scorename<-colnames(newTestData)[ncol(newTestData)]
            
            scores<-unique(round(newTestData[dec_id,scorename],3))
            length(scores)
            #tmpdata<-cbind(data$Xcorr,data$TorD)
            
            Results_FDR<- -1000
            score_threshold<-1000
            CIDS<-0
            
            for(k in 1:length(scores))
            {
              
              TP<-which(newTestData[,scorename]>=scores[k] & newTestData[,input$y]=="T")
              FP<-which(newTestData[,scorename]>=scores[k] & newTestData[,input$y]=="D")
              
              if(length(TP)>0)
              {
                FDR=length(FP)/length(TP)
                
                if(0.01-FDR>=0 & abs(0.01-FDR)<=abs(0.01-Results_FDR))
                {
                  
                  if(CIDS<-length(which(newTestData[,scorename]>=scores[k] & newTestData[,input$y]=="T")))
                  {
                    Results_FDR<-FDR
                    score_threshold<-scores[k]
                    CIDS<-length(which(newTestData[,scorename]>=scores[k] & newTestData[,input$y]=="T"))
                  }
                }
              }
              
            }
            
            Results_FDR
            score_threshold
            
            tar_candidate_id<-which(newTestData[,scorename]>=score_threshold & newTestData[,input$y]=="T")
            print(paste(paste("iteration ",iteration),": ", length(tar_candidate_id)))
            IDs[(iteration+1)]<-length(tar_candidate_id)
            IDs_name[(iteration+1)]<-paste0("Iter",iteration)
            
          }
          
          
          output$IDplot<-renderPlot({
            barplot(IDs, main="Identified PSMs", names=IDs_name,xlab="Model", ylab="# of identified PSMs")
          })
          
        })
        
        output$log<-renderPrint({
          fitting()
        })
        return("Done!")
      }
    }) 
  })
}




