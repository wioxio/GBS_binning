
data_name="" #your data name 
data.total<-read.table(data_name,header=FALSE,sep="\t")

data.total=as.matrix(data.total)
ids<-data.total[,1]
chromosome=ids
chromosome=strsplit(as.character(chromosome),":")
distance.total=sapply(chromosome,function(x) x[2])
chromosome=sapply(chromosome,function(x) x[1])
distance.total=as.numeric(as.character(distance.total))
chromosome=as.numeric(as.character(chromosome))
data.total=(data.total[,-1])
data.total[which(data.total=="A")]=1
data.total[which(data.total=="B")]=0
data.total[which(data.total=="H")]=0.5
data.total[which(data.total=="-")]=NA
class(data.total) <- "numeric"
break_points=matrix(list(),8)
for(i in 1:8){break_points[[i]]=0}


bin_starts=integer()
NofislandsR_total=integer()
binC=1
#x <- scan(break_points_file_name, what="", sep="\n")
#pre_break_points <- strsplit(x, ",[[:space:]]+")
#pre_joinmap_choice=read.table("joinmap_choice",header=FALSE)
#pre_joinmap_choice=as.numeric(as.character(as.matrix(pre_joinmap_choice)))
#pdf(paste("pdf_file_name.pdf",sep=""))
for(chr in min(chromosome):max(chromosome)){
  
  
  
  print(paste("chromosome",chr))
  data=data.total[which(chromosome==chr),]
  data=as.matrix(data)
  mode(data)<-"numeric"
  distance=distance.total[which(chromosome==chr)]
  mode(distance)<-"numeric"
  ids.chr=ids[which(chromosome==chr)]
  
  
  aftersmoothingA=matrix(0,length(data[1,]),length(data[,1])) #this will be the final result
  lastsnp=length(data[,1]) #number of snps
  Nofislands=array(0,dim=c(binC,length(data[1,])))#for the distribution of the number of contigous regions
  
  
 
  minD=50000 
  maxD=5000000
  Dseq=10000000*1/5
  #block when using a single value# Dseq=seq.int(minD,maxD,length.out=binC) #a sequence of bin sizes
  distanceIndex=array(0,dim=c(length(data[,1]),binC,2))
  
  
  
  
  for(x in 1:length(data[,1])){ #this loop is for assigning proper flanking snps according to the given bin size
    
    
    k=length(Dseq)
    
    
    temp=intersect(which(distance[x]-distance<=Dseq[k]),which(distance[x]-distance>=0))
	if(length(temp)==0){begin=distanceIndex[x-1,k,1]
	}else{    begin=min(temp)
	}
	
	temp=intersect(which(distance-distance[x]<=Dseq[k]),which(distance-distance[x]>=0))
	if(length(temp)==0){end=distanceIndex[x-1,k,2]
	}else{    end=max(temp)
	}
	
    distanceIndex[x,k,]=c(begin,end)
    
      k=k-1
    while(k>0){
      
      #begin=ifelse(begin+1>x,x,begin+1) missing the case that begin stays
      distanceIndex[x,k,1]=begin
      while(distance[x]-distance[begin]>Dseq[k]){
        begin=ifelse(begin+1>x,x,begin+1)
      }
      k=k-1
      
    }
    
    k=length(Dseq)-1
    while(k>0){
      distanceIndex[x,k,2]=end
      while(distance[end]-distance[x]>Dseq[k]){
        end=ifelse(end-1<x,x,end-1)
      }
      k=k-1
      
    }
    
    
    
  }
  
  
  
  
  
  ##the main part
  
  
  for(j in 1:length(data[1,])){#for each individual
    {
     # print(paste("at",j))
      for(t in 1:binC){# distanceIndex's number of testing threshold
       # print(paste("index",t))
        newdata=data
        
        #deal with the first region and the last region(to prevent the case that the first bin and the last bin are NAs)
        
        newdata[distanceIndex[1,t,1]:distanceIndex[1,t,2],j]=mean(data[distanceIndex[1,t,1]:distanceIndex[1,t,2],j],na.rm=TRUE)
        newdata[distanceIndex[lastsnp,t,1]:distanceIndex[lastsnp,t,2],j]=mean(data[distanceIndex[lastsnp,t,1]:distanceIndex[lastsnp,t,2],j],na.rm=TRUE) 	
        
        difference=matrix(NA,1,length(data[,1]))#difference of the left and right
        
        for( i in 1:(length(data[,1])-1)){
          #compare right and left
          difference[i]=mean(newdata[distanceIndex[i,t,1]:i,j],na.rm=TRUE)-mean(newdata[(i):distanceIndex[i,t,2],j],na.rm=TRUE)#changed from i+1 to i in the second term
          
          
        }
       # print(paste("index sub1",t))
        #pick peaks
        tempD=abs(difference)
        peaks=which(tempD==max(tempD,na.rm=TRUE))[1]
        peakV=max(tempD,na.rm=TRUE)[1]
        if(peakV<0.3){peaks=NA #if it is less than 0.3, it's regarded as small
        }else{
          tempD[distanceIndex[peaks,t,1]:distanceIndex[peaks,t,2]]=0
        }
        while(peakV>=0.3){# while having meaningful peaks
          peakstemp=which(tempD==max(tempD,na.rm=TRUE))[1]
          peakV=max(tempD,na.rm=TRUE)[1]
          if(peakV>=0.3){
            peaks=c(peaks,peakstemp)
            tempD[distanceIndex[peakstemp,t,1]:distanceIndex[peakstemp,t,2]]=0
          }
        }
      #  print(paste("index sub2",t))
        #smoothing
        aftersmoothing=data[,j]
        
        #if(!is.na(peaks)&length(data[,1])-peaks[1]<300){peaks=NA}	#short one peak?
        
        if(is.na(peaks)){#if there is no peak
        #  print(paste("index sub3_1",t))
          meantemp=mean(newdata[,j],na.rm=TRUE)
          if(is.na(meantemp)){ next }#if the is no data
          dis=c(abs(meantemp),abs(meantemp-1),abs(meantemp-0.5))
          resultvalue=integer(0)
          if(dis[1]<=dis[2]){resultvalue=0
          if(dis[3]<=dis[1]){resultvalue=0.5}
          }else if(dis[2]<=dis[3]){resultvalue=1
          }else{resultvalue=0.5}
          
       #   print(paste("index sub3",t))
          
          aftersmoothing=rep(resultvalue,length(aftersmoothing))
        }else{#if there are peaks
       #   print(paste("index sub4",t))
          peaks=sort(peaks,decreasing=FALSE)
          
          if(length(peaks)>=2){ 	
            l=length(peaks)
            for(s in 2:l){
              
              if(distance[peaks[s]]-distance[peaks[s-1]]<(Dseq[t]*2)&abs(mean(newdata[ifelse((s-2)>=1,peaks[(s-2)],1):peaks[s-1],j],na.rm=TRUE)-mean(newdata[peaks[s]:ifelse(s+1<=length(peaks),peaks[s+1],length(data[,1])),j],na.rm=TRUE))<=0.3){ peaks=peaks[-s]} #remove short region & flanking regions with small difference
              if((s+1)>length(peaks)){break}
              
            }
            
          }
          
          
          meantemp=mean(newdata[1:peaks[1],j],na.rm=TRUE)
          dis=c(abs(meantemp),abs(meantemp-1),abs(meantemp-0.5))
          resultvalue=integer(0)
          if(dis[1]<=dis[2]){resultvalue=0
          if(dis[3]<=dis[1]){resultvalue=0.5}
          }else if(dis[2]<=dis[3]){resultvalue=1
          }else{resultvalue=0.5}
          
          
          
          aftersmoothing[1:peaks[1]]=resultvalue #the first bin
          
          if(length(peaks)>=2){
            for(s in 2:length(peaks)){#form the second bin
              
              meantemp=mean(newdata[peaks[s-1]:peaks[s],j],na.rm=TRUE)
              dis=c(abs(meantemp),abs(meantemp-1),abs(meantemp-0.5))
              resultvalue=integer(0)
              if(dis[1]<=dis[2]){resultvalue=0
              if(dis[3]<=dis[1]){resultvalue=0.5}
              }else if(dis[2]<=dis[3]){resultvalue=1
              }else{resultvalue=0.5}
              
              
              
              aftersmoothing[(peaks[s-1]+1):peaks[s]]=resultvalue
              
            }
          }
          meantemp=mean(newdata[peaks[length(peaks)]:length(data[,1]),j],na.rm=TRUE)
          dis=c(abs(meantemp),abs(meantemp-1),abs(meantemp-0.5))
          resultvalue=integer(0)
          if(dis[1]<=dis[2]){resultvalue=0
          if(dis[3]<=dis[1]){resultvalue=0.5}
          }else if(dis[2]<=dis[3]){resultvalue=1
          }else{resultvalue=0.5}
          
          
      ##    print(paste("index sub5",t))
          aftersmoothing[(peaks[length(peaks)]+1):length(data[,1])]=resultvalue
        }#not NA
        
        # par(mfrow=c(2,1))
        
      #  print(paste("index sub6",t))
        aftersmoothingA[j,]=aftersmoothing
        
        #   plot(aftersmoothing,ylim=c(0,1))
        # plot(data[,j],col="green",xlab="", ylab="",lwd=1)
        
        
        
        #count the number of islands
        bin_starts_temp=integer()
        for(q in 2:lastsnp){
          if(aftersmoothing[q-1]!=aftersmoothing[q]){ Nofislands[t,j]=Nofislands[t,j]+1 ;
          bin_starts_temp=c(bin_starts_temp,distance[q])
          break_points[[chr]]=c(break_points[[chr]],distance[q])
          }
          #   
        }
        bin_starts=rbind(bin_starts,c(chr,j,t,bin_starts_temp))
        # print(paste("index sub7",t))
        if(chr==0){#closing this section because there are too many plots
          plot(as.numeric(as.matrix(aftersmoothingA[j,])),col="green",pch="*",ylim=c(0,1),main=paste("binning index:",t,"chromsome:",chr,"individual:",j));
          #abline(v=which(distance%in%as.numeric(pre_break_points[[chr]])),col="grey")
          #abline(v=pre_joinmap_choice,col="grey")
          }
          # if(length(bin_starts_temp)>0){ plot(as.numeric(as.matrix(aftersmoothingA[j,])),col="green",pch="*",ylim=c(0,1),main=paste("binning index:",t,"chromsome:",chr,"individual:",j));abline(v=which(distance%in%as.numeric(pre_break_points[[chr]])),col="grey")
        # }
        abc=1
       }#t iteration:distanceIndex's number of testing threshold
      if(chr==0){
       plot(as.numeric(as.matrix(data[,j])),col="blue",pch="*",ylim=c(0,1),main=paste("binning index:",t,"chromsome:",chr,"individual:",j)); #plots the original data
        #abline(v=which(distance%in%as.numeric(pre_break_points[[chr]])),col="grey")
        #abline(v=pre_joinmap_choice,col="grey")
        }
    }#j iteration/ individual
  }
  #dev.off()
  
  
  #for the bin distributoin
  NofislandsR=apply(Nofislands,1,sum) #to check the number of contiguous regions based on each bin size, plot this
  NofislandsR_total=cbind(NofislandsR_total,NofislandsR)
  ##plot(NofislandsR,main=chr)
  
  #below should be applied applying the optimal t
  
  break_points[[chr]]=sort(unique(break_points[[chr]]))
  break_points[[chr]]=  break_points[[chr]][-1]
  
  bp_start= 1
  bp_end= which(distance==break_points[[chr]][1])
  middle=round((bp_start+bp_end)/2)
  joinmap_choice=middle
  
  for(bps in 2:length(break_points[[chr]])){
    bp_start= which(distance==break_points[[chr]][bps-1])
    bp_end= which(distance==break_points[[chr]][bps])
    middle=round((bp_start+bp_end)/2)
    joinmap_choice=c(joinmap_choice,middle)
  }
  bp_start= which(distance==break_points[[chr]][length(break_points[[chr]])])
  bp_end= length(break_points[[chr]])
  middle=round((bp_start+bp_end)/2)
  joinmap_choice=c(joinmap_choice,middle)
  joinmap_choice=sort(unique(joinmap_choice))
  
  break_points.smoothed=t(aftersmoothingA[,joinmap_choice])
  break_points.data=data[joinmap_choice,]
  break_points.ids=ids.chr[joinmap_choice]
  mixed=integer()
  for(s in 1:length(joinmap_choice)){
    mixed=cbind(mixed,break_points.smoothed[s,],break_points.data[s,])
    colnames(mixed)[(2*(s-1)+1):(2*s)]=c(ids.chr[s],ids.chr[s])
  }
  break_points.smoothed[which(break_points.smoothed==1)]="A"
  break_points.smoothed[which(break_points.smoothed==0)]="B"
  break_points.smoothed[which(break_points.smoothed==0.5)]="H"
  break_points.smoothed[which(is.na(break_points.smoothed)==TRUE)]="-"
  if(chr==0){write.table((joinmap_choice),file="joinmap_choice",quote=FALSE,row.names=FALSE,col.names=FALSE)}
  write.table(mixed,file=paste(data_name,"chosen_smoothed_original_comparison.",chr,sep=""),col.names=colnames(mixed),row.names=FALSE,quote=FALSE)
  
  
  write.table(paste("NAME=Cr\nPOPT=F2\nNLOC=",length(break_points.smoothed[,1]),"\nNIND=",length(break_points.smoothed[1,]),sep=""),file=paste(data_name,".chosen_smoothed.chr",chr,sep=""),row.names=  FALSE,col.names=FALSE,quote=FALSE)
  write.table((break_points.smoothed),file=paste(data_name,".chosen_smoothed.chr",chr,sep=""),row.names=  break_points.ids,col.names=FALSE,quote=FALSE,append=TRUE,sep="\t")
  #write.table(break_points.data,file=paste(data_name,".chosen_smoothed.chr",chr,sep=""),col.names=  break_points.ids,row.names=FALSE,quote=FALSE)
  write.table(break_points.ids,file=paste(data_name,".chosen_markers",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
}#chromsome

NofislandsR_total

#dev.off()

