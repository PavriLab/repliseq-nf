#!/usr/bin/env Rscript

# MIT License
#
# Copyright (c) 2020 Tobias Neumann
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

library(getopt)

spec = matrix(c(
  'help' , 'h', 0, "logical","print the usage of the command",
  'RT', "r", 2,"character","txt table of RT values",
  'span', "s", 2,"integer","loess smoothing span"
),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==1 ) {
  #get the script name
  cmd = commandArgs(FALSE)
  self = strsplit(cmd[grep("--file",cmd)],"=")[[1]][2]
  cat(basename(self),": Create mismatch plots from rate tabs.\n\n")
  #print a friendly message and exit with a non-zero error code
  cat(getopt(spec,command = self,usage=T))
  q(status=1)
}


if ( is.null(opt$RT) ) stop("arg RT must be specified")
if ( is.null(opt$span) ) { opt$span = 300000 }


library(preprocessCore)

################################################
## Prepare counts                             ##
################################################

merge<-read.table("merge_RT.txt" , header=TRUE)

#################################################
## Loess smoothing on raw data                 ##
#################################################

chrs=grep(levels(merge$chr),pattern= "[_YM]",invert=TRUE,value=TRUE)

AllLoess=list()

for(i in 1:(ncol(merge)-3)){
  AllLoess[[i]]=data.frame();
  cat("Current dataset:", colnames(merge)[i+3], "\n");
  for(Chr in chrs){
    RTb=subset(merge, merge$chr==Chr);
    lspan=opt$span/(max(RTb$start)-min(RTb$start));
    cat("Current chrom:" , Chr, "\n");
    RTla=loess(RTb[,i+3] ~ RTb$start, span=lspan);
    RTl=data.frame(c(rep(Chr,times=RTla$n)), RTla$x, merge[which( merge$chr==Chr & merge$start %in% RTla$x),3],RTla$fitted)
    colnames(RTl)=c("chr" , "start" , "end" ,colnames(RTb)[i+3])
    if(length(AllLoess[[i]])!=0){
      AllLoess[[i]]=rbind(AllLoess[[i]],RTl)
    }
    if(length(AllLoess[[i]])==0){
      AllLoess[[i]] = RTl
    }
  }
}

for(i in 1:length(AllLoess)){
  write.table(AllLoess[[i]][complete.cases(AllLoess[[i]]),],
    gsub(".bg" , "loess.bg" , colnames(AllLoess[[i]])) [4],
    sep= "\t" , row.names=FALSE, quote=FALSE, col.names = FALSE
  )
}

################################################
## Quantile normalization                     ##
################################################

merge_values<-as.matrix(merge[,4:ncol(merge)])

ad<-stack(merge[,4:ncol(merge)])$values

norm_data<-normalize.quantiles.use.target(merge_values,ad)
merge_norm<-data.frame(merge[,1:3],norm_data)
colnames(merge_norm)<-colnames(merge)

for(i in 4:ncol(merge_norm)){
  write.table(merge_norm[complete.cases(merge_norm[,i]), c(1,2,3,i)],
    gsub(".bg" , "qnorm.bg", colnames(merge_norm)[i]),
    sep= "\t" ,row.names=FALSE, quote=FALSE, col.names = FALSE
  )
}

#################################################
## Loess smoothing on quantile-normalized data ##
#################################################

chrs=grep(levels(merge_norm$chr),pattern= "[_YM]",invert=TRUE,value=TRUE)

AllLoess=list()

for(i in 1:(ncol(merge_norm)-3)){
  AllLoess[[i]]=data.frame();
  cat("Current dataset:", colnames(merge_norm)[i+3], "\n");
  for(Chr in chrs){
    RTb=subset(merge_norm, merge_norm$chr==Chr);
    lspan=opt$span/(max(RTb$start)-min(RTb$start));
    cat("Current chrom:" , Chr, "\n");
    RTla=loess(RTb[,i+3] ~ RTb$start, span=lspan);
    RTl=data.frame(c(rep(Chr,times=RTla$n)), RTla$x, merge_norm[which( merge_norm$chr==Chr & merge_norm$start %in% RTla$x),3],RTla$fitted)
    colnames(RTl)=c("chr" , "start" , "end" ,colnames(RTb)[i+3])
    if(length(AllLoess[[i]])!=0){
      AllLoess[[i]]=rbind(AllLoess[[i]],RTl)
    }
    if(length(AllLoess[[i]])==0){
      AllLoess[[i]] = RTl
    }
  }
}

for(i in 1:length(AllLoess)){
  write.table(AllLoess[[i]][complete.cases(AllLoess[[i]]),],
    gsub(".bg" , "qnorm.loess.bg" , colnames(AllLoess[[i]])) [4],
    sep= "\t" , row.names=FALSE, quote=FALSE, col.names = FALSE
  )
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"
if (file.exists(RLogFile) == FALSE) {
    sink(RLogFile)
    a <- sessionInfo()
    print(a)
    sink()
}

################################################
################################################
################################################
##############################################
