#!/usr/bin/env Rscript

#######################
#Cross Validation step
#######################

# This script permits to compute the p-values from the adjusted UMI count
# After running the script generates a matrix file withe the p-values for each Fusion and sample

# Format example of the count file:
#     Fusion;Count_Full;Count_UMI;Adj_UMI
#     BRCA2E22G-BRCA2E23D;2019;669;522
#     BRCA2E2G-BRCA2E4D;62;20;17
#     BRCA1E22G-BRCA1E24D;10;5;2
#     BRCA1E20G-BRCA1E24D;8;3;2
#     BRCA1E12G-BRCA1E3D;1;1;0
#     BRCA1E2G-BRCA1E6D;5;3;2
#     BRCA1ES1613WTG-BRCA1ES1613GD;238;83;56
#     BRCA2E22G-BRCA2E24D;5;3;2
#     BRCA1E20G-BRCA1E15D;7;2;2

messageHelp = "Usage: appTrainedModel.r\n
    Mandatory \n
        -I, --input /path/to/countFileDirectory\t\tList of count file (.count)
		-O, --output /path/to/outputFile.txt\t\tOutput file tab-delimated (.txt)\n
    Options \n
        -s, --nSamp Integer\t\tNumber of ramdom samples used to perform cross validation [Default=10]
        -i, --nIter Integer\t\tNumber of iteration perfom during cross validation [Default=10]
        -p, --plot /path/to/plot.pdf\t\tPrint plot of count distribution
    -h, --help\t\tPrint this help message and exit\n

   You could : Rscript appTrainedModel.r -I /path/to/countFileDirectory/ -O /path/to/output.txt"

argsFull <- commandArgs()
Rscript <- argsFull[1]

scriptPath=dirname(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])))
args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]

if (length(which(argsFull=="--args"))==0){message(messageHelp);q(save = "no")}
if (length(args)<3){message(messageHelp);stop("Not enought arguments")}

i=1
nIter = 10
nSamp = 10
plotPath = NULL

while (i < length(args)){
    if(args[i]=="-I"|args[i]=="--input"){
        pathCount=args[i+1];i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        pathOutput=args[i+1];i = i+2
    }else if(args[i]=="-s"|args[i]=="--nSamp"){
        nSamp=as.numeric(args[i+1]);i = i+2
    }else if(args[i]=="-i"|args[i]=="--nIter"){
        nIter=as.numeric(args[i+1]);i = i+2
    }else if(args[i]=="-p"|args[i]=="--plot"){
        plotPath=args[i+1];i = i+2
    }else{
        message(paste("********Unknown option:",args[i],"\n"));message(messageHelp);stop()
    }
}

countFiles = list.files(pathCount, pattern=".count$", full.names=TRUE)
if(length(countFiles)<1){message(paste("No count file detected\n"));message(messageHelp);stop()}
if(length(countFiles)<=nSamp){message(paste("Not enought count file for cross-validation\nAt least ",nSamp," samples required\n"));message(messageHelp);stop()}

negbinom.n=10
separator = ";"
matrixCount = NULL

nbTimes <- function(x){length(x[x!=0 & !is.na(x)])}

removeNA <- function(x){
	if(length(x)!=length(x[!is.na(x)])){
		x[is.na(x)]=0
	}
	return(as.numeric(x))
}

getMatrixFromCountFile <- function(countFiles){
    for (i in countFiles ){
        count = read.table(i,sep=separator,header=T)
    	count = count[,c(1,3)]
    	colnames(count) = c("Fusion",basename(i))
        if(is.null(matrixCount)){
            matrixCount=count
        }else{
            matrixCount <- merge(matrixCount, count, by.x="Fusion", by.y="Fusion", all.x=TRUE, all.y=TRUE)
        }
    }
    nbSamp = apply(matrixCount[,2:ncol(matrixCount)],1,nbTimes)
    matrixCount = matrixCount[nbSamp>0,]
    tmp = as.data.frame(apply(matrixCount[,2:ncol(matrixCount)],2,removeNA))
    matrixCount = cbind(Fusion = matrixCount[,1],tmp)
    return(matrixCount)
}

#adjusted count
getMeanGeom <- function(x){
	m = length(x)
	gm = (prod(x[x!=0]))^(1/m)
	return(gm)
}

getAdjustedCount <- function(matrixCount){
    gm = apply(matrixCount[,2:ncol(matrixCount)],1,getMeanGeom)

    tmp= matrixCount[,2:ncol(matrixCount)]*(1/gm)
    matrixCountAdj = matrixCount
    for (i in 1:ncol(tmp)){
    	x = tmp[,i]
    	m = median(x[x!=0])
    	matrixCountAdj[,i+1] = matrixCountAdj[,i+1]/m
    }
    return(matrixCountAdj)
}

# functions to fit model

count<-function(x, values){
  res<-rep(NA,length(values))
  for(i in 1:length(values)) res[i]<-sum(x==values[i])
  return(res)
}

format.p<-function(p){
  if(p<1e-16) return("<1e-16")
  if(p<1e-6) return("<1e-6")
  if(p<1e-3) return("<1e-3")
  return(format(p,digits=3,nsmall=3))
}

cuter<-function (x, cuts = NULL, g = 2, mod = "[[", minmax = FALSE){
  mod <- match.arg(mod, c("[[", "]]"))
  if (is.null(cuts)) {
    cuts <- quantile(x, seq(0, 1, 1/g), na.rm = T)
    cuts <- cuts[-c(1, length(cuts))]
  }
  res <- x
  ncut <- length(cuts)
  levels <- rep(NA, ncut + 1)
  if (mod == "[[") {
    res[x < cuts[1]] <- 1
    res[x >= cuts[ncut]] <- ncut + 1
    levels[1] <- paste("<", cuts[1], sep = "")
    levels[ncut + 1] <- paste(">=", cuts[ncut], sep = "")
    if (minmax) {
      m1 <- min(x)
      m2 <- max(x)
      if (m1 < cuts[1]) {
        levels[1] <- paste("[", format.p(m1), ",", cuts[1],
                           "[", sep = "")
      }
      if (m2 >= cuts[ncut]) {
        levels[ncut + 1] <- paste("[", cuts[ncut], ",",
                                  format.p(m2), "]", sep = "")
      }
    }
    if (ncut > 1) {
      for (i in 1:(ncut - 1)) {
        res[x >= cuts[i] & x < cuts[i + 1]] <- i + 1
        levels[i + 1] <- paste("[", cuts[i], ",", cuts[i +
                                                         1], "[", sep = "")
      }
    }
  }
  else if (mod == "]]") {
    res[x <= cuts[1]] <- 1
    res[x > cuts[ncut]] <- ncut + 1
    levels[1] <- paste("<=", cuts[1], sep = "")
    levels[ncut + 1] <- paste(">", cuts[ncut], sep = "")
    if (minmax) {
      m1 <- min(x)
      m2 <- max(x)
      if (m1 <= cuts[1]) {
        levels[1] <- paste("[", format.p(m1), ",", cuts[1],
                           "]", sep = "")
      }
      if (m2 > cuts[ncut]) {
        levels[ncut + 1] <- paste("]", cuts[ncut], ",",
                                  format.p(m2), "]", sep = "")
      }
    }
    if (ncut > 1) {
      for (i in 1:(ncut - 1)) {
        res[x > cuts[i] & x <= cuts[i + 1]] <- i + 1
        levels[i + 1] <- paste("]", cuts[i], ",", cuts[i +
                                                         1], "]", sep = "")
      }
    }
  }
  res <- factor(res, levels = as.character(1:(ncut + 1)))
  levels(res) <- levels
  return(res)
}

vh.like <- function(pattern, x){
    w<-regexpr(pattern,x)
    w<-which(na(w==1) & na(nchar(x)==attr(w,"match.length")))
    res<-x ; res[]<-NA ; class(res)<-"logical" ; res[!is.na(x)]<-FALSE ; res[w]<-TRUE
    return(res)
  }

na<-function(x, to = FALSE){
    x[is.na(x)]<-to
    x
}

#Function to transform a table jid*sid in sid*jid
# a : initial table
# sid : samples ID (names of columns from a)
# jid.column : name of column with junctions ID

js2sj<-function(a, sid, jid.column,nbSamp){
  if(!all(c(jid.column,sid)%in%names(a))){
    stop("***** The 'jid' and 'sid' don't find in the columns of table")
  }else if(length(unique(a[,jid.column]))!=nrow(a)){
    stop("***** The junctions were not identify in unique way")
  }
  a = a[nbSamp>=5,]
  data = as.data.frame (t(a[,sid]))
  names(data) = a[,jid.column]
  row.names(data) = sid
  return(data)
}

# Function to create Gamma / NegBinomial models from dataset
# data : table of data sid*jid
# jid : Junctions ID to modelisate
# negbinom.n : the min number of classes to create for negative binomial ajustment )

fit.gamma.negbinomial<-function(x, negbinom.n=10){
    model<-NULL
    xp<-x
    bx<-boxplot.stats(x)
    nto<-0
    if(length(bx$out>0)){
      nto<-sum(x%in%bx$out)
      x<-x[!x%in%bx$out]
    }
	if(length(x[x!=0 & !is.na(x)])>2){
      # Modelisation
      m<-mean(x,na.rm=T)
      s<-sd(x,na.rm=T)

      # Assay gamma distribution
      scale<-s**2/m
      shape<-m/scale
      tryCatch({
          test<-ks.test(x, pgamma, shape=shape, scale=scale)$p.value
      },
      error=function(cond) {
          test = 0
      })
      if(is.na(test)) test=0
      if(test>=0.01){ # acceptable adequation to the gamma distribution
        model <- 1
        Param1 <- scale
        Param2 <- shape
        mean <- m
        sd <- s
      }

      # Assay negative binomiale distribution with intervalls of negbinom.n of same size
      if(is.null(model)){
          step<-max(x)/negbinom.n
          y<-as.numeric(cuter(x, cuts=c(seq(step, max(x), step)), mod="]]"))-1

          my<-mean(y, na.rm=T)
          vy<-var(y, na.rm=T)
          prob<-my/vy  ;  size<-my*prob/(1-prob)

          obs<-count(y, 0:max(y))
          th<-dnbinom(0:(max(y)-1), size=size, prob=prob) ; th<-c(th,1-sum(th))
          tryCatch({
            test<-chisq.test(obs, p=th)$p.value
          },
          error=function(cond) {
            test = 0
          })
          if(is.na(test)) test=0
          if(test>=0.01){ # acceptable adequation to the negative binomiale distribution
            model <- 2
            Param1 <- prob
            Param2 <- size
            mean <- m
            sd <- s
          }
      }
    }
    if(is.null(model)){mean<-sd<-model<-Param1<-Param2<-NA}
  result <<- c(model,Param1,Param2,mean,sd)
}

getCrossValResult <- function(x){
    model = x[seq(from=1,to=length(x),by=5)]
    Param1 = x[seq(from=2,to=length(x),by=5)]
    Param2 = x[seq(from=3,to=length(x),by=5)]
    mean = x[seq(from=4,to=length(x),by=5)]
    sd = x[seq(from=5,to=length(x),by=5)]
    # get final model
    nNA = length(model[is.na(model)])
    nGamma = length(model[model==1 & !is.na(model)])
    nNegBinom = length(model[model==2 & !is.na(model)])
    if(nNA == length(model)){
        modelFinal = NA
    }else if(nGamma>=nNegBinom){
        modelFinal = 1
    }else{
        modelFinal = 2
    }
    #print(as.numeric(Param1[model==modelFinal]))
    # get final parameters
    Param1f = median(Param1[model==modelFinal],na.rm = TRUE)
    Param2f = median(Param2[model==modelFinal],na.rm = TRUE)
    meanf = median(mean[model==modelFinal],na.rm = TRUE)
    sdf = median(sd[model==modelFinal],na.rm = TRUE)
    result <<- c(modelFinal,Param1f,Param2f,meanf,sdf)
}

#Function to apply the previously model on a serie of values
# model : is the model got by 'fit.gamma.negbinomial'
# x : a serie of values
# jid : junction ID designating the modele to apply

predict.gamma.negbinomial<-function(m, x){
  if(!is.na(m$model)){
    if(m$model=="Gamma"){
	  val<-1-abs(pgamma(x, scale=m$Param1, shape=m$Param2, lower.tail=T)-pgamma(x, scale=m$Param1, shape=m$Param2, lower.tail=F))
      val[na(val==0)]<-1e-324
    }else if(substr(m["model"],1,17)=="Negative binomial"){
      val<-1-abs(pnbinom(x, size=m$Param2, prob=m$Param1, lower.tail=T)-pnbinom(x, size=m$Param2, prob=m$Param1, lower.tail=F))
      val[na(val==0)]<-1e-324
    }
  }else{
    val<-rep(NA,length(x))
  }
  val[x==0]=1
  val<-p.adjust(val, "hochberg", n=nb.tests)
  return(val)
}

mergedResult = NULL
message("   Start Cross validation...")
for (i in 1:nIter){
    matrixCount = NULL
    message(paste(i,"/",nIter,"iteration\n"))
    # Sample selection
    countFilesRand = sample(countFiles,nSamp,replace=FALSE)
    matrixCount = getMatrixFromCountFile(countFilesRand)
    matrixCountAdj = getAdjustedCount(matrixCount)

    # Create table of data analysis
    sid<-names(matrixCountAdj)[2:ncol(matrixCountAdj)]
    jid.column<-names(matrixCountAdj)[1]
    a<-matrixCountAdj
    nbSamp = apply(matrixCountAdj[,2:ncol(matrixCountAdj)],1,nbTimes)
    data<-js2sj(a, sid, jid.column,nbSamp)
    jid<-names(data)

    # Model Fitting
    result <- apply(data,2,fit.gamma.negbinomial,negbinom.n=negbinom.n)
    result <- as.data.frame(t(result))
    names(result) = c("model","Param1","Param2","mean","sd")
    if(is.null(mergedResult)){
        mergedResult = result
    }else{
        mergedResult = merge(mergedResult,result,by = 0,all=TRUE,suffixes = paste0(".",c(i-1,i)))
        row.names(mergedResult)=mergedResult[,"Row.names"]
        mergedResult = mergedResult[,-1]
    }
}
message("   Finish Cross validation.")
CrossResult <- as.data.frame(t(apply(mergedResult,1,getCrossValResult)))
row.names(CrossResult) = row.names(mergedResult)
names(CrossResult) = c("model","Param1","Param2","mean","sd")
CrossResult$model[CrossResult$model==1 & !is.na(CrossResult$model)] = "Gamma"
CrossResult$model[CrossResult$model==2 & !is.na(CrossResult$model)] = paste("Negative binomial ",negbinom.n,"c",sep="")

write.table(cbind(row.names(CrossResult),CrossResult),sub(".txt","_model.txt",pathOutput),sep="\t",dec=".",quote=FALSE,row.names=FALSE,na="")

if(!is.null(plotPath)){
    pdf(plotPath)
    for(j in row.names(CrossResult)){
    	m=CrossResult[j,]
        if(is.na(m$model)){
        }else if(m$model=="Gamma"){
    		x = c(0:(m$mean*5))
    		p = dgamma(x,scale=m$Param1, shape=m$Param2)
    		plot(p~x,type="l",col="red",lwd=2,main=paste(m$model,j,sep="_"), sub = paste("scale =",round(m$Param1,3),", shape =", round(m$Param2,3)))
    	}else if(substr(m$model,1,17)=="Negative binomial"){
    		x = c(0:(m$mean*5))
    		p = dnbinom(x,size=m$Param2, prob=m$Param1)
    		plot(p~x,type="l",col="red",lwd=2,main=paste(m$model,j,sep="_"), sub = paste("prob =",round(m$Param1,3),", size =", round(m$Param2,3)))
    	}
    }
    dev.off()
}

#########################
# 4 - Application of model on the dataset
message("   Model apliccation...")
matrixCount = NULL
matrixCount = getMatrixFromCountFile(countFiles)
matrixCountAdj = getAdjustedCount(matrixCount)
write.table(matrixCountAdj, file = sub(".txt","_AdjCount.txt",pathOutput), quote = FALSE, sep = "\t", na = "", dec = ".", row.names = FALSE,)

# Create table of data analysis
sid<-names(matrixCountAdj)[2:ncol(matrixCountAdj)]
jid.column<-names(matrixCountAdj)[1]
nbSamp = apply(matrixCountAdj[,2:ncol(matrixCountAdj)],1,nbTimes)
newdata<-js2sj(matrixCountAdj, sid, jid.column,nbSamp)

# Count of the number of tests carried out
testable = row.names(CrossResult[!is.na(CrossResult$model),])
nb.tests <- sum(apply(newdata[,which(names(newdata)%in%testable)],2,function(x) sum(!is.na(x))))

matrixPval = as.data.frame(matrix(NA,nrow(matrixCountAdj),ncol(matrixCountAdj)-1))
colnames(matrixPval) = colnames(matrixCountAdj)[-1]
matrixPval = cbind(Fusion = matrixCountAdj[,1],matrixPval)

for(j in row.names(CrossResult)){
    m<-CrossResult[j,]
    val<-predict.gamma.negbinomial(m = m, x=newdata[,j])
	matrixPval[matrixPval[,1]==j,2:ncol(matrixPval)] = val
}

write.table(matrixPval, file = pathOutput, quote = FALSE, sep = "\t", na = "", dec = ".", row.names = FALSE,)
