#!/usr/bin/env Rscript

#######################
#APPLICATION STEP
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
        -m, --model /path/to/Model_XXXX.RData\t\tThe Model file (.RData)
		-O, --output /path/to/outputFile.txt\t\tOutput file tab-delimated (.txt)\n
    -h, --help\t\tPrint this help message and exit\n

   You could : Rscript appTrainedModel.r -I /path/to/countFileDirectory/ -m path/to/Model_XXXX.RData -O /path/to/output.txt"

argsFull <- commandArgs()
Rscript <- argsFull[1]

scriptPath=dirname(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])))
args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]

if (length(which(argsFull=="--args"))==0){message(messageHelp);q(save = "no")}
if (length(args)<4){message(messageHelp);stop("Not enought arguments")}

i=1
while (i < length(args)){
    if(args[i]=="-I"|args[i]=="--input"){
        pathCount=args[i+1];i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        pathOutput=args[i+1];i = i+2
    }else if(args[i]=="-m"|args[i]=="--model"){
        pathModel=args[i+1];i = i+2
	}else{
        message(paste("********Unknown option:",args[i],"\n"));message(messageHelp);stop()
    }
}
   
# pathModel = "//s-grp/grp/Biomol/ARN/raphael/analyse_RNAseq/RT-MLPA/BRCA2/Train/Model_1581515489.RData"
# pathCount = "//s-grp/grp/Biomol/ARN/raphael/analyse_RNAseq/RT-MLPA/BRCA2/App/"
# pathOutput = "//s-grp/grp/Biomol/ARN/raphael/analyse_RNAseq/RT-MLPA/BRCA2/result.txt"

countFiles = list.files(pathCount, pattern=".count$", full.names=TRUE)
if(length(countFiles)<1){message(paste("No count file detected\n"));message(messageHelp);stop()}
negbinom.n=10
separator = ";"
matrixCount = NULL
load(pathModel)

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
nbTimes <- function(x){length(x[x!=0 & !is.na(x)])}
nbSamp = apply(matrixCount[,2:ncol(matrixCount)],1,nbTimes)
matrixCount = matrixCount[nbSamp>0,]

removeNA <- function(x){
	if(length(x)!=length(x[!is.na(x)])){
		x[is.na(x)]=0
	}
	return(as.numeric(x))
}

tmp = as.data.frame(apply(matrixCount[,2:ncol(matrixCount)],2,removeNA))
matrixCount = cbind(Fusion = matrixCount[,1],tmp)

#adjusted count
getMeanGeom = function(x){
	m = length(x)
	gm = (prod(x[x!=0]))^(1/m)
	return(gm)
}
gm = apply(matrixCount[,2:ncol(matrixCount)],1,getMeanGeom)

tmp= matrixCount[,2:ncol(matrixCount)]*(1/gm)
matrixCountAdj = matrixCount
for (i in 1:ncol(tmp)){
	x = tmp[,i]
	m = median(x[x!=0])
	matrixCountAdj[,i+1] = matrixCountAdj[,i+1]/m
}

count<-function(x, values){
  res<-rep(NA,length(values))
  for(i in 1:length(values)) res[i]<-sum(x==values[i])
  return(res)
}

msd<-function(x){
  m<-mean(x,na.rm=T) ; sd<-sd(x,na.rm=T)
  if(m<1e-3){ m<-"<1e-3" }else{ m<-format(m,digits=3,nsmall=3) }
  if(sd<1e-3){ sd<-"<1e-3" }else{ sd<-format(sd,digits=3,nsmall=3) }
  return(paste(m,sd,sep=";"))
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

js2sj<-function(a, sid, jid.column){
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

sid<-names(matrixCountAdj)[2:ncol(matrixCountAdj)]
jid.column<-names(matrixCountAdj)[1]
a<-matrixCountAdj

# 2 - Fromating to an input ind*evenement : data
message("   Creating matrix...")
# Create table of data analysis
data<-js2sj(a, sid, jid.column)
jid<-names(data)

#Function to apply the previously model on a serie of values
# model : is the model got by 'fit.gamma.negbinomial'
# x : a serie of values
# jid : junction ID designating the modele to apply

predict.gamma.negbinomial<-function(m, x, jid){
  if(!is.null(m)){
    if(m$model=="Gamma"){
	  val<-1-abs(pgamma(x, shape=m$shape, scale=m$scale, lower.tail=T)-pgamma(x, shape=m$shape, scale=m$scale, lower.tail=F))
      val[na(val==0)]<-1e-324
    }else if(substr(m$model,1,17)=="Negative binomial"){
      val<-1-abs(pnbinom(x, size=m$size, prob=m$prob, lower.tail=T)-pnbinom(x, size=m$size, prob=m$prob, lower.tail=F))
      val[na(val==0)]<-1e-324
    }
  }else{
    val<-rep(NA,length(x))
  }
  val[x==0]=1
  return(val)
}

#########################
# 4 - Application of model on the dataset
message("   Model apliccation...")

newdata<-data # forme sid*jid

# Count of the number of tests carried out
nb.tests<-0
for(j in names(model)){
    if(j%in%names(newdata) & vh.like("(Gamma)|(Negative binomial.*)",model[[j]][["model"]])){
        nb.tests<-nb.tests+sum(!is.na(newdata[,j]))
    }
}

matrixPval = as.data.frame(matrix(NA,nrow(matrixCountAdj),ncol(matrixCountAdj)-1))
colnames(matrixPval) = colnames(matrixCountAdj)[-1]
matrixPval = cbind(Fusion = matrixCountAdj[,1],matrixPval)

for(j in names(model)){
  m<-model[[j]]
  if(j%in%names(newdata) & vh.like("(Gamma)|(Negative binomial.*)",m$model)){
    val<-predict.gamma.negbinomial(m = m, x=newdata[,j], jid=j)
    val<-p.adjust(val, "hochberg", n=nb.tests)
	matrixPval[matrixPval[,1]==j,2:ncol(matrixPval)] = val
  }
}

write.table(matrixPval, file = pathOutput, quote = FALSE, sep = "\t", na = "", dec = ".", row.names = FALSE,)

