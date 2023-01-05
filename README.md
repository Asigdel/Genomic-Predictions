# Genomic-Predictions
# Running Parallel jobs for Genomic Predictions
# Author: Anil Sigdel
# Jan 5, 2023

## Codes

> run.R.jobs.paral.txt

#!/bin/bash

SEED=(123 1001 1010 1100 2001 101 313 3002 231 4440)            #10 replicates
FOLD=(1 2 3 4 5 6 7 8 9 10)                                                                                     #10 folds crossvalidation

#SEED=(101)
#FOLD=(6)

for i in "${SEED[@]}"
do
        for j in "${FOLD[@]}"
        do
                sbatch R.jobs.paral.txt $i $j   #submit i replicates times j jobs
                sleep 1                                             #wait 1 secod to submit next job
        done
done


#############################################################################################################################################
> R.jobs.paral.txt

#!/bin/bash
#SBATCH --job-name=paralG_Holstein_Genic
#SBATCH --mail-user=hendyelpacheco@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50gb                      #memory per each job
#SBATCH --time=90:00:00                                 #time per each job
#SBATCH --qos=penagaricano-b

module load R

Rscript Rcode.Full.Pred $1 $2



###########################################################################################################################################

> Rcode.Full.Pred
#Prediction code using BGLR (under paralel job submission)
#RKHS linear kernel to perform GBLUP; RANDOM sampling to get TRAIN and TEST subsets; WITHOUT BTx and BTy
#July,2019

args = commandArgs(trailingOnly = TRUE)
SEED = as.numeric(args[1]); print(SEED)
FOLD = as.numeric(args[2]); print(FOLD)

library(BGLR)

### Read data (pheno + geno)
data = read.table("DataGen.txt", colClasses=c(rep("factor",4),rep("numeric",3),"factor",rep("numeric",290989)), header=TRUE)            

data = data[order(data$REL),]           #ATTENTION animal order must match with the order of the G matrix built on Rcode.Pred.methods.compare

### Read phenotypic data
Y = data [,1:8]
y = data[,5]                                            #SCR records
rel = data[,6]                                          #reliability

### Read genotypic data
X = as.matrix(data[,-c(1:8)])

### Number of markers
p = ncol(X)

### Number of animals
n = nrow(data)

### Computing K relationship matrix assuming Gaussian kernel
#S = scale(X, center =  TRUE, scale = TRUE)                     #S = matrix of centered in zero and standarized SNPs (not mandatory!)
#load(file = "S.rda")
#D = (as.matrix(dist(S,method='euclidean'))^2/p)
#h = 0.5
#K = exp(-h*D)
#save (K, file = 'K.rda')

load(file = "/ufrc/penagaricano/hendyelpacheco/BTAXY/Predictions/Gaussian/PredG_Auto/K.rda")


#Create equally size n folds for crossvalidation
nfolds = 10

#Random sampling data
set.seed(SEED)
sample = sample (1:n,replace=FALSE)                                                     #sample at random from 1 to number of animals
list = split(sample, sample(rep(1:nfolds), replace=FALSE))              #split at random samples into n folds
str(list)

tst = as.vector(unlist(list[FOLD], use.names=FALSE))            #convert list into vector and save the elements of each list
head(tst)

### Creating a testing set
yNA = y
yNA[tst] = NA                                                                                   #assing NA to rows that were listed in tst

### Setting predictors
ETA = list( list (~ factor(Evaluation), data = Y, model = 'FIXED'),
                   list (K = K, model = 'RKHS')                                 #Gaussian Kernel with scaled X
                   )

### Fitting the model
fm = BGLR (y=yNA, ETA = ETA, nIter = 100000, burnIn = 30000, thin = 5,  weights = rel, saveAt = 'RKHS_')
#save (fm, file = 'fm.rda')
str(fm)

### Assesment of predictive ability in TRN and TST data sets
corP_TRN = cor(fm$yHat[-tst],y[-tst])
corP_TST = cor(fm$yHat[tst],y[tst])
corS_TRN = cor(fm$yHat[-tst],y[-tst], method ="spearman")
corS_TST = cor(fm$yHat[tst],y[tst], method ="spearman")
MSE_TRN = mean((y[-tst]-fm$yHat[-tst])^2)
MSE_TST = mean((y[tst]-fm$yHat[tst])^2)
#fm$fit

result = cbind(corP_TRN,corP_TST,corS_TRN,corS_TST,MSE_TRN,MSE_TST)
#yhat = fm$yHat
file = paste(paste("Out",SEED,sep=""),FOLD,sep=".")
save(result,file=file)

#####################################################################################################################################################

> Rcode_Compile_AllOutputs.txt


#SUBMITTED 100 jobs (10 fold * 10 replicates) at the same time: run.R.jobs.paral.txt, R.jobs.paral.txt and Rcode.Full.Pred
                                #bash run.R.jobs.paral.txt (define specific seed, possible to reproduce the same results!)
                                #OutSEED.FOLD is a R file (only opens in R using load("Out")

#COMPILE ALL 100 jobs results into one file
files = list.files(pattern = "Out")
Outs = array()

for (i in 1:length(files)){
load(files[i]);
Outs = rbind(Outs,result)                              #result is the object inside the Out*
}
Outs=Outs[-1,]
colMeans(Outs)

write.table(Outs, file="Outs_LinearKernel", sep=" ",row.names=FALSE,col.names=TRUE,quote=FALSE)

######################################################################################################################################################








