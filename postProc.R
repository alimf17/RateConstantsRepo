library(rstan)
library(bayestestR)

args = commandArgs(trailingOnly=T)

# /anvil/scratch/x-farhat/bufferRates/ZiyuanRun1epe1D_final_parameters.rds

path = args[1]
name = args[2]

numChains = as.integer(args[3])

print(paste(path, name,"A", "_final_parameters.rds", sep = ""))

sampName =paste(path, name,"A", "_final_parameters.rds", sep = "")
sample = readRDS(sampName)

rates = array(dim = c(dim(sample[[1]]$rates), length(sample), numChains))


for(k in 1:numChains){

    chainName = paste(path, name, LETTERS[k],"_final_parameters.rds", sep = "")
    chain = readRDS(chainName)

    for(m in 1:length(chain)){
        rates[,,m,k] = chain[[m]]$rates
    }
}

rhats = array(dim = dim(sample[[1]]$rates))

#print(rates[1,2,,])

L = length(sample)

for(i in 1:nrow(rhats)){
    for(j in 1:ncol(rhats)){
        rhats[i,j] = Rhat(rates[i,j,,])
    }
}

print(rhats)

print(apply(rates[,,(floor(L/2):L),], c(1,2), mean))
M = (apply(rates[,,(floor(L/2):L),], c(1,2), ci, method = "HDI", ci = 0.95))
print(dim(M))

low = array(dim = dim(M))
high = array(dim = dim(M))

for(i in 1:nrow(low)){
    for(j in 1:ncol(low)){
        low[i,j] = M[i,j][[1]][1,2]
        high[i,j] = M[i,j][[1]][1,3]
    }
}

print(low)
print(high)

