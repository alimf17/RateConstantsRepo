library(doParallel)
library(foreach)

source("Eq.R")
source("Inference.R")

#Usage: Rscript runInference.R numberOfCores nameOfRun rdsFileArrayOfExperimentalData numberOfStepsForRun OutputFilesLocation OptionalBooleanForOnlyAdjacentTransitions OptionalInitialStateForMCMC

args = commandArgs(trailingOnly=TRUE)

cores = args[1]

prefix = args[2]

print(paste("This chain is ", prefix,".rds",sep = ""))
print(prefix)


experiment = readRDS(args[3])
equilibria = calcEQs(experiment)

steps = as.numeric(args[4])

output = args[5]

if(length(args)>5){
    adjacentOnly = as.logical(args[6])
} else{
    adjacentOnly = FALSE
}

if(is.na(adjacentOnly)){
    
    adjacentOnly = FALSE
    args[7] = args[6]

}

kill.rates = 1-diag(dim(experiment)[1])    

if(adjacentOnly){
    for(i in 1:nrow(kill.rates)){
        for(j in i:ncol(kill.rates)){
            if(abs(j-i)>1){
                kill.rates[i,j] =0
                kill.rates[j,i] = 0
            }
        }
    }    
}

if(length(args)>6){
    pars = readRDS(args[7])
    pars.init = pars[[length(pars)]]$rates
}



sink(paste(output,prefix,".out", sep = ""))

if(exists("pars.init")){

    print("program has an initial condition")

    print(pars.init)

    main(experiment,equilibria = equilibria, total.iter = steps,experiments.per.iter = 10000,cores = cores,name = paste(output,prefix,sep = ''),init.rates=pars.init, kill.rates = kill.rates)
} else{

    print("kill")
    print(kill.rates)
    main(experiment,equilibria = equilibria, total.iter = steps,experiments.per.iter = 10000,cores = cores,name = paste(output,prefix, sep =''), kill.rates = kill.rates)
}
Rprof(NULL)

print(paste("Time elapsed for run:",end-begin))

sink()
