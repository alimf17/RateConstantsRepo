library(doParallel)
library(foreach)

source("Eq.R")
source("Inference.R")

args = commandArgs(trailingOnly=TRUE)

nodes = as.numeric(args[1])

corespernode = as.numeric(args[2])

prefix = args[3]
print(paste("This chain is ", prefix,".rds",sep = ""))
print(prefix)

superexperiment = readRDS(args[4])

if(is.list(superexperiment)){
    experiment = superexperiment$data
    equilibria = superexperiment$eqs
} else{
    experiment = superexperiment
    equilibria = calcEQs(experiment)
}



steps = as.numeric(args[5])

output = args[6]

if(length(args)>6){
    adjacentOnly = as.logical(args[7])
} else{
    adjacentOnly = FALSE
}

if(is.na(adjacentOnly)){
    
    adjacentOnly = FALSE
    args[8] = args[7]

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

if(length(args)>7){
    pars = readRDS(args[8])
    pars.init = pars[[length(pars)]]$rates
}


cores = corespernode*nodes

sink(paste(output,prefix,".out", sep = ""))

if(exists("pars.init")){

    print("program has an initial condition")

    print(pars.init)

    main(experiment,equilibria = equilibria, total.iter = steps,experiments.per.iter = 10000,cores = corespernode*nodes,name = paste(output,prefix,sep = ''),init.rates=pars.init, kill.rates = kill.rates)
} else{

    main(experiment,equilibria = equilibria, total.iter = steps,experiments.per.iter = 10000,cores = corespernode*nodes,name = paste(output,prefix, sep =''), kill.rates = kill.rates)
}
Rprof(NULL)

print(paste("Time elapsed for run:",end-begin))

sink()
