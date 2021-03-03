library(doParallel)
library(foreach)

source("Eq.R")
source("InferenceEqMag20201117.R")

args = commandArgs(trailingOnly=TRUE)

print(args)

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


print("1")

steps = as.numeric(args[5])

output = args[6]

if(length(args)>6){
    adjacentOnly = as.logical(args[7])
} else{
    adjacentOnly = FALSE
}



print("2")

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
    pars.init = pars[[length(pars)]]$eqmag
    print("first pars init")
    print(args[8])
    print(pars.init)
}

print("should pars.init")

cores = corespernode*nodes

sink(paste(output,prefix,".out", sep = ""))
begin = Sys.time()

if(exists("pars.init")){
#    step = NULL
#    a = list.files(path = output, pattern = paste(prefix,"_[[:digit:]]*parameters.rds", sep = ""))
#    print(paste(prefix,"_[[:digit:]]*parameters.rds", sep = ""))
#    print(a)
#    if(length(a)>0){
#        last.file = a[length(a)]
#        print(last.file)
#        step = length(readRDS(paste(output,last.file, sep = "")))
#        if(grepl("\\d+",last.file, perl = TRUE)){
#            step = step + as.numeric(gsub(".*?([0-9]+).*", "\\1",last.file)) - 1
#        }
#    }

    print("program thinks it has an init")

    print(pars.init)

    main(experiment,equilibria = equilibria, total.iter = steps,experiments.per.iter = 10000,cores = corespernode*nodes,name = paste(output,prefix,sep = ''),init.eqmag=pars.init, kill.rates = kill.rates)
} else{

    main(experiment,equilibria = equilibria, total.iter = steps,experiments.per.iter = 10000,cores = corespernode*nodes,name = paste(output,prefix, sep =''), kill.rates = kill.rates)
}
Rprof(NULL)
end = Sys.time()
print(paste("Time elapsed for run:",end-begin))

sink()
