library(doParallel)
library(foreach)

source("Eq.R")
source("DivorcedSpecialInference20201225.R")

args = commandArgs(trailingOnly = T)

track = as.integer(args[1])-1

notate = numeric(4)

global.assign(readRDS("wild_type_transition_matrices.rds"), calcEQs(readRDS("wild_type_transition_matrices.rds")))

#### This block of code defines our master grid
P = c(0.51, 0.52, 0.53)
bigP = matrix(nrow = length(P), ncol = 4)
bigP[,1] = P
bigP[,2] = c(0.50, 0.505, 0.515)
bigP[,3] = c(0.525, 0.535, 0.54)
bigP[,4] = c(0.545, 0.55, 0.555)

Q = c(0.35, 0.40, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75)
bigQ = matrix(nrow = length(Q), ncol = 5)
for(i in 1:5){
    bigQ[,i] = Q+(i-1)*0.01
}

S = c(0.1, 0.15, 0.2, 0.25, 0.3)
bigS = matrix(nrow = length(S), ncol = 12)
for(i in 1:12){
        bigS[,i] = S+(i-1)*0.25
}

R1 = c(0.5, 1.0, 1.5, 2.0, 2.5, 5,10)
bigR =  matrix(nrow = length(R1), ncol = 5)
bigR[,1] = R1
bigR[,2] = c(R1[1:4]+2.5, 5.5, 6.0, 6.5)
bigR[,3] = c(7+(0:5)/2, 10.5)
bigR[,4] = 11+(0:6)/2
bigR[,5] = 0.15+(0:6)*0.05
####

#### This block of code converts our index (given by argument) to select which
#### elements of the grid will have their likelihood calculated   
notate[1] = (track/(ncol(bigQ)*ncol(bigS)*ncol(bigR)))-(track/(ncol(bigQ)*ncol(bigS)*ncol(bigR)))%%1 
track = track-notate[1]*(ncol(bigQ)*ncol(bigS)*ncol(bigR))
notate[2] = (track/(ncol(bigS)*ncol(bigR)))-(track/(ncol(bigS)*ncol(bigR)))%%1
track = track-notate[2]*(ncol(bigS)*ncol(bigR))
notate[3] = (track/(ncol(bigR)))-(track/(ncol(bigR)))%%1
track = track-notate[3]*(ncol(bigR))
notate[4] = track
inds = notate+1

P = bigP[,inds[1]]
Q = bigQ[,inds[2]]
S = bigS[,inds[3]]
R1 = bigR[,inds[4]]


#### This actually runs the likelihood calculations
cores = 28
cl <- makeForkCluster(cores) 
registerDoParallel(cores)

for(p in P){
    for(q in Q){
        for(s in S){
            for(r in R1){

                pars = matrix(c(p, q, s, r), nrow = 2)
                #We ran this likelihood with far fewer steps and linear splines
                likelihood = log.likelihood.main(pars, 30000)
                line = paste(p, q, s, r, likelihood, sep = '\t')

                #We originally ran this to output no matter what the likelihood was. 
                #We discovered empirically that the parameter space where ln likelihood
                #was finite had no holes, so this is a quality of life improvement
                if(likelihood > -Inf){
                    #We outputted all data to a single file. Strongly recommend that
                    #you do NOT do this: combined with our faster runs, outputs
                    #interfered with each other somewhat. There was enough to recover
                    #a fit with spline, but no need to cause hassle
                    write(line,file='Landscape.out',append=TRUE)
                }
            }
        }
    }
}

stopCluster(c1)

####



#We calculated the linear splines using the package polspline and the code:
  
#model = polymars(responses = <array of ln likelihoods of finite ln likelihood parameters>, predictors = <Nx4 matrix of parameters, where one row is one parameter set>, gcv = <integer between 1-500. No significant difference, which gave us confidence in our model>)

#Slices of the likelihood landscape can be extracted with
#predict.polymars(model, <array of c(p,q,R,S)>)

####


