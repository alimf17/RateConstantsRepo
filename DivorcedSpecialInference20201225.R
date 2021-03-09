library(abind)
library(foreach)
library(doParallel)
library(MASS)
library(fitdistrplus)
library(MGLM)

source("Eq.R")

#An integer indicating the number of single tracks in an experiment.
totalmols = 1500


#The total number of moves
total.move = 4


#Length of the experiment in seconds
tau = 0.04

#The number of time steps in a single simulation.
#Increase if your simulations fail to converge, but this comes with a compute time trade off.
steps = 100
dt = tau/steps


#Posterior Mean Four State Rates
super.rates = readRDS("DenseRates.rds")

#Posterior Mean clr4-delta Rates
sub.rates = readRDS("Clr4DeltRates.rds")

#Fixed rates denotes the rates that will be unaltered by the inference
#It's a double equal to the fixed rate if the rate is unchanged
#And NA for rates that will be changed as we infer our parameters
fixed.rates = matrix(nrow = 5, ncol = 5)
diag(fixed.rates) = 0
fixed.rates[c(1,2,5),1] = super.rates[c(1,2,4),1]
fixed.rates[c(1,2,5),2] = super.rates[c(1,2,4),2]
fixed.rates[c(1,2,5),5] = super.rates[c(1,2,4),4]
fixed.rates[3,1:3] = sub.rates[3,]
fixed.rates[1:3,3] = sub.rates[,3]
fixed.rates[4,1:2] = super.rates[3,1:2]-sub.rates[3,1:2]


eqs = theoreticalratesEq(super.rates)$eq


true.states = 5

trenchcoats = 3

extras = 1


# Summary: This method saves several global constants based on the experimental data
#
# Parameters:
#
#  raw.experiment: An array with dimensions ExExN, E = number of experimentally apparent states, N = number of experiments
#  equilibria: The equilibria of the experiments. This mainly matters to know how many molecules to bin in each state
#
# Returns: Nothing.

global.assign = function(raw.experiment, equilibria){

    #"trenchcoats" is an array of experimental states which are secretlu
    #multiple chemical states. We fix it to c(3) because this specific
    #parameterization only works if it's c(3)
    assign('trenchcoats',3,envir=.GlobalEnv)

    #"extras" is an array indicating how many additional states each trenchcoat
    #is hiding. We fix it to c(1) because this parameterization fails otherwise
    assign('extras',1,envir=.GlobalEnv)

    #How many types ,experiments.per.iterof small molecules there are at equilibirium
    #There can be a rounding defect where we lose one or two molecules because of the rounding
    mols = round(totalmols*equilibria)
    mode(mols) = "integer"


    #This helps generate the proportion of changed molecules
    matrix.of.oneovermols = matrix(rep(1/mols,nrow(raw.experiment)),ncol = ncol(raw.experiment),byrow = TRUE)


    experiment = array(dim = dim(raw.experiment))

    for(j in 1:ncol(raw.experiment)){
        experiment[,j,] = round(raw.experiment[,j,]*mols[j])
    }


    exprat = array(dim = dim(experiment))

    for(k in 1:dim(experiment)[3]){
        for(j in 1:ncol(experiment)){
            exprat[,j,k] = experiment[,j,k]/sum(experiment[,j,k])
        }
    }


    #We are not measuring concentrations with these: merely proportion of
    #the molecule of each original form. We can treat each initial
    #form as a pool distinct from the others
    initial.conditions = diag(nrow=nrow(experiment))

    #We want the initial conditions to be zero for everything except what's on the diagonal
    diag(initial.conditions) <- mols

    true.states = length(equilibria)+sum(extras)

    assign('experiment',experiment,envir=.GlobalEnv)
    assign('initial.conditions',initial.conditions,envir=.GlobalEnv)
    assign('matrix.of.oneovermols',matrix.of.oneovermols,envir=.GlobalEnv)
    assign('true.states',true.states,envir=.GlobalEnv)

}

# Summary: this method converts our parmeters to a set of rate constants
#
# Parameters:
#       pars: a 2x2 matrix, where p = pars[1,1], q = pars[2,1], R = pars[1,2], S = [2,2],
#			  as these parameters are defined in the supplementary methods
# Returns: a hollow true.statesXtrue.states nonnegative matrix, where the entry of the
#          of the jth column in the ith row represents the rate
#          of transition from chemical state j to chemical state i in units 1/sec.

convert.pars.to.rates= function(pars){


    rates = fixed.rates

    p = pars[1,1]

    rates[5,3] = pars[1,2]

    rates[c(1,2,5), 4] = (super.rates[c(1,2,4), 3]-rates[c(1,2,5), 3]*p)/(1-p)

    rates[3,5] = pars[2,1]*super.rates[3,4]
    rates[4,5] = (1-pars[2,1])*super.rates[3,4]


    expanded.eqs = c(eqs[1:2], p*eqs[3], (1-p)*eqs[3], eqs[4])

    #This should already be the case, but we HAVE to make sure
    diag(rates) = 0

    #This MUST be assigned BEFORE rates[3,4]
    rates[4,3] = pars[2,2]

    #This assignment is made for convenience purposes to clean up the actual assignment
    rates[3,4] = 0

    rates[3,4] = drop(rates[4,] %*% expanded.eqs)/expanded.eqs[4] - sum(rates[,4])

    return(rates)

}

# Summary: This method calculates a single simulation of the experiment
#          based one set of rate constants
#
# Parameters:
#       current.rates: a hollow true.statesXtrue.states nonnegative matrix, where the entry of the
#                      of the jth column in the ith row represents the rate
#                      of transition from chemical state j to chemical state i in units 1/sec.
#		real.init.conds: an S length vector indicating the current number of 
#						 molecules in each state at equilibrium. This has to 
#						 change for each calculation because pars[1,1] indicates
#						 how many of each beta state exists at equilibrium
#
# Returns: An true.statesXtrue.states matrix of integers, where the entry of the jth column and
#          the ith row is the number of molecules which were in state j at the
#          beginning of the simulated experiment and state i at the end

experiment.simulation = function(current.rates, real.initial.conds){
	
	time = 0
	amount = diag(real.initial.conds)

    stat = 1:true.states

	
    transition.rates = current.rates*dt

	damount = matrix(rep(0,true.states*true.states),nrow = true.states)
	while(time < tau){

    damount = damount-damount

    for(k in 1:true.states){
        notk = stat[stat!=k]
        outgoing.rates = sum(transition.rates[notk,k])
        probs = transition.rates[,k]/outgoing.rates
        for(i in 1:true.states){
				#This is calculated once so that the amount leaving is the same amount as the amount
				#going to the other pools of Swi6
				remaining = rbinom(1,amount[k,i],exp(-outgoing.rates))
				leaving = amount[k,i]-remaining

				#We need to account for issues with rounding. If we have more than two compartment that Swi6 can
				#transition into, then simply rounding has the potential to leave us one molecule short, if
				#everything gets rounded down. We correct that here, randomly selecting one of the pools to
				#get an additional molecule, if necessary

				#for our multinomial distribution, we are SUPER relying on the fact that self transitions
				#are counted as having probability zero
				probs = transition.rates[,k]/outgoing.rates
				move.vector = rmultinom(1,leaving,probs)
				move.vector[k] = -leaving
				
				if(sum(move.vector)!=0){
					stop("You're losing or gaining particles, wiseguy.")
				}
				
				damount[,i] = damount[,i]+move.vector

			}
		}
		amount = amount+damount
		time = time + dt
	}
	
	
	return(amount)
	
}

# Summary: This method calculates many simulations of the experimental data
#          based one set of rate constants
#
# Parameters:
#       current.rates: a hollow true.statesXtrue.states nonnegative matrix, where the entry of the
#                      of the jth column in the ith row represents the rate
#                      of transition from state j to state i in units 1/sec.
#       experiments.per.iter: an integer indicating how many simulations should
#                             be run by the calculator
#
# Returns: An true.statesXtrue.statesxexperiments.per.iter matrix of integers, where the entry
#          in the kth layer of the jth column and the ith row is the number
#          of molecules which were in chemical state j at the beginning of the kth
#          simulated experiment and chemical state i at the end

simulations.calculator = function(current.rates,experiments.per.iter){


    real.initial.conds = numeric(length(diag(initial.conditions))+sum(extras))


	#This code block assigns all the molecules in the experimental state equilbirium
	#to an equilibrium population which applies to the actual chemical states
    index = 1
    for(i in 1:length(diag(initial.conditions))){

        if(i %in% trenchcoats){

            much = extras[which(trenchcoats==i)]

            percent = current.rates[index+(0:much), index+(0:much)]%*%diag(1/(colSums(current.rates[index+(0:much), index+(0:much)])))

            hiddenEq = calcEQs(percent)

            curtain = floor(diag(initial.conditions)[i]*hiddenEq)
			
			#While curtain should contain a reasonable equilbrium mostly, the floor function means that 
			#we can be off by a molecule. This line helps us assign the molecule we lost
            curtain = curtain+rmultinom(1,diag(initial.conditions)[i]-sum(curtain),hiddenEq
)
            real.initial.conds[index+0:much] = curtain
            index = index + much
        }
        else{
            real.initial.conds[index] = diag(initial.conditions)[i]
        }
        index = index + 1
    }


	transition.big.array <- foreach(j=1:experiments.per.iter, .combine=special.abind, .export=c("experiment.simulation","real.initial.conds","dt","tau","experiment"),.packages= c('MASS','MGLM')) %dopar% {
   		transition.matrix = experiment.simulation(current.rates, real.initial.conds)  #calling a function


		#implements a pseudocount, which matches experimental data better, as no value ends up being less than 0.002 in the raw data
		commentOut = 'for(i in 1:true.states){
			zers = which((transition.matrix[,i]<1e-12))
            amts = which(!(transition.matrix[,i]<1e-12))
            znum = length(zers)
            pseudo = 3
			if(znum > 0){
				print("pseudocount correction")
			}
            picker = as.numeric(amts)
            picker = picker/sum(picker)
            picked = rmultinom(1,pseudo*znum,picker)
			transition.matrix[,i][amts] = transition.matrix[,i][amts] 
			transition.matrix[,i][zers] = pseudo
		}'

		transition.matrix
}


	
	return(transition.big.array)
}



# Summary: This method calculates the apparent experimental data given an experimental simulation
#
# Parameters:
#		transition.big.array: An true.statesXtrue.statesxexperiments.per.iter matrix of integers, where the entry
#	          in the kth layer of the jth column and the ith row is the number
#	          of molecules which were in chemical state j at the beginning of the kth
#	          simulated experiment and chemical state i at the end
#
# Returns: An (true.states-1)X(true.states-1)xexperiments.per.iter matrix of integers, where the entry
#          in the kth layer of the jth column and the ith row is the number
#          of molecules which were in experimental state j at the beginning of the kth
#          simulated experiment and experimental state i at the end

translate.sim.to.data = function(transition.big.array){

    data = array(dim=c(nrow(experiment),ncol(experiment),dim(transition.big.array)[3]))

    index = 1
    for(i in 1:nrow(experiment)){

        index2 = 1

        for(j in 1:nrow(experiment)){

            data[i,j,] = transition.big.array[index, index2,]

            index2 = index2+1

            if(j %in% trenchcoats){

                for(k in 1:extras[which(trenchcoats == j)]){
                    data[i,j,] = data[i,j,]+transition.big.array[index, index2,]
                    index2 = index2+1
                }
            }

        }

        index = index+1
        
        if(i %in% trenchcoats) {

            for(m in 1:extras[which(trenchcoats == i)]){

                index2 = 1

                for(j in 1:nrow(experiment)){

                    data[i,j,] = data[i,j,]+transition.big.array[index, index2,]

                    index2 = index2+1

                    if(j %in% trenchcoats){

                        for(k in 1:extras[which(trenchcoats == j)]){
                            data[i,j,] = data[i,j,]+transition.big.array[index, index2,]
                            index2 = index2+1
                        }
                    }

                }
          
                index = index + 1

            }
       
                 
        
        }

    

    }

    return(data)

}


# Summary: This method calculates K+ln density of the parameter prior distribution, 
#		   where K is an unknown constant
#
# Parameters:
#       pars: a 2x2 matrix, where p = pars[1,1], q = pars[2,1], R = pars[1,2], S = [2,2],
#			  as these parameters are defined in the supplementary methods
#
# Returns: the ln density of the joint rates prior

prior.pars = function(pars){

    #Per our new parameterization, we only care if our parameters give possible rates
    good = all(convert.pars.to.rates(pars) >= 0)



    if(good){
        #The log density of the prior is extremely unlikely to be zero
        #But we don't care: all possible rates have the same density, which gets 
        #normalized out, so 0 is as good a stand-in as any
        return(0)
    } else {
        return(-Inf)
    }
}


# Summary: This method calculates the ln likelihood of one set of rates with Approximate
#          Bayesian Computation
#
# Parameters:
#       current.rates: a hollow true.statesXtrue.states nonnegative matrix, where the entry of the
#                      of the jth column in the ith row represents the rate
#                      of transition from state j to state i in units 1/sec.
#       experiments.per.iter: an integer indicating how many simulations should
#                             be run by the calculator
#
# Returns: a float of the approximate ln likelihood of the rates

log.likelihood.experiment = function(current.rates,experiments.per.iter){

	sim = simulations.calculator(current.rates,experiments.per.iter)

    data = translate.sim.to.data(sim)


    means = apply(experiment,c(1,2),mean)

    per.mol = matrix(numeric(length(means)),nrow = nrow(means))

    per.mol[initial.conditions != 0] = 1/initial.conditions[initial.conditions != 0]

    mean.transitions = means%*%per.mol

    alphas.matrix = matrix(nrow=dim(data)[1], ncol = dim(data)[2])

    for(j in 1:dim(data)[2]){

        #This is transposed to fit MGLMfit, where COLUMNs are categories
        tofit = t(data[,j,])

        alphas.matrix[,j] = MGLMfit(tofit, dist="MN")@estimate
    }

	logliks = numeric(nrow(experiment))
    
    for(j in 1:nrow(experiment)){
        for(k in 1:dim(experiment)[3]){
            logliks[j] = logliks[j]+dmultinom(experiment[,j,k], prob = alphas.matrix[,j], log = TRUE)

        }
    }
	
    return(sum(logliks))

}


# Summary: this method calculates the ln density of one set of rates with approximate
#          bayesian computation
#
# Parameters:
#       pars: a 2x2 matrix, where p = pars[1,1], q = pars[2,1], R = pars[1,2], S = [2,2],
#			  as these parameters are defined in the supplementary methods
#       experiments.per.iter: an integer indicating how many simulations should
#                             be run by the calculator
#
# Returns: a float of the approximate ln posterior density of the rates

log.likelihood.main = function(current.pars, experiments.per.iter){

	prior.like = prior.pars(current.pars)
	if(prior.like > -Inf){
        exp.like = log.likelihood.experiment(convert.pars.to.rates(current.pars),experiments.per.iter)
        return(exp.like+prior.like)
    } else {
        return(prior.like)
    }
}

#this is simply a wrapper for the abind function to make it easier to implement for the parallelizer
special.abind = function(arr1,arr2){
	
	return(abind(arr1,arr2,along = 3))
	
}



