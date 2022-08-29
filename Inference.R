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
total.move = 5

#Length of the experiment in seconds
tau = 0.04

#The number of time steps in a single simulation.
#Increase if your simulations fail to converge, but this comes with a compute time trade off.
steps = 100
dt = tau/steps

#Move step size: rate constant moves are multiplications by the exponential of a variable selected
#from a t distrubution with sd ln.sd and df degrees of freedom. 
#This has the advantage of not altering the "rates" that indicate a "movement" 
#from a Swi6 state to itself: for the sake of these simulations, these are coded to stay at 0
ln.sd = 0.02
df = 10

#This defines the bounds of a log uniform prior
low.mag = -2
up.mag = 10



#We track the number of moves accepted by the inference: we want accepted/total moves ~ 0.1-0.2
#Higher, and the moves aren't as big as they could be to move the inference quickly
#Lower, and the algorithm will get stuck too often
accepted = 0

#We track the number of moves each type of move proposed by the inference
#we want accepted/total moves ~ 0.1-0.2 for each move type, for the same reasons as above
proposed.moves = rep(0,total.move)
accepted.moves = rep(0,total.move)

#Number of times the program saves its current chain. Saving every step is
#too computationally costly, but if the program quits in the middle, we want to
#have our progress saved. Tune this for your optimal trade off
save.runs = 40


# Summary: This is the actual method that needs to be executed for the simulation
#
# Parameters:
#   
#  raw.experiment: An array with dimensions SxSxN, S = number of chemical states, N = number of experiments
#  equilibria: The equilibria of the experiments. This mainly matters to know how many molecules to bin in each state
#  total.iter: The number of steps for the simulation to run
#  experiments.per.iter: The number of simulations to run for the 
#                        Approximate Bayesian Computation of the likelihood
#  cores: The number of CPUs available for the program to run
#  name: A name for the monte carlo chain. Unique names prevent output files from clobbering each other
#  init.rates: Some initial state for the MCMC to begin inference from
#  kill.rates: A binary hollow symmetric matrix with dimensions SxS. 
#              0 entries mean the corresponding transition is chemically forbidden
#              1 entries mean the corresponding transition is allowed
#              If rates are decreasing too much and not converging, it can mean
#              that the transition does not actually occur
#
# Returns: Nothing. Instead, it saves several files, including an RDS of the chain.
               

main = function(raw.experiment, equilibria=calcEQs(raw.experiment), total.iter = 20000, experiments.per.iter = 3000,cores=1, name = "",init.rates = NULL, kill.rates = (1-diag(length(equilibria)))) {

    print("In Method")

    #Picks which steps to save the inference chain at
    step.save = max(floor(total.iter/save.runs),1)

    #We forbid transition possibilities from being nonsymmetric
    if(!isSymmetric(kill.rates) || any(diag(kill.rates) != 0) || any(kill.rates !=0 & kill.rates != 1)){
        print("kill")
        print(kill.rates)
        stop("Non symmetric transition allowances are forbidden, and we require self transition rates to be set at zero")
    }

    
    global.assign(raw.experiment, equilibria, kill.rates)
        

    cl <- makeForkCluster(cores) 
    registerDoParallel(cl)

    #This is a flat out semi decent guess on rate constants
    if(is.null(init.rates)){
        init.rates = list(rates = init.cond.rates, likelihood = log.likelihood.main(init.cond.rates,experiments.per.iter),accepted = TRUE)
    } else{
	    init.rates = do.move(list(rates = init.rates,likelihood = log.likelihood.main(init.rates,experiments.per.iter),accepted = TRUE),experiments.per.iter)
    }

    #parameter.inference is the list containing all of the states of the monte carlo chain
    #It is a list of lists, where each sublist contains a matrix encoding the rate constants,
    #the log likelihood of these rate constants, and whether the preceeding move was accepted or not
    parameter.inference = list(init.rates)



    #With the stochastic way our likelihood is calculated, the inference chain 
    #can get stuck because the likelihood is calculated as higher than it truly is. 
    #This can't be "fixed" because there is a selection effect in play:
    #The number of likelihood calculations is high enough that eventually, 
    #our stochastic likelihood calculation randomly picks a likelihood much
    #higher than the actual one
    #This variable tracks the number of rejections in a row. When it's too high, 
    #the log likelihood is recalculated.
    reject.count = 0

    for(step in (1:total.iter)) {
	
	    print("Iter")
	    print(step)
	    
        #Each monte carlo step only depends on the state of the last
	    parameter.inference[[step+1]] = do.move(parameter.inference[[step]],experiments.per.iter)

	    #This calculates the number of moves proposed of each type
	    proposed.moves[parameter.inference[[step+1]]$move] = proposed.moves[parameter.inference[[step+1]]$move]+1
 
	    #This calculates the number of accepted moves (other than the initial state):
        #Tune move step sizes so that 20% of moves of each type are accepted
	    if(parameter.inference[[step+1]]$accepted){
		    accepted = accepted+1
            reject.count = 0
		    accepted.moves[parameter.inference[[step+1]]$move] = accepted.moves[parameter.inference[[step+1]]$move]+1
	    } else{

            reject.count = reject.count+1
            
            #As mentioned before, we recalculate likelihoods if there are too 
            #many rejections in a row. Change this if statement here to decide
            #what counts as "too many" for you
            if(reject.count > 20){
                parameter.inference[[step+1]]$likelihood = log.likelihood.main(parameter.inference[[step+1]]$rates, experiments.per.iter)
                reject.count = 0
            }

        }

	    if(((step %% step.save) == 0) || (step == total.iter)){
		    print(parameter.inference[[step+1]])
            print(paste(name,"_savestate.rds",sep=""))
		    
            #We save our saved state twice on purpose. 
            #The first says which step the inference stopped at
            #while the second gives the chain a convenient place for the chain
            #to pick up from
            saveRDS(parameter.inference,paste(name,"_",step,"parameters.rds",sep = ""))
            saveRDS(parameter.inference,paste(name,"_savestate.rds",sep=""))
        }

	    print("Acceptance by move type")
	    print(accepted.moves/proposed.moves)
    }

    saveRDS(parameter.inference,paste(name,"_final_parameters.rds",sep = ""))

    likelihoods = numeric(length(parameter.inference))

    #We also save a vector of likelihoods to track their evolution.
    #Convergent likelihoods look like a fuzzy caterpillar
    #To a first approximation
    for(i in 1:length(likelihoods)){
        tests = paramter.inference[[i]]
        likelihoods[i] = tests$likelihood
    }
    png(file=paste(name,"Likelihoods.png"))
    plot(likelihoods)
    dev.off()


    #This output is just a simple prediction of the experimental data as a sanity check.
    trans.predict = experiment.simulation(parameter.inference[[total.iter+1]][["rates"]])*matrix.of.oneovermols



    print("Parameters")
    print(parameter.inference[[total.iter+1]][["rates"]])
    print("Predicted transitions")
    print(trans.predict)
    print("Experimental transitions")
    print(experiment)
    print("Proportion of accepted moves")
    print(accepted/(total.iter))


    stopCluster(cl)
}

# Summary: This method saves several global constants based on the experimental data
#
# Parameters:
#   
#  raw.experiment: An array with dimensions SxSxN, S = number of chemical states, N = number of experiments
#  equilibria: The equilibria of the experiments. This mainly matters to know how many molecules to bin in each state
#  kill.rates: A binary hollow symmetric matrix with dimensions SxS. 
#              0 entries mean the corresponding transition is chemically forbidden
#              1 entries mean the corresponding transition is allowed
#              If rates are decreasing too much and not converging, it can mean
#              that the transition does not actually occur
#
# Returns: Nothing.

global.assign = function(raw.experiment, equilibria, kill.rates){
    
    #How many types of each type of molecule there are at equilibirium
    #There can be a rounding defect where we lose one or two molecules because of the rounding
    mols = round(totalmols*equilibria)
    mode(mols) = "integer"

    #This helps generate the proportion of molecules in each state from the numbers the simulations generate
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

    
    init.cond.rates  = -matrix(rep(log(diag(apply(raw.experiment,c(1,2),mean))),nrow(experiment)),ncol = ncol(experiment),byrow = TRUE)/tau

    #The program sets all transitions that have no reaction to 0  
    #since zeroes are easy to ignore in computations
    init.cond.rates = init.cond.rates * kill.rates
    assign('experiment',experiment,envir=.GlobalEnv)
    assign('initial.conditions',initial.conditions,envir=.GlobalEnv)
    assign('matrix.of.oneovermols',matrix.of.oneovermols,envir=.GlobalEnv)
    assign('init.cond.rates',init.cond.rates,envir=.GlobalEnv)
    assign('kill.rates',kill.rates,envir=.GlobalEnv)
}


# Summary: This method calculates a single simulation of the experimental data
#          based one set of rate constants
#
# Parameters:
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#
# Returns: An SxS matrix of integers, where the entry of the jth column and
#          the ith row is the number of molecules which were in state j at the 
#          beginning of the simulated experiment and state i at the end

experiment.simulation = function(current.rates){
	
	time = 0
	amount = initial.conditions

	transition.rates = current.rates*dt

    damount = matrix(rep(0,nrow(experiment)*ncol(experiment)),nrow = nrow(experiment))

	stat = 1:ncol(experiment)

	while(time < tau){

	damount = damount-damount 

        for(k in 1:ncol(experiment)){
            notk = stat[stat!=k]
            outgoing.rates = sum(transition.rates[notk,k])
            probs = transition.rates[,k]/outgoing.rates
            for(i in 1:ncol(experiment)){
				#This is calculated once so that the amount leaving is the same amount as the amount
				#going to the other pools of Swi6
				remaining = rbinom(1,amount[k,i],exp(-outgoing.rates))
				leaving = amount[k,i]-remaining

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
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#       experiments.per.iter: an integer indicating how many simulations should
#                             be run by the calculator
#
# Returns: An SxSxexperiments.per.iter matrix of integers, where the entry
#          in the kth layer of the jth column and the ith row is the number 
#          of molecules which were in state j at the beginning of the kth 
#          simulated experiment and state i at the end

simulations.calculator = function(current.rates,experiments.per.iter){


    print("init con")
    print(experiments.per.iter)

	transition.big.array <- foreach(j=1:experiments.per.iter, .combine=special.abind, .export=c("experiment.simulation","initial.conditions","dt","tau","experiment"),.packages= c('MASS','MGLM')) %dopar% {
   		
        transition.matrix = experiment.simulation(current.rates)  #calling a function


		#This block implements a pseudocount. It was not useful for our data, but could be for yours

		#for(i in 1:ncol(experiment)){
		#	zers = which((transition.matrix[,i]<1e-12))
        #    amts = which(!(transition.matrix[,i]<1e-12))
        #    znum = length(zers)
        #    pseudo = 1
        #    picker = as.numeric(amts)
        #    picker = picker/sum(picker)
        #    picked = rmultinom(1,pseudo*znum,picker)
		#	transition.matrix[,i][amts] = transition.matrix[,i][amts] 
		#	transition.matrix[,i][zers] = pseudo
		#}

		transition.matrix
    }


	
	return(transition.big.array)
}


# Summary: This method calculates the ln density of the rates prior distribution
#          assuming a set of independent Jeffries' priors for each positive rate
#
# Parameters:
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#
# Returns: the ln density of the joint rates prior

rate.prior = function(current.rates){
 	

    #We had originally intended a ln-uniform prior for the rates.
    #When we switched to a Jeffries prior, we forgot to remove the
    #constant log((up.mag-low.mag)) term. This had no ultimate effect 
    #on the simulation, as it was a constant that was normalized away from the posterior.
    #However, we are leaving it in for transparency, and will remove it with our next push.
	prior = -(log((current.rates))+log((up.mag-low.mag)))

    #This assumes that there's no way for active rates to reach zero with the MCMC moves.
    #This is true for our move set, but might not be for yours if you change them.
	prior[current.rates == 0] = 0

	return(sum(prior))
}

# Summary: This method calculates the ln likelihood of one set of rates with Approximate 
#          Bayesian Computation 
#
# Parameters:
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#       experiments.per.iter: an integer indicating how many simulations should
#                             be run by the calculator
#
# Returns: a float of the approximate ln likelihood of the rates

log.likelihood.experiment = function(current.rates,experiments.per.iter){

	data = simulations.calculator(current.rates,experiments.per.iter)

    alphas.matrix = matrix(nrow=dim(data)[1], ncol = dim(data)[2])

    for(j in 1:dim(data)[2]){

        #This is transposed for MGLMfit compatibility
        tofit = t(data[,j,]) 

        alphas.matrix[,j] = MGLMfit(tofit, dist="MN")@estimate

    }

	logliks = numeric(ncol(experiment))
    
    for(j in 1:ncol(experiment)){
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
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#       experiments.per.iter: an integer indicating how many simulations should
#                             be run by the calculator
#
# Returns: a float of the approximate ln posterior density of the rates

log.likelihood.main = function(current.rates, experiments.per.iter){

	prior.like = rate.prior(current.rates)
	exp.like = log.likelihood.experiment(current.rates,experiments.per.iter)
	print("prior likelihood")
	print(prior.like)
	print("experimental likelihood")
	print(exp.like)
	return(exp.like+prior.like)
		
}



# Summary: this method calculates whether we accept or reject a given move of the MCMC 
#
# Parameters:
#       lik.trial: the posterior density of the proposed new state in the inference
#       lik.cur: the posterior density of the prior state in the inference
#
# Returns: a boolean indicating whether the proposed move is accepted or not

accept.test = function(lik.trial, lik.cur) {
	lhr = lik.trial - lik.cur
	
	acceptance.r = log(runif(1))	
	print("lhr")
	print(lhr)
	print("lik.trial")
	print(lik.trial)
	print("lik.cur")
	print(lik.cur)
	print('accept')
	print(acceptance.r)
	if (is.na(lik.trial)) {
		return(FALSE)
	}
	
	if (lik.trial == -Inf & lik.cur == -Inf){
		print("Both states equally terribly bad (probability zero). Going with the new one for the sake of movement. Throw this step away for inference")
		return(TRUE)
	}	
	else{
	print(paste("Difference:",lhr,"Accepted:",(lhr > acceptance.r)))

	if (lhr > acceptance.r) {
		return(TRUE)
	}
	else{
		return(FALSE)
	}
	}
}


# Summary: this method proposes a move generated by multiplying all nonzero rates 
#          by ln-t distributed variables
#
# Parameters:
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#
# Returns: a hollow SxS nonnegative matrix, where the entry of the 
#          of the jth column in the ith row represents the changed rates 
#          of transition from state j to state i in units 1/sec.
move.all.rates = function(current.rates){

	new.rates = current.rates


	sdid = sample.int(length(ln.sd),1)	

	new.rates = current.rates*matrix(exp(rt(length(experiment[,,1]), df = df)* ln.sd[sdid]/(nrow(experiment)^2)),nrow = nrow(experiment))


	return(new.rates)

}

# Summary: this method proposes a move generated by multiplying one nonzero rate
#          by a ln-t distributed variable
#
# Parameters:
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#
# Returns: a hollow SxS nonnegative matrix, where the entry of the 
#          of the jth column in the ith row represents the changed rates 
#          of transition from state j to state i in units 1/sec.
move.one.rate = function(current.rates){
	
	new.rates = current.rates
	pick = 1:nrow(experiment)
	i = sample(pick,1,replace = TRUE)
    j = sample(pick[kill.rates[i,] != 0], 1,replace = TRUE)
	sdid = sample.int(length(ln.sd),1)	
	new.rates[i,j] = current.rates[i,j]*exp(rt(1, df = df)* ln.sd[sdid])
	print(paste("SD",ln.sd[sdid]))
	return(new.rates)

}

# Summary: this method proposes a move generated by multiplying two opposed nonzero rates 
#          by ln-t distributed variables
#
# Parameters:
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#
# Returns: a hollow SxS nonnegative matrix, where the entry of the 
#          of the jth column in the ith row represents the changed rates 
#          of transition from state j to state i in units 1/sec.
scale.opposite.rate = function(current.rates){

	new.rates = current.rates
    pick = 1:nrow(experiment)
	i = sample(pick,1,replace = TRUE)
    j = sample(pick[kill.rates[i,] != 0], 1,replace = TRUE)
    sdid = sample.int(length(ln.sd),1)	
	scale.exp = rt(1, df = df)* ln.sd[sdid]/2
	print(paste("SD",ln.sd[sdid]))
        new.rates[i,j] = current.rates[i,j]*exp(scale.exp)
        new.rates[j,i] = current.rates[j,i]*exp(scale.exp)

        return(new.rates)

}

# Summary: this method proposes a move generated by proposing a de novo replacement for a nonzero rate
#
# Parameters:
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#
# Returns: a hollow SxS nonnegative matrix, where the entry of the 
#          of the jth column in the ith row represents the changed rates 
#          of transition from state j to state i in units 1/sec.
change.one.rate = function(current.rates){

	new.rates = current.rates
    pick = 1:nrow(experiment)
	i = sample(pick,1,replace = TRUE)
    j = sample(pick[kill.rates[i,] != 0],1, replace = TRUE)
    changed.rate = c(i,j)
	new.rates[changed.rate[1],changed.rate[2]] = rgamma(1, shape = sqrt(init.cond.rates[changed.rate[1],changed.rate[2]]),scale = sqrt(init.cond.rates[changed.rate[1],changed.rate[2]])^3)	
	print("Small Aether")
	return(new.rates)

}

# Summary: this method proposes a move generated by proposing a de novo replacement for all nonzero rates
#
# Parameters:
#       current.rates: a hollow SxS nonnegative matrix, where the entry of the 
#                      of the jth column in the ith row represents the rate 
#                      of transition from state j to state i in units 1/sec.
#
# Returns: a hollow SxS nonnegative matrix, where the entry of the 
#          of the jth column in the ith row represents the changed rates 
#          of transition from state j to state i in units 1/sec.
change.all.rate = function(current.rates){

	new.rates = matrix(rgamma(nrow(experiment)*ncol(experiment), shape = sqrt(init.cond.rates),scale = sqrt(init.cond.rates)^3), nrow = nrow(experiment),byrow = TRUE)
	
    new.rates = new.rates*kill.rates

	print("Big Aether")
	return(new.rates)

}


# Summary: this method implements a single MCMC step
#
# Parameters:
#       full.current.rates: A list with the following components
#                      rates: a hollow SxS nonnegative matrix, where the entry of the 
#                       of the jth column in the ith row represents the rate 
#                       of transition from state j to state i in units 1/sec.
#                      likelihood: a double indicating the approximate ln density of the rates in the step
#                      accepted: a boolean which is TRUE if the rates were accepted from a new move
#                      move: an integer denoting which move was attempted at the last step
#
# Returns: A list with the same components as above, after a move was attempted on it

do.move = function(full.current.rates,experiments.per.iter){

	move = sample(1:total.move,size = 1,replace = TRUE)
	print("Move")
	print(move)

	new.rates = switch(move,move.all.rates(full.current.rates$rates),move.one.rate(full.current.rates$rates),scale.opposite.rate(full.current.rates$rates),change.one.rate(full.current.rates$rates),change.all.rate(full.current.rates$rates))

	print("old likelihood")
	print(full.current.rates$likelihood)

	

	new.like = log.likelihood.main(new.rates,experiments.per.iter)

    new.like.ast = new.like

    if(move %in% c(4,5)){
        old.pick = (dgamma(full.current.rates$rates, shape = sqrt(init.cond.rates),scale = sqrt(init.cond.rates)^3, log = T))
        new.pick = (dgamma(new.rates, shape = sqrt(init.cond.rates),scale = sqrt(init.cond.rates)^3, log = T))
        new.like.ast = new.like.ast+sum(old.pick[old.pick != -Inf])-sum(new.pick[new.pick != -Inf])
    }


	print("new like")
	print(new.like)

	if(accept.test(new.like.ast,full.current.rates$likelihood)){
		return(list(rates = new.rates, likelihood = new.like, accepted = TRUE,move = move))
	}
	else{
		return(list(rates = full.current.rates$rates,likelihood = full.current.rates$likelihood,accepted = FALSE,move = move))
	}

}

#this is simply a wrapper for the abind function to make it easier to implement for the parallelizer
special.abind = function(arr1,arr2){
	
	return(abind(arr1,arr2,along = 3))
	
}

