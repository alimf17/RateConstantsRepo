library(abind)
#library(dyplr)
library(foreach)
library(doParallel)
library(MASS)
library(fitdistrplus)
library(MGLM)

source("/home/alimf/Histones/Scripts/Eq.R")

#Column indicates the state where relevant Swi6 started, rows indicate where 
#relevant Swi6 ended at the end of the experiment. Note that I adjucated that <0.01
#is 0.005, and added 0.005 (largest plausible round errors) to these transitions 
#until each column summed to 1, starting first with the highest probability 
#transitions
#experiment = matrix(c(0.755,0.165,0.075,0.005,0.215,0.585,0.195,0.005,0.04,0.08,0.79,0.09,0.01,0.005,0.13,0.855), nrow = 4, ncol = 4)



#An integer indicating the number of molecules we're tracking in a trial. Chose 300 because it's in range and gives integer
#results with the equilibria

#Changed this to 1500 because, according to Saikat, experiments amalgamate single molecule tracks in individual cells, and total about 1000-2000 molecules
totalmols = 1500


#The total number of moves
total.move = 8

#Length of the experiment in seconds
tau = 0.04
steps = 100
dt = tau/steps

#Move step size: move is a multiplication by the exponential of a variable selected
#from a normal distrubution. This has the advantage of not altering the "rates" that
#indicate a "movement" from a Swi6 state to itself: for the sake of these simulations,
#these are coded to stay at 0
ln.sd = 0.2

df = 10

#This defines the bounds of a log uniform prior
low.mag = -2
up.mag = 10



#Begin the count of accepted moves
accepted = 0

#Begin the counts of each proposed move type
proposed.moves = rep(0,total.move)

#Begin the counts of each accepted move type
accepted.moves = rep(0,total.move)

#number of saves
save.runs = 40

#The executable main method. total.iter is the monte carlo chain length, experiments.per.iter
#is the number of simulations run per experiment
main = function(raw.experiment, equilibria=calcEQs(raw.experiment), total.iter = 100, experiments.per.iter = 3000,cores=1, name = "",init.eqmag = NULL, kill.rates = (1-diag(length(equilibria)))) {

    if(!is.null(init.eqmag)){
        kill.rates[!is.finite(init.eqmag)] = 0
    }

    print("In Method")

    step.save = max(floor(total.iter/save.runs),1)

    #We forbid transition possibilities from being nonsymmetric
    if(!isSymmetric(kill.rates) || any(diag(kill.rates) != 0)){
        stop("Non symmetric transition allowances are forbidden, and we require self transition rates to be set at zero")
    }

    
    global.assign(raw.experiment, equilibria, kill.rates)
        

    #setup parallel backend to use many processors
    cl <- makeForkCluster(cores) #not to overload your computer
    registerDoParallel(cores)

    #This is a flat out semi decent guess on rate constants
    if(is.null(init.eqmag)){
        init.eqmag = list(eqmag = rates.to.eqmag(init.cond.rates), likelihood = log.likelihood.main(rates.to.eqmag(init.cond.rates),experiments.per.iter),accepted = TRUE)
    } else{
	    init.eqmag = do.move(list(eqmag = init.eqmag,likelihood = log.likelihood.main(init.eqmag,experiments.per.iter),accepted = TRUE),experiments.per.iter)
    }

    #parameter.inference is the list containing all of the states of the monte carlo chain
    #It is a list of lists, where each sublist contains a matrix encoding the rate constants,
    #the log likelihood of these rate constants, and whether the preceeding move was accepted or not
    parameter.inference = list(init.eqmag)


    reject.count = 0

    #This steps through the monte carlo chain
    for(step in (1:total.iter)) {
	
	    print("Iter")
	    print(step)
	    #Each monte carlo step only depends on the state of the last
	    parameter.inference[[step+1]] = do.move(parameter.inference[[step]],experiments.per.iter)

	    #This calculates the number of moves proposed of each type
	    proposed.moves[parameter.inference[[step+1]]$move] = proposed.moves[parameter.inference[[step+1]]$move]+1
 
	    #This calculates the number of accepted moves (other than the initial state): I can tune my move steps according to this
	    if(parameter.inference[[step+1]]$accepted){
		    accepted = accepted+1
            reject.count = 0
		    accepted.moves[parameter.inference[[step+1]]$move] = accepted.moves[parameter.inference[[step+1]]$move]+1
	    } else{
            reject.count = reject.count+1
            if(reject.count > 20){
                parameter.inference[[step+1]]$likelihood = log.likelihood.main(parameter.inference[[step+1]]$eqmag, experiments.per.iter)
                reject.count = 0
            }

        }
	    if(((step %% step.save) == 0) || (step == total.iter)){
		    print(parameter.inference[[step+1]])
            print(paste(name,"_savestate.rds",sep=""))
		    saveRDS(parameter.inference,paste(name,"_",step,"parameters.rds",sep = ""))
            saveRDS(parameter.inference,paste(name,"_savestate.rds",sep=""))
        }

	    print("Acceptance by move type")
	    print(accepted.moves/proposed.moves)
    }

likelihoods = c()

#I want to know how the likelihood evolves	
for(tests in parameter.inference){
	likelihoods = c(likelihoods,tests$likelihood)
}


saveRDS(parameter.inference,paste(name,"_final_parameters.rds",sep = ""))

#The last output of the chain is a simple prediction of the rate constants.
#I also save an RDS of parameter.inference in case there is something more 
#sophisticated I want to try
trans.predict = experiment.simulation(eqmag.to.rates(parameter.inference[[total.iter+1]][["eqmag"]]))*matrix.of.oneovermols


png(file="Likelihoods.png")
plot(likelihoods)
dev.off()

print("Parameters")
print(eqmag.to.rates(parameter.inference[[total.iter+1]][["eqmag"]]))
print("Predicted transitions")
print(trans.predict)
print("Experimental transitions")
print(experiment)
print("Proportion of accepted moves")
print(accepted/(total.iter))


#stop cluster
stopCluster(cl)

}

global.assign = function(raw.experiment, equilibria, kill.rates){
    #How many types ,experiments.per.iterof small molecules there are at equilibirium
    mols = round(totalmols*equilibria)


    print("mols")
    print(mols)

    mode(mols) = "integer"

    print("rounding.defect")
    print(max(abs((totalmols*equilibria - mols)/(totalmols*equilibria))))

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



    print("raw")
    print(apply(raw.experiment, c(1,2), mean))
    print("exprat")
    print(apply(exprat, c(1,2), mean)) 
    print("Exp")
    print(apply(experiment, c(1,2), mean))
    
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
    assign('kill.eqmag',log(kill.rates),envir=.GlobalEnv)
}

rates.to.eqmag = function(current.rates){

    ro = nrow(current.rates)
    co = ncol(current.rates)

    current.eqmag = matrix(nrow = nrow(current.rates), ncol = ncol(current.rates))

    for(i in 1:(ro-1)){

        current.eqmag[i,i] = -Inf
        for(j in (i+1):co){
            if(current.rates[i,j] != 0 && current.rates[j,i] != 0){
                current.eqmag[j,i] = log(current.rates[j,i]) - log(current.rates[i,j])
            } else if(current.rates[i,j] == 0 & current.rates[j,i] == 0){
                current.eqmag[j,i] = -Inf
            } else{
                stop("Asymmetric zero rates!")
            }


            current.eqmag[i,j] = log(current.rates[i,j])
        }
    }

    current.eqmag[ro,ro] = -Inf
    
    return(current.eqmag)

}

eqmag.to.rates = function(current.eqmag){

    ro = nrow(current.eqmag)
    co = ncol(current.eqmag)

    current.rates = matrix(nrow = nrow(current.eqmag), ncol = ncol(current.eqmag))

    for(i in 1:(ro-1)){

        current.rates[i,i] = 0
        for(j in (i+1):co){
            current.rates[j,i] = exp(current.eqmag[i,j] + current.eqmag[j,i])
            current.rates[i,j] = exp(current.eqmag[i,j])
        }
    }

    current.rates[ro,ro] = 0
    
    return(current.rates)
}

#log likelihood calculations
experiment.simulation = function(current.rates){
	
	time = 0
	amount = initial.conditions

	transition.rates = current.rates*dt

    damount = matrix(rep(0,nrow(experiment)*ncol(experiment)),nrow = nrow(experiment))


    stat = 1:nrow(experiment)

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

simulations.calculator = function(current.rates,experiments.per.iter){


    print("init con")
    print(experiments.per.iter)

	transition.big.array <- foreach(j=1:experiments.per.iter, .combine=special.abind, .export=c("experiment.simulation","initial.conditions","dt","tau","experiment"),.packages= c('MASS','MGLM')) %dopar% {
   		transition.matrix = experiment.simulation(current.rates)  #calling a function


		#implements a pseudocount, which matches experimental data better, as no value ends up being less than 0.002 in the raw data
		#for(i in 1:ncol(experiment)){
		#	zers = which((transition.matrix[,i]<1e-12))
        #    amts = which(!(transition.matrix[,i]<1e-12))
        #    znum = length(zers)
        #    pseudo = 3
		#   if(znum > 0){
	#			print("pseudocount correction")
		#	}
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

eqmag.prior = function(current.eqmag){
 	
    current.rates = eqmag.to.rates(current.eqmag)

	prior = -(log((current.rates))+log((up.mag-low.mag)))

	prior[current.rates == 0] = 0

	return(sum(prior))
}

log.likelihood.experiment = function(current.rates,experiments.per.iter){

	data = simulations.calculator(current.rates,experiments.per.iter)

    means = apply(experiment,c(1,2),mean)

    per.mol = matrix(numeric(length(means)),nrow = nrow(means))

    per.mol[initial.conditions != 0] = 1/initial.conditions[initial.conditions != 0]

    mean.transitions = means%*%per.mol

    #alphas.matrix = sapply(apply(aperm(data,c(3,2,1)),2,MGLMfit,dist="MN"),slot, "estimate")

    alphas.matrix = matrix(nrow=dim(data)[1], ncol = dim(data)[2])

    for(j in 1:dim(data)[2]){

        #This is transposed to fit MGLMfit, where COLUMNs are categories
        tofit = t(data[,j,]) 

        alphas.matrix[,j] = MGLMfit(tofit, dist="MN")@estimate

    }

	logliks = numeric(ncol(experiment))
    
    for(j in 1:ncol(experiment)){
        for(k in 1:dim(experiment)[3]){
            logliks[j] = logliks[j]+dmultinom(experiment[,j,k], prob = alphas.matrix[,j], log = TRUE)
        }
    }


    #if(sample.int(500,1) == 500){
        print("Random Double Check of Rates Matrix and logliks")
        print(alphas.matrix)
        print(logliks)
        print("experimental rates")
        print(mean.transitions)
        
    #}
    print(which(logliks == -Inf,arr.ind = TRUE))
    print(experiment[which(logliks == -Inf,arr.ind = TRUE)])

    return(sum(logliks))

}

log.likelihood.main = function(current.eqmag, experiments.per.iter){

    current.rates = eqmag.to.rates(current.eqmag)
	prior.like = eqmag.prior(current.eqmag)
	exp.like = log.likelihood.experiment(current.rates,experiments.per.iter)
	print("prior likelihood")
	print(prior.like)
	print("experimental likelihood")
	print(exp.like)
	return(exp.like+prior.like)
		
}


pick.mag = function(current.eqmag){

    #The number of entries in the off diagonal upper triangular portion of an nxn square matrix is the n-1 triangle number
    entries = nrow(current.eqmag)*(nrow(current.eqmag)-1)/2

    i = 1
    j = 1

    #This while loop protects us if we try to change rates which shouldn't
    while(kill.eqmag[i,j] == -Inf){
        #We have to pick out which entry to sample this way, or else we decide what rate to change unevenly
        pick = sample.int(entries, size = 1)
        

        #Our indices correspond to the mags in this way:
        #   C1 C2 C3 C4 ...
        #R1 -- 1  2  4
        #R2 -- -- 3  5
        #R3 -- -- -- 6
        #R4 -- -- -- --
        #...

        #Etc. The column number is equal to the index of the smallest triangle number
        #larger than our pick, or the index of the largest triangle number less than
        #or equal to our pick plus one. We choose the latter parameterization to not 
        #break the edge case that we pick the largest possible pick 

        triangle = (0:(nrow(current.eqmag)-1))*(1:(nrow(current.eqmag)))/2

        j =  max(which(triangle < pick))+1

        i = pick-triangle[j-1]
        

        if(j <= i){
            stop("picked an eq or diag when we should pick a mag")
        }

    }
    
    return(c(i,j))

}

pick.eq = function(current.eqmag){

    #The number of entries in the off diagonal lower triangular portion of an nxn square matrix is the n-1 triangle number
    entries = nrow(current.eqmag)*(nrow(current.eqmag)-1)/2

    i = 1
    j = 1

    #This while loop protects us if we try to change rates which shouldn't
    while(kill.eqmag[i,j] == -Inf){
        #We have to pick out which entry to sample this way, or else we decide what rate to change unevenly
        pick = sample.int(entries, size = 1)
        

        #Our indices correspond to the mags in this way:
        #   C1 C2 C3 C4 ...
        #R1 -- -- -- --
        #R2 1  -- -- --
        #R3 2  3  -- --
        #R4 4  5  6  --
        #...

        #Etc. The column number is equal to the index of the smallest triangle number
        #larger than our pick, or the index of the largest triangle number less than
        #or equal to our pick plus one. We choose the latter parameterization to not 
        #break the edge case that we pick the largest possible pick 

        triangle = (0:(nrow(current.eqmag)-1))*(1:(nrow(current.eqmag)))/2

        i =  max(which(triangle < pick))+1

        j = pick-triangle[i-1]

        

        if(i <= j){
            stop("picked a mag or diag when we should pick an eq")
        }

    }
    
    return(c(i,j))


}

#moves


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


move.all.mags = function(current.eqmag){

	new.eqmag = current.eqmag

	sdid = sample.int(length(ln.sd),1)	

    magfin = upper.tri(new.eqmag) & is.finite(kill.rates)

	new.eqmag[magfin] = current.eqmag[magfin]+matrix(rt(length(experiment[,,1]), df = df)* ln.sd[sdid]/(nrow(experiment)^2),nrow = nrow(current.eqmag))[magfin]


	return(new.eqmag)

}

move.one.mag = function(current.eqmag){
	
	new.eqmag = current.eqmag
	
    pick = pick.mag(current.eqmag)

    i = pick[1]
    j = pick[2]

    sdid = sample.int(length(ln.sd),1)	
	new.eqmag[i,j] = current.eqmag[i,j]+(rt(1, df = df)* ln.sd[sdid])
	print(paste("SD",ln.sd[sdid]))
	return(new.eqmag)

}

move.all.eqs = function(current.eqmag){

	new.eqmag = current.eqmag

	sdid = sample.int(length(ln.sd),1)	

    eqfin = lower.tri(new.eqmag) & is.finite(kill.rates)

	new.eqmag[eqfin] = current.eqmag[eqfin]+matrix(rt(length(experiment[,,1]), df = df)* ln.sd[sdid]/(nrow(experiment)^2),nrow = nrow(current.eqmag))[eqfin]


	return(new.eqmag)

}

move.one.eq = function(current.eqmag){
	
	new.eqmag = current.eqmag
	
    pick = pick.eq(current.eqmag)

    i = pick[1]
    j = pick[2]

    sdid = sample.int(length(ln.sd),1)	
	new.eqmag[i,j] = current.eqmag[i,j]+(rt(1, df = df)* ln.sd[sdid])
	print(paste("SD",ln.sd[sdid]))
	return(new.eqmag)

}

conjure.all.mags = function(current.eqmag){

	new.eqmag = current.eqmag
    
    magfin = upper.tri(new.eqmag) & is.finite(kill.rates)

	new.eqmag[magfin] = matrix(rt(length(experiment[,,1]), df = df),nrow = nrow(current.eqmag))[magfin]

	return(new.eqmag)

}

conjure.one.mag = function(current.eqmag){
	
	new.eqmag = current.eqmag
	
    pick = pick.mag(current.eqmag)

    i = pick[1]
    j = pick[2]

	new.eqmag[i,j] = (rt(1, df = df))
	return(new.eqmag)

}

conjure.all.eqs = function(current.eqmag){

	new.eqmag = current.eqmag

    eqfin = lower.tri(new.eqmag) & is.finite(kill.rates)

	new.eqmag[eqfin] = matrix(rt(length(experiment[,,1]), df = df),nrow = nrow(current.eqmag))[eqfin]

	return(new.eqmag)

}

conjure.one.eq = function(current.eqmag){
	
	new.eqmag = current.eqmag
	
    pick = pick.eq(current.eqmag)

    i = pick[1]
    j = pick[2]

	new.eqmag[i,j] = (rt(1, df = df))
	return(new.eqmag)

}


do.move = function(full.current.rates,experiments.per.iter){

	move = sample(1:total.move,size = 1,replace = TRUE)
	print("Move")
	print(move)

	new.eqmag = switch(move,move.all.mags(full.current.rates$eqmag),move.one.mag(full.current.rates$eqmag),move.all.eqs(full.current.rates$eqmag),move.one.eq(full.current.rates$eqmag),conjure.one.mag(full.current.rates$eqmag),conjure.all.mags(full.current.rates$eqmag), conjure.one.eq(full.current.rates$eqmag),conjure.all.eqs(full.current.rates$eqmag))
	#	new.rates = move.all.rates(full.current.rates$eqmag)
	#	new.rates = move.one.rate(full.current.rates$eqmag)

	print("old likelihood")
	print(full.current.rates$likelihood)

	

	new.like = log.likelihood.main(new.eqmag,experiments.per.iter)

    new.like.ast = new.like

    if(move %in% 5:8){
        old.pick = dt(full.current.rates$eqmag, df = df, log = T)
        new.pick = dt(new.eqmag,df = df, log = T)
        new.like.ast = new.like.ast+sum(old.pick[is.finite(old.pick)])-sum(new.pick[is.finite(new.pick)])
    }


	print("new like")
	print(new.like)

	print(paste("Distance",sqrt(sum((new.eqmag-full.current.rates$eqmag)^2))))

	if(accept.test(new.like.ast,full.current.rates$likelihood)){
		return(list(eqmag = new.eqmag, likelihood = new.like, accepted = TRUE,move = move))
	}
	else{
		return(list(eqmag = full.current.rates$eqmag,likelihood = full.current.rates$likelihood,accepted = FALSE,move = move))
	}

}

special.abind = function(arr1,arr2){
	
	return(abind(arr1,arr2,along = 3))
	
}

