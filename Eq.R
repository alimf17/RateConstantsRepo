library(pracma)


# Summary: this method calculates the theoretical equilibria of a set of transition PROBABILITY matrices
#
# Parameters:
#       transitions: an SxS matrix M of floats between 0 and 1 where M[i,j] is the 
#                    probabilty of transitioning from state j to state i over
#                    the course of the experiment, OR an array of T such matrices
#
# Returns: An Sx1 matrix of the theoretical equilibiria of the transitions
calcEQs = function(transitions) {

    dimen = dim(transitions)

    mats = 1

    if(length(dimen) == 3) {mats = dimen[3]}

    print("dfsa")
    print(transitions)
    print((dim(transitions)))
    print(mats)

    eqs = matrix(nrow =(dim(transitions)[1]),ncol = mats)

    modtrans = array(transitions,c(dim(transitions)[1:2],mats))

    for(k in 1:mats){

        cleanup = modtrans[,,k]
        for(j in 1:ncol(cleanup)){
            cleanup[,j] = cleanup[,j]/sum(cleanup[,j])
        }

        eigs = eigen(cleanup)

        eqs[,k] = (1/sum(eigs$vector[,abs(eigs$values-1)<1e-6]))*eigs$vector[,abs(eigs$values-1)<1e-6]

    }

    return(apply(eqs,1,mean))

}

# Summary: this method calculates the theoretical equilibria of a set of transition RATE matrices
#
# Parameters:
#       transitions: an SxS matrix M of floats where M[i,j] is the 
#                    rate constant of the reaction from state j to state i
#                    in units 1/sec, OR an array of T such matrices
#
# Returns: An Sx1 matrix of the theoretical equilibiria of the transition rates
theoreticalratesEq = function(rates){

    dimen = dim(rates)

    mats = 1

    if(length(dimen) == 3) {mats = dimen[3]}


    eqs = matrix(nrow =dimen[1],ncol = mats)

    modtrans = array(rates,c(dimen[1:2],mats))


    for(k in 1:mats){

        cleanup = modtrans[,,k]
        for(j in 1:ncol(cleanup)){
            cleanup[j,j] = -sum(cleanup[,j])
        }

        cleanup = cbind(cleanup, rep(0, dimen[1]))

        cleanup[sample.int(dimen[1],1),] = rep(1,dimen[2]+1)

        eqs[,k] = rref(cleanup)[,dimen[2]+1] 

        print(rref(cleanup))


    }

        return(list(dat=eqs, eq=rowMeans(eqs)))
}



# Summary: this method calculates the chemical rates of transition if the 
#          experimental probabilities of transition over the experimental 
#          time frame truly correspond to single transitions over that time
#          frame
#
# Parameters:
#       transitions: an array of T SxS matrices M of floats between 0 and 1 where M[i,j] is the 
#                    probabilty of transitioning from state j to state i over
#                    the course of the experiment
#
# Returns: An SxSxT array R, where R[i,j,k] = the 1/sec rate constant of going from state j to state i
naiveRates = function(transitions){

    tau = .04
    naives = array(dim=dim(transitions))

    for(k in 1:dim(transitions)[3]){

        for(j in 1:dim(transitions)[2]){

            bigL = -log(transitions[j,j,k])/tau
            naives[,j,k] = bigL*transitions[,j,k]/(1-exp(-bigL*tau))
            naives[j,j,k] = 0


        }

    }

    return(naives)

}



# Summary: this method calculates the theoretical experimental probabilities if each
#          transition reaction can occur at most once over the course of the experiment
#
# Parameters:
#       transitions: an array of T SxS matrices M of floats where M[i,j] is the 
#                     1/sec rate constant of going from state j to state i
#                    the course of the experiment
#
# Returns: An SxSxT array P, where P[i,j,k] = the probability of a single transition
#          from state j to state i from the kth set of rate constants
naiveTrans = function(rates){

    tau = .04
    naives = array(dim=dim(rates))


    for(j in 1:dim(rates)[2]){

        bigL = sum(rates[,j])
        naives[,j] = (rates[,j]/bigL)*(1-exp(-bigL*tau))
        naives[j,j] = exp(-bigL*tau)


    }


    return(naives)

}

