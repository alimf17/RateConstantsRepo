library(pracma)

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



