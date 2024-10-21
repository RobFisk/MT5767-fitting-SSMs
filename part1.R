# simulating a BAS stochastic model by creating functions for each sub-process

# survival
survival <- function(pop, phi){
    newpop <- rep(0, length(pop))
    for (i in 1:length(pop)){
        newpop[i] <- rbinom(1, pop[i], phi)
    }
    return(newpop)
}

# ageing
ageing <- function(pop){
    newpop <- rep(0, length(pop))
    for (i in 1:(length(pop)-1)){
        newpop[i+1] <- pop[i]
    }
    return(newpop)
}

# reproduction
reproduction <- function(pop, rho){
    births <- 0
    for (i in 1:length(pop)){
        births <- births + rpois(1, pop[i]*rho[i])
    }
    pop[1] <- births
    return(pop)
}

# additional function for detection
detection <- function(pop, p){
    detected <- rep(0, length(pop))
    for (i in 1:length(pop)){
        detected[i] <- rbinom(1, pop[i], p)
    }
    return(detected)
}


# putting it all together to observe dynamics
model <- function(N0, phi, rho, p, nyears=25){
    pop <- matrix(, nrow=nyears, ncol=length(N0))
    obs <- matrix(, nrow=nyears, ncol=length(N0))
    pop[1,] <- N0
    obs[1,] <- detection(pop[1,], p)
    for (i in 2:nyears){
        pop[i,] <- pop[i-1,] |> 
            survival(phi) |> 
            ageing() |> 
            reproduction(rho)
        obs[i,] <- pop[i,] |> detection(p)
    }

    # return results as dataframe
    dat <- data.frame(
        year=1:nyears,
        pop=pop,
        obs=obs
    )
    return(dat)
}

# testing
n0 <- c(120, 70, 50, 30)
phi <- c(0.45, 0.7, 0.7, 0.0)
rho <- c(0.0, 0.9, 1.9, 0.0)
p <- 0.5

sim <- model(n0, phi, rho, p, 25)
print(sim)