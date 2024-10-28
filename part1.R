library(ggplot2)
library(gridExtra)
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
        pop[i,] <- pop[i-1,] |> survival(phi) |> ageing() |> reproduction(rho)
        obs[i,] <- pop[i,] |> detection(p)
    }

    # return results as data frame
    dat <- data.frame(
        year = 1:nyears,
        pop_1 = pop[,1], pop_2 = pop[,2], pop_3 = pop[,3], pop_4 = pop[,4],
        pop_total = rowSums(pop),
        obs_1 = obs[,1], obs_2 = obs[,2], obs_3 = obs[,3], obs_4 = obs[,4],
        obs_total = rowSums(obs)
    )
    return(dat)
}

# testing
n0 <- c(120, 70, 50, 30)
phi <- c(0.45, 0.7, 0.7, 0.0)
rho <- c(0.0, 0.9, 1.9, 0.0)
p <- 0.5

sim <- model(n0, phi, rho, p, 25)

# Plotting first year individuals and observation
first_year <- ggplot(sim, aes(x=year)) +
    geom_line(aes(y=pop_1), color="dodgerblue") +
    geom_line(aes(y=obs_1), color="deepskyblue", linetype="dashed") +
    theme_minimal() + 
    labs(title="First year individuals and observations", x="Year", y="Number of individuals")

# Similar for 2nd/3rd/4th year
second_year <- ggplot(sim, aes(x=year)) +
    geom_line(aes(y=pop_2), color="dodgerblue") +
    geom_line(aes(y=obs_2), color="deepskyblue", linetype="dashed") +
    theme_minimal() + 
    labs(title="Second year individuals and observations", x="Year", y="Number of individuals") 

third_year <- ggplot(sim, aes(x=year)) +
    geom_line(aes(y=pop_3), color="dodgerblue") +
    geom_line(aes(y=obs_3), color="deepskyblue", linetype="dashed") +
    theme_minimal() + 
    labs(title="Third year individuals and observations", x="Year", y="Number of individuals")

fourth_year <- ggplot(sim, aes(x=year)) +
    geom_line(aes(y=pop_4), color="dodgerblue") +
    geom_line(aes(y=obs_4), color="deepskyblue", linetype="dashed") +
    theme_minimal() + 
    labs(title="Fourth year individuals and observations", x="Year", y="Number of individuals")

# put plots in a list to produce the output png, this is a bit clunky but it works
plot_list <- list(first_year, second_year, third_year, fourth_year)
ggsave(file="plots/yearlybreakdown.png", arrangeGrob(grobs=plot_list, ncol=2, nrow=2), width=16, height=9)