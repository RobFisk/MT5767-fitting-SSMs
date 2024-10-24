source("part1.r")
library(ggplot2)
library(gridExtra)

# parameters
n0 <- c(120, 70, 50, 30)
phi <- c(0.45, 0.7, 0.7, 0.0)
rho <- c(0.0, 0.9, 1.9, 0.0)
p <- 0.5

# simulation many times to get ranges
nsim <- 1000
nyears <- 25
totals <- matrix(nrow=nsim, ncol=nyears)
esttotals <- matrix(nrow=nsim, ncol=nyears)
for (i in 1:nsim){
    dat <- model(n0, phi, rho, p, nyears)
    totals[i,] <- dat$pop_total
    esttotals[i,] <- dat$obs_total
}

# summarise
poptotals <- data.frame(
    year = 1:nyears,
    median = apply(totals, 2, median),
    lower = apply(totals, 2, quantile, probs=0.025),
    upper = apply(totals, 2, quantile, probs=0.975)
)
estimates <- data.frame(
    year = 1:nyears,
    median = apply(esttotals, 2, median),
    lower = apply(esttotals, 2, quantile, probs=0.025),
    upper = apply(esttotals, 2, quantile, probs=0.975)
)

# boxplots using ggplot
p <- ggplot(poptotals, aes(x=year)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey", alpha=0.5) +
    geom_line(aes(y=median), color="blue") +
    geom_ribbon(data=estimates, aes(ymin=lower, ymax=upper), fill="grey", alpha=0.5) +
    geom_line(data=estimates, aes(y=median), color="red") +
    labs(title="Population and estimated population over time", x="Year", y="Number of individuals") +
    theme_minimal() 

ggsave(file="plots/ul_quartiles.png", arrangeGrob(p), width=16, height=9)