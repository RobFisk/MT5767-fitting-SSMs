library(statsecol)
library(gridExtra)
library(ggplot2)
library(tidyverse)
head(wildebeest)

rain_rk <- function(pars, years, removals, Nhat, SEhat, rain, model = 1, type = "nll"){
  
  #parameter set up
  N0 <- exp(pars[1])
  N <- numeric(years)
  r <- numeric(years)
  k <- numeric(years)
  N[1] <- N0
  r[1] <- NA
  k[1] <- NA
  
  if(model == 1){ #constant r, k depending on rainfall in previous year
    r[2:years] <- rep(exp(pars[2]),years-1)
    k[2:years] <- exp(pars[3]+pars[4]*rain[1:years-1])
  }
  if(model == 2){ #constant r, k depending on rainfall in current year
    r[2:years] <- rep(exp(pars[2]),years-1)
    k[2:years] <- exp(pars[3]+pars[4]*rain[2:years])
  }
  if(model == 3){ #constant k, r depending on rainfall in current year
    r[2:years] <- exp(pars[2]+pars[3]*rain[2:years])
    k[2:years] <- rep(exp(pars[4]),years-1)
  }
  if(model == 4){ #constant k, r depending on rainfall in previous year
    r[2:years] <- exp(pars[2]+pars[3]*rain[1:years-1])
    k[2:years] <- rep(exp(pars[4]),years-1)
  }
  
  #generate population dynamics:
  for(i in 2:years){
    N[i]=N[i-1] + r[i] * N[i-1] * (1-N[i-1]/k[i]) - removals[i-1]
  }
  
  negloglik <- -sum(dnorm(Nhat,N,SEhat,log=TRUE), na.rm=TRUE)
  
  #what should be returned? 
  if(type=="nll"){  return(negloglik)}
  if(type=="proj"){ return(N)}
}

yrs <- nrow(wildebeest)
rmv <- wildebeest$Catch
Nhat <- wildebeest$Nhat
SEhat <- wildebeest$sehat
rain <- wildebeest$rain

#Optimization of parameters for each model####

fit_1 <- optim(par = c(log(0.1),log(0.25),log(1.5), 0), fn = rain_rk, years = yrs, 
               removals = rmv, Nhat = Nhat, SEhat = SEhat, rain = rain, model=1)

fit_2 <- optim(par = c(log(0.1), log(0.25), log(1.5), 0), fn = rain_rk, years = yrs, 
               removals = rmv, Nhat = Nhat, SEhat = SEhat, rain = rain, model=2)

fit_3 <- optim(par = c(log(0.1),log(0.25), 0, 0), fn = rain_rk, years = yrs, 
               removals = rmv, Nhat = Nhat, SEhat = SEhat, rain = rain, model=3)

fit_4 <- optim(par = c(log(0.1), log(0.25), 0, 0), fn = rain_rk, years = yrs, 
               removals = rmv, Nhat = Nhat, SEhat = SEhat, rain = rain, model=4)

#Alpha 1 and Beta 1 comparisons

Model <- c("Constant r. Current year k", "Constant r previous year k")
Model2 <- c("Constant k. Current year r", "Constant k previous year r")
#beta 1 results table
k_table <- tibble(
  Model = Model,
  beta_1 = c(fit_2$par[4], fit_1$par[4])
)

r_table <- tibble(
  Model = Model2,
  alpha_1 = c(fit_3$par[3], fit_4$par[3])
)

#Alpha 1 and beta 1 outputs for different modelling scenarios####
print(k_table)
print(r_table)

#AIC comparison####

aic1 <- 2*fit_1$value + 2*length(fit_1$par)
aic2 <- 2*fit_2$value + 2*length(fit_2$par)
aic3 <- 2*fit_3$value + 2*length(fit_3$par)
aic4 <- 2*fit_4$value + 2*length(fit_4$par)


#Table of dAIC and AIC values for each model####
aictab <- data.frame(
  Model = c("Constant r. Previous year k","Constant r. Current year k","Constant k. Current year r","Constant k. Previous year r"),
  AIC = c(aic1,aic2,aic3,aic4),
  dAIC = c(aic1,aic2,aic3,aic4)-min(c(aic1,aic2,aic3,aic4)))
aictab[order(aictab$dAIC),]

#population projection####

proj_1 <- rain_rk(fit_1$par, years = yrs, removals = rmv, Nhat = Nhat, 
                  SEhat = SEhat, rain = rain, model=1, type="proj")
proj_2 <- rain_rk(fit_2$par, years = yrs, removals = rmv, Nhat = Nhat, 
                  SEhat = SEhat, rain = rain, model=2, type="proj")
proj_3 <- rain_rk(fit_3$par, years = yrs, removals = rmv, Nhat = Nhat, 
                  SEhat = SEhat, rain = rain, model=3, type="proj")
proj_4 <- rain_rk(fit_4$par, years = yrs, removals = rmv, Nhat = Nhat, 
                  SEhat = SEhat, rain = rain, model=4, type="proj")

pred_df <- data.frame(years = rep(wildebeest$year,4),
                      N = c(proj_1,proj_2, proj_3, proj_4),
                      Model=rep(aictab$Model,each=nrow(wildebeest)))

#Table of predictions for when r varies by Rainfall
pred_df_r <- pred_df %>%
  filter(Model %in% c("Constant k. Current year r", "Constant k. Previous year r" ))

#Table of predictions for when k varies by Rainfall
pred_df_k <- pred_df %>%
  filter(Model %in% c("Constant r. Current year k", "Constant r. Previous year k" ))


#plots of changes ina abundance in different modelling scenarios####
pred_plot_r <- ggplot(wildebeest, aes(x=year, y=Nhat)) +
  geom_errorbar(aes(ymin=lci,ymax=uci), width=0) +
  geom_point(size=3) +
  geom_line(data=pred_df_r, aes(x=years,y=N,color=Model,group=Model),size=0.8) +
  ylim(0,2.1) + ylab("Abundance (millions)") + xlab("Year") 

pred_plot_k <- ggplot(wildebeest, aes(x=year, y=Nhat)) +
  geom_errorbar(aes(ymin=lci,ymax=uci), width=0) +
  geom_point(size=3) +
  geom_line(data=pred_df_k, aes(x=years,y=N,color=Model,group=Model),size=0.8) +
  ylim(0,2.1) + ylab("Abundance (millions)") + xlab("Year") 

grid.arrange(pred_plot_k, pred_plot_r, ncol = 2)
