library(vegan)
library (plotfunctions) ## alpha lines within plots
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#library(bayesplot)
#library(tidybayes)
library(tidyverse)        # ggplot, dplyr, %>%, and friends
#library(brms)             # Bayesian modeling through Stan
library(tidybayes)        # Manipulate Stan objects in a tidy way
#library(emmeans)          # Calculate marginal effects in even fancier ways
library(brms)
library(easyGgplot2)
library(gridGraphics)

## set working directory

setwd("~/Documentos/R_ejercicio/montana_tesis/chapter_1")

## loading insect assemblages matrix
## all matrices (Euclidean distance, distance over river network, assemblages) have sites (rows) in same order 
whole <- read.csv ("whole.csv") ## whole assemblage
whole <- whole[,-1]
non_sh <- read.csv("non_shredders.csv") ## non-shredder assemblage
non_sh <- non_sh[,-1]
non_sh<- non_sh[-8,]
sh <- read.csv("shredders.csv")
sh <- sh[,-1]
sh <- sh[-8,]

## sites
coord <- read.csv ("coordenadas.csv")
elevation<-(coord$Elevation)
elevation <-elevation[-8]


## shredder
(raremax_sh <- min(rowSums(sh)))
max(rowSums(sh))
Srare_sh <- rarefy(sh, raremax_sh)

sh_rep <- rep(list(sh), 100)

z<- list()
for (i in 1: 100){
  z[[i]]<- matrix(1:15)
}

for(i in seq_along(z)){
  z[[i]]<- rrarefy(sh_rep[[i]], raremax_sh,MARGIN = 1)
}

r_sh <- as.data.frame.list(z) ## 
colnames(r_sh) <- seq(1:100)
r_sh$FFG <- c("sh")
r_sh$elevation <- elevation

rich_r_sh <- pivot_longer(r_sh, cols=1:100, names_to = "run", values_to = "r_rich")
mean(rich_r_sh$r_rich[rich_r_sh$elevation == 156])
sd(rich_r_sh$r_rich [rich_r_sh$elevation == 156])

rrich_sh_mean <- rich_r_sh%>%
                  group_by(elevation)%>%
                  summarise(mean =mean(r_rich))
  
rrich_sh_sd <- rich_r_sh%>%
  group_by(elevation)%>%
  summarise(sd = sd(r_rich))


summarize (rich_r_sh, mean = mean(r_rich))

rich_r_sh$r_rich <- (round(rich_r_sh$r_rich, 0))

## non-shredder

S_sh <- specnumber(sh)
min(S_sh)
max(S_sh)
S_non <- specnumber(non_sh) # observed number of species

min(S_non)
max(S_non)

(raremax_non <- min(rowSums(non_sh)))
min(rowSums(non_sh))
max(rowSums(non_sh))
Srare_non <- rarefy(non_sh, raremax_non)

plot(S_non, Srare_non, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)

nosh_rep <- rep(list(non_sh), 100)

z<- list()
for (i in 1: 100){
  z[[i]]<- matrix(1:15)
}

for(i in seq_along(z)){
  z[[i]]<- rarefy(nosh_rep[[i]], raremax_non)
}


r_nosh <- as.data.frame.list(z) ## 
colnames(r_nosh) <- seq(1:100)
r_nosh$FFG <- c("nosh")
r_nosh$elevation <- elevation

rich_r_nosh <- pivot_longer(r_nosh, cols=1:100, names_to = "run", values_to = "r_rich")
#rich_r_nosh$r_rich <- (round(rich_r_nosh$r_rich, 0))


## put together as one matrix

rrich_matrix <- rbind (rich_r_sh, rich_r_nosh)



rrich_sites <- as.data.frame(Srare_non)
rrich_sites$elevation <- elevation
rrich_sites$FFG <- c("nonsh")
colnames(rrich_sites) <- c("r_rich", "elevation", "FFG")

rrich_sites_sh <- as.data.frame(Srare_sh)
rrich_sites_sh$elevation <- elevation
rrich_sites_sh$FFG <- c("sh")
colnames(rrich_sites_sh) <- c("r_rich", "elevation", "FFG")

rrich_sites <- rbind (rrich_sites, rrich_sites_sh)
rrich_sites$r_rich <- round(rrich_sites$r_rich, 0)


## brms models
priors <- c(set_prior("normal( 0.005, 0.005)", class = "b"))
weakly_prior <- prior(normal(0,3))
fit_rich <- brm(r_rich ~ elevation*FFG, data = rrich_sites, chains =4, iter = 6000, warmup = 5000, family = poisson, prior = weakly_prior )

fit_rich <- brm(r_rich ~ elevation*FFG, data = rrich_sites, chains =4, iter = 6000, warmup = 5000, family = poisson)

plot(fit_rich)
print(fit_rich)
plot(conditional_effects(fit_rich), points = TRUE)


###
posterior_b <- fit_rich %>%
  gather_draws(`b_.*`, regex = TRUE)

x <- seq(156, 2720, by= 5 ) ## distance vector

mu_mu_non <- mean(posterior_b$.value[posterior_b$.variable == "b_Intercept"]) + mean(posterior_b$.value[posterior_b$.variable == "b_elevation"]) * x

mu_mu_sh <- (mean(posterior_b$.value[posterior_b$.variable == "b_Intercept"]) +  mean(posterior_b$.value[posterior_b$.variable == "b_FFGsh"]))+ 
  (mean(posterior_b$.value[posterior_b$.variable == "b_elevation"]) +mean(posterior_b$.value[posterior_b$.variable == "b_elevation:FFGsh"]))* x


## mean and quantiles
slope_non_sh <- as.data.frame(c(mean(posterior_b$.value[posterior_b$.variable == "b_elevation"]),
                                quantile(posterior_b$.value[posterior_b$.variable == "b_elevation"], 
                                         probs= c(0.025,0.975))))
slope_non_sh <- t(slope_non_sh)  
slope_non_sh$FGG <- c("no-sh")
slope_non_sh <- as.data.frame(slope_non_sh)
colnames(slope_non_sh)<- c("mean", "lowCI","UpCI", "FFG" )

slope_sh <- as.data.frame(c((mean(posterior_b$.value[posterior_b$.variable == "b_elevation"])+mean(posterior_b$.value[posterior_b$.variable == "b_elevation:FFGsh"])),
                            quantile(posterior_b$.value[posterior_b$.variable == "b_elevation"]+ posterior_b$.value[posterior_b$.variable == "b_elevation:FFGsh"], 
                                     probs= c(0.025,0.975))))
slope_sh <- t(slope_sh)  
slope_sh$FFG <- c("sh")
slope_sh <- as.data.frame(slope_sh)
colnames(slope_sh)<- c("mean", "lowCI","UpCI", "FFG" )

slope_ddr_fit <- rbind(slope_sh, slope_non_sh)
slope_ddr_fit$Index <- c(1.5, 1)
plot(slope_ddr_fit$Index, slope_ddr_fit$mean)

plot(slope_ddr_fit$Index, slope_ddr_fit$mean, ylim = rev(c(0.001, -0.001)), xlim=c(0.75,1.75), xaxt= "n", 
     pch=16, cex=1.5, ylab= "Estimated slope", xlab="", cex.lab= 1.5)
mtext("Non-Shredder",cex = 1.5, side = 1, at=1)
mtext("Shredder",cex = 1.5, side = 1, at=1.5)
arrows(slope_ddr_fit$Index, x1= (slope_ddr_fit$Index), y0= (slope_ddr_fit$lowCI), y1= (slope_ddr_fit$UpCI) , code=0, angle=90, length=0.15, lwd =2, lty=1);

### plotting ggplot violins

## Slopes violins
slope_non_sh_v <- as.data.frame(c((posterior_b$.value[posterior_b$.variable == "b_elevation"])))
colnames(slope_non_sh_v) <- c("Value")
slope_non_sh_v$Group <- c("Non-shredder")

slope_sh_v <- as.data.frame(c((posterior_b$.value[posterior_b$.variable == "b_elevation"])
                              + (posterior_b$.value[posterior_b$.variable == "b_elevation:FFGsh"])))

colnames(slope_sh_v)<- c("Value")
slope_sh_v$Group <- c("Shredder")

slope_rich <- rbind (slope_non_sh_v, slope_sh_v)


b1_plot <- ggplot2.violinplot(data=slope_rich, xName='Group',yName="Value",
                              groupName='Group',mainTitle="Slope" ,groupColors=c('grey','grey'))+
  labs(y="Estimated parameter", x='')+
  #scale_x_discrete(labels=c('Los Musgos', 'Machacas'))+
  theme(legend.position = "none", axis.line = element_line(linewidth =1, colour="black"),
        #axis.text.x = element_blank(),
        axis.line.y.left = element_line(linewidth = 0.5, colour = "black"),
        axis.line.x.bottom = element_line(linewidth = 0.5, colour = "black"),
        panel.background = element_rect(fill = "white"),
        #strip.background = element_rect(colour = "white", fill = "white"),
        #strip.text = element_text(size = 18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face=""),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.title = element_text(hjust = 0.5, size= 18)
  )

b1_plot

## ggplot mean slope and 1000 lines from posterior distribution
##filtering data

rrich_nonsh <- filter (rrich_sites, FFG == "nonsh") 

rrich_sites$r_rich [rrich_sites$FFG== "nonsh"] ~ rrich_sites$elevation[rrich_sites$FFG== "nonsh"]

mu_non_data <- as.data.frame(cbind(x, mu_mu_non))
mu_non_data$expo <- exp(mu_non_data$mu_mu_non)

posterior_b_1000 <-  fit_rich%>%
  gather_draws(`b_.*`, regex = TRUE, ndraws = 1000)

mu_mu_non_1000 <- mean(posterior_b_1000$.value[posterior_b_1000$.variable == "b_Intercept"]) + mean(posterior_b_1000$.value[posterior_b_1000$.variable == "b_elevation"]) * x
mu_non_data_1000 <- as.data.frame(cbind(x, mu_mu_non_1000))
mu_non_data_1000$expo_1000 <- exp(mu_non_data_1000$mu_mu_non_1000)

for(i in 1:1000){
  mu<- posterior_b$.value[posterior_b$.variable == "b_Intercept"][i] + posterior_b$.value[posterior_b$.variable == "b_elevation"][i]*x
}

mu_100_non <- as.data.frame (cbind(x,mu))

ggplot(data = mu_100_non, aes(x, mu))+
  geom_line(size=1, color= "lightblue")
  geom_line(data= mu_non_data, aes(x, expo), color = "black")

non_sh_ggplot <- ggplot(data= mu_non_data, aes(x, expo))+
  ylim(0, 10)+
  xlim(156,2720)+
  labs(y="Rarefied richness", x='Elevation (m asl)')+
  geom_line(size=1.5)+
  geom_line(data= mu_non_data_1000 ,aes(x, expo_1000),
            alpha= 0.9, color= "lightblue")+
  geom_point(data= rrich_nonsh, aes(elevation, r_rich), size= 3)+
  theme(legend.position = "none", axis.line = element_line(linewidth =1, colour="black"),
        #axis.text.x = element_blank(),
        axis.line.y.left = element_line(linewidth = 0.8, colour = "black"),
        panel.background = element_rect(fill = "white"),
        #strip.background = element_rect(colour = "white", fill = "white"),
        #strip.text = element_text(size = 18),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=14),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5, size= 18)
  )
  
non_sh_ggplot





## Plotting

par(mfrow= c(1,1), mar =c(5.1,6,2,1.5))
plot(rrich_sites$r_rich [rrich_sites$FFG== "nonsh"] ~ rrich_sites$elevation[rrich_sites$FFG== "nonsh"], pch=16, cex= 1.5, ylim=c(0,10), xlim = c(156,2720), 
     cex.lab= 2, cex.axis = 1.5, xlab= "Elevation (m asl)"
     ,ylab="Rarefied richness",
    
)
mtext("Non-shredder", cex = 1.5, side = 3)

#text ( 65, 0.8, "Non-Shredder assemblage", cex = 1.5 )
for(i in 1:1000){
  mu<- posterior_b$.value[posterior_b$.variable == "b_Intercept"][i] + posterior_b$.value[posterior_b$.variable == "b_elevation"][i]*x
  lines(x, exp(mu), col = alpha("lightblue", 0.1) )
}
points(rrich_sites$r_rich [rrich_sites$FFG== "nonsh"] ~ rrich_sites$elevation[rrich_sites$FFG== "nonsh"], pch=16, cex= 1.5, col = "black")
lines(x, exp(mu_mu_non), col= "black", lwd=2)

#plot.new();
plot(rrich_sites$r_rich [rrich_sites$FFG== "sh"] ~ rrich_sites$elevation[rrich_sites$FFG== "sh"], pch=16, cex= 1.5, ylim=c(0,10), xlim = c(156,2720), cex.lab= 1.5, xlab= "Elevation (m asl)"
     ,ylab= "", cex.lab= 2, cex.axis = 1.5
     #,ylab="Rarefied richness"
)

mtext("Shredder", cex = 1.5, side = 3)
#text ( 65, 0.8, "Shredder assemblage", cex = 1.5 )
for(i in 1:1000){
  mu<- (posterior_b$.value[posterior_b$.variable == "b_Intercept"][i] +  posterior_b$.value[posterior_b$.variable == "b_FFGsh"][i])+ 
    (posterior_b$.value[posterior_b$.variable == "b_elevation"][i] + posterior_b$.value[posterior_b$.variable == "b_elevation:FFGsh"][i])* x 
    lines(x, exp(mu), col = alpha("lightblue", 0.1) )
}
points(rrich_sites$r_rich [rrich_sites$FFG== "sh"] ~ rrich_sites$elevation[rrich_sites$FFG== "sh"], pch=16, cex= 1.5, col= "black")
lines(x,exp(mu_mu_sh), col= "black", lwd = 2)

plot(slope_ddr_fit$Index, slope_ddr_fit$mean, ylim = rev(c(0.001, -0.001)), xlim=c(0.75,1.75), xaxt= "n", 
     pch=16, cex=1.5, ylab= "Estimated parameter", xlab="", cex.lab= 2, cex.axis = 1.5, )
mtext("Slope", cex = 1.5, side = 3)
mtext("Non-shredder",cex = 1.5, side = 1, at=1, line=3)
mtext("Shredder",cex = 1.5, side = 1, at=1.5, line = 3)
arrows(slope_ddr_fit$Index, x1= (slope_ddr_fit$Index), y0= (slope_ddr_fit$lowCI), y1= (slope_ddr_fit$UpCI) , code=0, angle=90, length=0.15, lwd =2, lty=1);
abline(h=0, lwd=2, lty=2, col="red")


### Machacas communities

machacas <- read.csv("rrarefy/nonsh/machacas_non.csv", row.names = 1)
machacas<- t(machacas)
r_machacas<- rrarefy(machacas, sample = 10)

mach_1 <- specaccum(machacas, method="random")

sp4 <- specaccum(machacas, "rarefaction")

plot(sp4, xvar = "individuals", xlim= c(0,2000))
mach_mm<- fitspecaccum(sp4, method = "rarefaction", model = "michaelis-menten" )
fitted(mach_mm)

plot(sp4)
lines(mach_mm, col=2, lwd=2)
plot(mach_mm, col=2, lwd=2)

fit(sp4$individuals, sp4$richness)

rarecurve(machacas)



data(BCI)
S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI)))
data(BCI)

Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)


sp1 <- specaccum(BCI)
sp2 <- specaccum(BCI, "random")
sp2
sp3 <- specaccum(BCI, "rarefaction")
summary(sp2)

plot(sp3)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="yellow", add=TRUE, pch="+")
## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, "lomolino")
coef(mod1)
fitted(mod1)
plot(sp1)
## Add Lomolino model using argument 'add'
plot(mod1, add = TRUE, col=2, lwd=2)
## Fit Arrhenius models to all random accumulations
mods <- fitspecaccum(sp2, "arrh")
plot(mods, col="hotpink")
boxplot(sp2, col = "yellow", border = "blue", lty=1, cex=0.3, add= TRUE)
## Use nls() methods to the list of models
sapply(mods$models, AIC)



bci_rep <- rep(list(BCI), 100)





q<- list()
for (i in 1: 100){
  q[[i]]<- matrix(1:50)
}

for(i in 1:100){
  q[[i]]<- rarefy(bci_rep[[i]], 300,MARGIN = 1)
}

for(i in seq_along(z)){
  z[[i]]<- rarefy(sh_rep[[i]], raremax_sh,MARGIN = 1)
}

r_sh <- as.data.frame.list(z) ## 
colnames(r_sh) <- seq(1:100)
r_sh$FFG <- c("sh")
r_sh$elevation <- elevation

rich_r_sh <- pivot_longer(r_sh, cols=1:100, names_to = "run", values_to = "r_rich")
mean(rich_r_sh$r_rich[rich_r_sh$elevation == 156])
sd(rich_r_sh$r_rich [rich_r_sh$elevation == 156])

rrich_sh_mean <- rich_r_sh%>%
  group_by(elevation)%>%
  summarise(mean =mean(r_rich))

rrich_sh_sd <- rich_r_sh%>%
  group_by(elevation)%>%
  summarise(sd = sd(r_rich))


Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)
as.data.frame(Srare)

###

## set working directory

setwd("~/Documentos/R_ejercicio/montana_tesis/chapter_1")

kixpur <- read.csv("rrarefy/kixpur_nonsh.csv")
kixpur <- kixpur[,-1]
##Loading data
## sites
coord <- read.csv ("coordenadas.csv")

## loading insect assemblages matrix
## all matrices (Euclidean distance, distance over river network, assemblages) have sites (rows) in same order 
whole <- read.csv ("whole.csv") ## whole assemblage
whole <- whole[,-1]
non_sh <- read.csv("non_shredders.csv") ## non-shredder assemblage
non_sh <- non_sh[,-1]
non_sh<- non_sh[-8,]
sh <- read.csv("shredders.csv")
sh <- sh[,-1]
sh <- sh[-8,]

## non-shredder
S_non <- specnumber(non_sh) # observed number of species
(raremax_non <- min(rowSums(non_sh)))
min(rowSums(non_sh))
max(rowSums(non_sh))
Srare_non <- rarefy(non_sh, raremax_non)

plot(S_non, Srare_non, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)



## shredder

S_sh <- specnumber(sh) # observed number of species
(raremax_sh <- min(rowSums(sh)))
Srare_sh <- rarefy(sh, raremax_sh)

min(rowSums(sh))
max(rowSums(sh))
plot(S_sh, Srare_sh, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)

rarecurve(non_sh, step = 6, sample = raremax_non, col = "blue", cex = 0.6)
rarecurve(sh, step = 6, sample = raremax_sh, col = "blue", cex = 0.6)

elevation<-(coord$Elevation)
elevation <-elevation[-8]

non_sh_rrarefy <- as.data.frame(cbind(elevation, Srare_non, Srare_sh))

non_sh_rrarefy$Srare_non_r <- (round(non_sh_rrarefy$Srare_non, 0))

non_sh_rrarefy$Srare_sh_r <- (round(non_sh_rrarefy$Srare_sh, 0))

priors <- c(set_prior("normal(0.005, 1)", class = "b"))
fit_non<- brm(Srare_non_r~ elevation, data = non_sh_rrarefy, chains =4, iter = 4000, warmup = 3000, family = poisson, prior = priors )
plot(fit_non)
print(fit_non)
plot(conditional_effects(fit_non))

fit_sh <- brm(Srare_sh_r~ elevation, data = non_sh_rrarefy, chains =4, iter = 4000, warmup = 3000, family = poisson, prior = priors )
plot(fit_sh)
print(fit_sh)
plot(conditional_effects(fit_sh))


## reaches
S_kixpur <- specnumber(kixpur)
(raremax_kixpur <- min(rowSums(kixpur)))
Srare_kixpur <- rarefy(kixpur, raremax_kixpur)

plot(S_kixpur, Srare_kixpur, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
## shredder

S_sh <- specnumber(sh) # observed number of species
(raremax_sh <- min(rowSums(sh)))
Srare_sh <- rarefy(sh, raremax_sh)

min(rowSums(sh))
max(rowSums(sh))
plot(S_sh, Srare_sh, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)

rarecurve(non_sh, step = 6, sample = raremax_non, col = "blue", cex = 0.6)
rarecurve(sh, step = 6, sample = raremax_sh, col = "blue", cex = 0.6)

elevation<-(coord$Elevation)
elevation <-elevation[-8]

non_sh_rrarefy <- as.data.frame(cbind(elevation, Srare_non, Srare_sh))

non_sh_rrarefy$Srare_non_r <- (round(non_sh_rrarefy$Srare_non, 0))

non_sh_rrarefy$Srare_sh_r <- (round(non_sh_rrarefy$Srare_sh, 0))

priors <- c(set_prior("normal(0.005, 1)", class = "b"))
fit_non<- brm(Srare_non_r~ elevation, data = non_sh_rrarefy, chains =4, iter = 4000, warmup = 3000, family = poisson, prior = priors )
plot(fit_non)
print(fit_non)
plot(conditional_effects(fit_non))

fit_sh <- brm(Srare_sh_r~ elevation, data = non_sh_rrarefy, chains =4, iter = 4000, warmup = 3000, family = poisson, prior = priors )
plot(fit_sh)
print(fit_sh)
plot(conditional_effects(fit_sh))


plot(coord$Elevation, Srare_non, ylim= c(0,9), xlab = "Elevation (m asl)", ylab = "Rarefied No. of Species", pch= 16, main = "Non-shredder")
abline(lm(Srare_non~coord$Elevation))
brms::brm(Srare_non~coord$Elevation, chains =4, iter = 4000, warmup = 3000, family = poisson, prior = NULL )

plot(coord$Elevation, Srare_sh, ylim= c(0,3), xlab = "Elevation (m asl)", ylab = "Rarefied No. of Species", pch= 16, main = "Shredder")
abline(lm(Srare_sh~coord$Elevation))

## slope
get_variables(fit_non)

b_nonshredder <-fit_non %>%
  spread_draws(b_V1)

b_shredder <-fit_sh %>%
  spread_draws(b_V1)

mean((b_nonshredder$b_V1))
quantile(b_nonshredder$b_V1, probs= c(0.025,0.975))

mean((b_shredder$b_V1))
quantile(b_shredder$b_V1, probs= c(0.025,0.975))


par(mfrow= c(2,2), mar =c(4,5,1.5,1.5))
plot(non_sh_rrarefy$elevation, S_non, xlab = "", cex.lab= 2, cex.axis=1.5, cex= 1.5,
     ylab = "Observed richness (Taxa number)", pch= 16, main = "Non-shredder", cex.main=2)
plot(non_sh_rrarefy$elevation, S_sh, xlab = "", cex.lab= 2,cex.axis=1.5,cex= 1.5,
     ylab = "", 
     pch= 16, main = "Shredder", cex.main=2)
plot(non_sh_rrarefy$elevation, non_sh_rrarefy$Srare_non_r, cex.lab= 1.5,cex.axis=1.5, cex= 1.5,ylim=c(0,9),xlab = "Elevation (m asl)", ylab = "Rarefied richness (Taxa number)", pch= 16)
plot(non_sh_rrarefy$elevation, non_sh_rrarefy$Srare_sh_r, cex.lab= 1.5,cex.axis=1.5,cex= 1.5,ylim=c(0,9),xlab = "Elevation (m asl)", 
     ylab = "", 
     pch= 16)



## bayesian model in rstan format

sink("enviro_dist.stan")
cat("
    data {
    
    int <lower=1> N; //number of data points
    vector [N] S; // rarefied No. Species
    vector [N] D; // Elevation of sites
    
    }
    
    parameters {
    real a; // intercept
    real b; // slope
   	real <lower=0> sigma;
    }
    
    model { 
    
    //priors. 
    a ~ uniform (0 , 10);
    b ~ normal(0 , 10);
    
    //likelihood    	
    E ~ normal( a + b*D,sigma);
    }
    
    "
    ,fill=TRUE)
sink()

#make list for stan
env_Dist_data<- list(E= matrix_complete$eucl_envi , D=matrix_complete$Dist_km, N=length(matrix_complete$eucl_envi))

env_dist_fit<-stan(file='enviro_dist.stan', data = env_Dist_data, 
                   iter = 4000, chains = 4, warmup = 3000)
print(env_dist_fit, digits_summary = 6)
traceplot(env_dist_fit)
plot(env_dist_fit)

## Plotting fitted model

extract_env_dist_fit <- extract (env_dist_fit)
print(names(extract_env_dist_fit))

a <- (mean(extract_env_dist_fit$a))
b <- mean(extract_env_dist_fit$b)
x <- seq(0, 133, length.out=100000)
fitted_mean <- a + b*x

bsteps<- rstan::extract(env_dist_fit, pars = "b")
asteps<- rstan::extract(env_dist_fit, pars = "a")

plot(matrix_complete$eucl_envi ~ matrix_complete$Dist_km)
par(mai=c(0.7,0.7,0.1,0.1), mgp=c(2,1,0))
plot(matrix_complete$eucl_envi~ matrix_complete$Dist_km, xlab= "Distance between sites (km)",
     ylab= "Enviromental distance (Euclidean distance)", pch = 16, cex= 1.5, col="blue",
     cex.lab = 1.5, xlim =c(0,133), ylim=c(0,9)
)

