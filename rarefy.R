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
