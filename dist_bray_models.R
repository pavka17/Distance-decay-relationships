## Script to estimate euclidean distances between points
## Author: Pavel Garc?a
## April 12, 2023

## Loading packages

library(topoDistance)
library(raster)
library (riverdist)
library(vegan)
library (tidyr)
#library(dplyr)## summarizing and reshaping data
library (ggplot2)
#library (brms)
library(gridExtra)
#library(broom)
#library(coda)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#library(bayesplot)
#library(tidybayes)
library(tidyverse)        # ggplot, dplyr, %>%, and friends
#library(brms)             # Bayesian modeling through Stan
library(tidybayes)        # Manipulate Stan objects in a tidy way
library(broom)            # Convert model objects to data frames
library(broom.mixed)      # Convert brms model objects to data frames
#library(vdemdata)         # Use data from the Varieties of Democracy (V-Dem) project
#library(betareg)          # Run beta regression models
#library(extraDistr)       # Use extra distributions like dprop()
#library(ggdist)           # Special geoms for posterior distributions
#library(gghalves)         # Special half geoms
#library(ggbeeswarm)       # Special distribution-shaped point jittering
library(ggrepel)          # Automatically position labels
library(patchwork)        # Combine ggplot objects
#library(marginaleffects)  # Calculate marginal effects for frequentist models
library(emmeans)          # Calculate marginal effects in even fancier ways
library(modelsummary)     # Create side-by-side regression tables
library (plotfunctions) ## alpha lines within plots
library(otuSummary) ## transforming all triangle matrices to long matrices
library(bayesplot)
library(gridExtra)
library(ggplotify)
library(patchwork)
library(cowplot)
library(easyGgplot2)
library(gridGraphics)


#library(rasterVis)
#library(rgdal)
#library(rgl)
#library(dismo)

## set working directory

#setwd("~/Documentos/R_ejercicio/montana_tesis/chapter_1")
setwd ("C:/Users/Pavka/Documents/R_ejercicio/montana_tesis/chapter_1")


## loading file with geographic coordinates
## Data with geographic coordinates systems in decimal degrees
## Datum: WGS1984 

coord <- read.csv ("coordenadas.csv")
coord_s <- read.csv ("coordenadas_solo.csv", header = T)
## defining spatial points
sites <- SpatialPoints(coord_s)

## Calculate straight line distance between points

strtDists <- pointDistance (sites, lonlat=TRUE)
strtDists <- as.dist(strtDists)

strtDists
### Exporting distance matrix as csv file
#write.csv(strtDists, file= "straight_dist.csv") ## export as csv for further results
#strtDists <- read.csv ("straight_dist.csv")

### Distances over a river network

## loading coordinate in utm
site_utm <- read.csv ("sitios_utm.csv")
sites_bucq <- site_utm[c(12,13,14,15),]
site_utm <- site_utm[-c(8,12,13,14,15),]
## removing Kixpur, machacas, caoba para probar
site_utm2 <- site_utm[-c(3),]
#Use estos para poder encontrar el numero de segmento y vertice para el cleanup()
##site_utm <- read.csv ("utm_chixoy.csv")
##site_utm <- site_utm [site_utm$Site== "Chixoy",]
## reading and loading the river network shape

#salinas <- line2network(layer = "rios_salinas2")
#salinas <- cleanup (salinas)
#y
#y
#y
#50 # meters
#60 ## mouth
#1 ## vertice
## plotting the river network and points
## all have to be on a projected geographyc system
## I am using UTM WGS_1984, 15 N

#plot (salinas)
#points(site_utm$Long_utm, site_utm$Lat_utm, pch=15, col=4)

#riverdistancemat (seg= salinas , vert = site_utm, rivers = salinas)

## checking

#topologydots(rivers= salinas)

## snapping distances
#sites_riv <- xy2segvert(x=site_utm$Long_utm, y=site_utm$Lat_utm, rivers= salinas)
#sites_riv2 <- xy2segvert(x=site_utm2$Long_utm, y=site_utm2$Lat_utm, rivers= salinas)

#head(sites_riv)  # a look at the first few rows of the output

#hist(sites_riv$snapdist, main="snapping distance (m)")

## Displaying point data in river locations using riverpoints()

#zoomtoseg(seg=c(4,39), rivers= salinas)
#plot (salinas)
#points(site_utm$Long_utm, y= site_utm$Lat_utm, pch=16, col="red")
#points(site_utm$Long_utm, y= site_utm$Lat_utm, pch=16, col="red")

## which point is which site
#riverpoints (seg = 180, vert = 123, rivers= salinas, pch= 15, col = "blue")
#riverpoints (seg = 1, vert = 122, rivers= salinas, pch= 15, col = "black")
#riverpoints (seg = 4, vert = 39, rivers= salinas, pch= 15, col = "red")

#riverpoints(seg=sites_riv$seg, vert=sites_riv$vert, rivers=salinas, pch=15, col="blue")

## distances between pair of sites

#riverdistance(startseg = 1, startvert = 122, 
#             endseg = 4, endvert = 39, rivers= salinas,
#            map=T)


## Sacmoc1- Sacmoc2
#riverdistance(startseg = 265, startvert = 43, 
#             endseg = 265, endvert = 41, rivers= salinas,
#            map=T)
## Sacmoc 1 - Kixpur
#riverdistance(startseg = 265, startvert = 43, 
#             endseg = 180, endvert = 123, rivers= salinas,
#            map=F)

#dist_stream_Sal <- riverdistancemat(sites_riv$seg, sites_riv$vert, salinas)

#write.csv(dist_stream_Sal, file= "river_dist.csv") ## export as csv for further results
#dist_stream_Sal <- read.csv("river_dist.csv")
#dist_stream_Sal <- as.dist (dist_stream_Sal[lower.tri(dist_stream_Sal)])
#dist_stream_Sal

#distantcia entre arroyos biotopo

## Loading Polochic network
#bucq <- line2network(layer = "rios_polochic_bucq")
#plot (bucq)

## loading utm coordinates after snapping points for BUCQ

#sites_bucq <- read.csv("biotopo_utm.csv")
#points(sites_bucq$Long_UTM, y= sites_bucq$Lat_UTM, pch=16, col="red")

## snapping points to streams
#sites_bucq_riv <- xy2segvert(x=sites_bucq$Long_UTM, y=sites_bucq$Lat_UTM, rivers= bucq)
#riverpoints (seg = 13, vert = 2, rivers= bucq, pch= 15, col = "blue")
#riverpoints (seg = 24, vert = 2, rivers= bucq, pch= 15, col = "blue")
#riverpoints (seg = 26, vert = 2, rivers= bucq, pch= 15, col = "blue")
#riverpoints (seg = 29, vert = 2, rivers= bucq, pch= 15, col = "blue")

## calculating distance matrix along streams network between BUCQ points
#dist_stream_bucq <- riverdistancemat(sites_bucq_riv$seg, sites_bucq_riv$vert, bucq)


#write.csv(dist_stream_bucq, file= "river_dist_bucq.csv") ## export as csv for further results
#dist_stream_bucq <- read.csv("river_dist_bucq.csv")

#dist_stream_Sal<- as.dist(dist_strea_Sal[lower.tri(dist_strea_Sal)])
#dist_stream_bucq <- as.dist(dist_strea_bucq[lower.tri(dist_strea_bucq)])
#dist_stream_bucq 

### loading the unique file with old the stream distances
stream_dist <- read.csv("river_dist_whole.csv")
stream_dist <- stream_dist[,-1]
colnames(stream_dist) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
stream_dist
stream_dist <- as.dist (stream_dist)
stream_dist


## loading insect assemblages matrix
## all matrices (Euclidean distance, distance over river network, assemblages) have sites (rows) in same order 
whole <- read.csv ("whole.csv") ## whole assemblage
whole <- whole[,-1]
non_sh <- read.csv("non_shredders.csv") ## non-shredder assemblage
non_sh <- non_sh[,-1]
sh <- read.csv("shredders.csv")
sh <- sh[,-1]

## getting bray-curtis simmilarity between sites
whole_dist <- vegdist(whole, method = "bray", binary= F, upper=F)
non_sh_dist <- vegdist(non_sh, method = "bray", binary= F, upper=F)
sh_dist <- vegdist(sh, method = "bray", binary= F, upper=F)


## transforming all triangle matrices to long matrices

## Euclidean distance
euclidean_sitios <- matrixConvert(strtDists, colname = c("s1", "s2", "Dist"))

## River distance
dist_river <- matrixConvert(stream_dist, colname= c("s1", "s2", "River_Dist"))

## Environmental distances
envi_sitios <- read.csv ("environment.csv")
envi_sitios_b <- envi_sitios[,-1]
## Standarizing variables

envi_sitios_b <- envi_sitios_b %>% mutate_at (c("Elevation", "Canopy_cover",
                                                "z.cm.", "w.m.", "v.m.s.", "Q.m..s.", "T_mean",
                                                "T_min", "T_max", "DO_mean", "DO_min", "DO_max", 
                                                "EC", "pH", "NO3", "NH4", "SRP"), 
                                              ~(scale(.)%>% as.vector))

eucl_envi <- vegdist(envi_sitios_b, method = "euclidean", na.rm = T, upper=F)
eucl_envi_dist <- matrixConvert(eucl_envi, colname = c("s1", "s2", "eucl_envi"))

## assemblage
whole_dist <- matrixConvert(whole_dist, colname =c("s1", "s2", "bray_whole"))
whole_dist$bray_whole <- 1- whole_dist$bray_whole ## converint to 1 = 100% similarity
non_sh_dist <- matrixConvert(non_sh_dist, colname =c("s1", "s2", "bray_non_sh") )
non_sh_dist$bray_non_sh <- 1- non_sh_dist$bray_non_sh
sh_dist  <- matrixConvert (sh_dist, colname =c("s1", "s2", "bray_sh") )
sh_dist$bray_sh <- 1- sh_dist$bray_sh

## combining matrixes 
matrix_complete <- cbind(euclidean_sitios, dist_river, eucl_envi_dist, whole_dist, non_sh_dist, sh_dist)

### bayesian model

## First re scaling Distances from meters to kilometers
matrix_complete$Dist_km <- matrix_complete$Dist/1000

matrix_complete$River_Dist_km <- matrix_complete$River_Dist/1000

## Plotting environment ~ Distances between sites
#par(mai=c(0.7,0.7,0.1,0.1), mgp=c(2,1,0))

par(mai=c(0.7,0.7,0.1,0.1), mgp=c(2,1,0))
plot(matrix_complete$eucl_envi~ matrix_complete$River_Dist_km, xlab= "Distance between sites (km)",
     ylab= "Enviromental distance (Euclidean distance)", pch = 16, cex= 1.5,
     cex.lab = 1.5)
     
## bayesian model in rstan format

sink("enviro_dist.stan")
cat("
    data {
    
    int <lower=1> N; //number of data points
    vector [N] E; // Environmental distances between pair of sites
    vector [N] D; // Distances on Km between pair of sites
    
    }
    
    parameters {
    real a; // intercept
    real b; // slope
   	real <lower=0> sigma;
    }
    
    model { 
    
    //priors. 
    a ~ uniform (0 , 10);
    b ~ normal(-0.06 , 0.05);
    
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
par(mfrow= c(1,1))
par(mai=c(1,1,0.1,0.1), mgp=c(2,1,0))
plot(matrix_complete$eucl_envi~ matrix_complete$Dist_km, xlab= "",
     ylab= "", pch = 16, cex= 1.5, col="blue",
     cex.lab = 2, cex.axis = 1.5, xlim =c(0,133), ylim=c(0,9)
)
mtext("Enviromental dissimilarity", cex = 2, side = 2, line = 3.5)
mtext("(Euclidean distance)", cex = 1.5, side=2, line= 2)
mtext("Distance between sites (km)", cex = 2, side = 1, line = 3)
text(25,0.5, expression(italic('E~N')~"(4.7 + 0.006 "%*%""~italic('D')~",1.4)"), cex = 1.5)
#points(x, fitted_mean, type="l", col="black", lwd = 3)

#plot(matrix_complete$eucl_envi ~ matrix_complete$Dist_km, xlab= "Distance between sites (km)",
#     ylab= "Enviromental distance (Euclidean distance)", pch=16, col="blue" , xlim=c(0,500)   )

for(i in 1:1000){
  lines(x, asteps$a[i] + bsteps$b[i]*as.numeric(x), col = alpha("lightblue", 0.1) )
}
points(matrix_complete$eucl_envi ~ matrix_complete$Dist_km, pch=16, cex= 1.5)
lines(x, median(asteps$a) + median(bsteps$b)*as.numeric(x))
#text(20,2, expression(italic('E~N(4.7~+~0.006~D')), cex = 1.5)

### Bayesian model to test Assemblage similarity ~ Distance

library(brms)

## Assemblabe similarity ~ Euclidean Distance, shortest path between two points
## Bayesian zero-inflated beta regression

hist(matrix_complete$bray_sh, col = "dark grey", border = FALSE, xlim = c(0,1))
hist(matrix_complete$bray_non_sh, col = "dark grey", border = FALSE, xlim = c(0,1))
hist(matrix_complete$bray_whole, col = "dark grey", border = FALSE, xlim = c(0,1))

priors <- c(set_prior("normal(-0.16, 0.14)", class = "b")) ## values from Brown and Swam 2010
priors2 <- c(set_prior("normal(-0.00746, 0.022)", class = "b"))

Sh_dist_bray <- as.data.frame (cbind(matrix_complete$Dist_km, matrix_complete$River_Dist_km, matrix_complete$eucl_envi, matrix_complete$bray_sh))
colnames(Sh_dist_bray) <- c("Dist_km", "River_dist_km", "eucl_envi", "bray")
Sh_dist_bray$FFG <- c("sh")

no_Sh_dist_bray <- as.data.frame (cbind(matrix_complete$Dist_km, matrix_complete$River_Dist_km, matrix_complete$eucl_envi, matrix_complete$bray_non_sh))
colnames(no_Sh_dist_bray) <- c("Dist_km", "River_dist_km", "eucl_envi", "bray")
no_Sh_dist_bray$FFG <- c("no-sh")

matrix_dist_bray <- rbind(Sh_dist_bray, no_Sh_dist_bray)

## Assemblabe similarity ~ in-network Distance, shortest over watercourse path between two points


ddr_fit <- brm (bf(bray ~ Dist_km*FFG, phi~Dist_km*FFG, zi ~ 1), data = matrix_dist_bray, family = zero_inflated_beta(), chains = 4, iter = 6000, warmup = 5000,
                         cores = 4,  prior = priors2)

print(summary(ddr_fit),digits =4)
plot(conditional_effects(ddr_fit), points = TRUE)

## getting posterior draws for parameters from fit model

posterior_b <- ddr_fit %>%
  gather_draws(`b_.*`, regex = TRUE)

x <- seq(0, 133, length.out=100) ## distance vector

mean(posterior_b$.value[posterior_b$.variable == "b_FFGsh"])
quantile(posterior_b$.value[posterior_b$.variable == "b_FFGsh"], probs= c(0.025,0.975))

mean(posterior_b$.value[posterior_b$.variable == "b_Dist_km:FFGsh"])
quantile(posterior_b$.value[posterior_b$.variable == "b_Dist_km:FFGsh"], probs= c(0.025,0.975))

mu_mu_non <- mean(posterior_b$.value[posterior_b$.variable == "b_Intercept"]) + mean(posterior_b$.value[posterior_b$.variable == "b_Dist_km"]) * x

mu_mu_sh <- (mean(posterior_b$.value[posterior_b$.variable == "b_Intercept"])+ mean(posterior_b$.value[posterior_b$.variable == "b_FFGsh"])) + 
  (mean(posterior_b$.value[posterior_b$.variable == "b_Dist_km"])+mean(posterior_b$.value[posterior_b$.variable == "b_Dist_km:FFGsh"])) * x

## mean and quantiles
slope_non_sh <- as.data.frame(c(mean(posterior_b$.value[posterior_b$.variable == "b_Dist_km"]),
                      quantile(posterior_b$.value[posterior_b$.variable == "b_Dist_km"], 
                               probs= c(0.025,0.975))))
slope_non_sh <- t(slope_non_sh)  
slope_non_sh$FGG <- c("no-sh")
slope_non_sh <- as.data.frame(slope_non_sh)
colnames(slope_non_sh)<- c("mean", "lowCI","UpCI", "FFG" )

slope_sh <- as.data.frame(c((mean(posterior_b$.value[posterior_b$.variable == "b_Dist_km"])+mean(posterior_b$.value[posterior_b$.variable == "b_Dist_km:FFGsh"])),
                                quantile(posterior_b$.value[posterior_b$.variable == "b_Dist_km"]+ posterior_b$.value[posterior_b$.variable == "b_Dist_km:FFGsh"], 
                                         probs= c(0.025,0.975))))
slope_sh <- t(slope_sh)  
slope_sh$FFG <- c("sh")
slope_sh <- as.data.frame(slope_sh)
colnames(slope_sh)<- c("mean", "lowCI","UpCI", "FFG" )

slope_ddr_fit <- rbind(slope_sh, slope_non_sh)
slope_ddr_fit$Index <- c(1.5, 1)
plot(slope_ddr_fit$Index, slope_ddr_fit$mean)

plot(slope_ddr_fit$Index, slope_ddr_fit$mean, ylim = rev(c(0, -0.01)), xlim=c(0.75,1.75), xaxt= "n", 
     pch=16, cex=1.5, ylab= "Estimated slope", xlab="", cex.lab= 1.5)
mtext("Non-Shredder",cex = 1.5, side = 1, at=1)
mtext("Shredder",cex = 1.5, side = 1, at=1.5)
arrows(slope_ddr_fit$Index, x1= (slope_ddr_fit$Index), y0= (slope_ddr_fit$lowCI), y1= (slope_ddr_fit$UpCI) , code=0, angle=90, length=0.15, lwd =2, lty=1);

b1_nonsh <- as.data.frame(c(posterior_b$.value[posterior_b$.variable == "b_Dist_km"]))
colnames(b1_nonsh) <- c("Slope") 
b1_nonsh$FFG <- c("Non-shredder")

b1_sh <- as.data.frame(c(posterior_b$.value[posterior_b$.variable == "b_Dist_km"] + posterior_b$.value[posterior_b$.variable == "b_Dist_km:FFGsh"]))
colnames(b1_sh) <- c("Slope") 
b1_sh$FFG <- c("Shredder")

b1_ddr_fit <- rbind(b1_nonsh, b1_sh)    

b1_ddr_fit_plot <- ggplot2.violinplot(data=b1_ddr_fit, xName='FFG',yName="Slope", mainTitle="Slope")+
                      #        groupColors=c('black','black'))+
  geom_violin(aes(fill= FFG))+
  geom_hline(yintercept= 0, color= "red", linetype="dashed")+
  scale_fill_manual(values=c("grey", "grey"))+
  labs(y="Estimated parameter", x='')+
  scale_x_discrete(labels=c('', ''))+
  theme(legend.position = "none", axis.line = element_line(linewidth =1, colour="black"),
        #axis.text.x = element_blank(),
        axis.line.y.left = element_line(linewidth = 0.5, colour = "black"),
        axis.line.x.bottom = element_line(linewidth = 0.5, colour = "black"),
        panel.background = element_rect(fill = "white"),
        #strip.background = element_rect(colour = "white", fill = "white"),
        #strip.text = element_text(size = 18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5, size= 18)
  )

b1_ddr_fit_plot


### Assemblage ~ river dist

ddr_fit_r <- brm (bf(bray ~ River_dist_km*FFG, phi~ River_dist_km*FFG, zi ~ 1), data = matrix_dist_bray, family = zero_inflated_beta(), chains = 4, iter = 6000, warmup = 5000,
                cores = 4,  prior = priors2)

print(summary(ddr_fit_r),digits =5)
plot(conditional_effects(ddr_fit_r), points = TRUE)

## getting posterior draws for parameters from fit model

posterior_b_r <- ddr_fit_r %>%
  gather_draws(`b_.*`, regex = TRUE)

x_r <-seq(0:415) ## distance vector over watercourse

mean(posterior_b_r$.value[posterior_b_r$.variable == "b_FFGsh"])
quantile(posterior_b_r$.value[posterior_b_r$.variable == "b_FFGsh"], probs= c(0.025,0.975))

mean(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"])
quantile(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"], probs= c(0.025,0.975))

mu_mu_non_r <- mean(posterior_b_r$.value[posterior_b_r$.variable == "b_Intercept"]) + mean(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"]) * x_r

mu_mu_sh_r <- (mean(posterior_b_r$.value[posterior_b$.variable == "b_Intercept"])+ mean(posterior_b_r$.value[posterior_b_r$.variable == "b_FFGsh"])) + 
  (mean(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"])+mean(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"])) * x_r

## mean and quantiles
slope_non_sh_r <- as.data.frame(c(mean(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"]),
                                quantile(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"], 
                                         probs= c(0.025,0.975))))
slope_non_sh_r <- t(slope_non_sh_r)  
slope_non_sh_r$FGG <- c("no-sh")
slope_non_sh_r <- as.data.frame(slope_non_sh_r)
colnames(slope_non_sh_r)<- c("mean", "lowCI","UpCI", "FFG" )

slope_sh_r <- as.data.frame(c((mean(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"])+mean(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"])),
                            quantile(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"]+ posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"], 
                                     probs= c(0.025,0.975))))
slope_sh_r <- t(slope_sh_r)  
slope_sh_r$FFG <- c("sh")
slope_sh_r <- as.data.frame(slope_sh_r)
colnames(slope_sh_r)<- c("mean", "lowCI","UpCI", "FFG" )

slope_ddr_fit_r <- rbind(slope_sh_r, slope_non_sh_r)
slope_ddr_fit_r$Index <- c(1.5, 1)

b1_nonsh_r <- as.data.frame(c(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"]))
colnames(b1_nonsh_r) <- c("Slope") 
b1_nonsh_r$FFG <- c("Non-shredder")

b1_sh_r <- as.data.frame(c(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"] + 
                             posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"]))
colnames(b1_sh_r) <- c("Slope") 
b1_sh_r$FFG <- c("Shredder")

b1_ddr_fit_r <- rbind(b1_nonsh_r, b1_sh_r)    

b1_ddr_fit_r_plot <- ggplot2.violinplot(data=b1_ddr_fit_r, xName='FFG',yName="Slope", mainTitle="Slope")+
  #        groupColors=c('black','black'))+
  geom_violin(aes(fill= FFG))+
  geom_hline(yintercept= 0, color= "red", linetype="dashed")+
  scale_fill_manual(values=c("grey", "grey"))+
  labs(y="Estimated parameter", x='')+
  scale_x_discrete(labels=c('', ''))+
  theme(legend.position = "none", axis.line = element_line(linewidth =1, colour="black"),
        #axis.text.x = element_blank(),
        axis.line.y.left = element_line(linewidth = 0.5, colour = "black"),
        axis.line.x.bottom =  element_line(linewidth = 0.5, colour = "black"),
        panel.background = element_rect(fill = "white"),
        #strip.background = element_rect(colour = "white", fill = "white"),
        #strip.text = element_text(size = 18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5, size= 18)
  )

b1_ddr_fit_r_plot

### Assemblage ~ environment

priors_env <- c(set_prior("normal(-0.06, 0.05)", class = "b"))

ddr_fit_e <- brm (bf(bray ~ eucl_envi*FFG, phi~ eucl_envi*FFG, zi ~ 1), data = matrix_dist_bray, family = zero_inflated_beta(), chains = 4, iter = 6000, warmup = 5000,
                  cores = 4,  prior = priors_env)


print(summary(ddr_fit_e),digits =4)
plot(conditional_effects(ddr_fit_e), points = TRUE)

## getting posterior draws for parameters from fit model

posterior_b_e <- ddr_fit_e %>%
  gather_draws(`b_.*`, regex = TRUE)


x_e <-seq(1,9, by = 0.01) ## distance vector over watercourse

#mean(posterior_b_r$.value[posterior_b_r$.variable == "b_FFGsh"])
#quantile(posterior_b_r$.value[posterior_b_r$.variable == "b_FFGsh"], probs= c(0.025,0.975))

#mean(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"])
#quantile(posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"], probs= c(0.025,0.975))

mu_mu_non_e <- mean(posterior_b_e$.value[posterior_b_e$.variable == "b_Intercept"]) + mean(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"]) * x_e

mu_mu_sh_e <- (mean(posterior_b_e$.value[posterior_b_e$.variable == "b_Intercept"])+ mean(posterior_b_e$.value[posterior_b_e$.variable == "b_FFGsh"])) + 
  (mean(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"])+mean(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi:FFGsh"])) * x_e

## mean and quantiles
slope_non_sh_e <- as.data.frame(c(mean(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"]),
                                  quantile(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"], 
                                           probs= c(0.025,0.975))))
slope_non_sh_e <- t(slope_non_sh_e)  
slope_non_sh_e$FGG <- c("no-sh")
slope_non_sh_e <- as.data.frame(slope_non_sh_e)
colnames(slope_non_sh_e)<- c("mean", "lowCI","UpCI", "FFG" )

slope_sh_e <- as.data.frame(c((mean(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"])+mean(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi:FFGsh"])),
                              quantile(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"]+ posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi:FFGsh"], 
                                       probs= c(0.025,0.975))))
slope_sh_e <- t(slope_sh_e)  
slope_sh_e$FFG <- c("sh")
slope_sh_e <- as.data.frame(slope_sh_e)
colnames(slope_sh_e)<- c("mean", "lowCI","UpCI", "FFG" )

slope_ddr_fit_e <- rbind(slope_sh_e, slope_non_sh_e)
slope_ddr_fit_e$Index <- c(1.5, 1)

b1_nonsh_e <- as.data.frame(c(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"]))
colnames(b1_nonsh_e) <- c("Slope") 
b1_nonsh_e$FFG <- c("Non-shredder")

b1_sh_e <- as.data.frame(c(posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"] + 
                             posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi:FFGsh"]))
colnames(b1_sh_e) <- c("Slope") 
b1_sh_e$FFG <- c("Shredder")

b1_ddr_fit_e <- rbind(b1_nonsh_e, b1_sh_e)    

b1_ddr_fit_e_plot <- ggplot2.violinplot(data=b1_ddr_fit_e, xName='FFG',yName="Slope", mainTitle="Slope")+
  #        groupColors=c('black','black'))+
  geom_violin(aes(fill= FFG))+
  geom_hline(yintercept= 0, color= "red", linetype="dashed")+
  scale_fill_manual(values=c("grey", "grey"))+
  labs(y="Estimated parameter", x="")+
  #scale_x_discrete(labels=c('Non-shredder', 'Shredder'))+
  theme(legend.position = "none", axis.line = element_line(linewidth =1, colour="black"),
        #axis.text.x = element_blank(),
        axis.line.y.left = element_line(linewidth = 0.5, colour = "black"),
        axis.line.x.bottom = element_line(linewidth = 0.5, colour = "black"),
        panel.background = element_rect(fill = "white"),
        #strip.background = element_rect(colour = "white", fill = "white"),
        #strip.text = element_text(size = 18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5, size= 18)
  )

b1_ddr_fit_e_plot


### assenblages dist ggplot

matrix_complete$bray_non_sh~ matrix_complete$Dist_km

matrix_complete2 <- matrix_complete[,c(-4,-5,-7,-8,-10,-11,-13,-14,-16,-17)]
mu_mu_non_p <- plogis(mu_mu_non)

ggplot(matrix_complete, aes(x = Dist_km, y= bray_non_sh))

ggplot(aes(x = x, y= mu_mu_non_p))
nonsh_ggplot <- ggplot2.scatterplot(data = matrix_complete2, xName='Dist_km', yName='bray_non_sh', mainTitle="Non-shredder assemblage")+
  ylim(0, 0.8)+
  xlim(0,133)+
  #geom_abline(intercept = (posterior_b_1000$.value[posterior_b_1000$.variable== "b_Intercept"] ), 
  #            slope = (posterior_b_1000$.value[posterior_b_1000$.variable=="b_initial_size"]), color="lightblue"
   #           , alpha = 0.05, linetype="solid", size=1.5)+
  #geom_abline(intercept = (posterior_b_1000$.value[posterior_b_1000$.variable== "b_Intercept"] + posterior_b_1000$.value[posterior_b_1000$.variable== "b_TreatmentDialiumguianense"] ), 
  #            slope = (posterior_b_1000$.value[posterior_b_1000$.variable=="b_initial_size"] + posterior_b_1000$.value[posterior_b_1000$.variable== "b_initial_size:TreatmentDialiumguianense"]), color="indianred1", 
  #            linetype="solid", size=1.5, alpha = 0.05)+  

  #mu_mu_non <- mean(posterior_b$.value[posterior_b$.variable == "b_Intercept"]) + mean(posterior_b$.value[posterior_b$.variable == "b_Dist_km"]) * x
  
  geom_line(aes(x = x, y= mu_mu_non_p))+
  #  geom_line(aes(x, mu_mu_non_p), color="black", linetype="solid", linewidth =1.5)+
  #geom_abline(intercept = ((mean(b_i$b_Intercept) + mean(b_D$b_TreatmentDialiumguianense))), 
  #            slope = (mean(b_I$b_initial_size + posterior_b$.value[posterior_b$.variable== "b_initial_size:TreatmentDialiumguianense"])), color="red", 
  #            linetype="solid", size=1.5)+  
  geom_point(aes(size=3))+
  labs(y= "", x='Distance between sites (Km)')+
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

nonsh_ggplot

y.grob <- textGrob("Assamblage similarity (Bray-Curtis)", 
                   gp=gpar(fontface="bold", col="black", fontsize=20), rot=90)

grid.arrange(arrangeGrob(plot, left = y.grob))

### Plotting

#plot.new();
#par(mfrow= c(1,1), mar =c(5,6,2,1.5))
par(mfrow= c(3,3), mar =c(5.1,6,2,1.5))
plot(matrix_complete$bray_non_sh~ matrix_complete$Dist_km, pch=16, cex= 1.5, cex.axis=2,
     ylim=c(0,0.8), xlim = c(0,133), cex.lab= 1.5, xlab= "Distance between sites (km)"
     #,ylab= "Assemblage similarity (Bray-Curtis)"
     ,ylab=""
)
mtext("Non-shredder assemblage", cex = 1.5, side = 3)
#text ( 65, 0.8, "Non-Shredder assemblage", cex = 1.5 )
for(i in 1:1000){
 mu<- posterior_b$.value[posterior_b$.variable == "b_Intercept"][i] + posterior_b$.value[posterior_b$.variable == "b_Dist_km"][i]*x
  lines(x, plogis(mu), col = alpha("lightblue", 0.1) )
}
points(matrix_complete$bray_non_sh~ matrix_complete$Dist_km, pch=16, cex= 1.5)
lines(x,plogis(mu_mu_non))

#plot.new();
plot(x,plogis(mu_mu_sh), type="l", ylim=c(0,0.8), xlim = c(0,133), cex.axis=2,
     cex.lab= 1.5, xlab= "Distance between sites (km)"
     #,ylab= "Assemblage similarity (Bray-Curtis)"
     ,ylab="")

mtext("Shredder assemblage", cex = 1.5, side = 3)
#text ( 65, 0.8, "Shredder assemblage", cex = 1.5 )
for(i in 1:1000){
  mu<- (posterior_b$.value[posterior_b$.variable == "b_Intercept"][i]+ posterior_b$.value[posterior_b$.variable == "b_FFGsh"][i]) + 
    (posterior_b$.value[posterior_b$.variable == "b_Dist_km"][i]+ posterior_b$.value[posterior_b$.variable == "b_Dist_km:FFGsh"][i]) *x
  lines(x, plogis(mu), col = alpha("lightblue", 0.1) )
}
points(matrix_complete$bray_sh~ matrix_complete$Dist_km, pch=16, cex= 1.5)
lines(x,plogis(mu_mu_sh))

plot(slope_ddr_fit$Index, slope_ddr_fit$mean, ylim = rev(c(0.001, -0.01)), xlim=c(0.75,1.75), 
     xaxt= "n", cex.axis=2,
     pch=16, cex=2, ylab= "Estimated parameter", xlab="", cex.lab= 1.5)
mtext("Slope", cex = 1.5, side = 3)
#mtext("Non-Shredder",cex = 1, side = 1, at=1, line=3)
#mtext("Shredder",cex = 1, side = 1, at=1.5, line = 3)
arrows(slope_ddr_fit$Index, x1= (slope_ddr_fit$Index), y0= (slope_ddr_fit$lowCI), y1= (slope_ddr_fit$UpCI) , code=0, angle=90, length=0.15, lwd =2, lty=1);
abline(h= 0, lty=2, lwd= 2, col="red")

## River distances
#plot.new();
plot(matrix_complete$bray_non_sh~ matrix_complete$River_Dist_km, pch=16, cex= 1.5, cex.axis=2,
     ylim=c(0,0.8), xlim = c(0,420), cex.lab= 1.5, xlab= "Distance in-network between sites (km)"
     ,ylab= "Assemblage similarity (Bray-Curtis)"
     #,ylab= ""
     
)
#mtext("Slope", cex = 1.5, side = 3)
#mtext("Assemblage similarity (Bray-Curtis)", cex = 2, side = 2, line= 3)
for(i in 1:1000){
  mu<-posterior_b_r$.value[posterior_b_r$.variable == "b_Intercept"][i] + 
  (posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"][i]) *x_r
  lines(x_r, plogis(mu), col = alpha("lightblue", 0.1) )
}
points(matrix_complete$bray_non_sh~ matrix_complete$River_Dist_km, pch=16, cex= 1.5)
lines(x_r,plogis(mu_mu_non_r))


#plot.new();
plot(x_r,plogis(mu_mu_sh_r), type="l", ylim=c(0,0.8), xlim = c(0,420), cex.lab= 1.5, cex.axis=2,
     xlab= "Distance in-network between sites (km)"
     #,ylab= "Assemblage similarity (Bray-Curtis)"
     , ylab=""
)
#text ( 65, 1, "Shredder assemblage", cex = 1.5 )
for(i in 1:1000){
 mu <- (posterior_b_r$.value[posterior_b_r$.variable == "b_Intercept"][i]+ posterior_b_r$.value[posterior_b_r$.variable == "b_FFGsh"][i]) + 
  (posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km"][i]+ posterior_b_r$.value[posterior_b_r$.variable == "b_River_dist_km:FFGsh"][i]) *x_r
  lines(x_r, plogis(mu), col = alpha("lightblue", 0.1) )
}
  
points(matrix_complete$bray_sh~ matrix_complete$River_Dist_km, pch=16, cex= 1.5)
lines(x_r,plogis(mu_mu_sh_r))


plot(slope_ddr_fit_r$Index, slope_ddr_fit_r$mean, ylim = rev(c(0.001, -0.01)), cex.axis=2,
     xlim=c(0.75,1.75), xaxt= "n", 
     pch=16, cex=2, ylab= "Estimated parameter", xlab="", cex.lab= 1.5)
#points(slope_ddr_fit_r$Index, slope_ddr_fit_r$mean, cex=2, pch=16)
#mtext("Slope", cex = 1.5, side = 3)
mtext("Non-Shredder",cex = 1, side = 1, at=1, line=3)
mtext("Shredder",cex = 1, side = 1, at=1.5, line = 3)
arrows(slope_ddr_fit_r$Index, x1= (slope_ddr_fit_r$Index), y0= (slope_ddr_fit_r$lowCI), y1= (slope_ddr_fit_r$UpCI) , code=0, angle=90, length=0.15, lwd =2, lty=1);
abline(h= 0, lty=2, lwd= 2, col="red")


## env distances
#plot.new();
plot(matrix_complete$bray_non_sh~ matrix_complete$eucl_envi, pch=16, cex= 1.5, cex.axis=2, 
     ylim=c(0,0.8), xlim = c(1.7, 8.2), cex.lab= 2, 
     xlab= "Environmental dissimilarity",
     #ylab= "Assemblage similarity (Bray-Curtis)"
     ,ylab=""
)
mtext("(Euclidean distance)", side = 1, cex= 1, line= 4.2)
#text ( 65, 1, "Shredder assemblage", cex = 1.5 )
for(i in 1:1000){
  mu<- posterior_b_e$.value[posterior_b_e$.variable == "b_Intercept"][i] + (posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"][i]) * x_e
      lines(x_e, plogis(mu), col = alpha("lightblue", 0.1) )
}
points(matrix_complete$bray_non_sh~ matrix_complete$eucl_envi, pch=16, cex= 1.5)
lines(x_e ,plogis(mu_mu_non_e))

#plot.new();
plot(matrix_complete$bray_sh~ matrix_complete$eucl_envi, pch=16, cex= 1.5, ylim=c(0,0.8), xlim = c(1.7, 8.2), 
     cex.lab= 2, cex.axis=2, 
     xlab= "Environmental dissimilarity"
     #,ylab= "Assemblage similarity (Bray-Curtis)"
     , ylab=""
)
mtext("(Euclidean distance)", side = 1, cex= 1, line= 4.2)
#text ( 65, 1, "Shredder assemblage", cex = 1.5 )
for(i in 1:1000){
  mu<- (posterior_b_e$.value[posterior_b_e$.variable == "b_Intercept"][i]+ posterior_b_e$.value[posterior_b_e$.variable == "b_FFGsh"][i]) + 
  (posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi"][i]+ posterior_b_e$.value[posterior_b_e$.variable == "b_eucl_envi:FFGsh"][i])*x_e
  lines(x_e, plogis(mu), col = alpha("lightblue", 0.1) )
}
points(matrix_complete$bray_sh~ matrix_complete$eucl_envi, pch=16, cex= 1.5)
lines(x_e,plogis(mu_mu_sh_e))

plot(slope_ddr_fit_e$Index, slope_ddr_fit_e$mean, ylim = rev(c(0.01, -0.14)), 
     xlim=c(0.75,1.75), xaxt= "n", cex.axis=2,
     pch=16, cex=2, ylab= "Estimated parameter", xlab="", cex.lab= 1.5)
#mtext("Slope", cex = 1.5, side = 3)
mtext("Non-Shredder",cex = 1, side = 1, at=1, line=3)
mtext("Shredder",cex = 1, side = 1, at=1.5, line = 3)
arrows(slope_ddr_fit_e$Index, x1= (slope_ddr_fit_e$Index), y0= (slope_ddr_fit_e$lowCI), y1= (slope_ddr_fit_e$UpCI) , code=0, angle=90, length=0.15, lwd =2, lty=1)
abline(h= 0, lty=2, lwd= 2, col="red")


### plotting ggplot violins

## Slopes violins
## spatial distances



