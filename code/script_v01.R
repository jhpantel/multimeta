#################################################################################################
##  R script for analyzing results given in Pantel et al. 202x, Ecological Monographs          ##
##  Metapopulation dynamics of multiple species in a heterogeneous landscape                   ##
##                                                                                             ##
#################################################################################################

# libraries: lme4, rjags, data.table, tidyr, gridExtra, ggplot2, ggrepel, mcmcplots, adegenet, car, plotly, RColorBrewer

##### Step 1.  Input directory with results -----------------------------------
dir <- "/Users/jhpantel/OneDrive/Documents/CNRS/pantel_work/manuscript/1_metapopulation_manuscript/Accepted_Version"

##### Step 2.  Assemble globally used variables -------------------------------
spec_correct <- c("G. radiata","D. surinamense","D. depressissimum","D. cimex","D. aeruginosum","A. marmorata","Ph. acuta","Ps. columella","G. cubensis","Pl. guadeloupensis","B. kuhniana","B. glabrata","B. schrammi","Po. glauca","Py. parvulus","Ma. cornuarietis","Me. tuberculata PAP", "Me. tuberculata GOS","Me. tuberculata FAL","Me. tuberculata MAD","Me. tuberculata CPF", "Me. tuberculata SEN", "T. granifera", "N. virginea", "E. viridans", "I. exustus", "H. duryi" )

logit_to_prob <- function(x){
  ret <- exp(x) / (1+exp(x))
  return(ret)
}

###Data 1. Assemble the probability a site was wet over course of data collection, estimated from a GLM--------
##var name: Pr_W = a number that gives the global probability a site was wet. It is estimated as the intercept from a GLM with site and year as random effects.
##Step 1: Read in dry data - will use dry_reduced_modified.txt. This is derived from dat.array[,,"etat"], which has the raw record of observed dry and wet values, but where I also filled in a few embedded NA values on the advice of Patrice and Philippe. This is the 'final' wet/dry data before the filling in of embedded 0 values.
dry <- matrix(scan(paste(dir,"/data/dry_reduced_modified.txt",sep=""),n = 285*15,skip=1), 285, 15, byrow = TRUE)
dry[dry == 0] <- NA   #Replace all 0 with NA
dry[dry > 1] <- 0   #Replace all 2/dry with 0

miss <- c(19,24,58,72,235)  #Sites 19, 24, 58, 73, 246
dry <- dry[-miss,]  #Exclude sites with too few data points
also_miss <- c(28,129)
dry <- dry[-also_miss,]  #Remove additional 2 sites

##Step 2: Arrange data for glm to estimate site and year effects
island <- read.table(paste(dir,"/data/island_site_id.txt",sep=""),header=T)
island <- island[-miss,]
island <- island[-also_miss,]

#Convert site-by-year matrix to list of observations
dry_list <- as.vector(t(dry))
dry_list <- data.frame(dry_list)
dry_list$site <-  rep(island$site, each=15)
dry_list$year <- rep(1:15,278)

##Step 3: Run glm with site and year as random effects of P(wet)
library(lme4)

mod1 <- glmer(dry_list ~ (1|site) + (1|year), data=dry_list, family=binomial())

re_site <- ranef(mod1, condVar=T, whichel="site")
re_year <- ranef(mod1, condVar=T, whichel="year")

glm.fitted <- fitted(mod1)

dry_list$fit<- NA
dry_list$fit[is.na(dry_list$dry_list)==F] <- glm.fitted

dry_list$re_site <- rep(re_site$site[1:278,], each=15)
dry_list$re_year <- rep(re_year$year[,], times=278)

dry_list$estimate <- dry_list$re_site + dry_list$re_year + summary(mod1)$coefficients[1]

dry_list$prob <- exp(dry_list$estimate) / (1 + exp(dry_list$estimate))

Pr_W <- mean(dry_list$prob)
wet_prob <- matrix(dry_list$prob, nrow=278, ncol=15, byrow=T)
wet_prob_year <- colMeans(wet_prob)

#just fyi
obs_Pr_W <- (sum(dry[dry!=2],na.rm=T)) / (sum(is.na(dry)==F))

#########################
##### NO COVARIATES #####
#########################
###Data 2. Assemble results of Bayesian model output - no covariates.  ----
##var name: full = a 27x12x3 array = number of species x number of parameters estimated in the Bayesian model x mean, 2.5%, 97.5%
#Pre-allocate a variable to hold results
spec <- c("1_Gundlachia_radiata", "2_Drepanotrema_lucidum", "3_Drepanotrema_depressissimum", "4_Drepanotrema_cimex", "5_Drepanotrema_aeruginosum", "6_Aplexa_marmorata", "7_Physa_acuta", "8_Lymnaea_columella", "9_Lymnaea_cubensis", "10_Plesiophysa_guadeloupensis", "11_Biomphalaria_straminea", "12_Biomphalaria_glabrata", "13_Biomphalaria_schrammi", "14_Pomacea_glauca", "15_Pyrgophorus_coronatus", "16_Marisa_cornuarietis", "17_Melanoides_tuberculata_PAP", "18_Melanoides_tuberculata_GOS", "19_Melanoides_tuberculata_FAL", "20_Melanoides_tuberculata_MAD", "21_Melanoides_tuberculata_CPF", "22_Melanoides_tuberculata_SEN", "23_Tarebia_granifera", "24_Neritina_virginea", "25_Eupera_viridans", "26_Indoplanorbis_exustus", "27_Helisoma_duryi")
num_spec <- length(spec)
full <- array(NA, dim=c(num_spec,12,3), dimnames=list(spec, c("phiW", "levins_EW", "cW", "phiD", "levins_ED", "cD", "pW","Ew_cW","EW_ED_cW","log_Ew_cW","log_EW_ED_cW","exp_cW"), c("Mean","2_5","97_5")))
library(rjags)
for(z in 1:num_spec){
  #Load results file
  load(paste(dir,"/results/no_covar/", spec[z], "/jsample.RData", sep=""))

  #Eliminate z estimates
  jsample.full[[1]] <- jsample.full[[1]][,-(6:4175)]
  jsample.full[[2]] <- jsample.full[[2]][,-(6:4175)]
  jsample.full[[3]] <- jsample.full[[3]][,-(6:4175)]  
  
  #Record mean and 2.5% - 97.5%
  full[z,c(1,3,4,6:7),1] <- summary(jsample.full)$stat[1:5,1]
  full[z,c(1,3,4,6:7),2] <- summary(jsample.full)$quantile[1:5,1]
  full[z,c(1,3,4,6:7),3] <- summary(jsample.full)$quantile[1:5,5]
  
  df=as.data.frame(rbind(jsample.full[[1]][,1:3],jsample.full[[2]][,1:3],jsample.full[[3]][,1:3]))
  
  df$levins_EW <- -log(df$`beta[1]`)
  df$levins_ED <- -log(df$`beta[3]`)
  
  full[z,2,1] <- mean(df$levins_EW)
  full[z,2,2] <- quantile(df$levins_EW,prob=0.025)
  full[z,2,3] <- quantile(df$levins_EW,prob=0.975)
  
  full[z,5,1] <- mean(df$levins_ED)
  full[z,5,2] <- quantile(df$levins_ED,prob=0.025)
  full[z,5,3] <- quantile(df$levins_ED,prob=0.975)
  
  df$EW_cW <- df$levins_EW / df$`beta[2]`
  df$EW_ED_cW <- ((df$levins_EW)*Pr_W + (df$levins_ED)*(1-Pr_W))/(df$`beta[2]`*Pr_W)
  
  full[z,8,1] <- mean(df$EW_cW)
  full[z,9,1] <- mean(df$EW_ED_cW)
  
  full[z,8,2:3] <- quantile(df$EW_cW, probs=c(.025,.975))
  full[z,9,2:3] <- quantile(df$EW_ED_cW, probs=c(.025,.975))
  
  full[z,10,1] <- mean(log(df$EW_cW))
  full[z,11,1] <- mean(log(df$EW_ED_cW))
  
  full[z,10,2:3] <- quantile(log(df$EW_cW), probs=c(.025,.975))
  full[z,11,2:3] <- quantile(log(df$EW_ED_cW), probs=c(.025,.975))
}

##### Step 3. Format results for manuscript -------------------------------
### Table 1, no covar. Formatted table of Bayesian posterior estimates. ---------------
##Number of occurrences of species in dataset
##var name: detect = a 27x3 array = number of species x (1) the number of appearances / detections of a species in the dataset, (2) the total number of sites a species ever appears in, and (3) the number of sites a species ever appears in that were observed to be dry at least once
##var name: ordered = an list of the number of appearances / detections of each species in the dataset, in decreasing order. All results will be sorted by this.

detect <- array(NA, dim=c(num_spec,3), dimnames=list(spec, c("detect","num_site_detect","num_site_detect_dry")))
##How often was each species detected?
for (i in 1:length(spec)){
  #load species presence info
  # 0 = non visited site or missing data
  # 1 = species not detected (i.e., either dry or wet site)
  # 2 = species detected (i.e., necessarily a wet site)
  S <- matrix(scan(paste(dir,"/results/no_covar/",spec[i],"/state.txt", sep=""),n = 278*15,skip=1),278,15,byrow = TRUE)
  
  #How often was each species detected?
  detect[i,1] <- sum(sum(S==2))
  
  #Sites where the species was detected at least once
  temp <- apply(S,1,function(x) any(x==2))
  detect[i,2] <- sum(temp==T)
  
  #Sites where the species was detected at least once and the site was dry at least once
  temp2 <- apply(dry,1,function(x) any(x==2))
  detect[i,3] <- sum(temp==T & temp2==T)
}

ordered <- order(detect[,1],decreasing=T)

##Parameter estimates
#gather data
tab2 <- full[ordered,c(1,4,3,7),1]
tab2_25 <- full[ordered,c(1,4,3,7),2]
tab2_975 <- full[ordered,c(1,4,3,7),3]

#format
tab2_format <- paste(round(tab2,2), " (", round(tab2_25,2), ", ", round(tab2_975,2), ")", sep="")
tab2_final <- data.frame(phiW=tab2_format[1:27], phiD=tab2_format[28:54], cW=tab2_format[55:81], pW=tab2_format[82:108], row.names=spec_correct[ordered])
write.csv(tab2_final,paste(dir,"/raw_output/Table_1/Table1_no_cov.csv",sep=""))

### Results misc. Information about parameter estimates given in the results ---------------
# detection prob: mean and sd with >= 150 appearances in dataset
mean(tab2[detect[ordered,1] >= 150,4])
sd(tab2[detect[ordered,1] >= 150,4])

# width of 95% CI with >= 150 apprearances in dataset
mean(tab2_975[,4] - tab2_25[,4])
sd(tab2_975[,4] - tab2_25[,4])
mean(tab2_975[detect[ordered,1] >= 150,4] - tab2_25[detect[ordered,1] >= 150,4])
sd(tab2_975[detect[ordered,1] >= 150,4] - tab2_25[detect[ordered,1] >= 150,4])

# persistence
mean(tab2[,1])
sd(tab2[,1])
mean(tab2[detect[ordered,1] >= 150,1])
sd(tab2[detect[ordered,1] >= 150,1])

# colonization
mean(tab2[,3])
sd(tab2[,3])
mean(tab2[detect[ordered,1] >= 150,3])
sd(tab2[detect[ordered,1] >= 150,3])

# persistence dry
mean(tab2[,2])
sd(tab2[,2])
mean(tab2[detect[ordered,1] >= 150,2])
sd(tab2[detect[ordered,1] >= 150,2])

# persistence dry: width of 95% CI with >= 150 apprearances in dataset
mean(tab2_975[,2] - tab2_25[,2])
sd(tab2_975[,2] - tab2_25[,2])
mean(tab2_975[detect[ordered,1] >= 150,2] - tab2_25[detect[ordered,1] >= 150,2])
sd(tab2_975[detect[ordered,1] >= 150,2] - tab2_25[detect[ordered,1] >= 150,2])

### Figure 1. Figure S2-S1. Observed and model-estimated occupancy over time. -------------
## function to take the column means of a matrix, which means the mean across sites for each year
mat_mean <- function(dat){
  blah <- matrix(dat, nrow=278, ncol=15, byrow=F)
  blah <- as.data.frame(blah)
  blah <- data.matrix(blah)
  result <- colMeans(blah)
  return(result)
}

## function to take the column means of a matrix, which means the mean across sites for each year, corrected for whether the species was observed at that site and whether the site was wet that year
mat_mean_present <- function(dat,pres,dry){
  blah <- matrix(dat, nrow=278, ncol=15, byrow=F)
  blah <- as.data.frame(blah)
  blah <- data.matrix(blah)
  
  #record of whether a site was visited that year or not
  blah2 <- matrix(NA, nrow=278, ncol=15)
  blah2[pres==0] <- 0  #site was not visited that year
  blah2[pres !=0] <- 1 #site was visited that year
  
  #modify record to only included sites that were wet that year
  blah2[dry != 1] <- 0  #all sites that are not wet have detected data removed
  
  #filter z for only information where site was visited that year and wet
  new_z <- matrix(NA, nrow=278, ncol=15)
  new_z <- as.data.frame(new_z)
  new_z <- data.matrix(new_z)
  
  new_z[blah2 == 1] <- blah[blah2 == 1]
  
  result <- apply(new_z, 2, function(x) sum(x==2,na.rm=T)/(sum(is.na(x)==F)))
  return(result)
}

##Assemble the model-estimated proportion of sites occupied each year.
## var name: pr_occupied: a 27*15*4*3 array that gives different measures of occupancy of each species (27) for each year (15). Those (4) measures are: 1. recorded pr. sites occupied each year, 2. Bayesian posterior estimated pr. sites occupied each year modified by the detectability probability, 3. Bayesian posterior estimated pr. sites occupied each year, and 4. the model-predicted occupancy pr. for sites that were actually visited that year, that were observed as wet that year, multipled by the detectability probability. And for each we give (3) the mean and the 2.5%-97.5% estimate (if taken from a posterior distribution).
pr_occupied <- array(NA, dim=c(num_spec,15,4,3), dimnames=list(spec,NULL,c("observed","model_predict_detect","model_predict","model_predict_observed_sites"),c("Mean","2_5","97_5")))

for(z in 1:num_spec){
  print(z)

  #1. observed
  ##Import detected data
  # state matrix
  # 0 = non visited site or missing data
  # 1 = species not detected (i.e., either dry or wet site)
  # 2 = species detected (i.e., necessarily a wet site)
  S <- matrix(scan(paste(dir,"/results/no_covar/",spec[z],"/state.txt", sep=""),n = 278*15,skip=1), 278, 15, byrow = TRUE)
  
  #The number of sites the species is observed in each year
  pr_occupied[z,,1,1] <- apply(S, 2, function(x) (sum(x==2)/sum(x==1 | x==2)))
  
  #using model data (2-4)
  #Load results file
  load(paste(dir,"/results/no_covar/",spec[z],"/jsample.RData", sep=""))
  
  #Get z estimates
  df=as.data.frame(rbind(jsample.full[[1]],jsample.full[[2]],jsample.full[[3]]))
  dp <- df[,5]
  df <- df[,6:4175]
  
  #Extract mean
  occupy <- apply(df, 1, mat_mean)
  
  #3. model_predict
  pr_occupied[z,,3,1] <- rowMeans(occupy)-1
  pr_occupied[z,,3,2:3] <- t(apply(occupy, 1, function(x) quantile(x, probs=c(.025,.975))))-1
  
  #2. model_predict_detect
  occupy_dp <- apply(occupy,1,function(x) (x-1)*dp)
  occupy_dp_wet <- apply(occupy_dp,1,function(x) x*wet_prob_year)
  
  pr_occupied[z,,2,1] <- rowMeans(occupy_dp_wet)
  pr_occupied[z,,2,2:3] <- t(apply(occupy_dp_wet,1,function(x) quantile(x,probs=c(.025,.975))))
  
  #4. model_predict_observed_sites
  dry <- matrix(scan(paste(dir,"/results/no_covar/",spec[z],"/dry_final.txt",sep=""), n = 278*15,skip=1), 278, 15, byrow = TRUE)
  
  #Extract mean
  occupy <- apply(df, 1, function(x)  mat_mean_present(x,S,dry))
  occupy_dp <- apply(occupy,1,function(x) x*dp)
  
  pr_occupied[z,,4,1] <- colMeans(occupy_dp)
  pr_occupied[z,,4,2:3] <- t(apply(occupy_dp, 2, function(x) quantile(x, probs=c(.025,.975))))
}

# Re-arrange data, create via ggplot
library(data.table)
blah <- as.data.table(pr_occupied)
colnames(blah) <- c("spec","year","model","val","P")
blah <- blah[blah$model != "model_predict_detect",]
library(tidyr)
library(gridExtra)
library(ggplot2)
blah.wide <- pivot_wider(blah, names_from = val, values_from = P)
colnames(blah.wide)[c(4,5)] <- c("v25","v975")
sub <- blah.wide[blah.wide$spec == "10_Plesiophysa_guadeloupensis",]
pd <- position_dodge(1) # move them .05 to the left and right
group_col <- c("model_predict" = "gray", "model_predict_observed_sites" = "skyblue", "observed" ="orange")

# For one species
ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(),axis.title=element_blank(), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)),text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"),axis.line = element_line(colour = "black"),panel.border = element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim)

## Appendix S2, Figure S1
# Arranged in a grid
plot_list = list()
for(z in 1:27){
  if(z == 6|z==3|z==11|z==2|z==7|z==8|z==17|z==5){
    ylim <- c(0,1)
    nb <- 5
  } else if(z == 16|z==14|z==25|z==1){
    ylim <- c(0,0.5)
    nb <- 5
  } else if(z == 19|z==20|z==10|z==24|z == 21|z==26|z==22|z==27|z==23){
    ylim <- c(0,0.1)
    nb <- 3
  }
  else{
    ylim <- c(0,0.25)
    nb <- 5
  }
  
  sub <- blah.wide[blah.wide$spec == spec[z],]
  
  if (z==24|z==21|z==26|z==22|z==27){
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none", panel.grid=element_blank(), axis.title.y=element_blank(), axis.title.x=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)),text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"),axis.line = element_line(colour = "black",size=.1),panel.border = element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim)
  } else if (z==6|z==16|z==1) {
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank(), axis.text.x=element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied")
  } else if (z==10) {
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied")
  } else {
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(),axis.title=element_blank(), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)),text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"),axis.line = element_line(colour = "black",size=.1),panel.border = element_blank(),axis.text.x=element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim)
  }
  plot_list[[z]] <- p
}

ordered_plot_list <- plot_list[ordered]
g <- grid.arrange(grobs=ordered_plot_list,legend,ncol=7,nrow=4)

ggsave(file=paste(dir,"/raw_output/S2_Figure_S1.svg",sep=""), g, units="in", width=8.75, height=5.25, dpi=600)

## Figure 1
# Arranged in a grid
plot_list = list()
plot_count <- 1
for(z in c(6,7,12,16:18)){
  if(z == 6|z==3|z==11|z==2|z==7|z==8|z==17|z==5){
    ylim <- c(0,1)
    nb <- 5
  } else if(z == 16|z==14|z==25|z==1){
    ylim <- c(0,0.5)
    nb <- 5
  } else if(z == 19|z==20|z==10|z==24|z == 21|z==26|z==22|z==27|z==23){
    ylim <- c(0,0.1)
    nb <- 3
  }
  else{
    ylim <- c(0,0.25)
    nb <- 5
  }
  sub <- blah.wide[blah.wide$spec == spec[z],]
  if (z==18) {
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd, size=1.5,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.5, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied") + scale_x_continuous(breaks = round(seq(min(sub$year), max(sub$year), by = 1),1))
  } else if (z==6){
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd, size=1.5,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.5, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank(), axis.text.x=element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied") + scale_x_continuous(breaks = round(seq(min(sub$year), max(sub$year), by = 1),1))
  }
  else{
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd, size=1.5,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.5, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank(), axis.text.x=element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied") + scale_x_continuous(breaks = round(seq(min(sub$year), max(sub$year), by = 1),1))
  }
  plot_list[[plot_count]] <- p
  plot_count <- plot_count + 1
}

g <- grid.arrange(grobs=plot_list,legend,ncol=1,nrow=6)

ggsave(file=paste(dir,"/raw_output/Figure_1/Figure_1_no_cov.svg",sep=""), g, units="in", width=3, height=8, dpi=600)

### Figure 2. Plot extinction / colonization for each species. ---------------

## Figure 2
##Re-arrange data in desired order of plot
new_full <- array(NA,dim=c(1,6))
bg <- array(NA,dim=c(27,2))
bg[,1] <- "black"
bg[,2] <- "white"
for (z in ordered){
  for (i in 10:11){
    new_full <- rbind(new_full,c(full[z,i,],dimnames(full)[1][[1]][z],dimnames(full)[2][[1]][i],bg[z,i-9]))
  }
}
new_full <- new_full[2:55,]

svg(file=paste(dir,"/raw_output/Figure_2.svg", sep=""), width=4.608, height=4.858)
par(mar=c(1,5,1,1),mgp=c(3,1,0))
plot(new_full[,1],1:54,type="n",xlim=c(-2.5,2.5),yaxt="n",ylab="",xlab="",cex.axis=.5)
axis(2, at=seq(1.5,54.5,2), labels=rev(spec_correct[ordered]), las=1,cex.axis=.7,font=3)
points(rev(new_full[,1]),1:54,pch=21,col="black",bg=rev(new_full[,6]),cex=.9)

abline(v=0,lty=2)
abline(h=seq(1.5,54.5,2),col=rgb(190,190,190,200,maxColorValue=255),lty=2,lwd=.5)
abline(v=log(0.1),lty=2,col=rgb(0,0,0,100,maxColorValue=255),lwd=.5)
abline(v=log(0.25),lty=2,col=rgb(0,0,0,100,maxColorValue=255),lwd=.5)
abline(v=log(0.5),lty=2,col=rgb(0,0,0,100,maxColorValue=255),lwd=.5)
abline(v=log(0.75),lty=2,col=rgb(0,0,0,100,maxColorValue=255),lwd=.5)

arrows(as.numeric(rev(new_full[,1])), 1:54, as.numeric(rev(new_full[,2])), 1:54, length=0.025, angle=90, code=2, lwd=.3, col="black")

arrows(as.numeric(rev(new_full[,1])), 1:54, as.numeric(rev(new_full[,3])), 1:54, length=0.025, angle=90, code=2, lwd=.3, col="black")
dev.off()


### Figure 3. Observed vs. model-predicted equilibrium occupancy----
## Figure 3a
# x-axis: p* = [1-(e/c)]*Pr(detection)
# y-axis: observed detection frequency (averaged across years)
# compile Levins expectation for persistence posterior distribution
#Pre-allocate a variable to hold results
num_spec <- length(spec)
lev_detect <- array(NA, dim=c(num_spec,3), dimnames=list(spec, c("Mean","2_5","97_5")))
occ <- array(NA,dim=length(spec))
library(rjags)
library(ggrepel)

for(z in 1:num_spec){
  # Load results file
  load(paste(dir,"/results/no_covar/", spec[z], "/jsample.RData", sep=""))
  
  # Eliminate z estimates
  jsample.full[[1]] <- jsample.full[[1]][,-(6:4175)]
  jsample.full[[2]] <- jsample.full[[2]][,-(6:4175)]
  jsample.full[[3]] <- jsample.full[[3]][,-(6:4175)]  
  
  # Compile p* data, modified by detection probability
  df=as.data.frame(rbind(jsample.full[[1]],jsample.full[[2]],jsample.full[[3]]))
  df$levins_EW <- -log(df$`beta[1]`)
  df$levins_ED <- -log(df$`beta[3]`)
  df$p_star <- (1-(((df$levins_EW)*Pr_W + (df$levins_ED)*(1-Pr_W))/(df$`beta[2]`*Pr_W)))*(df$pM)
  df$p_star_01 <- df$p_star
  df$p_star_01[df$p_star < 0] <- 0
  
  lev_detect[z,1] <- mean(df$p_star)
  lev_detect[z,2:3] <- quantile(df$p_star, probs=c(.025,.975))
  
  # Compile pr_occupied_observed, which is the proportion of occupied sites each year, then take the average ONLY including years from beginning to end of that species appearances in the dataset.
  
  ##Import detected data
  # state matrix
  # 0 = non visited site or missing data
  # 1 = species not detected (i.e., either dry or wet site)
  # 2 = species detected (i.e., necessarily a wet site)
  S <- matrix(scan(paste(dir,"/results/no_covar/",spec[z],"/state.txt", sep=""),n = 278*15,skip=1), 278, 15, byrow = TRUE)
  
  #The number of sites the species is observed in each year
  num_occ <- apply(S, 2, function(x) (sum(x==2)/sum(x==1 | x==2)))
  
  ############################################################
  # compute the year a site is first visited                 #
  ############################################################
  nsite<-dim(S)[1]
  nyear<-15
  e <- NULL
  for (i in 1:nsite)
  {
    temp <- 1:nyear
    tempo <- min(temp[S[i,]==1],temp[S[i,]==2],temp[S[i,]==3])
    e <- c(e,tempo)
  }
  
  ############################################################
  # compute the year a site is last visited                  #
  ############################################################
  fin <- NULL
  for (i in 1:nsite)
  {
    temp <- 1:nyear
    tempo <- max(temp[S[i,]==1],temp[S[i,]==2],temp[S[i,]==3])
    fin <- c(fin,tempo)
  }
  
  start <- min(e)
  end <- max(fin)
  
  occ[z] <- mean(num_occ[start:end],na.rm=T)
}
lev_detect_01 <- lev_detect[,1]
lev_detect_01[lev_detect[,1] < 0] <- 0
cor(lev_detect_01,occ)
gdat <- as.data.frame(cbind(lev_detect_01,occ))
rownames(gdat) <- spec_correct
rownames(gdat)[17:22] <- substr(spec_correct[17:22],17,25)
gdat$inv <- rep("N",27)
gdat$inv[c(11,7,8,17,16,18,23,19,20,21,26,27,22)] <- "I"
group_col <- c("I" = "black", "N" = "white")
dt.triangle <- data.table(group = c(1,1,1), polygon.x = c(0,0.6,0.6), polygon.y = c(0,0,0.6))

# Scatterplot colored by native or invasive
p <- ggplot(data=gdat,aes(lev_detect_01,occ,fill=inv)) + geom_point(shape = 21,colour = "black",size=1) + scale_fill_manual(values=group_col) + theme_bw()+ theme(legend.position="none",axis.text=element_text(size=7),axis.title=element_blank()) + geom_segment(aes(x = 0, y = 0, xend = 0.6, yend = 0.6)) + geom_polygon(data=dt.triangle,aes(x=polygon.x,y=polygon.y),inherit.aes=FALSE,alpha=0.1, color='black')
#p <- ggplot(data=gdat,aes(lev_detect_01,occ,fill=inv)) + geom_point(shape = 21,colour = "black",size=1) + scale_fill_manual(values=group_col) + theme_bw()+ theme(legend.position="none",axis.title=element_text(size=7,margin=margin(0,0,0,0)),axis.text=element_text(size=7),plot.title = element_text(size=7,margin=margin(0,0,0,0))) + xlab("expected detection frequency") + ylab("observed detection frequency") + geom_segment(aes(x = 0, y = 0, xend = 0.6, yend = 0.6)) + geom_polygon(data=dt.triangle,aes(x=polygon.x,y=polygon.y),inherit.aes=FALSE,alpha=0.1, color='black') + ggtitle("b) without covariates")
# Add labels for closely-spaced values
p <- p + geom_text_repel(data=subset(gdat,lev_detect_01 < 0.05 & occ < 0.2),aes(label=rownames(gdat)[gdat$lev_detect_01 < 0.05 & gdat$occ < 0.2]),max.overlaps=Inf,direction="y",hjust=1,nudge_x=-0.05,xlim = c(-Inf, Inf),min.segment.length = 0,size=2.2,segment.size=0.2) + scale_x_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6),limits=c(-.1,.6),expand = expansion(mult = 0.05),labels=c(0,.1,.2,.3,.4,.5,.6)) + scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6),limits=c(0,0.61),labels=c(0,.1,.2,.3,.4,.5,.6))
# Add labels for well-spaced values
p <- p + geom_text_repel(data=subset(gdat,lev_detect_01 > 0.05),aes(label=rownames(gdat)[gdat$lev_detect_01 > 0.05]),max.overlaps=Inf,size=2.2)
# Save to external file
svg(file=paste(dir,"/raw_output/Figure_3/Figure_3a.svg", sep=""), width=6, height=3)
par(mar=c(0,0,0,0),mgp=c(0,0,0))
p
dev.off()

### Figure 4a. Plot extinction / colonization for each species. ---------------
# x-axis: p* = [1-(e/c)]
# y-axis: e (eD and eW, weighed by Pr(wet))
p_star_post_mean <- array(NA, num_spec)
ext <- array(NA, num_spec)

for(z in 1:num_spec){
  # Load results file
  load(paste("/Users/jhpantel/OneDrive/Documents/CNRS/pantel_work/2015_model_run/for_manuscript/1_metapopulation_manuscript/no_covar/all_sites/results/", spec[z], "/jsample.RData", sep=""))
  
  # Eliminate z estimates
  jsample.full[[1]] <- jsample.full[[1]][,-(6:4175)]
  jsample.full[[2]] <- jsample.full[[2]][,-(6:4175)]
  jsample.full[[3]] <- jsample.full[[3]][,-(6:4175)]  
  
  # Compile p* data, modified by detection probability
  df=as.data.frame(rbind(jsample.full[[1]],jsample.full[[2]],jsample.full[[3]]))
  df$levins_EW <- -log(df$`beta[1]`)
  df$levins_ED <- -log(df$`beta[3]`)
  df$p_star <- 1-(((df$levins_EW)*Pr_W + (df$levins_ED)*(1-Pr_W))/(df$`beta[2]`*Pr_W))
  df$ext <- (df$levins_EW)*Pr_W + (df$levins_ED)*(1-Pr_W)
  
  p_star_post_mean[z] <- mean(df$p_star)
  ext[z] <- mean(df$ext)
}

p_star_post_mean_01 <- p_star_post_mean
p_star_post_mean_01[p_star_post_mean < 0] <- 0

# Scatterplot
gdat <- as.data.frame(cbind(p_star_post_mean_01,ext))
rownames(gdat) <- spec_correct
rownames(gdat)[17:22] <- substr(spec_correct[17:22],17,25)

p <- ggplot(data=gdat,aes(p_star_post_mean_01,ext)) + geom_point(shape = 19,colour = "black",size=1.5) + theme_bw()+ theme(legend.position="none",axis.text=element_text(size=7),axis.title=element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
# Add labels for closely-spaced values
p <- p + geom_text_repel(data=subset(gdat,p_star_post_mean_01 < 0.05 & ext < 0.3),aes(label=rownames(gdat)[gdat$p_star_post_mean_01 < 0.05 & gdat$ext < 0.3]),max.overlaps=Inf,direction="y",hjust=1,nudge_x=-0.05,xlim = c(-Inf, Inf),min.segment.length = 0,size=3,segment.size=0.6,segment.color="#B3B3B3") + scale_x_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7),limits=c(-.16,.75),expand = expansion(mult = 0.05),labels=c(0,.1,.2,.3,.4,.5,.6,.7)) + scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6),limits=c(0,0.6),labels=c(0,.1,.2,.3,.4,.5,.6))
# Add labels for well-spaced values
p <- p + geom_text_repel(data=subset(gdat,p_star_post_mean_01 > 0.05 | ext > 0.3),aes(label=rownames(gdat)[gdat$p_star_post_mean_01 > 0.05 | gdat$ext > 0.3]),max.overlaps=Inf,size=3)
# Save to external file
svg(file=paste(dir,"/raw_output/Figure_3/Fig_4a.svg", sep=""), width=5.333, height=4)
par(mar=c(0,0,0,0),mgp=c(0,0,0))
p
dev.off()

## Simulation model input and results.----

##varname: wet, a list of all sites for which a species was ever observed.
wet <- matrix(dry_list$prob, nrow=278, ncol=15, byrow=T)
write.csv(wet,file=paste(dir,"/data/simulation/parameters/site_table.csv",sep=""),row.names=F)
write.csv(wet,file=paste(dir,"/data/simulation/no_covar/parameters/site_table.csv",sep=""),row.names=F)

##varname: psi = a 278*27 record of average occupancy probabilities
psi <- array(NA, dim=c(278,num_spec), dimnames=list(NULL,spec))
for(z in 1:num_spec){
  #Load results file
  load(paste(dir,"/results/no_covar/",spec[z],"/jsample.RData", sep=""))
  #Get z estimates
  df=as.data.frame(rbind(jsample.full[[1]],jsample.full[[2]],jsample.full[[3]]))
  df <- df[,6:4175]
  #Extract mean
  occupy <- apply(df, 1, mat_mean)
  psi[,z] <- rowMeans(occupy)
}

psi <- psi-1  #Gives the percentage of years (across all JAGS iterations) that each site is occupied
write.csv(psi,file=paste(dir,"/data/simulation/no_covar/parameters/param_psi_bysite.csv",sep=""),row.names=F)

##varname: site = a 27*278*15*3 array with the site- and year-specific cW, EW, and ED values for the model with no covariates
site <- array(NA, dim=c(num_spec,278,15,3), dimnames=list(spec,NULL,c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015"),c("cW","EW","ED")))

for(z in 1:num_spec){
  #cW
  site[z,,,1] <- array(full[z,3,1],dim=c(278,15))
  #EW
  site[z,,,2] <- array(full[z,2,1],dim=c(278,15))
  #ED
  site[z,,,3] <- array(full[z,5,1],dim=c(278,15))
}

##Write eW, cW, and eD to .csv file for Mathematica
for(z in c(1:27)){
  #Scenario = intercept (no covariates) (model V05)
  write.csv(site[z,,,1],file=paste(dir,"/data/simulation/no_covar/parameters/",spec[z],"_site_intercept_cW.csv",sep=""),row.names=F)
  write.csv(site[z,,,2],file=paste(dir,"/data/simulation/no_covar/parameters/",spec[z],"_site_intercept_eW.csv",sep=""),row.names=F)
  write.csv(site[z,,,3],file=paste(dir,"/data/simulation/no_covar/parameters/",spec[z],"_site_intercept_eD.csv",sep=""),row.names=F)
}

#########################
##### COVARIATES ########
#########################
### Data 3. Assemble results of Bayesian model output - with covariates -------------------------------

##var name: full = a 27*33*4 array = number of species * number of parameters estimated x mean, 2.5%, 97.5%, n (number of iterations that incuded this parameter, used for SVSS or beta parameters)
##var name: things = a 27*9*3 array = number of species * persistence ("ln_Ew_cW","ln_EW_ED_cW","Ew_cW","EW_ED_cW","logit_phiW","logit_phiD","exp_cW","EW", "ED") x mean, 2.5%, 97.5%

#Names of all parameters estimated in Bayesian model
param <- c("phiW", "phiD", "cW", "cD", "b_phiW_size", "b_phiW_veg", "b_phiW_stability", "b_phiW_lrs", "b_phiW_mangrove", "b_phiW_river", "b_cW_size", "b_cW_veg", "b_cW_stability", "b_cW_connec", "b_cW_rs", "b_cW_colsource", "b_cW_mangrove", "b_cW_river", "g_phiW_size", "g_phiW_veg", "g_phiW_stability", "g_phiW_lrs", "g_phiW_mangrove", "g_phiW_river", "g_cW_size", "g_cW_veg", "g_cW_stability", "g_cW_connec", "g_cW_rs", "g_cW_colsource", "g_cW_mangrove", "g_cW_river", "pW")

#Pre-allocate a variable to hold results of Bayesian model
num_spec <- length(spec)
num_param <- length(param)
full <- array(NA, dim=c(num_spec,num_param,4), dimnames=list(spec, param, c("Mean","2_5","97_5","n")))
things <- array(NA, dim=c(num_spec,9,3), dimnames=list(spec, c("ln_EW_cW","ln_EW_ED_cW","EW_cW","EW_ED_cW","logit_phiW","logit_phiD","exp_cW","EW","ED"), c("Mean","2_5","97_5")))

for(z in 1:num_spec){
  #Load results file
  load(paste(dir,"/results/covar/",spec[z],"/jsample.RData", sep=""))
  
  df=as.data.frame(rbind(jsample.full[[1]][,1:33],jsample.full[[2]][,1:33],jsample.full[[3]][,1:33]))
  
  #Mean, CI, and n when beta values when gamma = 1
  for(i in 5:18){
    sub <- df[df[,i+14] == 1,i]
    full[z,i,1] <- mean(sub)
    full[z,i,2:3] <- quantile(sub, probs=c(.025,.975))
    full[z,i,4] <- length(sub)  
  }
  
  #Mean and CI values for alpha and pW parameters
  sub <- df[,c(1:4,33)]
  sub$logit_phiW <- logit_to_prob(sub$`alpha[1]`)
  sub$logit_phiD <- logit_to_prob(sub$`alpha[2]`)
  sub$exp_cW <- exp(sub$`alpha[3]`)
  
  sub$EW <- -log(sub$logit_phiW)
  sub$ED <- -log(sub$logit_phiD)
  
  sub$ln_EW_cW <- log((sub$EW) / (sub$exp_cW))
  sub$ln_EW_ED_cW <- log(((sub$EW)*Pr_W + (sub$ED)*(1-Pr_W))/(sub$exp_cW*Pr_W))
  sub$EW_cW <- sub$EW / sub$exp_cW
  sub$EW_ED_cW <- ((sub$EW*Pr_W) + (sub$ED)*(1-Pr_W))/(sub$exp_cW*Pr_W)
  
  things[z,1,1] <- mean(sub$ln_EW_cW)
  things[z,1,2:3] <- quantile(sub$ln_EW_cW,probs=c(.025,.975))
  
  things[z,2,1] <- mean(sub$ln_EW_ED_cW)
  things[z,2,2:3] <- quantile(sub$ln_EW_ED_cW,probs=c(.025,.975))
  
  things[z,3,1] <- mean(sub$EW_cW)
  things[z,3,2:3] <- quantile(sub$EW_cW,probs=c(.025,.975))
  
  things[z,4,1] <- mean(sub$EW_ED_cW)
  things[z,4,2:3] <- quantile(sub$EW_ED_cW,probs=c(.025,.975))
  
  things[z,5,1] <- mean(sub$logit_phiW)
  things[z,5,2:3] <- quantile(sub$logit_phiW,probs=c(.025,.975))
  
  things[z,6,1] <- mean(sub$logit_phiD)
  things[z,6,2:3] <- quantile(sub$logit_phiD,probs=c(.025,.975))
  
  things[z,7,1] <- mean(sub$exp_cW)
  things[z,7,2:3] <- quantile(sub$exp_cW,probs=c(.025,.975))
  
  things[z,8,1] <- mean(sub$EW)
  things[z,8,2:3] <- quantile(sub$EW,probs=c(.025,.975))
  
  things[z,9,1] <- mean(sub$ED)
  things[z,9,2:3] <- quantile(sub$ED,probs=c(.025,.975))
  
  sub <- sub[,1:5]
  
  full[z,c(1:4,33),1] <- apply(sub,2,mean)
  full[z,c(1:4,33),2:3] <- t(apply(sub,2,function(x) quantile(x,probs=c(.025,.975))))
  
  #Extract proportion of models that included covariates via SVSS
  #Remove non-gamma parameters from jsample MCMC object
  sub <- df[,19:32]
  nbMCMC <- nrow(sub)
  full[z,19:32,1] <- apply(sub,2,function(x) sum(x==1)/nbMCMC)
}

# To quickly look at convergence
library(mcmcplots)
mcmcplot(jsample.full,"alpha")
mcmcplot(jsample.full,"beta")

### Table 1, covar. Formatted table of Bayesian posterior estimates for parameter intercepts. ---------------
#gather data
tab2 <- c(things[ordered,c(5,7,6),1],full[ordered,33,1])
tab2_25 <- c(things[ordered,c(5,7,6),2],full[ordered,33,2])
tab2_975 <- c(things[ordered,c(5,7,6),3],full[ordered,33,3])
#format
tab2_format <- paste(round(tab2,2), " (", round(tab2_25,2), ", ", round(tab2_975,2), ")", sep="")
tab2_final <- data.frame(phiW=tab2_format[1:27], cW=tab2_format[28:54], phiD=tab2_format[55:81], pW=tab2_format[82:108], row.names=spec_correct[ordered])
write.csv(tab2_final,paste(dir,"/raw_output/Table_1/Table1_cov.csv",sep=""))

# mean and sd, all and with >= 150 appearances in dataset
# phiW
mean(things[,5,1])
sd(things[,5,1])
mean(things[detect[ordered,1] >= 150,5,1])
sd(things[detect[ordered,1] >= 150,5,1])
# cW
mean(things[,7,1])
sd(things[,7,1])
mean(things[detect[ordered,1] >= 150,7,1])
sd(things[detect[ordered,1] >= 150,7,1])
# phiD
mean(things[,6,1])
sd(things[,6,1])
mean(things[detect[ordered,1] >= 150,6,1])
sd(things[detect[ordered,1] >= 150,6,1])
# detection prob
mean(full[,33,1])
sd(full[,33,1])
mean(full[detect[ordered,1] >= 150,33,1])
sd(full[detect[ordered,1] >= 150,33,1])

### Table S1-S5 and S1-S6. Formatted table of Bayesian posterior estimates for significant beta coefficients --------

##Beta estimates
#gather data
beta <- full[ordered,5:10,1]
beta_25 <- full[ordered,5:10,2]
beta_975 <- full[ordered,5:10,3]

sig <- full[ordered,19:24,1]
beta_sig <- array(NA,dim=c(27,6))
beta_sig[which(sig>=.6)] <- beta[which(sig>=.6)]

#format
beta_format <- paste(round(beta,2)," (", round(beta_25,2), ", ", round(beta_975,2), ")", sep="")

#Add an asterick if parameter is significant.
beta_format_sig <- beta_format
beta_format_sig[as.vector(sig) >= .6 & is.na(as.vector(sig))==F] <- paste(beta_format_sig[as.vector(sig) >= .6 & is.na(as.vector(sig))==F],"*",sep="")

beta_final <- data.frame(size=beta_format_sig[1:27], veg=beta_format_sig[28:54], stability=beta_format_sig[55:81], lrs=beta_format_sig[82:108], mangrove=beta_format_sig[109:135], river=beta_format_sig[136:162], row.names=spec[ordered])

write.csv(beta_final,paste(dir,"/raw_output/Table_S5/Table_S5_phiW.csv",sep=""))

## Version that gives SSVS value
beta_format <- paste(round(beta,2),"    (", round(beta_25,2), ", ", round(beta_975,2), ")    ",round(sig,2), sep="")

beta_final <- data.frame(size=beta_format[1:27], veg=beta_format[28:54], stability=beta_format[55:81], lrs=beta_format[82:108], mangrove=beta_format[109:135], river=beta_format[136:162], row.names=spec_correct[ordered])

write.csv(beta_final,paste(dir,"/raw_output/Table_S5/Table_S5_phiW_SSVS.csv",sep=""))

#gather data
beta <- full[ordered,11:18,1]
beta_25 <- full[ordered,11:18,2]
beta_975 <- full[ordered,11:18,3]

sig <- full[ordered,25:32,1]
beta_sig <- array(NA,dim=c(27,8))
beta_sig[which(sig>=.6)] <- beta[which(sig>=.6)]

#format
beta_format <- paste(round(beta,2), " (", round(beta_25,2), ", ", round(beta_975,2), ")", sep="")

#Add an asterick if parameter is significant.
beta_format_sig <- beta_format
beta_format_sig[as.vector(sig) >= .6 & is.na(as.vector(sig))==F] <- paste(beta_format_sig[as.vector(sig) >= .6 & is.na(as.vector(sig))==F],"*",sep="")

beta_final <- data.frame(size=beta_format_sig[1:27], veg=beta_format_sig[28:54], stability=beta_format_sig[55:81], connec=beta_format_sig[82:108], rs=beta_format_sig[109:135], colsource=beta_format_sig[136:162], mangrove=beta_format_sig[163:189], river=beta_format_sig[190:216], row.names=spec[ordered])

write.csv(beta_final,paste(dir,"/raw_output/Table_S6/Table_S6_cW.csv",sep=""))

## Version that gives SSVS value
beta_format <- paste(round(beta,2),"    (", round(beta_25,2), ", ", round(beta_975,2), ")    ",round(sig,2), sep="")

beta_final <- data.frame(size=beta_format[1:27], veg=beta_format[28:54], stability=beta_format[55:81], connec=beta_format[82:108], rs=beta_format[109:135], colsource=beta_format[136:162], mangrove=beta_format[163:189], river=beta_format[190:216], row.names=spec_correct[ordered])

write.csv(beta_final,paste(dir,"/raw_output/Table_S6/Table_S6_cW_SSVS.csv",sep=""))

### Figure 1. Figure S2-S2. Observed and model-estimated occupancy over --------
pr_occupied <- array(NA, dim=c(num_spec,15,4,3), dimnames=list(spec,NULL,c("observed","model_predict_detect","model_predict","model_predict_observed_sites"),c("Mean","2_5","97_5")))

for(z in 1:num_spec){
  print(z)
  #1. observed
  ##Import detected data
  # state matrix
  # 0 = non visited site or missing data
  # 1 = species not detected (i.e., either dry or wet site)
  # 2 = species detected (i.e., necessarily a wet site)
  S <- matrix(scan(paste(dir,"/results/no_covar/",spec[z],"/state.txt", sep=""),n = 278*15,skip=1), 278, 15, byrow = TRUE)
  
  #The number of sites the species is observed in each year
  pr_occupied[z,,1,1] <- apply(S, 2, function(x) (sum(x==2)/sum(x==1 | x==2)))
  
  #using model data (2-4)
  #Load results file
  load(paste(dir,"/results/covar/",spec[z],"/jsample.RData", sep=""))
  
  #Get z estimates
  df=as.data.frame(rbind(jsample.full[[1]],jsample.full[[2]],jsample.full[[3]]))
  
  dp <- df[,33]
  df <- df[,34:4203]
  
  #Extract mean
  occupy <- apply(df, 1, mat_mean)
  
  #3. model_predict
  pr_occupied[z,,3,1] <- rowMeans(occupy)-1
  pr_occupied[z,,3,2:3] <- t(apply(occupy, 1, function(x) quantile(x, probs=c(.025,.975))))-1
  
  #2. model_predict_detect
  occupy_dp <- apply(occupy,1,function(x) (x-1)*dp)
  occupy_dp_wet <- apply(occupy_dp,1,function(x) x*wet_prob_year)
  
  pr_occupied[z,,2,1] <- rowMeans(occupy_dp_wet)
  pr_occupied[z,,2,2:3] <- t(apply(occupy_dp_wet,1,function(x) quantile(x,probs=c(.025,.975))))
  
  #4. model_predict_observed_sites
  dry <- matrix(scan(paste(dir,"/results/no_covar/",spec[z],"/dry_final.txt",sep=""), n = 278*15,skip=1), 278, 15, byrow = TRUE)  
  #Extract mean
  occupy <- apply(df, 1, function(x)  mat_mean_present(x,S,dry))
  occupy_dp <- apply(occupy,1,function(x) x*dp)
  
  pr_occupied[z,,4,1] <- colMeans(occupy_dp)
  pr_occupied[z,,4,2:3] <- t(apply(occupy_dp, 2, function(x) quantile(x, probs=c(.025,.975))))
}

# Re-arrange data, create via ggplot
blah <- as.data.table(pr_occupied)
colnames(blah) <- c("spec","year","model","val","P")
blah <- blah[blah$model != "model_predict_detect",]
blah.wide <- pivot_wider(blah, names_from = val, values_from = P)
colnames(blah.wide)[c(4,5)] <- c("v25","v975")
sub <- blah.wide[blah.wide$spec == "10_Plesiophysa_guadeloupensis",]
pd <- position_dodge(1) # move them .05 to the left and right
group_col <- c("model_predict" = "gray", "model_predict_observed_sites" = "skyblue", "observed" ="orange")

# For one species
ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(),axis.title=element_blank(), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)),text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"),axis.line = element_line(colour = "black"),panel.border = element_blank()) + ggtitle(as.character(spec_correct[10])) + scale_y_continuous(n.breaks=nb,limits=ylim)

## Appendix S2, Figure S2
# Arranged in a grid
plot_list = list()
for(z in 1:27){
  if(z == 6|z==3|z==11|z==2|z==7|z==8|z==17|z==5){
    ylim <- c(0,1)
    nb <- 5
  } else if(z == 16|z==14|z==25|z==1){
    ylim <- c(0,0.5)
    nb <- 5
  } else if(z == 19|z==20|z==10|z==24|z == 21|z==26|z==22|z==27|z==23){
    ylim <- c(0,0.1)
    nb <- 3
  }
  else{
    ylim <- c(0,0.25)
    nb <- 5
  }
  
  sub <- blah.wide[blah.wide$spec == spec[z],]
  
  if (z==24|z==21|z==26|z==22|z==27){
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none", panel.grid=element_blank(), axis.title.y=element_blank(), axis.title.x=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)),text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"),axis.line = element_line(colour = "black",size=.1),panel.border = element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim)
  } else if (z==6|z==16|z==1) {
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank(), axis.text.x=element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied")
  } else if (z==10) {
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied")
  } else {
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd,size=1.2,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.7, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(),axis.title=element_blank(), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)),text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"),axis.line = element_line(colour = "black",size=.1),panel.border = element_blank(),axis.text.x=element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim)
  }
  plot_list[[z]] <- p
}

ordered_plot_list <- plot_list[ordered]
g <- grid.arrange(grobs=ordered_plot_list,legend,ncol=7,nrow=4)

ggsave(file=paste(dir,"/raw_output/S2_Figure_S2.svg",sep=""), g, units="in", width=8.75, height=5.25, dpi=600)

## Figure 1
# Arranged in a grid
plot_list = list()
plot_count <- 1
for(z in c(6,7,12,16:18)){
  if(z == 6|z==3|z==11|z==2|z==7|z==8|z==17|z==5){
    ylim <- c(0,1)
    nb <- 5
  } else if(z == 16|z==14|z==25|z==1){
    ylim <- c(0,0.5)
    nb <- 5
  } else if(z == 19|z==20|z==10|z==24|z == 21|z==26|z==22|z==27|z==23){
    ylim <- c(0,0.1)
    nb <- 3
  }
  else{
    ylim <- c(0,0.25)
    nb <- 5
  }
  sub <- blah.wide[blah.wide$spec == spec[z],]
  if (z==18) {
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd, size=1.5,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.5, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied") + scale_x_continuous(breaks = round(seq(min(sub$year), max(sub$year), by = 1),1))
  } else if (z==6){
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd, size=1.5,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.5, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank(), axis.text.x=element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied") + scale_x_continuous(breaks = round(seq(min(sub$year), max(sub$year), by = 1),1))
  }
  else{
    p <- ggplot(sub,aes(x=year,y=Mean,fill=model,group=model)) + geom_point(position=pd, size=1.5,shape = 21,colour = "black") + geom_errorbar(aes(ymin=v25, ymax=v975), colour="black", width=.5, size=.2, position=pd) + scale_fill_manual(values=group_col) + theme_bw() + theme(legend.position="none",panel.grid=element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=5), plot.title=element_text(face='italic',size=6,margin=margin(0,0,0,0)), text=element_text(size=6),plot.margin=unit(c(5,3,5,3),"pt"), axis.line = element_line(colour = "black",size=.1), panel.border = element_blank(), axis.text.x=element_blank()) + ggtitle(as.character(spec_correct[z])) + scale_y_continuous(n.breaks=nb,limits=ylim) + labs(y= "proportion of sites occupied") + scale_x_continuous(breaks = round(seq(min(sub$year), max(sub$year), by = 1),1))
  }
  plot_list[[plot_count]] <- p
  plot_count <- plot_count + 1
}

g <- grid.arrange(grobs=plot_list,legend,ncol=1,nrow=6)

ggsave(file=paste(dir,"/raw_output/Figure_1/Figure_1_cov.svg",sep=""), g, units="in", width=3, height=8, dpi=600)


### Results misc. Information about parameter estimates given in the results --------
sig <- full[ordered,19:32,1]
sig_vals <- apply(sig,1,function(x) sum(x >= 0.6))
mean(sig_vals)
sd(sig_vals)

### Figure S1-S1 and S1-S2, covariate results -------------------------------
##beta plot function
plot_beta <- function(x,y,lab){
  ##plot outline
  plot(x[,1],1:27,type="n",xlim=c(-2.6,2.6),yaxt="n",xlab="",ylab="",main=lab,cex.axis=0.75)
  axis(2, at=1:27,labels=spec_correct[order(x[,1])],las=1,cex.axis=0.75,font=3)
  #significant highlighted in red
  order <- rep(NA,length(which(y >= 0.6)))
  for(i in 1:length(which(y >= 0.6))){
    order[i] <- which(order(x[,1]) == which(y >= 0.6)[i])
  }
  xleft <- -2.6
  ybottom <- order-.4
  xright <- 2.8
  ytop <- order+.4
  rect(xleft,ybottom,xright,ytop,col=rgb(169,169,169,100,maxColorValue=255),border=NA)
  #points
  points(sort(x[,1]),1:27,pch=1,col="black",cex=1)
  #95% C.I.
  arrows(sort(x[,1]), 1:27, x[(order(x[,1])),2], 1:27, length=0.025, angle=90, code=2, lwd=.75, col="black")
  arrows(sort(x[,1]), 1:27, x[(order(x[,1])),3], 1:27, length=0.025, angle=90, code=2, lwd=.75, col="black")
  abline(v=0,lty=2)
}

spec_lab <- c("Gundlachia radiata", "Drepanotrema lucidum", "Drepanotrema depressissimum", "Drepanotrema cimex", "Drepanotrema aeruginosum", "Aplexa marmorata", "Physa acuta", "Lymnaea columella", "Lymnaea cubensis", "Plesiophysa guadeloupensis", "Biomphalaria straminea", "Biomphalaria glabrata", "Biomphalaria schrammi", "Pomacea glauca", "Pyrgophorus coronatus", "Marisa cornuarietis", "Melanoides tuberculata PAP", "Melanoides tuberculata GOS", "Melanoides tuberculata FAL", "Melanoides tuberculata MAD", "Melanoides tuberculata CPF", "Melanoides tuberculata SEN", "Tarebia granifera", "Neritina virginea", "Eupera viridans", "Indoplanorbis exustus", "Helisoma duryi")

##beta phiW
#pdf(file=paste(dir,"/raw_output/S1_Figure_S1.svg", sep=""), width=8.75, height=5.25)
svg(file=paste(dir,"/raw_output/S1_Figure_S1.svg", sep=""), width=8.75, height=5.25)
par(mar=c(2.5,6.5,1,1),mfrow=c(2,3),mgp=c(3,.75,0))

#beta 5_size
plot_beta(full[,5,],full[,19,1],lab=expression(paste(phi," size",sep="")))
#beta 6_veg
plot_beta(full[,6,],full[,20,1],lab=expression(paste(phi," vegetation",sep="")))
#beta 7_stability
plot_beta(full[,7,],full[,21,1],lab=expression(paste(phi," stability",sep="")))
#beta 8_lrs
plot_beta(full[,8,],full[,22,1],lab=expression(paste(phi," LRS",sep="")))
#beta 9_mangrove
plot_beta(full[,9,],full[,23,1],lab=expression(paste(phi," back-mangrove",sep="")))
#beta 10_river
plot_beta(full[,10,],full[,24,1],lab=expression(paste(phi," river",sep="")))

dev.off()

##beta cW
#pdf(file=paste(dir,"/raw_output/S1_Figure_S2.svg", sep=""), width=16, height=12)
svg(file=paste(dir,"/raw_output/S1_Figure_S2.svg", sep=""), width=8.75, height=5.25)
par(mar=c(2.5,6.5,1,1),mfrow=c(2,4))

#beta 11_size
plot_beta(full[,11,],full[,25,1],lab="c size")
#beta 12_veg
plot_beta(full[,12,],full[,26,1],lab="c vegetation")
#beta 13_stability
plot_beta(full[,13,],full[,27,1],lab="c stability")
#beta 14_connec
plot_beta(full[,14,],full[,28,1],lab="c connectivity")
#beta 15_rs
plot_beta(full[,15,],full[,29,1],lab="c RS")
#beta 16_colsource
plot_beta(full[,16,],full[,30,1],lab="c propagule pressure")
#beta 17_mangrove
plot_beta(full[,17,],full[,31,1],lab="c back-mangrove")
#beta 18_river
plot_beta(full[,18,],full[,32,1],lab="c river")

dev.off()

### Figure 4b, 5. Snail covariate PCA nd DFA results -----------------------------------
##########Conduct PCA
sub <- full[,5:18,1]
cat <- c(rep("pul",13),rep("caeno",10),rep("other",2),rep("pul",2))
cat[rowSums(is.na(sub))>1] <- "missing"
pc_sub <- prcomp(sub,scale=T)
blah <- pc_sub

spec_correct_short <- c("G.rad", "D.sur", "D.dep", "D.cim", "D.aer", "A.mar", "P.acu", "Ps.col", "Ga.cub", "Pl.gua", "B.kuh", "B.gla", "B.sch", "Po.gla", "Py.par", "Ma.cor", "PAP", "GOS", "FAL", "MAD", "CPF", "SEN", "T.gra", "N.vir","E.vir","I.exu", "H.dur")
param_short <-c("φsize","φveg","φstab","φlrs","φman","φriv","Csize","Cveg","Cstab","Cconnec","Crs","Ccol","Cman","Criv")

################## Plot PCA, Figure 4b
lambda <- blah$sdev * sqrt(nrow(blah$x))
svg(file=paste(dir,"/raw_output/Figure_4/Fig_4b.svg",sep=""), width=5.333, height=4)
par(mar=c(3.25,1,0,0),mgp=c(3,1,0))
plot (t(t(blah$x)/lambda),type = "n", asp = 1, axes = F, cex = 2, cex.lab = 2, xlab = "", ylab = "")
#points(t((t(blah$x)/lambda)),pch=19,cex=2,col=c(rep("red",13),rep("blue",9),rep("red",2)))
#text(t((t(blah$x)/lambda)),labels=spec_short,cex = 1, col=c(rep("red",13),rep("blue",9),rep("red",2)), lwd = 3, font=3)
text(t((t(blah$x)/lambda)),labels=spec_correct_short,cex = 0.625, col=c(rep("red",13),rep("blue",9),rep("grey47",2),rep("red",2)), lwd = 3, font=3)
par (new=T) 
Rot <- t(t(blah$rotation)*lambda) 
XLIM <- c(-max(abs(Rot[,1])),max(abs(Rot[,1]))) 
XLIM <- XLIM+(XLIM*0.1) 
#YLIM <- XLIM
#XLIM[1] <- XLIM[1] - 1
plot(Rot,col=4,axes=FALSE,xlim=c(-5,XLIM[2]),ylim=XLIM,pch="", xlab = "", ylab = "") 
arrows (rep(0,nrow(blah$rotation)),rep(0,nrow(blah$rotation)),Rot[,1],Rot[,2],col="black",lwd = 1,length=.1) 
text (Rot[,1:2],param_short,col="black", pos = 2, offset = 1, cex = 0.625, font=2) 
#abline(h = -4.4, lwd = 3)
#abline(v = -4.4, lwd = 3)
axis(1, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4), labels = T,  pos = -5, las = 0, lwd = 1.5, cex.axis = 0.667)
#title(xlab = "PC1", ylab = "PC2", line = 2, cex.lab = 2)
axis(2, at = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4), labels = T, pos = -5, las = 2, lwd = 1.5, cex.axis = 0.667)
dev.off()

###########DFA on the PC axes
sub_analyze <- sub[cat != "other",]
blah <- as.character(cat[cat != "other"])
spec_lab <- spec_correct_short[cat != "other"]
library(adegenet)
#Retain 2 PC axes, 1 LD axis (be definition, only 2 groups), center (to mean 0), and scale (to unit variance)
df <- dapc(sub_analyze, grp = blah, n.pca=2, n.da=1, scale=T, var.loadings = T)

########### DFA plot with modifications, for one dimension (Modified from scatter.dapc, library adegenet). Figure 5

############Establish function inputs
x <- df; xax = 1; yax = 2; grp = x$grp; col = seasun(length(levels(grp))) 
pch = 20; bg = "white"; solid = 0.7; scree.da = TRUE; scree.pca = FALSE 
posi.da = "bottomright"; posi.pca = "bottomleft"; bg.inset = "white"
ratio.da = 0.25; ratio.pca = 0.25; inset.da = 0.02; inset.pca = 0.02
inset.solid = 0.5; onedim.filled = TRUE; mstree = FALSE
lwd = 1; lty = 1; segcol = "black"; legend = FALSE; posi.leg = "topright"
cleg = 1; txt.leg = levels(grp); cstar = 1; cellipse = 1.5
axesell = FALSE; label = levels(grp); clabel = 1; xlim = NULL
ylim = NULL; grid = FALSE; addaxes = TRUE;origin = c(0,0); include.origin = TRUE; sub = "";
csub = 1; possub = "bottomleft"; cgrid = 1; pixmap = NULL; contour = NULL; area = NULL

############Preliminary plotting needs
ONEDIM <- xax == yax | ncol(x$ind.coord) == 1
col <- rep(col, length(levels(grp)))
pch <- rep(pch, length(levels(grp)))
col <- transp(col, solid)
bg.inset <- transp(bg.inset, inset.solid)

############Plotting for single dimension (DF)
scree.da <- FALSE
if (ncol(x$ind.coord) == 1) {
  pcLab <- 1
} else {
  pcLab <- xax
}
ldens <- tapply(x$ind.coord[, pcLab], grp, density)
allx <- unlist(lapply(ldens, function(e) e$x))
ally <- unlist(lapply(ldens, function(e) e$y))
#### STOPPED
svg(file=paste(dir,"/raw_output/Figure_5/Figure_5_part1.svg",sep=""),width=6,height=4)
par(bg = bg,mar=c(0,0,0,0),mgp=c(3,1,0))
plot(allx, ally, type = "n", xlab = paste(""), ylab = "", axes=F,ylim=c(-.1,.5))
axis(1, at = c(-3, -2, -1, 0, 1, 2, 3), labels = T,  pos = -.075, las = 0, lwd = 2, cex.axis = 0.883)
axis(2, at = c(0, .1, .2, .3, .4, .5), labels = T,  pos = -3.6, las = 0, lwd = 2, cex.axis = 0.883)

############Add on species labels: create a new vector with locations along DF for text labels
loc <- x$ind.coord
lab <- spec_lab
#Space out PAP (17), Po. glauca (14), and MAD (20). And Ma. cor (16), Py. par (15), FAL (19).
##Caeno
loc[17] <- -1.8
loc[15] <- -0.9
loc[16] <- -0.68
loc[19] <- -0.5

#Pulmo
loc[11] <- 0.2
loc[4] <- 0.4
loc[25] <- 0.7 # H. duryi
loc[10] <- 0.95

loc[6] <- 1.2
loc[3] <- 1.4
loc[24] <- 1.6
loc[12] <- 1.8
loc[13] <- 2
loc[2] <- 2.2

loc[5] <- 2.6

for (i in 1:length(ldens)) {
  if (!onedim.filled) {
    lines(ldens[[i]]$x, ldens[[i]]$y, col = col[i], 
          lwd = 2)
  }
  else {
    polygon(c(ldens[[i]]$x, rev(ldens[[i]]$x)), c(ldens[[i]]$y, rep(0, length(ldens[[i]]$x))), col = col[i], lwd = 2, border = col[i])
  }
  points(x = x$ind.coord[grp == levels(grp)[i], pcLab], y = rep(0, sum(grp == levels(grp)[i])), pch = "|", col = col[i])
}

#Species labels
i <- 1
text(x = loc[grp == levels(grp)[i], pcLab], y = rep(-.04, sum(grp == levels(grp)[i])), labels = lab[grp == levels(grp)[i]], cex=0.833, col="black", srt=45)

i <- 2
text(x = loc[grp == levels(grp)[i], pcLab], y = rep(0.03, sum(grp == levels(grp)[i])), labels = lab[grp == levels(grp)[i]], cex=0.833, col="black", srt=45)

############Add cut-out lines for variables in offset plot
segments(x0=-.4,  y0=-.15, x1=-.4, y1=-.1, lwd=2)
segments(x0=.4,  y0=-.15, x1=.4, y1=-.1, lwd=2)
dev.off()

############Now create inset plot of environmental variables along DF1
#Plot DF1
svg(file=paste(dir,"/raw_output/Figure_5/Figure_5_part2.svg",sep=""),width=6,height=4)
plot(df$var.load, rep(0,length(df$var.load)), type="n", axes=F, xlab="", ylab="",xlim=c(-.4,.4))
axis(1, at = c(-.4, -.2, 0, .2, .4), labels = T,  pos = 0, las = 0, lwd = 2, cex.axis = 0.883)
#Add variables (beta coefficients)
param <-c("φsize","φveg","φstab","φlrs","φman","φriv","Csize","Cveg","Cstab","Cconnec","Crs","Ccol","Cman","Criv")
df_shift <- df$var.load
df_shift[c(11,6,7,10,12,4),] <- c(-.32,-.229,-.2,-.17,-.14,.13)
points(df$var.load, rep(0,length(df$var.load)), pch=19, cex=1, col="black")
text(x=df_shift, y=rep(.22,length(df$var.load)), labels=param, cex=0.883, col="black", srt=45)
dev.off()

### Simulation model input and results. Needs preparation of some files for QGIS maps in 'Figure 6' below --------------------------------------
##varname: psi = a 278*27 record of average occupancy probabilities
####Get psi as average occupancy per site
mat_mean <- function(dat){
  blah <- matrix(dat, nrow=278, ncol=15, byrow=F)
  blah <- as.data.frame(blah)
  blah <- data.matrix(blah)
  result <- rowMeans(blah)
  return(result)
}
psi <- array(NA, dim=c(278,num_spec), dimnames=list(NULL,spec))

for(z in 1:num_spec){
  #Load results file
  load(paste(dir,"/results/covar/",spec[z],"/jsample.RData", sep=""))
  #Get z estimates
  df=as.data.frame(rbind(jsample.full[[1]],jsample.full[[2]],jsample.full[[3]]))
  df <- df[,34:4203]
  #Extract mean
  occupy <- apply(df, 1, mat_mean)
  psi[,z] <- rowMeans(occupy)
}
psi <- psi-1  #Gives the percentage of years (across all JAGS iterations) that each site is occupied
write.csv(psi,file=paste(dir,"/data/simulation/parameters/param_psi_bysite.csv",sep=""),row.names=F)

##varname: site_obs, a list of all sites for which a species was ever observed.
# #Read in data for whether or not a species was ever detected at a site
obs <- read.csv(paste(dir,"/data/simulation/species_obs.csv",sep=""), header=F)
site_obs <- vector("list", num_spec)
i <- 1
while(i < 28) {
  site_obs[[i]] <- which(obs[,i]==1) +1
  #Add 1 to all elements since Mathematica has a header on the first row (and thus sites are 2-279)
  i <- i + 1
}
lapply(site_obs,write,file=paste(dir,"/data/simulation/parameters/sites_observed_in.txt",sep=""),append=TRUE,ncolumns=1000)
#Note if for some reason I need to re-create this file, I have to delete the old version (otherwise it just appends the new values to the bottom of the text file).


### Figure 6. QGIS maps for species 3, 6, 7, 12, and 16 showing site e/c --------
##varname: site = a 27*5*278*15*5 array with the site- and year-specific eW, eD, cW, EW, and ED values for a few scenarios: (1) using the model intercept only, (2) using all model covariates, (3) using signigicant model covariates only, (4) using all covariates with number of sites limited to only where the species was ever observed, and (5) using significant covariates with number of sites limited to only where the species was ever observed.

##varname: obs = a 278*27 record of whether each species was ever observed at each site

##varname: obs_site = a 278*27 record of the # of occurrences / the # of site visits total

##varname: site_coord_Ec, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site e/c value)

##varname: site_coord_pstar, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site p* value)

##varname: site_coord_pstar_dev, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site deviation of observed from p* value)

##varname: site_coord_pstar_sim, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site ci*p / ci*p + Ei)

##varname: site_coord_pstar_sim_detect, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site ci*p / ci*p + Ei multipled by detection probability)

##varname: site_coord_pstar_sim_detect_dev, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site deviation of site ci*p / ci*p + Ei multipled by detection probability from p* value)

site <- array(NA, dim=c(num_spec,5,278,15,5), dimnames=list(spec,c("site_intercept","site_all_covar","site_sig_covar","site_reduced_all_covar","site_reduced_sig_covar"),NULL,c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015"),c("eW","eD","cW","EW","ED")))

obs <- array(NA, dim=c(278,27), dimnames=list(NULL, spec))
obs_site <- array(NA, dim=c(278,27), dimnames=list(NULL, spec))

##Read in covariates
ndir <- paste(dir,"/results/covar/1_Gundlachia_radiata/",sep="")

# covariates site-specific
cov <- read.table(paste(ndir,"cov.txt",sep=""),h=T)

size = cov$size
size <- (size-mean(size))/sd(size)
connec = cov$connec
connec <- (connec-mean(connec))/sd(connec)
veg = cov$veg
veg <- (veg-mean(veg))/sd(veg)
stability = cov$stability
stability <- (stability-mean(stability))/sd(stability)

habit <- read.table(paste(ndir,"habit.txt",sep=""),h=T)

mangrove = habit$mangrove
mangrove <- (mangrove-mean(mangrove))/sd(mangrove)
river = habit$river
river <- (river-mean(river))/sd(river)

# covariates year-specific
# little rainy season (lrs): cumulative rainfall during the little rainy season from March 1 to May 31, 92 days, e.g. during the period of putative extinction
# rainy season (rs): cumulative rainfall during the rainy season from (i) July 1 - Dec. 31 of year preceding sampling campaign (184 days, lrs_July_Dec) and (ii) 60 days preceding the first date of the sampling campaign (60 days, rs_60_days)

pluvio <- read.table(paste(ndir,"pluvio.txt",sep=""),h=T)
lrs=pluvio$lrs
rs_1=pluvio$rs_July_Dec
#rs_2=pluvio$rs_60_days

lrs <- (lrs-mean(lrs))/sd(lrs)
rs_1 <- (rs_1-mean(rs_1))/sd(rs_1)
#rs_2 <- (rs_2-mean(rs_2))/sd(rs_2)

for(z in 1:num_spec){
  #Read in covariates year and site-specific
  colsource <- matrix(scan(paste(dir,"/results/covar/",spec[z],"/colsource.txt",sep=""),n = 278*15), 278,15,byrow = TRUE)
  colsource <- (colsource-mean(colsource))/sqrt(mean((colsource-mean(colsource))*(colsource-mean(colsource))))
  
  ##site_intercept
  #eW
  site[z,1,,,1] <- 1-array(things[z,5,1],dim=c(278,15))
  #eD
  site[z,1,,,2] <- 1-array(things[z,6,1],dim=c(278,15))
  #cW
  site[z,1,,,3] <- array(things[z,7,1],dim=c(278,15))
  #EW
  site[z,1,,,4] <- array(things[z,8,1],dim=c(278,15))
  #ED
  site[z,1,,,5] <- array(things[z,9,1],dim=c(278,15))
  
  ##site_all_covar
  include <- rep(1,33)
  #eW
  site[z,2,,,1] <- 1-logit_to_prob(matrix(full[z,1,1],nrow=278,ncol=15) + matrix(include[5]*full[z,5,1]*size,nrow=278,ncol=15,byrow=F) + matrix(include[6]*full[z,6,1]*veg,nrow=278,ncol=15,byrow=F) + matrix(include[7]*full[z,7,1]*stability,nrow=278,ncol=15,byrow=F) + matrix(include[8]*full[z,8,1]*lrs,nrow=278,ncol=15,byrow=T) + matrix(include[9]*full[z,9,1]*mangrove,nrow=278,ncol=15,byrow=F) + matrix(include[10]*full[z,10,1]*river,nrow=278,ncol=15,byrow=F))
  #eD
  site[z,2,,,2] <- 1-logit_to_prob(matrix(full[z,2,1],nrow=278,ncol=15))
  #cW
  site[z,2,,,3] <- exp(matrix(full[z,3,1],nrow=278,ncol=15) + matrix(include[11]*full[z,11,1]*size,nrow=278,ncol=15,byrow=F) + matrix(include[12]*full[z,12,1]*veg,nrow=278,ncol=15,byrow=F) + matrix(include[13]*full[z,13,1]*stability,nrow=278,ncol=15,byrow=F) + matrix(include[14]*full[z,14,1]*connec,nrow=278,ncol=15,byrow=F) + matrix(include[15]*full[z,15,1]*rs_1,nrow=278,ncol=15,byrow=T) + include[16]*full[z,16,1]*colsource + matrix(include[17]*full[z,17,1]*mangrove,nrow=278,ncol=15,byrow=F) + matrix(include[18]*full[z,18,1]*river,nrow=278,ncol=15,byrow=F))
  #EW
  site[z,2,,,4] <- -log(logit_to_prob(matrix(full[z,1,1],nrow=278,ncol=15) + matrix(include[5]*full[z,5,1]*size,nrow=278,ncol=15,byrow=F) + matrix(include[6]*full[z,6,1]*veg,nrow=278,ncol=15,byrow=F) + matrix(include[7]*full[z,7,1]*stability,nrow=278,ncol=15,byrow=F) + matrix(include[8]*full[z,8,1]*lrs,nrow=278,ncol=15,byrow=T) + matrix(include[9]*full[z,9,1]*mangrove,nrow=278,ncol=15,byrow=F) + matrix(include[10]*full[z,10,1]*river,nrow=278,ncol=15,byrow=F)))
  #ED
  site[z,2,,,5] <- -log(logit_to_prob(matrix(full[z,2,1],nrow=278,ncol=15)))
  
  ##site_sig_covar
  include <- rep(0,33)
  for(i in 5:18){
    if (full[z,i+14,1] >= 0.6) {include[i] <- 1}
  }
  #eW
  site[z,3,,,1] <- 1-logit_to_prob(matrix(full[z,1,1],nrow=278,ncol=15) + matrix(include[5]*full[z,5,1]*size,nrow=278,ncol=15,byrow=F) + matrix(include[6]*full[z,6,1]*veg,nrow=278,ncol=15,byrow=F) + matrix(include[7]*full[z,7,1]*stability,nrow=278,ncol=15,byrow=F) + matrix(include[8]*full[z,8,1]*lrs,nrow=278,ncol=15,byrow=T) + matrix(include[9]*full[z,9,1]*mangrove,nrow=278,ncol=15,byrow=F) + matrix(include[10]*full[z,10,1]*river,nrow=278,ncol=15,byrow=F))
  #eD
  site[z,3,,,2] <- 1-logit_to_prob(matrix(full[z,2,1],nrow=278,ncol=15))
  #cW
  site[z,3,,,3] <- exp(matrix(full[z,3,1],nrow=278,ncol=15) +  matrix(include[11]*full[z,11,1]*size,nrow=278,ncol=15,byrow=F) + matrix(include[12]*full[z,12,1]*veg,nrow=278,ncol=15,byrow=F) + matrix(include[13]*full[z,13,1]*stability,nrow=278,ncol=15,byrow=F) + matrix(include[14]*full[z,14,1]*connec,nrow=278,ncol=15,byrow=F) + matrix(include[15]*full[z,15,1]*rs_1,nrow=278,ncol=15,byrow=T) + include[16]*full[z,16,1]*colsource + matrix(include[17]*full[z,17,1]*mangrove,nrow=278,ncol=15,byrow=F) + matrix(include[18]*full[z,18,1]*river,nrow=278,ncol=15,byrow=F))
  #EW
  site[z,3,,,4] <- -log(logit_to_prob(matrix(full[z,1,1],nrow=278,ncol=15) + matrix(include[5]*full[z,5,1]*size,nrow=278,ncol=15,byrow=F) + matrix(include[6]*full[z,6,1]*veg,nrow=278,ncol=15,byrow=F) + matrix(include[7]*full[z,7,1]*stability,nrow=278,ncol=15,byrow=F) + matrix(include[8]*full[z,8,1]*lrs,nrow=278,ncol=15,byrow=T) + matrix(include[9]*full[z,9,1]*mangrove,nrow=278,ncol=15,byrow=F) + matrix(include[10]*full[z,10,1]*river,nrow=278,ncol=15,byrow=F)))
  #ED
  site[z,3,,,5] <- -log(logit_to_prob(matrix(full[z,2,1],nrow=278,ncol=15)))
  
  ##site_reduced_all_covar
  S <- matrix(scan(paste(dir,"/results/covar/",spec[z],"/state.txt", sep=""), n = 278*15,skip=1), 278, 15, byrow = TRUE)
  
  S[S==0] <- NA
  S[S==1] <- 0
  S[S==2] <- 1
  
  obs[,z] <- rowSums(S,na.rm=TRUE) > 0
  obs <- obs*1
  
  obs_site[,z] <- (rowSums(S,na.rm=T)) / (apply(S,1,function(x) sum(is.na(x)==F)))
  
  #eW
  site[z,4,obs[,z]==1,,1] <- site[z,2,obs[,z]==1,,1]
  #eD
  site[z,4,obs[,z]==1,,2] <- site[z,2,obs[,z]==1,,2]
  #cW
  site[z,4,obs[,z]==1,,3] <- site[z,2,obs[,z]==1,,3]
  #EW
  site[z,4,obs[,z]==1,,4] <- site[z,2,obs[,z]==1,,4]
  #ED
  site[z,4,obs[,z]==1,,5] <- site[z,2,obs[,z]==1,,5]
  
  ##site_reduced_sig_covar
  #eW
  site[z,5,obs[,z]==1,,1] <- site[z,3,obs[,z]==1,,1]
  #eD
  site[z,5,obs[,z]==1,,2] <- site[z,3,obs[,z]==1,,2]
  #cW
  site[z,5,obs[,z]==1,,3] <- site[z,3,obs[,z]==1,,3]
  #EW
  site[z,5,obs[,z]==1,,4] <- site[z,3,obs[,z]==1,,4]
  #ED
  site[z,5,obs[,z]==1,,5] <- site[z,3,obs[,z]==1,,5]
}

#write files to simulation folder for later use
for(z in 1:num_spec){
  for(d in 1:dim(site)[2]){
    write.csv(site[z,d,,,1], paste(dir,"/data/simulation/parameters/",spec[z],"_",dimnames(site)[2][[1]][d],"_","eW",".csv",sep=""),row.names=F)
    
    write.csv(site[z,d,,,2], paste(dir,"/data/simulation/parameters/",spec[z],"_",dimnames(site)[2][[1]][d],"_","eD",".csv",sep=""),row.names=F)
    
    write.csv(site[z,d,,,3], paste(dir,"/data/simulation/parameters/",spec[z],"_",dimnames(site)[2][[1]][d],"_","cW",".csv",sep=""),row.names=F)
    
    write.csv(site[z,d,,,4], paste(dir,"/data/simulation/parameters/",spec[z],"_",dimnames(site)[2][[1]][d],"_","Levins_EW",".csv",sep=""),row.names=F)
    
    write.csv(site[z,d,,,5], paste(dir,"/data/simulation/parameters/",spec[z],"_",dimnames(site)[2][[1]][d],"_","Levins_ED",".csv",sep=""),row.names=F)
  }
}
#Build site_coord_Ec, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site E/c value)
site_coord_Ec <- data.frame(array(NA,dim=c(278,32)))

#Build site_coord_pstar, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site p* value)
site_coord_pstar <- data.frame(array(NA,dim=c(278,32)))

#Build site_coord_pstar_dev, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site deviation from p* value)
site_coord_pstar_dev <- data.frame(array(NA,dim=c(278,32)))

#Build site_coord_pstar_sim, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site ci*p / ci*p + Ei)
site_coord_pstar_sim <- data.frame(array(NA,dim=c(278,32)))

#Build site_coord_pstar_sim, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site ci*p / ci*p + Ei multiplied by detection probability)
site_coord_pstar_sim_detect <- data.frame(array(NA,dim=c(278,32)))

#Build site_coord_pstar_sim_dev, a 278*32 file with information needed to create GIS layers for each species (site name, lat-long coordinates, site observed deviation from ci*p / ci*p + Ei)
site_coord_pstar_sim_detect_dev <- data.frame(array(NA,dim=c(278,32)))

#Import site coordinant and names
site_list <- read.csv(paste(dir,"/data/final_site_list_coordinates.csv",sep=""))
#Import list of final sites used in order of all data (278 sites)
site_used <- read.table(paste(dir,"/data/final_island_list.txt",sep=""), header=T)
#Trim excluded sites
dat <- subset(site_list, Numero.du.site %in% site_used$site)
#Re-order to final 1-278 order
dat_order <- dat[match(site_used$site, dat$Numero.du.site),]
#Read in R-friendly site names
site_format_name <- read.table(paste(dir,"/data/Visite_offline.txt",sep=""), colClasses=c(NA,NA,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL"),header=T, blank.lines.skip=F, sep="\t")
#Trim R-friendly site names to sites used
site_format_name <- site_format_name[duplicated(site_format_name$site)==F,]
site_format_name <- subset(site_format_name, site %in% dat_order$Numero.du.site)
site_format_name <- site_format_name[match(dat_order$Numero.du.site,site_format_name$site),]
dat_order <- cbind(dat_order,site_format_name)
#Paste R-friendly site names into site_coord_ec
site_coord_Ec[,1:5] <- dat_order[,c(6,7,5,3,4)]
colnames(site_coord_Ec)[1:5] <- colnames(dat_order)[c(6,7,5,3,4)]

site_coord_pstar[,1:5] <- dat_order[,c(6,7,5,3,4)]
colnames(site_coord_pstar)[1:5] <- colnames(dat_order)[c(6,7,5,3,4)]

site_coord_pstar_dev[,1:5] <- dat_order[,c(6,7,5,3,4)]
colnames(site_coord_pstar_dev)[1:5] <- colnames(dat_order)[c(6,7,5,3,4)]

site_coord_pstar_sim[,1:5] <- dat_order[,c(6,7,5,3,4)]
colnames(site_coord_pstar_sim)[1:5] <- colnames(dat_order)[c(6,7,5,3,4)]

site_coord_pstar_sim_detect[,1:5] <- dat_order[,c(6,7,5,3,4)]
colnames(site_coord_pstar_sim_detect)[1:5] <- colnames(dat_order)[c(6,7,5,3,4)]

site_coord_pstar_sim_detect_dev[,1:5] <- dat_order[,c(6,7,5,3,4)]
colnames(site_coord_pstar_sim_detect_dev)[1:5] <- colnames(dat_order)[c(6,7,5,3,4)]

##Manually change coordinates for 2 sites that were measured incorrectly in this version of the Access database. As per an e-mail from JHP to PJ and PD on October 2, 2015, the Access values for two sites were measured incorrectly. Patrice re-measured them to be: Ducos = 1758.25 N, 683.4 E and Menard_Cambrai = 1770 N, 685.1 E. To put these in the WSG 84 UTM Zone 20N scale (as explained in an e-mail to Ana Lozano-Campo on March 29, 2016), I convert these values from km to m (multiply by 1000) and subtract 425m from the Easting direction and subtracting 290m from the Northing direction.
site_coord_Ec$E[site_coord_Ec$nom == "Ducos"] <- (683.4*1000)-425
site_coord_Ec$N[site_coord_Ec$nom == "Ducos"] <- (1758.25*1000)-290
site_coord_Ec$E[site_coord_Ec$nom == "Menard_Cambrai"] <- (685.1*1000)-425
site_coord_Ec$N[site_coord_Ec$nom == "Menard_Cambrai"] <- (1770*1000)-290

site_coord_pstar$E[site_coord_pstar$nom == "Ducos"] <- (683.4*1000)-425
site_coord_pstar$N[site_coord_pstar$nom == "Ducos"] <- (1758.25*1000)-290
site_coord_pstar$E[site_coord_pstar$nom == "Menard_Cambrai"] <- (685.1*1000)-425
site_coord_pstar$N[site_coord_pstar$nom == "Menard_Cambrai"] <- (1770*1000)-290

site_coord_pstar_dev$E[site_coord_pstar_dev$nom == "Ducos"] <- (683.4*1000)-425
site_coord_pstar_dev$N[site_coord_pstar_dev$nom == "Ducos"] <- (1758.25*1000)-290
site_coord_pstar_dev$E[site_coord_pstar_dev$nom == "Menard_Cambrai"] <- (685.1*1000)-425
site_coord_pstar_dev$N[site_coord_pstar_dev$nom == "Menard_Cambrai"] <- (1770*1000)-290

site_coord_pstar_sim$E[site_coord_pstar$nom == "Ducos"] <- (683.4*1000)-425
site_coord_pstar_sim$N[site_coord_pstar$nom == "Ducos"] <- (1758.25*1000)-290
site_coord_pstar_sim$E[site_coord_pstar$nom == "Menard_Cambrai"] <- (685.1*1000)-425
site_coord_pstar_sim$N[site_coord_pstar$nom == "Menard_Cambrai"] <- (1770*1000)-290

site_coord_pstar_sim_detect$E[site_coord_pstar$nom == "Ducos"] <- (683.4*1000)-425
site_coord_pstar_sim_detect$N[site_coord_pstar$nom == "Ducos"] <- (1758.25*1000)-290
site_coord_pstar_sim_detect$E[site_coord_pstar$nom == "Menard_Cambrai"] <- (685.1*1000)-425
site_coord_pstar_sim_detect$N[site_coord_pstar$nom == "Menard_Cambrai"] <- (1770*1000)-290

site_coord_pstar_sim_detect_dev$E[site_coord_pstar$nom == "Ducos"] <- (683.4*1000)-425
site_coord_pstar_sim_detect_dev$N[site_coord_pstar$nom == "Ducos"] <- (1758.25*1000)-290
site_coord_pstar_sim_detect_dev$E[site_coord_pstar$nom == "Menard_Cambrai"] <- (685.1*1000)-425
site_coord_pstar_sim_detect_dev$N[site_coord_pstar$nom == "Menard_Cambrai"] <- (1770*1000)-290

#Add Ei/ci, p*, and deviation from p* for each species, weighted by Pr(wet) and averaged across years
for(z in 1:num_spec){
  site_coord_Ec[,(z+5)] <- rowMeans((site[z,2,,,4]*wet_prob + site[z,2,,,5]*(1-wet_prob))/(site[z,2,,,3]*wet_prob))
  site_coord_pstar[,(z+5)] <- (1 - (rowMeans((site[z,2,,,4]*wet_prob + site[z,2,,,5]*(1-wet_prob))/(site[z,2,,,3]*wet_prob))))*full[z,33,1]
  #There can never be observed values less than 0, so I change all the deviations where site_coord_pstar is < 0 to 0
  site_coord_pstar[site_coord_pstar[,(z+5)] < 0,(z+5)] <- 0
  site_coord_pstar_dev[,(z+5)] <- obs_site[,z] - site_coord_pstar[,(z+5)]
  
  #Get p* from simulation and use it to calculate ci*p / ci*p + Ei and the deviation of observed value from that
  #This is from simulation model run using all covariates
  dat <- read.csv(paste(dir,"/raw_output/simulation/simulmetapop_V04.01/",spec[z],"_result.csv",sep=""),header=F)
  pstar <- colMeans(dat)[999]
  
  #predicted occupancy
  site_coord_pstar_sim[,(z+5)] <- rowMeans(((site[z,2,,,3]*wet_prob)*pstar)/(((site[z,2,,,3]*wet_prob)*pstar) + (site[z,2,,,4]*wet_prob + site[z,2,,,5]*(1-wet_prob))))
  #predicted detection
  site_coord_pstar_sim_detect[,(z+5)] <- site_coord_pstar_sim[,(z+5)]*full[z,33,1]
  #difference between predicted and observed detection
  site_coord_pstar_sim_detect_dev[,(z+5)] <- obs_site[,z] - site_coord_pstar_sim_detect[,(z+5)]
}

colnames(site_coord_Ec)[6:32] <- spec
colnames(site_coord_pstar)[6:32] <- spec
colnames(site_coord_pstar_dev)[6:32] <- spec
colnames(site_coord_pstar_sim)[6:32] <- spec
colnames(site_coord_pstar_sim_detect)[6:32] <- spec
colnames(site_coord_pstar_sim_detect_dev)[6:32] <- spec

write.csv(site_coord_Ec, paste(dir,"/raw_output/maps/site_coord_Ec.csv",sep=""),col.names=F)
write.csv(site_coord_pstar, paste(dir,"/raw_output/maps/site_coord_pstar.csv",sep=""),col.names=F)
write.csv(site_coord_pstar_dev, paste(dir,"/raw_output/maps/site_coord_pstardev.csv",sep=""),col.names=F)
write.csv(site_coord_pstar_sim, paste(dir,"/raw_output/maps/site_coord_pstar_sim.csv",sep=""),col.names=F)
write.csv(site_coord_pstar_sim_detect, paste(dir,"/raw_output/maps/site_coord_pstar_sim_detect.csv",sep=""),col.names=F)
write.csv(site_coord_pstar_sim_detect_dev, paste(dir,"/raw_output/maps/site_coord_pstar_sim_detect_dev.csv",sep=""),col.names=F)

### These files are used to create maps of habitat suitability in QGIS. The .svg files for each species are located in paste(dir,"/raw_output/maps/Ec_maps",sep=""). The maps in Figure 6 and Appendix 2, Figure S3 are created from the file site_coord_Ec.csv. The maps in Appendix 2, Figure S4 are created from the file site_coord_pstar_sim_detect_dev.csv.

### Figure 7, Figure S1.S3, Figure S1.S4. Plot e x c data ---------------------------------------------------------

site_ce <- function(x,lim,lab){
  #Function to plot data and confidence ellipses for 5 scenarios
  col=c("black","red","forestgreen","blue","orange")
  
  if(lab == "x"){
    plot(x[,1], x[,2],xlim=lim,ylim=lim,axes=FALSE,type="n",ylab="",xlab="ln(e)",bty="n")
    title(spec_correct[z],line=-1,font.main=3)
    axis(1,at=c(-6,-4,-2,0,2))
    axis(2,at=c(-6,-4,-2,0,2),labels=FALSE)
  
    #Add points
    points(x[,1], x[,2],pch=rep(1,dim(x)[1]),col=rep("gray",dim(x)[1]))
  
    #Add ellipse for each group
    x[,3] <- droplevels(x[,3])  #Drop unused factor levels
  
    if (any(x[,3] == 1)){ #If scenario 1 was included
      points(mean(x[x[,3]==1,1]),mean(x[x[,3]==1,2]),pch=19,cex=2,col="black")
      x <- x[x[,3] != 1,]  #Drop scenario 1 from the data
      x[,3] <- droplevels(x[,3])  #Drop unused factor levels
    }
  
    if (any(x[,3] ==4) | any(x[,3] ==5)){ #If scenario 3 and 4 are included
      x <- x[is.na(x[,1]) == F,]  #Trim NA values
    }
  
    for(i in 1:length(levels(x[,3]))){  #For each scenario...
    
      #Draw data ellipse
      try(dataEllipse(x[x[,3]==levels(x[,3])[i],1], x[x[,3]==levels(x[,3])[i],2],add=T,center.cex=2,levels=.95,ellipse.label="",lty=2,fill=TRUE,fill.alpha=0,grid=F,col=col[as.numeric(levels(x[,3])[i])],plot.points = F),silent=T)
  }
  
    #Add 1:1 line
    lines(c(-6,2),c(-6,2))
    }
  else if (lab == "y"){
    plot(x[,1], x[,2],xlim=lim,ylim=lim,axes=FALSE,type="n",ylab="ln(c)",xlab="")
    title(spec_correct[z],line=-1,font.main=3)
    axis(1,at=c(-6,-4,-2,0,2),labels=FALSE)
    axis(2,at=c(-6,-4,-2,0,2))
    
    #Add points
    points(x[,1], x[,2],pch=rep(1,dim(x)[1]),col=rep("gray",dim(x)[1]))
    
    #Add ellipse for each group
    x[,3] <- droplevels(x[,3])  #Drop unused factor levels
    
    if (any(x[,3] == 1)){ #If scenario 1 was included
      points(mean(x[x[,3]==1,1]),mean(x[x[,3]==1,2]),pch=19,cex=2,col="black")
      x <- x[x[,3] != 1,]  #Drop scenario 1 from the data
      x[,3] <- droplevels(x[,3])  #Drop unused factor levels
    }
    
    if (any(x[,3] ==4) | any(x[,3] ==5)){ #If scenario 3 and 4 are included
      x <- x[is.na(x[,1]) == F,]  #Trim NA values
    }
    
    for(i in 1:length(levels(x[,3]))){  #For each scenario...
      
      #Draw data ellipse
      try(dataEllipse(x[x[,3]==levels(x[,3])[i],1], x[x[,3]==levels(x[,3])[i],2],add=T,center.cex=2,levels=.95,ellipse.label="",lty=2,fill=TRUE,fill.alpha=0,grid=F,col=col[as.numeric(levels(x[,3])[i])],plot.points = F),silent=T)
    }
    
    #Add 1:1 line
    lines(c(-6,2),c(-6,2))
  }
  else if (lab == "xy"){
    plot(x[,1], x[,2],xlim=lim,ylim=lim,axes=FALSE,type="n",ylab="ln(c)",xlab="ln(e)")
    title(spec_correct[z],line=-1,font.main=3)
    axis(1,at=c(-6,-4,-2,0,2))
    axis(2,at=c(-6,-4,-2,0,2))
    
    #Add points
    points(x[,1], x[,2],pch=rep(1,dim(x)[1]),col=rep("gray",dim(x)[1]))
    
    #Add ellipse for each group
    x[,3] <- droplevels(x[,3])  #Drop unused factor levels
    
    if (any(x[,3] == 1)){ #If scenario 1 was included
      points(mean(x[x[,3]==1,1]),mean(x[x[,3]==1,2]),pch=19,cex=2,col="black")
      x <- x[x[,3] != 1,]  #Drop scenario 1 from the data
      x[,3] <- droplevels(x[,3])  #Drop unused factor levels
    }
    
    if (any(x[,3] ==4) | any(x[,3] ==5)){ #If scenario 3 and 4 are included
      x <- x[is.na(x[,1]) == F,]  #Trim NA values
    }
    
    for(i in 1:length(levels(x[,3]))){  #For each scenario...
      
      #Draw data ellipse
      try(dataEllipse(x[x[,3]==levels(x[,3])[i],1], x[x[,3]==levels(x[,3])[i],2],add=T,center.cex=2,levels=.95,ellipse.label="",lty=2,fill=TRUE,fill.alpha=0,grid=F,col=col[as.numeric(levels(x[,3])[i])],plot.points = F),silent=T)
    }
    
    #Add 1:1 line
    lines(c(-6,2),c(-6,2))
  } else {
    plot(x[,1], x[,2],xlim=lim,ylim=lim,axes=FALSE,type="n",ylab="",xlab="")
    title(spec_correct[z],line=-1,font.main=3)
    axis(1,at=c(-6,-4,-2,0,2),labels = FALSE)
    axis(2,at=c(-6,-4,-2,0,2),labels = FALSE)
    
    #Add points
    points(x[,1], x[,2],pch=rep(1,dim(x)[1]),col=rep("gray",dim(x)[1]))
    
    #Add ellipse for each group
    x[,3] <- droplevels(x[,3])  #Drop unused factor levels
    
    if (any(x[,3] == 1)){ #If scenario 1 was included
      points(mean(x[x[,3]==1,1]),mean(x[x[,3]==1,2]),pch=19,cex=2,col="black")
      x <- x[x[,3] != 1,]  #Drop scenario 1 from the data
      x[,3] <- droplevels(x[,3])  #Drop unused factor levels
    }
    
    if (any(x[,3] ==4) | any(x[,3] ==5)){ #If scenario 3 and 4 are included
      x <- x[is.na(x[,1]) == F,]  #Trim NA values
    }
    
    for(i in 1:length(levels(x[,3]))){  #For each scenario...
      
      #Draw data ellipse
      try(dataEllipse(x[x[,3]==levels(x[,3])[i],1], x[x[,3]==levels(x[,3])[i],2],add=T,center.cex=2,levels=.95,ellipse.label="",lty=2,fill=TRUE,fill.alpha=0,grid=F,col=col[as.numeric(levels(x[,3])[i])],plot.points = F),silent=T)
    }
    
    #Add 1:1 line
    lines(c(-6,2),c(-6,2))
  }
}

##To plot all species with enough observations
library(car)
lab <- c("y","none","none","xy","x","x")
count <- 1

## Figure 7
svg(file=paste(dir,"/raw_output/Figure_7.svg", sep=""), width=6, height=4)
par(mfrow = c(2,3),mar=c(3,3,0,0))

for (z in c(6,12,17,7,16,18)){
  combine_scenario <- data.frame()
  #Ei, ci (e.g. weighted by the probability a site was wet)
  combine_scenario <- as.data.frame(cbind(c(rowMeans(site[z,1,,,4]*wet_prob + site[z,1,,,5]*(1-wet_prob)),rowMeans(site[z,2,,,4]*wet_prob + site[z,2,,,5]*(1-wet_prob)),rowMeans(site[z,3,,,4]*wet_prob + site[z,3,,,5]*(1-wet_prob)),rowMeans(site[z,4,,,4]*wet_prob + site[z,4,,,5]*(1-wet_prob)),rowMeans(site[z,5,,,4]*wet_prob + site[z,5,,,5]*(1-wet_prob))),c(rowMeans(site[z,1,,,3]*wet_prob),rowMeans(site[z,2,,,3]*wet_prob),rowMeans(site[z,3,,,3]*wet_prob),rowMeans(site[z,4,,,3]*wet_prob),rowMeans(site[z,5,,,3]*wet_prob))))
  
  ##eW, cW
  combine_scenario$scenario <- as.factor(c(rep(1,278),rep(2,278),rep(3,278),rep(4,278),rep(5,278)))
  ln_combine_scenario <- combine_scenario
  ln_combine_scenario[,1:2] <- log(ln_combine_scenario[,1:2])
  lim <- c(-6,3) # for subset
  site_ce(ln_combine_scenario[ln_combine_scenario[,3] == 1 | ln_combine_scenario[,3] == 2,],lim,lab=lab[count])
  #site_ce(ln_combine_scenario[ln_combine_scenario[,3] == 1 | ln_combine_scenario[,3] == 2 | ln_combine_scenario[,3] == 4,],lim,lab=lab[count])
  count <- count + 1
}
dev.off()

## Figure S1_S3
site_ce_full <- function(x,lim,lab){
  #Function to plot data and confidence ellipses for 5 scenarios
  col=c("black","red","forestgreen","blue","orange")
  plot(x[,1], x[,2],xlim=lim,ylim=lim,axes=FALSE,type="n",bty="n",xaxt="n",xlab="",ylab="")
  title(spec_correct[z],line=-0.4,font.main=3,cex.main=0.883,adj=0)
  axis(1,at=c(-6,-4,-2,0,2),labels=FALSE)
  axis(2,at=c(-6,-4,-2,0,2),labels=FALSE)
    
  #Add points
  points(x[,1], x[,2],pch=rep(1,dim(x)[1]),col=rep("gray",dim(x)[1]),cex=1)
    
  #Add ellipse for each group
  x[,3] <- droplevels(x[,3])  #Drop unused factor levels
    
  if (any(x[,3] == 1)){ #If scenario 1 was included
    points(mean(x[x[,3]==1,1]),mean(x[x[,3]==1,2]),pch=19,cex=.8,col="black")
    x <- x[x[,3] != 1,]  #Drop scenario 1 from the data
    x[,3] <- droplevels(x[,3])  #Drop unused factor levels
  }
    
  if (any(x[,3] ==4) | any(x[,3] ==5)){ #If scenario 3 and 4 are included
    x <- x[is.na(x[,1]) == F,]  #Trim NA values
  }
    
  for(i in 1:length(levels(x[,3]))){  #For each scenario...
    #Draw data ellipse
    try(dataEllipse(x[x[,3]==levels(x[,3])[i],1], x[x[,3]==levels(x[,3])[i],2],add=T,center.cex=1,levels=.95,ellipse.label="",lty=2,fill=TRUE,fill.alpha=0,grid=F,col=col[as.numeric(levels(x[,3])[i])],plot.points = F),silent=T)
  }
  #Add 1:1 line
  lines(c(-6,2),c(-6,2))
}

lab <- c("y","none","none","none","none","y","none","none","none","none","y","none","none","none","none","y","none","none","none","none","y","none","x","x","x","xy","x")
count <- 1
svg(file=paste(dir,"/raw_output/S1_Figure_S3.svg", sep=""), width=8.75, height=5.25)
par(mfrow = c(6,5),mar=c(2,1,0.25,0),mgp=c(0,0,0))

for (z in ordered){
  combine_scenario <- data.frame()
  #Ei, ci (e.g. weighted by the probability a site was wet)
  combine_scenario <- as.data.frame(cbind(c(rowMeans(site[z,1,,,4]*wet_prob + site[z,1,,,5]*(1-wet_prob)),rowMeans(site[z,2,,,4]*wet_prob + site[z,2,,,5]*(1-wet_prob)),rowMeans(site[z,3,,,4]*wet_prob + site[z,3,,,5]*(1-wet_prob)),rowMeans(site[z,4,,,4]*wet_prob + site[z,4,,,5]*(1-wet_prob)),rowMeans(site[z,5,,,4]*wet_prob + site[z,5,,,5]*(1-wet_prob))),c(rowMeans(site[z,1,,,3]*wet_prob),rowMeans(site[z,2,,,3]*wet_prob),rowMeans(site[z,3,,,3]*wet_prob),rowMeans(site[z,4,,,3]*wet_prob),rowMeans(site[z,5,,,3]*wet_prob))))
  
  ##eW, cW
  combine_scenario$scenario <- as.factor(c(rep(1,278),rep(2,278),rep(3,278),rep(4,278),rep(5,278)))
  ln_combine_scenario <- combine_scenario
  ln_combine_scenario[,1:2] <- log(ln_combine_scenario[,1:2])
  ifelse(z %in% ordered[1:10], lim <- c(-6,3), lim <- c(-7.3,3.5))
  site_ce_full(ln_combine_scenario[ln_combine_scenario[,3] == 1 | ln_combine_scenario[,3] == 2 | ln_combine_scenario[,3] == 4,],lim,lab=lab[count])
  count <- count + 1
}
dev.off()

## Figure S1_S4
lab <- c("y","none","none","none","none","y","none","none","none","none","y","none","none","none","none","y","none","none","none","none","y","none","x","x","x","xy","x")
count <- 1
svg(file=paste(dir,"/raw_output/S1_Figure_S4.svg", sep=""), width=8.75, height=5.25)
par(mfrow = c(6,5),mar=c(2,1,0.25,0),mgp=c(0,0,0))

for (z in ordered){
  combine_scenario <- data.frame()
  #Ei, ci (e.g. weighted by the probability a site was wet)
  combine_scenario <- as.data.frame(cbind(c(rowMeans(site[z,1,,,4]*wet_prob + site[z,1,,,5]*(1-wet_prob)),rowMeans(site[z,2,,,4]*wet_prob + site[z,2,,,5]*(1-wet_prob)),rowMeans(site[z,3,,,4]*wet_prob + site[z,3,,,5]*(1-wet_prob)),rowMeans(site[z,4,,,4]*wet_prob + site[z,4,,,5]*(1-wet_prob)),rowMeans(site[z,5,,,4]*wet_prob + site[z,5,,,5]*(1-wet_prob))),c(rowMeans(site[z,1,,,3]*wet_prob),rowMeans(site[z,2,,,3]*wet_prob),rowMeans(site[z,3,,,3]*wet_prob),rowMeans(site[z,4,,,3]*wet_prob),rowMeans(site[z,5,,,3]*wet_prob))))
  
  ##eW, cW
  combine_scenario$scenario <- as.factor(c(rep(1,278),rep(2,278),rep(3,278),rep(4,278),rep(5,278)))
  ln_combine_scenario <- combine_scenario
  ln_combine_scenario[,1:2] <- log(ln_combine_scenario[,1:2])
  ifelse(z %in% ordered[1:10], lim <- c(-6,3), lim <- c(-7.3,3.5))
  site_ce_full(ln_combine_scenario[ln_combine_scenario[,3] == 1 | ln_combine_scenario[,3] == 3 | ln_combine_scenario[,3] == 5,],lim,lab=lab[count])
  count <- count + 1
}
dev.off()


### Table S1.S7, Figure S1.S5. Correlation, covariance of ei and ci. ----

## Table S1.S7, Figure S1.S5. Correlation of ei and ci.
# Ei and ci variance and covariance for all covariates
#var_covar <- array(NA,c(num_spec,7),dimnames=list(spec,c("var(Ei)","var(ci)","covar(Ei,ci)","cor(Ei,ci)","t","df","p")))
var_covar <- array(NA,c(num_spec,4),dimnames=list(spec_correct,c("var(Ei)","var(ci)","covar(Ei,ci)","cor(Ei,ci)")))
for(z in 1:num_spec){
  recEi <- rowMeans(site[z,2,,,4]*wet_prob + site[z,2,,,5]*(1-wet_prob))
  recci <- rowMeans(site[z,2,,,3]*wet_prob)
  
  var_covar[z,1] <- var(recEi)
  var_covar[z,2] <- var(recci)
  var_covar[z,3] <- cov(recEi,recci)
  var_covar[z,4] <- cor(recEi,recci)
  #var_covar[z,5] <- cor.test(recEi,recci)$statistic
  #var_covar[z,6] <- cor.test(recEi,recci)$parameter
  #var_covar[z,7] <- cor.test(recEi,recci)$p.value
}
write.csv(var_covar[ordered,], paste(dir,"/raw_output/S1_Table_S7.csv",sep=""))

## Supplement S1, Figure S5. Plot of correlation coefficients for all species
spec_proper <- c("G. radiata", "D. lucidum", "D. depressissimum", "D. cimex", "D. aeruginosum", "A. marmorata", "Ph. acuta", "L. columella", "L. cubensis", "Pl. guadeloupensis", "B. straminea", "B. glabrata", "B. schrammi", "Po. glauca", "Py. coronatus", "Ma. cornuarietis", "Me. tuberculata PAP", "Me. tuberculata GOS", "Me. tuberculata FAL", "Me. tuberculata MAD", "Me. tuberculata CPF", "Me. tuberculata SEN", "T. granifera", "N. virginea","E. viridans","I .exustus", "H. duryi")
spec_proper[var_covar[,7] <= 0.05] <- paste(spec_proper[var_covar[,7] <= 0.05],"*",sep="") #add asterick symbol above significant correlation values

svg(paste(dir,"/raw_output/S1_Figure_S5.svg",sep=""),width=6,h=4.1)
par(mar=c(7,4.1,0,0))
plot(sort(var_covar[,4]),ylim=c(-1,1),type="n",axes=F,xlab="",ylab="ρ",cex.lab=1.5)
points(sort(var_covar[,4]),pch=19,col="black",cex=1.5)
axis(2,at=c(-1,-.5,0,.5,1),cex.axis=.833,las=1)
labs <- spec_correct[order(var_covar[,4])]
axis(1, at=1:27, labels = FALSE)
text(x=1.2:27.2, y=-1.2, labs, xpd=TRUE, srt=60, adj=1,cex=.833)
dev.off()


### Figure 8, Figure S1.S6, Table S1.S8. Simulation model results ------

### T1000 density plot ###
#Names of simulation models
mod <- c("V04","V04.01","V04.02", "V04.03","V04.04","V05")

#gather data from models
sim <- array(NA,dim=c(27,length(mod),3,999),dimnames=list(spec,mod,c("mean","2_5","97_5"),NULL))
sim_p1000 <- array(NA,dim=c(27,length(mod),999),dimnames=list(spec,mod,NULL))

for(k in 1:6){
  for(z in 1:27){
    #Import the results of the simulation model
    dat <- read.csv(paste(dir,"/raw_output/simulation/simulmetapop_",mod[k],"/",spec[z],"_result.csv",sep=""),header=F)
    #Get mean, 2.5% CI, and 97.5% CI
    sim[z,k,1,] <- colMeans(dat)
    sim[z,k,2,] <- apply(dat,2,function(x) quantile(x,probs=.025))
    sim[z,k,3,] <- apply(dat,2,function(x) quantile(x,probs=.975))
    
    #Get p1000 data for all 999 iterations
    sim_p1000[z,k,] <- dat[,999]
  }
}

#Change model V04.03 and V04.04 to be Pr(occupied) out of the total number of sites the species was ever observed in (not 278)
#Number of sites each species was ever observed in
site_obs <- readLines(paste(dir,"/data/simulation/parameters/sites_observed_in.txt",sep=""))
for(k in 4:5){
  for(z in 1:27){
    sim[z,k,1,] <- sim[z,k,1,]*((length(strsplit(site_obs[z]," ")[[1]]))/278)
    sim[z,k,2,] <- sim[z,k,2,]*((length(strsplit(site_obs[z]," ")[[1]]))/278)
    sim[z,k,3,] <- sim[z,k,3,]*((length(strsplit(site_obs[z]," ")[[1]]))/278)
    sim_p1000[z,k,] <- sim_p1000[z,k,]*((length(strsplit(site_obs[z]," ")[[1]]))/278)
  }
}
### Plot results of 5 scenarios
library(plotly)
library(RColorBrewer)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 5
col_use = gg_color_hue(n)
col_use <- c("black","red","blue","forestgreen","gray")

plot_list = list()
for(z in 1:27){
  unif <- rep(0,5)
  # Scenario 1: V05 (no covariates)
  scen1 <- data.frame(length = sim_p1000[z,6,])
  if(length(unique(scen1$length)) == 1){unif[1] <- 1}
  # Scenario 2: V04 (covariates, intercept only)
  scen2 <- data.frame(length = sim_p1000[z,1,])
  if(length(unique(scen2$length)) == 1){unif[2] <- 1}
  # Scenario 3: V04.01 (all covariates, all sites)
  scen3 <- data.frame(length = sim_p1000[z,2,])
  if(length(unique(scen3$length)) == 1){unif[3] <- 1}
  # Scenario 4: V04.03 (all covariates, reduced sites)
  scen4 <- data.frame(length = sim_p1000[z,4,])
  if(length(unique(scen4$length)) == 1){unif[4] <- 1}
  # Scenario 5: V04.02 (sig covariates, all sites)
  scen5 <- data.frame(length = sim_p1000[z,3,])
  if(length(unique(scen5$length)) == 1){unif[5] <- 1}
  
  #Now, combine your two dataframes into one.  First make a new column in each.
  scen1$scen <- "s1"
  scen2$scen <- "s2"
  scen3$scen <- "s3"
  scen4$scen <- "s4"
  scen5$scen <- "s5"
  scen1$unif <- unif[1]
  scen2$unif <- unif[2]
  scen3$unif <- unif[3]
  scen4$unif <- unif[4]
  scen5$unif <- unif[5]
  scen1$col <- col_use[1]
  scen2$col <- col_use[2]
  scen3$col <- col_use[3]
  scen4$col <- col_use[4]
  scen5$col <- col_use[5]
  
  #and combine into your new data frame scenLengths
  scenLengths <- rbind(scen1,scen2,scen3,scen4,scen5)
  use <- scenLengths[scenLengths$unif == 0,]
  
  #now make your lovely plot
  p <- ggplot(use, aes(length,color = scen))
  
  p <- p + geom_density(aes(x=length,y=..scaled..,color=scen),alpha=0.5,size=0.5) + 
    scale_color_manual(values=col_use[unif == 0]) + theme_bw() + 
    theme(panel.grid=element_blank(),legend.position="none",axis.title=element_blank(),plot.title=element_text(face='italic',size=7.25),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) + xlim(0,1) + ggtitle(as.character(spec_correct[z]))
  
  plot_list[[z]] <- p
}

ordered_plot_list <- plot_list[ordered]
g <- grid.arrange(grobs=ordered_plot_list,legend,ncol=7,nrow=4)

## S1, Figure S6
ggsave(file=paste(dir,"/raw_output/S1_Figure_S6.svg",sep=""), g, width=8.75, height=5.75)\

## Figure 8, subset of 8 species
plot_list = list()
count <- 1
for(z in c(6,12,17,7,16,18)){
  unif <- rep(0,5)
  # Scenario 1: V05 (no covariates)
  scen1 <- data.frame(length = sim_p1000[z,6,])
  if(length(unique(scen1$length)) == 1){unif[1] <- 1}
  # Scenario 2: V04 (covariates, intercept only)
  scen2 <- data.frame(length = sim_p1000[z,1,])
  if(length(unique(scen2$length)) == 1){unif[2] <- 1}
  # Scenario 3: V04.01 (all covariates, all sites)
  scen3 <- data.frame(length = sim_p1000[z,2,])
  if(length(unique(scen3$length)) == 1){unif[3] <- 1}
  # Scenario 4: V04.03 (all covariates, reduced sites)
  scen4 <- data.frame(length = sim_p1000[z,4,])
  if(length(unique(scen4$length)) == 1){unif[4] <- 1}
  # Scenario 5: V04.02 (sig covariates, all sites)
  scen5 <- data.frame(length = sim_p1000[z,3,])
  if(length(unique(scen5$length)) == 1){unif[5] <- 1}
  
  #Now, combine your two dataframes into one.  First make a new column in each.
  scen1$scen <- "s1"
  scen2$scen <- "s2"
  scen3$scen <- "s3"
  scen4$scen <- "s4"
  scen5$scen <- "s5"
  scen1$unif <- unif[1]
  scen2$unif <- unif[2]
  scen3$unif <- unif[3]
  scen4$unif <- unif[4]
  scen5$unif <- unif[5]
  scen1$col <- col_use[1]
  scen2$col <- col_use[2]
  scen3$col <- col_use[3]
  scen4$col <- col_use[4]
  scen5$col <- col_use[5]
  
  #and combine into your new data frame scenLengths
  scenLengths <- rbind(scen1,scen2,scen3,scen4,scen5)
  use <- scenLengths[scenLengths$unif == 0,]
  
  #now make your lovely plot
  p <- ggplot(use, aes(length,color = scen))
  
  p <- p + geom_density(aes(x=length,y=..scaled..,color=scen),alpha=0.5,size=0.75) + 
    scale_color_manual(values=col_use[unif == 0]) + theme_bw() + 
    theme(panel.grid=element_blank(),legend.position="none",axis.title=element_blank(),plot.title=element_text(face='italic',size=10),axis.text.x=element_text(size=8),axis.text.y=element_text(size=6)) + xlim(0,1) + ggtitle(as.character(spec_correct[z]))
  
  plot_list[[count]] <- p
  count <- count + 1
}

g <- grid.arrange(grobs=plot_list,ncol=3,nrow=2)
ggsave(file=paste(dir,"/raw_output/Figure_8.svg",sep=""), g, width=6, height=3.5)

## S1, Table S8: Simulation model results
#Assemble results into information for table
res <- array(NA,dim=c(27,6,3,2),dimnames=list(spec_correct,mod,c("mean","25","975"),c("p1000","textinct")))

for(z in 1:num_spec){
  for(i in 1:6){
    #p1000
    res[z,i,1,1] <- sim[z,i,1,999] #mean
    res[z,i,2,1] <- sim[z,i,2,999] #2.5%
    res[z,i,3,1] <- sim[z,i,3,999] #97.5%
    #textinct
    res[z,i,1,2] <- try(min(which(sim[z,i,1,]<0.001)),silent=T)
    res[z,i,2,2] <- try(min(which(sim[z,i,2,]<0.001)),silent=T)
    res[z,i,3,2] <- try(min(which(sim[z,i,3,]<0.001)),silent=T)
  }
}

res[is.infinite(res)==T] <- NA

# Write summary of results to a Table S9
tab <- cbind(res[ordered,6,2,1],res[ordered,6,1,1],res[ordered,6,3,1],res[ordered,6,1,2])  # Scenario 1,V05
tab <- cbind(tab,res[ordered,1,2,1],res[ordered,1,1,1],res[ordered,1,3,1],res[ordered,1,1,2])  # Scenario 2,V04
tab <- cbind(tab,res[ordered,2,2,1],res[ordered,2,1,1],res[ordered,2,3,1],res[ordered,2,1,2])  # Scenario 3,V04.01
tab <- cbind(tab,res[ordered,4,2,1],res[ordered,4,1,1],res[ordered,4,3,1],res[ordered,4,1,2])  # Scenario 4,V04.03
tab <- cbind(tab,res[ordered,3,2,1],res[ordered,3,1,1],res[ordered,3,3,1],res[ordered,3,1,2])  # Scenario 5,V04.02
tab_format <- round(tab,2)

write.csv(tab_format, file=paste(dir,"/raw_output/S1_Table_S8.csv",sep=""))

##### Time series plot
#Plot function
cols <- matrix(c(255,165,0 ,255,0,0 ,34,139,34,  0,0,255, 0,0,0),nrow=5, ncol=3,byrow=T)

mod_plot <- function(obs, mod, spe, cols){
  #background plot
  k <- c(2,1,3,4,5)
  k <- k[1:length(mod)]
  for (i in 1:length(spe)){
    plot(1:999,obs[1,1,1,],type='n',ylim=c(0,1), main=spe[i])
    for (j in k){
      polygon(c(1:999, rev(1:999), 1), c(obs[i,j,2,], rev(obs[i,j,3,]), obs[i,j,2,1]),col=rgb(cols[j,1],cols[j,2],cols[j,3],100,maxColorValue = 255), border=NA)
    }
    
    for (j in k){
      lines(1:999,obs[i,j,1,],col=rgb(cols[j,1],cols[j,2],cols[j,3],255,maxColorValue=255))
      lines(1:999,obs[i,j,2,],col=rgb(cols[j,1],cols[j,2],cols[j,3],255,maxColorValue = 255))
      lines(1:999,obs[i,j,3,],col=rgb(cols[j,1],cols[j,2],cols[j,3],255,maxColorValue = 255))
    }
  }
}

### Figure 3b, Appendix S1.Table_S9. Simulation equilibrium compared to observed results----
# Mean, standard deviation of deviations between p* and p_obs
mean(unlist(c(site_coord_pstar_sim_detect_dev[,6:32])))
sd(unlist(c(site_coord_pstar_sim_detect_dev[,6:32])))

## Appendix S1, Table S9
# Column 1.
colMeans(site_coord_pstar_sim_detect_dev[,6:32])
# Column 2.
apply(site_coord_pstar_sim_detect_dev[,6:32],2,sd)
# Column 3.
p_for_corr <- site_coord_pstar_sim_detect[,6:32]
cor_val <- array(NA,27)
for(i in 1:27){
  cor_val[i] <- cor(obs_site[,i],p_for_corr[,i])
}
cbind(colMeans(site_coord_pstar_sim_detect_dev[,6:32]),apply(site_coord_pstar_sim_detect_dev[,6:32],2,sd),cor_val)

#gather data
tab3 <- cbind(colMeans(site_coord_pstar_sim_detect_dev[,6:32])[ordered],apply(site_coord_pstar_sim_detect_dev[,6:32],2,sd)[ordered],cor_val[ordered])
write.csv(tab3,paste(dir,"/raw_output/S1_Table_S9.csv",sep=""))

## Figure 3b
# x-axis: p* = [(c_i x p*) /  (c_i x p* + e_i)]*Pr(detection)
# y-axis: observed detection frequency (averaged across years)
# x-axis - model-estimated p* (corrected for detection)
p1000_detect <- colMeans(site_coord_pstar_sim_detect[,6:32])

gdat <- as.data.frame(cbind(p1000_detect,occ))
rownames(gdat) <- spec_correct
rownames(gdat)[17:22] <- substr(spec_correct[17:22],17,25)
gdat$inv <- rep("N",27)
gdat$inv[c(11,7,8,17,16,18,23,19,20,21,26,27,22)] <- "I"
group_col <- c("I" = "black", "N" = "white")
dt.triangle <- data.table(group = c(1,1,1), polygon.x = c(0,0.6,0.6), polygon.y = c(0,0,0.6))

# Scatterplot colored by native or invasive
p <- ggplot(data=gdat,aes(p1000_detect,occ,fill=inv)) + geom_point(shape = 21,colour = "black",size=1) + scale_fill_manual(values=group_col) + theme_bw()+ theme(legend.position="none",axis.text=element_text(size=7),axis.title=element_blank()) + geom_segment(aes(x = 0, y = 0, xend = 0.6, yend = 0.6)) + geom_polygon(data=dt.triangle,aes(x=polygon.x,y=polygon.y),inherit.aes=FALSE,alpha=0.1, color='black')
#p <- ggplot(data=gdat,aes(lev_detect_01,occ,fill=inv)) + geom_point(shape = 21,colour = "black",size=1) + scale_fill_manual(values=group_col) + theme_bw()+ theme(legend.position="none",axis.title=element_text(size=7,margin=margin(0,0,0,0)),axis.text=element_text(size=7),plot.title = element_text(size=7,margin=margin(0,0,0,0))) + xlab("expected detection frequency") + ylab("observed detection frequency") + geom_segment(aes(x = 0, y = 0, xend = 0.6, yend = 0.6)) + geom_polygon(data=dt.triangle,aes(x=polygon.x,y=polygon.y),inherit.aes=FALSE,alpha=0.1, color='black') + ggtitle("b) without covariates")
# Add labels for closely-spaced values
p <- p + geom_text_repel(data=subset(gdat,p1000_detect < 0.01),aes(label=rownames(gdat)[gdat$p1000_detect < 0.01]),max.overlaps=Inf,direction="y",hjust=1,nudge_x=-0.05,xlim = c(-Inf, Inf),min.segment.length = 0,size=2.2,segment.size=0.2) + scale_x_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6),limits=c(-.1,.6),expand = expansion(mult = 0.05),labels=c(0,.1,.2,.3,.4,.5,.6)) + scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6),limits=c(0,0.61),labels=c(0,.1,.2,.3,.4,.5,.6))
# Add labels for well-spaced values
p <- p + geom_text_repel(data=subset(gdat,p1000_detect > 0.01),aes(label=rownames(gdat)[gdat$p1000_detect > 0.01]),max.overlaps=Inf,size=2.2)
# Save to external file
svg(file=paste(dir,"/raw_output/Figure_3/Figure_3b.svg", sep=""), width=6, height=3)
par(mar=c(0,0,0,0),mgp=c(0,0,0))
p
dev.off()