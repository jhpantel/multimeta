#######################
# data loading        #
#######################
# matrix of site state
# 0 = non visited site or missing data
# 1 = site in state W (Wet)
# 2 = site in state D (Dry)
dry <- matrix(scan("dry_final.txt", n = 278*15,skip=1), 278, 15, byrow = TRUE)

# state matrix
# 0 = non visited site or missing data
# 1 = species not detected (i.e., either dry or wet site)
# 2 = species detected (i.e., necessarily a wet site)
S <- matrix(scan("state.txt", n = 278*15,skip=1), 278, 15, byrow = TRUE)

# matrixes of visits and repeated visits within the same occasion 
# correspond  to two 280*15 matrices (sites by year)
# first eleven rows correspond to visits
# next eleven rows correspond to repeated visits within the same occasion
presence <- matrix(scan("presence.txt", n = 278*30,skip=1), 278, 30, byrow = TRUE)

y <- structure(presence,.Dim = c(278L,15L,2L))

####Read in covariates####
# covariates site-specific
cov <- read.table("cov.txt",h=T)

size = cov$size
size <- (size-mean(size))/sd(size)
connec = cov$connec
connec <- (connec-mean(connec))/sd(connec)
veg = cov$veg
veg <- (veg-mean(veg))/sd(veg)
stability = cov$stability
stability <- (stability-mean(stability))/sd(stability)

habit <- read.table("habit.txt",h=T)

mangrove = habit$mangrove
mangrove <- (mangrove-mean(mangrove))/sd(mangrove)
river = habit$river
river <- (river-mean(river))/sd(river)

# covariates year-specific
# little rainy season (lrs): cumulative rainfall during the little rainy season from March 1 to May 31, 92 days, e.g. during the period of putative extinction
# rainy season (rs): cumulative rainfall during the rainy season from (i) July 1 - Dec. 31 of year preceding sampling campaign (184 days, lrs_July_Dec) and (ii) 60 days preceding the first date of the sampling campaign (60 days, rs_60_days)

pluvio <- read.table("pluvio.txt",h=T)
lrs=pluvio$lrs
rs_1=pluvio$rs_July_Dec
#rs_2=pluvio$rs_60_days

lrs <- (lrs-mean(lrs))/sd(lrs)
rs_1 <- (rs_1-mean(rs_1))/sd(rs_1)
#rs_2 <- (rs_2-mean(rs_2))/sd(rs_2)


# covariates year and site-specific
colsource <- matrix(scan("colsource.txt", n = 278*15), 278, 15, byrow = TRUE)
colsource <- (colsource-mean(colsource))/sqrt(mean((colsource-mean(colsource))*(colsource-mean(colsource))))

# number of sites
nsite<-dim(y)[1]
# number de repetitions
nrep<-2
# number of visits to sites
nyear<-15

############################################################
# compute the year a site is first visited                 #
############################################################
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

##############################
# compute X(t)               #
##############################
# X(t) = 1 if a site is in state D
# X(t) = 0 if a site is in state W or not visited
#X <- matrix(0,nsite,ncol(S))
#X[dry==2]=1
X <- matrix(0,nsite,ncol(S))
X[dry==2]=1
#Include state uncertanity for some sites, give Pr(dry)
X[dry > 0 & dry < 1] <- 1 - dry[dry > 0 & dry < 1]

Zst<- S
# Winbugs g?re pas les 0, il faut les remplacer par des NA dans Zst et dans les donn?es
Zst[(Zst == 0)] <- NA
y[y==0]=NA

# donn?es 
mydatax <- list (y=y,nsite=nsite,nrep=nrep,nyear=nyear,e=e,X=X,fin=fin,size=size,veg=veg,stability=stability,connec=connec,colsource=colsource,lrs=lrs,rs_1=rs_1,mangrove=mangrove,river=river)

# nb d'it?rations
ni=20000
# nb burn-in
nb=10000
# nb thin
nt=1
# nb de chaines ? lancer en //
nc=3
# charger jags depuis R
library(rjags)
library(coda)


##############################
# run model FULL             #
##############################

# valeurs initiales
##Add at least one presence each year for initial Zst
GT_list <- apply(Zst[1:250,],2,function(x) max(unique(x),na.rm=T)==2)
MG_list <- apply(Zst[251:278,],2,function(x) max(unique(x),na.rm=T)==2)
list <- cbind(GT_list,MG_list)

for (i in 1:2){ #For each island...
  if (length(which(list[,i]==F)) > 0){ #If the species is missing in any year
    for (q in 1:sum(list[,i]==F)){ #For each year without that species present...
      if (i==1){Zst[sample(1:250,1),which(list[,i]==F)[q]] <- 2}
      else {Zst[sample(251:278,1),which(list[,i]==F)[q]] <- 2}
    }
  }
}

init1 <- list(pM=.1,psiGT = 0.1,psiMG = 0.1,beta=rep(-.2,14),alpha=c(-.2,-.2,-.2,NA),z=Zst)
init2 <- list(pM=.5,psiGT = 0.5,psiMG = 0.5,beta=rep(0,14),alpha=c(0,0,0,NA),z=Zst)
init3 <- list(pM=.9,psiGT = 0.9,psiMG = 0.9,beta=rep(.2,14),alpha=c(.2,.2,.2,NA),z=Zst)
inits <- list (init1,init2,init3)

parameters <- c("pM","beta","gam","alpha","z")

start2<-as.POSIXlt(Sys.time())
jmodel.full <- jags.model("model_2_cov.txt", mydatax, inits, n.chains = nc,n.adapt = nb)
jsample.full <- coda.samples(jmodel.full, parameters, n.iter=ni, thin = nt)
end2 <-as.POSIXlt(Sys.time())
duration.full = end2-start2
write.table(duration.full,"TempsCalculModeleFULL.txt")

# Sauvegarder les r?sultats
save(jsample.full, file='jsample.Rdata')

##############################
# # Plot des r?sultats         #
# ##############################
# # check convergence et sauvegarde dans un fichier pdf
# pdf("Graphic_trace1_modelfull.pdf")
# plot(jsample.full, trace = TRUE, density = FALSE,ask = dev.interactive())
# dev.off()
# pdf("Graphic_trace2_modelfull.pdf")
# xyplot(jsample.full)
# dev.off()
# 
# #gelman.diag.sans <- gelman.diag(jsample.sans)
# #write.csv(c(gelman.diag.sans[1],gelman.diag.sans[2]),"gelman.diag.sans.csv")
# 
# #gelman.diag.selec <- gelman.diag(jsample.selec)
# #write.csv(c(gelman.diag.selec[1],gelman.diag.selec[2]),"gelman.diag.selec.csv")
# 
# # numerical summaries and posterior distributions
# pdf("Graphic_posteriordistribution1_full.pdf")
# plot(jsample.full, trace = FALSE, density = TRUE,ask = dev.interactive()) 
# dev.off()
# pdf("Graphic_posteriordistribution2_full.pdf")
# densityplot(jsample.full,ylab="Density")
# dev.off()
# 
# # tableau des summary
# m.full <- summary(jsample.full)
# mu_sd=m.full$stat[,1:2]  ## make columns for mean and sd
# q=m.full$quantile[,c(3,1,5)] ## make columns for median and CI
# table.full=cbind(mu_sd,q)
# write.csv(table.full,file='table.full.csv')
# 
# ### Graph des proba ? post?riori combin?es des diff?rents mod?les
# df2=as.data.frame(rbind(jsample.full[[1]],jsample.full[[2]],jsample.full[[3]]))
#  
# dim(df2)[1]
# nr= 100        ### Pour s?lectionner le nombre d'it?ration ? garder
# 
# phiD <- df2$"beta[3]"[(dim(df2)[1]-nr):dim(df2)[1]]
# phiW <- df2$"beta[1]"[(dim(df2)[1]-nr):dim(df2)[1]]
# prior = runif(length(phiD))
# 
# pdf("Graphic_rates.pdf")
# par(mfrow=c(2,2))
# plot(density(phiD), xlab='Persistence rate',ylab='Posterior probability', xaxt="n",yaxt="n", xlim = c(0, 1), ylim = c(0, 24), col=rgb(15, 5, 107, maxColorValue = 255),lwd=2,lty="solid", main="")
# par(new=T)
# plot(density(phiW), xlab='',ylab='', xaxt="n",yaxt="n", xlim = c(0, 1), ylim = c(0, 24), col=rgb(38, 196, 236, maxColorValue = 255),lwd=2,lty="solid",main="")
# par(new=T)
# plot(density(prior), xlab='',ylab='', xaxt="n",yaxt="n", xlim = c(0, 1), ylim = c(0, 24), col=rgb(96, 96, 96, maxColorValue = 255),lwd=2, lty="dashed",main="")
# 
# axis(1,labels=c(0,0.2,0.4,0.6,0.8,1),at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.01)
# axis(2,labels=c(0,6,12,18,24),at=c(0,6,12,18,24), las=1,tck=-0.01)
# 
# colD <- df2$"beta[4]"[(dim(df2)[1]-nr):dim(df2)[1]]
# colW <- df2$"beta[2]"[(dim(df2)[1]-nr):dim(df2)[1]]
# prior = runif(length(phiD), min=0, max=10)
# 
# plot(density(colD), xlab='Colonization rate',ylab='Posterior probability', xaxt="n",yaxt="n", xlim = c(0, 10), ylim = c(0, 24), col=rgb(15, 5, 107, maxColorValue = 255),lwd=2,lty="solid", main="")
# par(new=T)
# plot(density(colW), xlab='',ylab='', xaxt="n",yaxt="n", xlim = c(0, 1), ylim = c(0, 24), col=rgb(38, 196, 236, maxColorValue = 255),lwd=2,lty="solid",main="")
# par(new=T)
# plot(density(prior), xlab='',ylab='', xaxt="n",yaxt="n", xlim = c(0, 1), ylim = c(0, 24), col=rgb(96, 96, 96, maxColorValue = 255),lwd=2, lty="dashed",main="")
# 
# axis(1,labels=c(0,.2,.4,.6,.8,1),at=c(0,.2,.4,.6,.8,1), tck=-0.01)
# axis(2,labels=c(0,6,12,18,24),at=c(0,6,12,18,24), las=1,tck=-0.01)
# dev.off()

