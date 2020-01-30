#setwd("F:\\R_bubeng_seedling_CNDD drought\\bayes")
library(rstan)
library(StanHeaders)
library(magrittr)
library(tibble)
#install.packages("StanHeaders")
#load("bayes_analysis_PC1_20200105")

options(mc.cores = parallel::detectCores())

seedling=read.csv("seedling_for_drought.csv")
wp=read.csv("wp.csv")

colnames(seedling)
colnames(wp)


#选取和water potential相同的物种
seedling_wp=seedling[seedling$sp%in%wp$Species,]
shared_wp=wp[wp$Species%in%seedling$sp,]
nrow(seedling)
length(unique(seedling_wp$sp))
nrow(shared_wp)

############################### year level
X <- cbind(rep(1, nrow(seedling_wp)),
           seedling_wp[,c("S_scon", "S_acon","S_ahet","S_shet","S_soilpc1","S_soilpc2","S_soilpc3","S_log_h1")])

colnames(X)[1] <- "Int"

dim(X)
U <- cbind(rep(1,length(unique(seedling_wp$sp))),
           shared_wp$avg_spwp)

U[, 2] <- scale(U[,2])

dim(U)

n_sp <- length(unique(seedling_wp$sp))
n_para <- ncol(X)
n_plot <- length(unique(seedling_wp$quadrat))
n_census <- length(unique(seedling_wp$census))
list_dat <- list(N = nrow(seedling_wp),
                 J = n_sp,
                 K = n_para,
                 S = n_plot,
                 T = n_census,
                 L = ncol(U),
                 suv = seedling_wp$surv,
                 plot = seedling_wp$quadrat %>%
                   as.character %>% as.factor %>% as.integer,
                 census = seedling_wp$census %>%
                   as.character %>% as.factor %>% as.integer,
                 sp = seedling_wp$sp %>%
                   as.character %>% as.factor %>% as.integer,
                 x = X %>% as.matrix,
                # x = t(X),
                 u = U)
str(list_dat)


fit3 <- stan(file = "model.stan",
            data = list_dat,
            verbose = TRUE,
            iter = 2000,
            warmup = 1000,
            thin = 1,
            chains =  4,
            refresh = 200,
            control = list(adapt_delta = 0.9, max_treedepth = 20))
 print(fit3, pars = c("gamma", "sigma_phi", "sigma_tau"))





###########dry season
seedling_wpd=seedling_wp[seedling_wp$season=="dry",]
X <- cbind(rep(1, nrow(seedling_wpd)),
           seedling_wpd[,c("S_scon", "S_acon","S_ahet","S_shet","S_soilpc1","S_soilpc2","S_soilpc3","S_log_h1")])

colnames(X)[1] <- "Int"
shared_wpd=shared_wp[shared_wp$Species%in%seedling_wpd$sp,]
dim(X)
U <- cbind(rep(1,length(unique(seedling_wpd$sp))),
           shared_wpd$avg_spwp)

U[, 2] <- scale(U[,2])

dim(U)
n_sp <- length(unique(seedling_wpd$sp))
n_para <- ncol(X)
n_plot <- length(unique(seedling_wpd$quadrat))
n_census <- length(unique(seedling_wpd$census))
list_dat <- list(N = nrow(seedling_wpd),
                 J = n_sp,
                 K = n_para,
                 S = n_plot,
                 T = n_census,
                 L = ncol(U),
                 suv = seedling_wpd$surv,
                 plot = seedling_wpd$quadrat %>%
                   as.character %>% as.factor %>% as.integer,
                 census = seedling_wpd$census %>%
                   as.character %>% as.factor %>% as.integer,
                 sp = seedling_wpd$sp %>%
                   as.character %>% as.factor %>% as.integer,
                 x = X %>% as.matrix,
                # x = t(X),
                 u = U)
str(list_dat)

#检验U的物种顺序是否和X的物种顺序一致
unique(cbind(sp,as.character(seedling_wpd$sp)))


fit <- stan(file = "model.stan",
            data = list_dat,
            verbose = TRUE,
            iter = 3000,
            warmup = 2000,
            thin = 1,
            chains =  4,
            refresh = 200,
            control = list(adapt_delta = 0.99, max_treedepth = 20))
 print(fit, pars = c("gamma", "sigma_phi", "sigma_tau"))

################ rainy season

seedling_wpr=seedling_wp[seedling_wp$season=="rainy",]
X <- cbind(rep(1, nrow(seedling_wpr)),
           seedling_wpr[,c("S_scon", "S_acon","S_ahet","S_shet","S_soilpc1","S_soilpc2","S_soilpc3","S_log_h1")])

colnames(X)[1] <- "Int"

shared_wpr=shared_wp[shared_wp$Species%in%seedling_wpr$sp,]

dim(X)
U <- cbind(rep(1,length(unique(seedling_wpr$sp))),
           shared_wpr$avg_spwp)

U[, 2] <- scale(U[,2])

dim(U)
n_sp <- length(unique(seedling_wpr$sp))
n_para <- ncol(X)
n_plot <- length(unique(seedling_wpr$quadrat))
n_census <- length(unique(seedling_wpr$census))
list_dat <- list(N = nrow(seedling_wpr),
                 J = n_sp,
                 K = n_para,
                 S = n_plot,
                 T = n_census,
                 L = ncol(U),
                 suv = seedling_wpr$surv,
                 plot = seedling_wpr$quadrat %>%
                   as.character %>% as.factor %>% as.integer,
                 census = seedling_wpr$census %>%
                   as.character %>% as.factor %>% as.integer,
                 sp = seedling_wpr$sp %>%
                   as.character %>% as.factor %>% as.integer,
                 x = X %>% as.matrix,
                # x = t(X),
                 u = U)
str(list_dat)


fit1 <- stan(file = "model.stan",
            data = list_dat,
            verbose = TRUE,
            iter = 3000,
            warmup = 2000,
            thin = 1,
            chains =  4,
            refresh = 200,
            control = list(adapt_delta = 0.99, max_treedepth = 20))
 print(fit1, pars = c("gamma", "sigma_phi", "sigma_tau"))



#fit  Osmometer in dry season
#fit1 Osmometer in rainy season
#fit3 year

############ 利用fit 和fit1进行计算
 print(fit, pars = c("gamma", "sigma_phi", "sigma_tau"))
 print(fit1, pars = c("gamma", "sigma_phi", "sigma_tau"))

 print(fit, pars = c("beta", "sigma_phi", "sigma_tau"))


###
#dry or rainy 每个参数的中文名
num=seedling_wpd$sp %>%as.character %>% as.factor %>% as.integer
sp=seedling_wpd$sp %>%as.character
spnum=cbind(sp,num)
spnum2=spnum[!duplicated(spnum),]
spnum3=spnum2[order(as.numeric(spnum2[,2])),]
name_dry=spnum3[,1]
#name_rainy=spnum3[,1]


 test1 <- rstan::extract(fit, "beta")
pr_dry=apply(test1$beta, c(2,3), mean)
 test2 <- rstan::extract(fit1, "beta")
pr_rainy=apply(test2$beta, c(2,3), mean)

pr_dry=as.data.frame(pr_dry)
colnames(pr_dry)=c("int","S_scon", "S_acon","S_ahet","S_shet","S_soilpc1","S_soilpc2","S_soilpc3","S_log_h1")
pr_rainy=as.data.frame(pr_rainy)
colnames(pr_rainy)=c("int","S_scon", "S_acon","S_ahet","S_shet","S_soilpc1","S_soilpc2","S_soilpc3","S_log_h1")

pr_dry$sp=as.character(name_dry)
pr_rainy$sp=as.character(name_rainy)


#test seasonal difference in conspeific adult
par(mfcol=c(2,2))
hist(pr_dry[,3],xlim=c(-3,2),xlab="CNDD of adult neighbour",ylab="number of species",main="dry season",col="red")
hist(pr_rainy[,3],xlim=c(-3,2),xlab="CNDD of adult neighbour",ylab="number of species",main="rainy season",col="blue")
t.test(pr_dry[,3],pr_rainy[,3])
#test seasonal difference in conspecific seedling

hist(pr_dry[,2],xlim=c(-4,2),xlab="CNDD of seedling neighbour",ylab="number of species",main="dry season",col="red")
hist(pr_rainy[,2],xlim=c(-4,2),xlab="CNDD of seedling neighbour",ylab="number of species",main="rainy season",col="blue")
t.test(pr_dry[,2],pr_rainy[,2])
#为啥季节又显著了呢？这是个问题，这个不是群落水平，而是物种水平，更多的物种受到了更严重的CNDD，需不需要和abundance联系呢？

#把wp和pr结合在一起。
res1=merge(pr_dry,wp,by.x="sp",by.y="Species",all.x=T)
pr_dry$sp%in%wp$Species
res2=merge(pr_rainy,wp,by.x="sp",by.y="Species",all.x=T)
colnames(res1)
par(mfcol=c(2,2))
#dry season
plot(res1$avg_spwp,res1$S_acon,xlab="osmotic pressure",pch=16,ylab="CNDD of adult neighbour",main="dry season",col="red")
lm_fit=lm(avg_spwp~S_acon,data=res1[res1$avg_spwp<4000,])
summary(lm_fit)

#rainy season
plot(res2$avg_spwp,res2$S_acon,xlab="osmotic pressure",pch=16,ylab="CNDD of adult neighbour",main="rainy season",col="blue")
summary(lm(res2$avg_spwp~res2$S_acon))
plot(res1$avg_spwp,res1$S_scon,xlab="osmotic pressure",ylab="CNDD of seedling neighbour",main="dry season",col="red")
summary(lm(res1$avg_spwp~res1$S_scon))
plot(res2$avg_spwp,res2$S_scon,xlab="osmotic pressure",ylab="CNDD of seedling neighbour",main="rainy season",col="blue")
summary(lm(res2$avg_spwp~res2$S_scon))

# how about Intercept? CNDD and Osmometer?
par(mfcol=c(2,2))
plot(res1$avg_spwp,res1$int,xlab="osmotic pressure",pch=16,ylab="Intercept (reference suvival rate)",main="dry season",col="red")
plot(res1$avg_spwp,res1$S_acon,xlab="osmotic pressure",pch=16,ylab="CNDD of adult neighbour",main="dry season",col="red")
plot(res1$int,res1$S_acon,xlab="Intercept (reference survival rate)",pch=16,ylab="CNDD of adult neighbour",main="dry season",col="red")
summary(lm(res1$int~res1$S_acon))
summary(lm(res1$avg_spwp~res1$int))
summary(lm(res1$avg_spwp~res1$S_acon))

plot(res1$int,res1$S_scon)
summary(lm(res1$int~res1$S_scon))


#save.image("")
