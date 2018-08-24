# Course: COMPG099 - Dissertation
# Institution: University College London
# Developer: Antonio Remiro Azocar

rm(list=ls())

setwd("~/Desktop/PhD/COMPG099/")  # set working directory here
source("functions.R")
source("main.R")
library("survHE")
library("hydroGOF")
results.dir <- "Results/Bladder2/"

load("Data_real/bladder.rda")
survival.data <- bladderTP[,c(1:4)]

file.identifier <- "bladder_cancer_study"

risk.every <- 6
fit.distr <- "exponential"
extrapolate.to <- 60
admin.cutoff <- 60
PSA.simulations <- 1000

seed <- 55
set.seed(seed)

lower.rho <- 0
upper.rho <- 1
no.rho.estimates <- 101
q.allowed <- rbind(c(1,1,1),c(0,1,1),c(0,0,0))

colnames(survival.data) <- c("tPFS","ePFS","tOS","eOS")

results <- partSA.trans(survival.data, risk.every, fit.distr, extrapolate.to,
                        admin.cutoff, PSA.simulations, seed, q.allowed,
                        results.dir, file.identifier)

save(results, file=paste(c(results.dir, "final_results/", file.identifier, ".rda"), collapse=""))

final.results.table <- data.frame(mean.rho=double(),
                                 stdev.rho=double(),
                                 true.lambda.11=double(),
                                 mae.11=double(),
                                 true.lambda.12=double(),
                                 mae.12=double(),
                                 true.lambda.13=double(),
                                 mae.13=double(),
                                 true.lambda.22=double(),
                                 mae.22=double(),
                                 true.lambda.23=double(),
                                 mae.23=double(),
                                 mae.true.pi1=double(),
                                 mae.partsa.pi1=double(),
                                 mae.true.pi2=double(),
                                 mae.partsa.pi2=double(),
                                 mae.true.pi3=double(),
                                 mae.partsa.pi3=double())

load(paste(c(results.dir, "final_results/", file.identifier, ".rda"), collapse=""))
plot_transitions(results$avg.lambda.ed[,1,1],
                 results$avg.lambda.lm[,1,1],
                 results$true.lambda.homog.avg[1,1],
                 results$true.lambda.homog.inf[1,1],
                 results$true.lambda.homog.sup[1,1],
                 y.lim=c(0.9,1),
                 title.name= expression(paste(lambda[11]^(t)," (PFS->PFS)")),
                 legend.position="bottomright",
                 extrapolate.to, admin.cutoff)
plot_transitions(results$avg.lambda.ed[,1,2],
                 results$avg.lambda.lm[,1,2],
                 results$true.lambda.homog.avg[1,2],
                 results$true.lambda.homog.inf[1,2],
                 results$true.lambda.homog.sup[1,2],
                 y.lim=c(0,0.1),
                 title.name= expression(paste(lambda[12]^(t)," (PFS->PD)")),
                 legend.position="topright",
                 extrapolate.to, admin.cutoff)
plot_transitions(results$avg.lambda.ed[,1,3],
                 results$avg.lambda.lm[,1,3],
                 results$true.lambda.homog.avg[1,3],
                 results$true.lambda.homog.inf[1,3],
                 results$true.lambda.homog.sup[1,3],
                 y.lim=c(0,0.1),
                 title.name= expression(paste(lambda[13]^(t)," (PFS->D)")),
                 legend.position="topright",
                 extrapolate.to, admin.cutoff)
plot_transitions(results$avg.lambda.ed[,2,2],
                 results$avg.lambda.lm[,2,2],
                 results$true.lambda.homog.avg[2,2],
                 results$true.lambda.homog.inf[2,2],
                 results$true.lambda.homog.sup[2,2],
                 y.lim=c(0.9,1),
                 title.name= expression(paste(lambda[22]^(t)," (PD->PD)")),
                 legend.position="bottomright",
                 extrapolate.to, admin.cutoff)
plot_transitions(results$avg.lambda.ed[,2,3],
                 results$avg.lambda.lm[,2,3],
                 results$true.lambda.homog.avg[2,3],
                 results$true.lambda.homog.inf[2,3],
                 results$true.lambda.homog.sup[2,3],
                 y.lim=c(0,0.1),
                 title.name= expression(paste(lambda[23]^(t)," (PD->D)")),
                 legend.position="topright",
                 extrapolate.to, admin.cutoff)
plot_proportions(est.props.ed=results$opt.MM.ed$est.props$mest[,1],
                 est.props.lm=results$opt.MM.ed$est.props$mest[,1],
                 props.partsa=results$opt.MM.ed$est.props$m[,1],
                 true.props.msm.avg=results$true.props.homog.avg[,1],
                 true.props.msm.inf=results$true.props.homog.inf[,1],
                 true.props.msm.sup=results$true.props.homog.sup[,1],
                 title.name=expression(paste(pi[1]^(t)," (PFS)")),
                 legend.position="topright",
                 extrapolate.to=extrapolate.to, 
                 admin.cutoff=admin.cutoff)
plot_proportions(est.props.ed=results$opt.MM.ed$est.props$mest[,2],
                 est.props.lm=results$opt.MM.ed$est.props$mest[,2],
                 props.partsa=results$opt.MM.ed$est.props$m[,2],
                 true.props.msm.avg=results$true.props.homog.avg[,2],
                 true.props.msm.inf=results$true.props.homog.inf[,2],
                 true.props.msm.sup=results$true.props.homog.sup[,2],
                 title.name=expression(paste(pi[2]^(t)," (PD)")),
                 legend.position="topright",
                 extrapolate.to=extrapolate.to, 
                 admin.cutoff=admin.cutoff)
plot_proportions(est.props.ed=results$opt.MM.ed$est.props$mest[,3],
                 est.props.lm=results$opt.MM.ed$est.props$mest[,3],
                 props.partsa=results$opt.MM.ed$est.props$m[,3],
                 true.props.msm.avg=results$true.props.homog.avg[,3],
                 true.props.msm.inf=results$true.props.homog.inf[,3],
                 true.props.msm.sup=results$true.props.homog.sup[,3],
                 title.name=expression(paste(pi[3]^(t)," (D)")),
                 legend.position="bottomright",
                 extrapolate.to=extrapolate.to, 
                 admin.cutoff=admin.cutoff)
mean.abs.error.PFS.true <- mae(results$opt.MM.ed$est.props$mest[,1],
                               results$true.props.homog.avg[,1])
mean.abs.error.PD.true <- mae(results$opt.MM.ed$est.props$mest[,2],
                              results$true.props.homog.avg[,2])
mean.abs.error.D.true <- mae(results$opt.MM.ed$est.props$mest[,3],
                             results$true.props.homog.avg[,3])
mean.abs.error.PFS.partsa <- mae(results$opt.MM.ed$est.props$mest[,1],
                                 results$opt.MM.ed$est.props$m[,1])
mean.abs.error.PD.partsa <- mae(results$opt.MM.ed$est.props$mest[,2],
                                results$opt.MM.ed$est.props$m[,2])
mean.abs.error.D.partsa <- mae(results$opt.MM.ed$est.props$mest[,3],
                               results$opt.MM.ed$est.props$m[,3])
mean.abs.error.trans11 <- mae(results$avg.lambda.ed[,1,1],
                              rep(results$true.lambda.homog.avg[1,1],extrapolate.to-1))
mean.abs.error.trans12 <- mae(results$avg.lambda.ed[,1,2],
                              rep(results$true.lambda.homog.avg[1,2],extrapolate.to-1))
mean.abs.error.trans13 <- mae(results$avg.lambda.ed[,1,3],
                              rep(results$true.lambda.homog.avg[1,3],extrapolate.to-1))
mean.abs.error.trans22 <- mae(results$avg.lambda.ed[,2,2],
                              rep(results$true.lambda.homog.avg[2,2],extrapolate.to-1))
mean.abs.error.trans23 <- mae(results$avg.lambda.ed[,2,3],
                              rep(results$true.lambda.homog.avg[2,3],extrapolate.to-1))
me.trans11 <- me(results$avg.lambda.ed[,1,1],rep(results$true.lambda.homog.avg[1,1],extrapolate.to-1))
me.trans12 <- me(results$avg.lambda.ed[,1,2],rep(results$true.lambda.homog.avg[1,2],extrapolate.to-1))
me.trans13 <- me(results$avg.lambda.ed[,1,3],rep(results$true.lambda.homog.avg[1,3],extrapolate.to-1))
me.trans22 <- me(results$avg.lambda.ed[,2,2],rep(results$true.lambda.homog.avg[2,2],extrapolate.to-1))
me.trans23 <- me(results$avg.lambda.ed[,2,3],rep(results$true.lambda.homog.avg[2,3],extrapolate.to-1))
relative.mad.trans11 <- mean.abs.error.trans11/results$true.lambda.homog.avg[1,1]
relative.mad.trans12 <- mean.abs.error.trans12/results$true.lambda.homog.avg[1,2]
relative.mad.trans13 <- mean.abs.error.trans13/results$true.lambda.homog.avg[1,3]
relative.mad.trans22 <- mean.abs.error.trans22/results$true.lambda.homog.avg[2,2]
relative.mad.trans23 <- mean.abs.error.trans23/results$true.lambda.homog.avg[2,3]
load(paste(c(results.dir, "optimal_rho_data/ED_", file.identifier, ".rda"), collapse=""))


final.results.table[1, "mean.rho"]  <- round(mean(optimal.rho.ed), digits=4)
final.results.table[1, "stdev.rho"] <- round(sd(optimal.rho.ed), digits=4)
final.results.table[1, "true.lambda.11"] <- round(results$true.lambda.homog.avg[1,1],4)        
final.results.table[1, "mae.11"] <- round(mean.abs.error.trans11,4)        
final.results.table[1, "true.lambda.12"] <- round(results$true.lambda.homog.avg[1,2],4)        
final.results.table[1, "mae.12"] <- round(mean.abs.error.trans12,4)        
final.results.table[1, "true.lambda.13"] <- round(results$true.lambda.homog.avg[1,3],4)        
final.results.table[1, "mae.13"] <- round(mean.abs.error.trans13,4)        
final.results.table[1, "true.lambda.22"] <- round(results$true.lambda.homog.avg[2,2],4)        
final.results.table[1, "mae.22"] <- round(mean.abs.error.trans22,4)        
final.results.table[1, "true.lambda.23"] <- round(results$true.lambda.homog.avg[2,3],4)        
final.results.table[1, "mae.23"] <- round(mean.abs.error.trans23,4)        
final.results.table[1, "mae.true.pi1"] <- round(mean.abs.error.PFS.true,4)
final.results.table[1, "mae.partsa.pi1"] <- round(mean.abs.error.PFS.partsa*10000, 4)
final.results.table[1, "mae.true.pi2"] <- round(mean.abs.error.PD.true,4)
final.results.table[1, "mae.partsa.pi2"] <- round(mean.abs.error.PD.partsa*10000, 4)
final.results.table[1, "mae.true.pi3"] <- round(mean.abs.error.D.true,4)
final.results.table[1, "mae.partsa.pi3"] <- round(mean.abs.error.D.partsa*10000, 4)

