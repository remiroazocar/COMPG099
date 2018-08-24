# Course: COMPG099 - Dissertation
# Institution: University College London
# Developer: Antonio Remiro Azocar

rm(list=ls())

setwd("~/Desktop/PhD/COMPG099/")  # set working directory here
source("functions.R")
source("main.R")
library("survHE")
library("hydroGOF")
results.dir <- "Results/Simulated1/"

seed <- 55
set.seed(seed)

admin.cutoff <- c(150,120)
prop <- c(0.5,0.65,0.8)
phi <- c(0,0.5,1)
no_subjects <- c(200,500,1000)
cens.par <- c(4,3)

int.censoring=TRUE
cens.model = "uniform"
theta1 <- 0.4
theta2 <- 1
theta <- 1
rescale.factor <- 120 # arbitrary rescaling factor
admin.censoring <- TRUE
risk.every <- 10
fit.distr <- "gengamma"
extrapolate.to <- c(500,500)
PSA.simulations <- 1000
lower.rho <- 0.01
upper.rho <- 0.99
no.rho.estimates <- 99
q.allowed <- rbind(c(1,1,1),c(0,1,1),c(0,0,0))
count <- 1

print("Simulation study. Data generated from FGM copula.")


# for(a in 1:length(no_subjects)) {
# for(b in 1:length(phi)) {        
# for(c in 1:length(prop)) { 
# for(d in 1:length(cens.par)) {
#         print(paste(c("Study ", count, " of ", length(no_subjects)*length(phi)*length(prop)*length(cens.par)),
#                     collapse=""))
#         print(paste(c("Number of subjects: ", no_subjects[a]), collapse=""))
#         print(paste(c("phi (association parameter) = ", phi[b]), collapse=""))
#         print(paste(c("Proportion of patients passing through transient state: ", 
#                       prop[c]),collapse=""))
#         print(paste(c("Upper bound of uniform intermittent censoring distribution: ",
#                       cens.par[d]), collapse=""))
#         file.identifier <- paste(c("N_", no_subjects[a], "_prop_", prop[c],
#                                  "_phi_", phi[b], "_censpar_", cens.par[d]),
#                                  collapse="")
#         survival.data <- sim.FGM.copula(N=no_subjects[a], prop=prop[c], phi=phi[b],
#                                         theta1=theta1, theta2=theta2,
#                                         theta=theta, int.censoring=TRUE,
#                                         cens.model="uniform", cens.par=cens.par[d])
#         survival.data <- process.data(survival.data, rescale.factor,
#                                       admin.censoring=TRUE, admin.cutoff[d])
#         results <- partSA.trans(survival.data, risk.every, fit.distr, extrapolate.to,
#                                 admin.cutoff[d], PSA.simulations, seed, q.allowed,
#                                 results.dir, file.identifier)
#         save(results, file=paste(c(results.dir, "final_results/", 
#                                    file.identifier, ".rda"), collapse=""))
#         count=count+1
# }
# }
# }
# }

final.results.table <- data.frame(N=integer(),
                                  xi=double(),
                                  phi=double(),
                                  cens.setting=integer(),
                                  mean.rho=double(),
                                  stdev.rho=double(),
                                  true.lambda.11=double(),
                                  mae.11=double(),
                                  me.11=double(),
                                  true.lambda.12=double(),
                                  mae.12=double(),
                                  me.12=double(),
                                  true.lambda.13=double(),
                                  mae.13=double(),
                                  me.13=double(),
                                  true.lambda.22=double(),
                                  mae.22=double(),
                                  me.22=double(),                                  
                                  true.lambda.23=double(),
                                  mae.23=double(),
                                  me.23=double(),                                   
                                  mae.true.pi1=double(),
                                  mae.partsa.pi1=double(),
                                  mae.true.pi2=double(),
                                  mae.partsa.pi2=double(),
                                  mae.true.pi3=double(),
                                  mae.partsa.pi3=double())

for(a in 1:length(no_subjects)) {
       for(b in 1:length(phi)) {
                for(c in 1:length(prop)) {
                        for(d in 1:length(cens.par)) {
                                file.identifier <- paste(c("N_", no_subjects[a], "_prop_", prop[c],
                                                           "_phi_", phi[b], "_censpar_", cens.par[d]),
                                                         collapse="")
                                load(paste(c(results.dir, "final_results/", file.identifier, ".rda"), collapse=""))
                                plot_transitions(results$avg.lambda.ed[,1,1],
                                                 results$avg.lambda.lm[,1,1],
                                                 results$true.lambda.homog.avg[1,1],
                                                 results$true.lambda.homog.inf[1,1],
                                                 results$true.lambda.homog.sup[1,1],
                                                 y.lim=c(0.95,1),
                                                 title.name= expression(paste(lambda[11]^(t)," (PFS->PFS)")),
                                                 legend.position="bottomright",
                                                 extrapolate.to[d], admin.cutoff[d])
                                plot_transitions(results$avg.lambda.ed[,1,2],
                                                 results$avg.lambda.lm[,1,2],
                                                 results$true.lambda.homog.avg[1,2],
                                                 results$true.lambda.homog.inf[1,2],
                                                 results$true.lambda.homog.sup[1,2],
                                                 y.lim=c(0,0.05),
                                                 title.name= expression(paste(lambda[12]^(t)," (PFS->PD)")),
                                                 legend.position="topright",
                                                 extrapolate.to[d], admin.cutoff[d])
                                plot_transitions(results$avg.lambda.ed[,1,3],
                                                 results$avg.lambda.lm[,1,3],
                                                 results$true.lambda.homog.avg[1,3],
                                                 results$true.lambda.homog.inf[1,3],
                                                 results$true.lambda.homog.sup[1,3],
                                                 y.lim=c(0,0.05),
                                                 title.name= expression(paste(lambda[13]^(t)," (PFS->D)")),
                                                 legend.position="topright",
                                                 extrapolate.to[d], admin.cutoff[d])
                                plot_transitions(results$avg.lambda.ed[,2,2],
                                                 results$avg.lambda.lm[,2,2],
                                                 results$true.lambda.homog.avg[2,2],
                                                 results$true.lambda.homog.inf[2,2],
                                                 results$true.lambda.homog.sup[2,2],
                                                 y.lim=c(0.95,1),
                                                 title.name= expression(paste(lambda[22]^(t)," (PD->PD)")),
                                                 legend.position="bottomright",
                                                 extrapolate.to[d], admin.cutoff[d])
                                plot_transitions(results$avg.lambda.ed[,2,3],
                                                 results$avg.lambda.lm[,2,3],
                                                 results$true.lambda.homog.avg[2,3],
                                                 results$true.lambda.homog.inf[2,3],
                                                 results$true.lambda.homog.sup[2,3],
                                                 y.lim=c(0,0.05),
                                                 title.name= expression(paste(lambda[23]^(t)," (PD->D)")),
                                                 legend.position="topright",
                                                 extrapolate.to[d], admin.cutoff[d])
                                plot_proportions(est.props.ed=results$opt.MM.ed$est.props$mest[,1],
                                                 est.props.lm=results$opt.MM.lm$est.props$mest[,1],
                                                 props.partsa=results$opt.MM.ed$est.props$m[,1],
                                                 true.props.msm.avg=results$true.props.homog.avg[,1],
                                                 true.props.msm.inf=results$true.props.homog.inf[,1],
                                                 true.props.msm.sup=results$true.props.homog.sup[,1],
                                                 title.name=expression(paste(pi[1]^(t)," (PFS)")),
                                                 legend.position="topright",
                                                 extrapolate.to=extrapolate.to[d], 
                                                 admin.cutoff=admin.cutoff[d])
                                plot_proportions(est.props.ed=results$opt.MM.ed$est.props$mest[,2],
                                                 est.props.lm=results$opt.MM.lm$est.props$mest[,2],
                                                 props.partsa=results$opt.MM.ed$est.props$m[,2],
                                                 true.props.msm.avg=results$true.props.homog.avg[,2],
                                                 true.props.msm.inf=results$true.props.homog.inf[,2],
                                                 true.props.msm.sup=results$true.props.homog.sup[,2],
                                                 title.name=expression(paste(pi[2]^(t)," (PD)")),
                                                 legend.position="topright",
                                                 extrapolate.to=extrapolate.to[d], 
                                                 admin.cutoff=admin.cutoff[d])
                                plot_proportions(est.props.ed=results$opt.MM.ed$est.props$mest[,3],
                                                 est.props.lm=results$opt.MM.lm$est.props$mest[,3],
                                                 props.partsa=results$opt.MM.ed$est.props$m[,3],
                                                 true.props.msm.avg=results$true.props.homog.avg[,3],
                                                 true.props.msm.inf=results$true.props.homog.inf[,3],
                                                 true.props.msm.sup=results$true.props.homog.sup[,3],
                                                 title.name=expression(paste(pi[3]^(t)," (D)")),
                                                 legend.position="bottomright",
                                                 extrapolate.to=extrapolate.to[d], 
                                                 admin.cutoff=admin.cutoff[d])
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
                                                              rep(results$true.lambda.homog.avg[1,1],extrapolate.to[d]-1))
                                mean.abs.error.trans12 <- mae(results$avg.lambda.ed[,1,2],
                                                              rep(results$true.lambda.homog.avg[1,2],extrapolate.to[d]-1))
                                mean.abs.error.trans13 <- mae(results$avg.lambda.ed[,1,3],
                                                              rep(results$true.lambda.homog.avg[1,3],extrapolate.to[d]-1))
                                mean.abs.error.trans22 <- mae(results$avg.lambda.ed[,2,2],
                                                              rep(results$true.lambda.homog.avg[2,2],extrapolate.to[d]-1))
                                mean.abs.error.trans23 <- mae(results$avg.lambda.ed[,2,3],
                                                              rep(results$true.lambda.homog.avg[2,3],extrapolate.to[d]-1))
                                me.trans11 <- me(results$avg.lambda.ed[,1,1],rep(results$true.lambda.homog.avg[1,1],extrapolate.to[d]-1))
                                me.trans12 <- me(results$avg.lambda.ed[,1,2],rep(results$true.lambda.homog.avg[1,2],extrapolate.to[d]-1))
                                me.trans13 <- me(results$avg.lambda.ed[,1,3],rep(results$true.lambda.homog.avg[1,3],extrapolate.to[d]-1))
                                me.trans22 <- me(results$avg.lambda.ed[,2,2],rep(results$true.lambda.homog.avg[2,2],extrapolate.to[d]-1))
                                me.trans23 <- me(results$avg.lambda.ed[,2,3],rep(results$true.lambda.homog.avg[2,3],extrapolate.to[d]-1))
                                relative.mad.trans11 <- mean.abs.error.trans11/results$true.lambda.homog.avg[1,1]
                                relative.mad.trans12 <- mean.abs.error.trans12/results$true.lambda.homog.avg[1,2]
                                relative.mad.trans13 <- mean.abs.error.trans13/results$true.lambda.homog.avg[1,3]
                                relative.mad.trans22 <- mean.abs.error.trans22/results$true.lambda.homog.avg[2,2]
                                relative.mad.trans23 <- mean.abs.error.trans23/results$true.lambda.homog.avg[2,3]
                                load(paste(c(results.dir, "optimal_rho_data/ED_", file.identifier, ".rda"), collapse=""))
                                if (cens.par[d]==3) {  
                                        cens.setting <- 2        
                                } else if (cens.par[d]==4) {
                                        cens.setting <- 1        
                                }
                                final.results.table[count, "N"] <- no_subjects[a]
                                final.results.table[count, "xi"] <- prop[c]        
                                final.results.table[count, "phi"] <- phi[b]        
                                final.results.table[count, "cens.setting"] <- cens.setting        
                                final.results.table[count, "mean.rho"]  <- round(mean(optimal.rho.ed), digits=4)
                                final.results.table[count, "stdev.rho"] <- round(sd(optimal.rho.ed), digits=4)
                                final.results.table[count, "true.lambda.11"] <- round(results$true.lambda.homog.avg[1,1],4)        
                                final.results.table[count, "mae.11"] <- round(mean.abs.error.trans11,4)        
                                final.results.table[count, "me.11"]  <- round(me.trans11,4)
                                final.results.table[count, "true.lambda.12"] <- round(results$true.lambda.homog.avg[1,2],4)        
                                final.results.table[count, "mae.12"] <- round(mean.abs.error.trans12,4)        
                                final.results.table[count, "me.12"]  <- round(me.trans12,4)      
                                final.results.table[count, "true.lambda.13"] <- round(results$true.lambda.homog.avg[1,3],4)        
                                final.results.table[count, "mae.13"] <- round(mean.abs.error.trans13,4)        
                                final.results.table[count, "me.13"]  <- round(me.trans13,4)    
                                final.results.table[count, "true.lambda.22"] <- round(results$true.lambda.homog.avg[2,2],4)        
                                final.results.table[count, "mae.22"] <- round(mean.abs.error.trans22,4)        
                                final.results.table[count, "me.22"]  <- round(me.trans22,4)
                                final.results.table[count, "true.lambda.23"] <- round(results$true.lambda.homog.avg[2,3],4)        
                                final.results.table[count, "mae.23"] <- round(mean.abs.error.trans23,4)        
                                final.results.table[count, "me.23"]  <- round(me.trans23,4) 
                                final.results.table[count, "mae.true.pi1"] <- round(mean.abs.error.PFS.true,4)
                                final.results.table[count, "mae.partsa.pi1"] <- round(mean.abs.error.PFS.partsa*10000, 4)
                                final.results.table[count, "mae.true.pi2"] <- round(mean.abs.error.PD.true,4)
                                final.results.table[count, "mae.partsa.pi2"] <- round(mean.abs.error.PD.partsa*10000, 4)
                                final.results.table[count, "mae.true.pi3"] <- round(mean.abs.error.D.true,4)
                                final.results.table[count, "mae.partsa.pi3"] <- round(mean.abs.error.D.partsa*10000, 4)
                                count=count+1
                         }
                }
        }
}






