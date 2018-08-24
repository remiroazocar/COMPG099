library(rms)
library(msm)
library(TPmsm)
source("functions.R")

partSA.trans <- function(survival.data, risk.every, fit.distr, extrapolate.to,
                         admin.cutoff, PSA.simulations, seed, q.allowed,
                         results.dir, file.identifier) {
        KM.fits <- KM.curves(survival.data, results.dir=results.dir, 
                             file.identifier=file.identifier)
        digit.surv.input(KM.fits$PFSfit, KM.fits$OSfit, results.dir=results.dir,
                         file.identifier=file.identifier)
        make.at.risk(fit=KM.fits$PFSfit, admin.cutoff=admin.cutoff, risk.every=risk.every,
                     results.dir=results.dir, endpoint="PFS")
        make.at.risk(fit=KM.fits$OSfit, admin.cutoff=admin.cutoff, risk.every=risk.every,
                     results.dir=results.dir, endpoint="OS")
        reconstructed.IPD <- reconstruct.IPD(results.dir=results.dir,
                                             file.identifier=file.identifier)
        hmc.fit <-fit.models(formula=Surv(time,event)~as.factor(arm),
                             data = reconstructed.IPD,
                             distr = fit.distr, method = "hmc", seed=seed)
        png(filename=paste(c(results.dir, "surv_fit_plots/" , file.identifier, ".png"),
                           collapse=""))
        plot(hmc.fit, xlab="t", ylab="survival probability S(t)")
        dev.off()
        png(filename=paste(c(results.dir, "PSA_plots/", file.identifier, ".png"),
                           collapse=""))
        fit.PSA <- make.surv(hmc.fit, t=seq(0,extrapolate.to-1),
                             nsim=PSA.simulations)
        psa.plot(fit.PSA)
        dev.off()
        PFS.PSA <- matrix(nrow=PSA.simulations, ncol=extrapolate.to)
        OS.PSA <- matrix(nrow=PSA.simulations, ncol=extrapolate.to)
        for(i in 1:PSA.simulations) {
                PFS.PSA[i,] <- fit.PSA$S[[i]][[1]][,2]
                OS.PSA[i,] <- fit.PSA$S[[i]][[2]][,2]
        }
        save(PFS.PSA, file=paste(c(results.dir, "PSA_trace/PFS_", 
                                   file.identifier, ".rda"), collapse=""))
        save(OS.PSA, file=paste(c(results.dir, "PSA_trace/OS_", 
                                  file.identifier, ".rda"), collapse=""))
        load(paste(c(results.dir, "PSA_trace/PFS_", file.identifier, ".rda"), collapse=""))
        load(paste(c(results.dir, "PSA_trace/OS_", file.identifier, ".rda"), collapse=""))
        optimal.rho.ed <- rho.optimisation.homog(PFS.PSA=PFS.PSA, OS.PSA=OS.PSA,
                                                 lower.rho=lower.rho, upper.rho=upper.rho,
                                                 no.rho.estimates=no.rho.estimates,
                                                 method="euclidean_distance",
                                                 results.dir=results.dir,
                                                 file.identifier=file.identifier)
        # optimal.rho.lm <- rho.optimisation.homog(PFS.PSA=PFS.PSA, OS.PSA=OS.PSA,
        #                                          lower.rho=lower.rho, upper.rho=upper.rho,
        #                                          no.rho.estimates=no.rho.estimates,
        #                                          method="linear_model",
        #                                          results.dir=results.dir,
        #                                          file.identifier=file.identifier)
        save(optimal.rho.ed, file=paste(c(results.dir, "optimal_rho_data/ED_",
                                          file.identifier, ".rda"), collapse=""))
        # save(optimal.rho.lm, file=paste(c(results.dir, "optimal_rho_data/LM_",
        #                                   file.identifier, ".rda"), collapse=""))
        load(paste(c(results.dir, "optimal_rho_data/ED_", file.identifier, ".rda"), 
                   collapse=""))
        # load(paste(c(results.dir, "optimal_rho_data/LM_", file.identifier, ".rda"),
        #            collapse=""))
        # run MM with transitions estimated using computed optimal rho
        opt.MM.ed <- recover.opt.transitions(PFS.PSA=PFS.PSA, OS.PSA=OS.PSA, optimal.rho.ed) 
        # opt.MM.lm <- recover.opt.transitions(PFS.PSA=PFS.PSA, OS.PSA=OS.PSA, optimal.rho.lm)
        lambda.ed <- opt.MM.ed$est.lambda
        # lambda.lm <- opt.MM.lm$est.lambda
        # average estimated lambda over PSA iterations
        avg.lambda.ed <- apply(opt.MM.ed$est.lambda, c(2,3,4), mean)
        # avg.lambda.lm <- apply(opt.MM.lm$est.lambda, c(2,3,4), mean)
        panel.data <- panel_processing(survival.data)
        surv.msm <- msm(state~time, subject=ID, data=panel.data,
                        qmatrix=q.allowed,gen.inits=TRUE, obstype= 2)
        true.lambda.homog <- pmatrix.msm(surv.msm, t=1, ci="normal", cl=0.95, B=PSA.simulations)
        true.lambda.homog.avg <- true.lambda.homog$estimates
        true.lambda.homog.inf <- true.lambda.homog$L
        true.lambda.homog.sup <- true.lambda.homog$U
        true.props.homog.avg <- run.MM.true.homog(true.lambda.homog.avg, extrapolate.to)
        true.props.homog.inf <- run.MM.true.homog(true.lambda.homog.inf, extrapolate.to)
        true.props.homog.sup <- run.MM.true.homog(true.lambda.homog.sup, extrapolate.to)
        
        true.lambda.inhomog <- inhomog.true.probs(survival.data, extrapolate.to, PSA.simulations)
        true.lambda.inhomog.avg <- true.lambda.inhomog$avg
        true.lambda.inhomog.inf <- true.lambda.inhomog$inf
        true.lambda.inhomog.sup <- true.lambda.inhomog$sup
        true.props.inhomog.avg <- run.MM.true.inhomog(true.lambda.inhomog.avg, extrapolate.to)
        true.props.inhomog.inf <- run.MM.true.inhomog(true.lambda.inhomog.inf, extrapolate.to)
        true.props.inhomog.sup <- run.MM.true.inhomog(true.lambda.inhomog.sup, extrapolate.to)
        return(list(lambda.ed=lambda.ed, 
                    # lambda.lm=lambda.lm,
                    avg.lambda.ed=avg.lambda.ed,
                    #avg.lambda.lm=avg.lambda.lm,
                    opt.MM.ed=opt.MM.ed,
                    #opt.MM.lm=opt.MM.lm,
                    true.lambda.inhomog=true.lambda.inhomog,
                    true.lambda.homog.avg=true.lambda.homog.avg,
                    true.lambda.homog.inf=true.lambda.homog.inf,
                    true.lambda.homog.sup=true.lambda.homog.sup,
                    true.props.homog.avg=true.props.homog.avg,
                    true.props.homog.inf=true.props.homog.inf,
                    true.props.homog.sup=true.props.homog.sup,
                    true.lambda.inhomog.avg=true.lambda.inhomog.avg,
                    true.lambda.inhomog.inf=true.lambda.inhomog.inf,
                    true.lambda.inhomog.sup=true.lambda.inhomog.sup,
                    true.props.inhomog.avg=true.props.inhomog.avg,
                    true.props.inhomog.inf=true.props.inhomog.avg,
                    true.props.inhomog.sup=true.props.inhomog.sup))
}
        


