# No stable copula random number generator in R
# Copula random number generator implemented following steps in
# http://www.ce.utexas.edu/prof/bhat/ABSTRACTS/Supp_material.pdf
sim.FGM.copula<-function(N, prop, phi, theta1,
                         theta2, theta, int.censoring=TRUE,
                         cens.model="uniform", cens.par) {
        N1 <- ceiling(N*prop)
        N2 <- N-N1
        # 1 -> 2 -> 3
        U1 <- runif(N1, min=0, max=1)
        U2 <- runif(N1, min=0, max=1)
        A <- phi*(2*U1-1)-1
        B <- (1-phi*(2*U1-1))^2+(4*phi*U2*(2*U1-1))
        V <- 2*U2/(sqrt(B)-A)
        t12 <- -theta1*log(1-U1)
        t23 <- -theta2*log(1-V)
        tPFS <- t12
        tOS <- t12+t23
        # 1 -> 3
        t13 <- rexp(n=N2, rate=theta)
        tPFS <- c(tPFS, t13)
        tOS <- c(tOS, t13)
        eventPFS <- rep(1, N)
        eventOS <- rep(1, N)
        if (int.censoring) {
                if (cens.model=="uniform") {
                        u <- runif(N, min=0, max=cens.par)
                        eventPFS[tPFS>u]<-0
                        eventOS[tOS>u]<-0
                        tOS[eventPFS==0] <- tPFS[eventPFS==0]
                }
        }
        return(data.frame("tPFS"=tPFS, "ePFS"=eventPFS,
                          "tOS"=tOS, "eOS"=eventOS))
}

generate.sim2 <- function(N, prop, meanlog.transient,
                          sdlog.transient, meanlog.direct,
                          sdlog.direct,
                          int.censoring=TRUE,
                          cens.model, cens.par) {
        N1 <- ceiling(N*prop)
        N2 <- N-N1
        # 1 -> 2 -> 3
        t12 <- rlnorm(N1, meanlog=meanlog.transient, sdlog=sdlog.transient)
        R2 <- runif(N1, min=0, max=1)
        DT1 <- plnorm(t12, meanlog=meanlog.transient,sdlog=sdlog.transient)
        t123<- qlnorm((DT1+R2*(1-DT1)), meanlog=meanlog.transient,sdlog=sdlog.transient)
        tPFS <- t12
        tOS <- t123
        # 1 -> 3
        t13 <- rlnorm(N2, meanlog=meanlog.direct, sdlog=sdlog.direct)
        tPFS <- c(tPFS, t13)
        tOS <- c(tOS, t13)
        eventPFS <- rep(1, N)
        eventOS <- rep(1, N)
        if (int.censoring) {
                if (cens.model=="uniform") {
                        u <- runif(N, min=0, max=cens.par)
                        eventPFS[tPFS>u]<-0
                        eventOS[tOS>u]<-0
                }
        }
        return(data.frame("tPFS"=tPFS, "ePFS"=eventPFS,
                          "tOS"=tOS, "eOS"=eventOS))
}

process.data <- function(survival.data, rescale.factor,
                         admin.censoring=TRUE, admin.cutoff,
                         exact="both") {
        survival.data$tPFS <- survival.data$tPFS*rescale.factor
        survival.data$tOS <- survival.data$tOS*rescale.factor
        N <- nrow(survival.data)
        if (exact=="death_only") {
                survival.data$tPFS[survival.data$tPFS==survival.data$tOS] <- ceiling(survival.data$tPFS[survival.data$tPFS==survival.data$tOS])
                survival.data$tOS <- ceiling(survival.data$tOS)
        } else if (exact=="none") {
                survival.data$tPFS <- ceiling(survival.data$tPFS)
                survival.data$tOS <- ceiling(survival.data$tOS)
        }
        if (admin.censoring) {
                survival.data$ePFS[survival.data$tPFS>admin.cutoff]<-0
                survival.data$eOS[survival.data$tOS>admin.cutoff]<-0
                survival.data$tPFS[survival.data$tPFS>admin.cutoff]<-admin.cutoff
                survival.data$tOS[survival.data$tOS>admin.cutoff]<-admin.cutoff
        }
        return(survival.data)
}

KM.curves <- function(survival.data, results.dir, file.identifier) {
        PFSsurv <- Surv(time=survival.data$tPFS, event=survival.data$ePFS)
        PFSfit <- survfit(PFSsurv~1)
        OSsurv <- Surv(time=survival.data$tOS, event=survival.data$eOS)
        OSfit <- survfit(OSsurv~1)
        png(filename=paste(c(results.dir, 'KM_curves/', file.identifier, ".png"),
                           collapse=""))
        plot(OSfit, conf.int=TRUE, col="red", xlab="t",
             ylab="survival probability S(t)", mark.time=TRUE)
        lines(PFSfit, conf.int=TRUE, col="orange", mark.time=TRUE)
        title("KM curves")
        legend(x="topright", legend=c("OS", "PFS"), col=c("red","orange"),
               lty=c(1,1))
        dev.off()
        return(list(PFSfit=PFSfit,OSfit=OSfit))
}

digit.surv.input <- function(PFSfit, OSfit, results.dir, file.identifier) {
        PFS.surv.times <- summary(PFSfit)$time
        PFS.surv.probs <- summary(PFSfit)$surv
        OS.surv.times <- summary(OSfit)$time
        OS.surv.probs <- summary(OSfit)$surv
        first_line <- c(0,1) # first line of digitised file is time=0, S(t)=1
        PFS.digitised <- rbind(first_line, data.frame(PFS.surv.times,
                                                      PFS.surv.probs))
        OS.digitised <- rbind(first_line,data.frame(OS.surv.times,
                                                    OS.surv.probs))
        # survHE digitise function requires input as .txt files
        write.table(PFS.digitised,
                    file=paste(c(results.dir, "survHE_input/surv_times/PFS_", 
                                 file.identifier, ".txt"),
                               collapse=""))
        write.table(OS.digitised,
                    file=paste(c(results.dir, "survHE_input/surv_times/OS_", 
                                 file.identifier, ".txt"),
                               collapse=""))
}

# helper function to create numbers at risk
make.at.risk=function(fit, admin.cutoff, risk.every, results.dir, endpoint) {
        surv.times <- summary(fit)$time
        t.risk <- seq(0,admin.cutoff-1,risk.every)
        lower <- vector()
        lower <- c(lower,1)
        upper <- vector()
        n.risk <- vector()
        for (i in 1:length(t.risk)) {
                if (i==length(t.risk)) {
                        n.risk.covered <- summary(fit)$n.risk[((summary(fit)$time)>=t.risk[i])]
                        upper <- c(upper, length(summary(fit)$time)+1)
                } else {
                        n.risk.covered <- summary(fit)$n.risk[((summary(fit)$time)>t.risk[i])
                                                              &((summary(fit)$time)<=t.risk[i+1])]
                        upper <- c(upper, lower[i] + length(n.risk.covered)-1)
                        lower <- c(lower, upper[i]+1)
                }
                if (is.na(summary(fit)$n.risk[lower[i]])) {
                        t.risk <- t.risk[1:i-1]
                        lower <- lower[1:i-1]
                        upper <- upper[1:i-1]
                        upper[i-1] <- length(surv.times)+1
                        break()
                }
                else {
                        n.risk <- c(n.risk, summary(fit)$n.risk[lower[i]])
                }
        }
        at.risk.data <- data.frame(t.risk, lower, upper, n.risk)
        write.table(at.risk.data, file=paste(c(results.dir, "survHE_input/at_risk_data/",  
                                               endpoint, "_", file.identifier, ".txt"),
                                               collapse=""))
}

reconstruct.IPD <- function(results.dir, file.identifier) {
        PFS.digit.out <- survHE::digitise(surv_inp=paste(c(results.dir,
                                                           "survHE_input/surv_times/PFS_", 
                                                           file.identifier,".txt"),
                                                         collapse=""),
                                          nrisk_inp=paste(c(results.dir,
                                                            "survHE_input/at_risk_data/PFS_",
                                                            file.identifier, ".txt"),
                                                          collapse=""),
                                          km_output=paste(c(results.dir,
                                                            "survHE_output/reconstructed_KM/PFS_", 
                                                            file.identifier, ".txt"),
                                                          collapse=""),
                                          ipd_output=paste(c(results.dir,
                                                             "survHE_output/reconstructed_IPD/PFS_", 
                                                             file.identifier, ".txt"),
                                                           collapse=""))
        PFS.digit.out <- survHE::digitise(surv_inp=paste(c(results.dir,
                                                           "survHE_input/surv_times/OS_", 
                                                           file.identifier, ".txt"),
                                                         collapse=""),
                                          nrisk_inp=paste(c(results.dir,
                                                            "survHE_input/at_risk_data/OS_", 
                                                            file.identifier, ".txt"),
                                                          collapse=""),
                                          km_output=paste(c(results.dir,
                                                            "survHE_output/reconstructed_KM/OS_", 
                                                            file.identifier, ".txt"),
                                                          collapse=""),
                                          ipd_output=paste(c(results.dir,
                                                             "survHE_output/reconstructed_IPD/OS_",
                                                             file.identifier, ".txt"),
                                                           collapse=""))
        ### reconstruct IPD ###
        PFS.IPD <- make.ipd(ipd_files=paste(c(results.dir, "survHE_output/reconstructed_IPD/PFS_", 
                                              file.identifier, ".txt"),
                                            collapse=""))
        OS.IPD <- make.ipd(ipd_files=paste(c(results.dir, "survHE_output/reconstructed_IPD/OS_", 
                                             file.identifier, ".txt"),
                                           collapse=""))
        OS.IPD[,3]<-1
        full.IPD <-  rbind(PFS.IPD, OS.IPD)
        return(full.IPD)
}

# make.lambda=function(os,pfs,rho){
#         #Define number of post progression
#         pp<-os-pfs
#         #Divide by pfs so must only calculate for pfs>0
#         t.max<-length(which(pfs>1e-15))
#         #Calculate transition probabilities for full time horizen
#         t.length<-length(pfs)
#         lambda=array(0,c(t.length-1,3,3))
# 
#         #Transitions from pre-progression
#         lambda[1:(t.max-1),1,1]=pfs[2:t.max]/pfs[1:(t.max-1)]
#         lambda[1:(t.max-1),1,2]=(1-pfs[2:t.max]/pfs[1:(t.max-1)])*(1-rho)
#         lambda[1:(t.max-1),1,3]=(1-pfs[2:t.max]/pfs[1:(t.max-1)])*rho
# 
#         #Transition from progressed
#         #At time point 1 there is no one in post-progression so cannot divide
#         #(lambda[1,2,2] is set to 0)
#         lambda[2:(t.max-1),2,2]<-pmin(1,pmax(0,(pp[3:t.max]-pfs[2:(t.max-1)]*
#                                                         lambda[2:(t.max-1),1,2]))/
#                                               pp[2:(t.max-1)])
#         lambda[,2,3]<-1-lambda[,2,2]
# 
#         #Death is an absorbing state
#         lambda[,3,3]<-1
# 
#         #Calculation of state membership
#         #Truth
#         m=cbind(pfs,pp,1-os)
#         #Estimated
#         mest=matrix(0,t.length,3)
#         mest[1,]=c(1,0,0)
#         for (t in 2:nrow(mest)) {
#                 mest[t,]=mest[(t-1),]%*%lambda[t-1,,]
#         }
#         list(mest=mest,m=m,lambda=lambda)
# }

###Function that can do the optimisation
make.lambda=function(os,pfs,rho){
        #Define number of post progression
        pp<-os-pfs
        #Divide by pfs so must only calculate for pfs>0
        t.max<-length(which(pfs>1e-15))
        #Calculate transition probabilities for full time horizen
        t.length<-length(pfs)
        lambda=array(0,c(t.length-1,3,3))

        #Transitions from pre-progression
        lambda[1:(t.max-1),1,1]=pfs[2:t.max]/pfs[1:(t.max-1)]
        lambda[1:(t.max-1),1,2]=(1-pfs[2:t.max]/pfs[1:(t.max-1)])*(1-rho)
        lambda[1:(t.max-1),1,3]=(1-pfs[2:t.max]/pfs[1:(t.max-1)])*rho

        #Transition from progressed
        #At time point 1 there is no one in post-progression so cannot divide (lambda[1,2,2] is set to 0)
        lambda[2:(t.max-1),2,2]<-pmin(1,pmax(0,(pp[3:t.max]-pfs[2:(t.max-1)]*lambda[2:(t.max-1),1,2]))/
                                              pp[2:(t.max-1)])
        lambda[,2,3]<-1-lambda[,2,2]

        #Death is an absorbing state
        lambda[,3,3]<-1

        #Calculation of state membership
        #Truth
        m=cbind(pfs,pp,1-os)
        #Estimated
        mest=matrix(0,t.length,3)
        mest[1,]=c(1,0,0)
        for (t in 2:nrow(mest)) {
                mest[t,]=mest[(t-1),]%*%lambda[t-1,,]
        }
        list(mest=mest,m=m,lambda=lambda)
}

check.rho.ed=function(os,pfs,rho){
        x=make.lambda(os=os,pfs=pfs,rho=rho)
        diff<-sum((x$m-x$mest)^2)
        return(diff)
}

## Computes Linear Model Difference
check.rho.lm=function(os, pfs, rho){
        x=make.lambda(os=os,pfs=pfs,rho=rho)
        sdy<-sd(x$mest[,2])
        sdx<-sd(x$m[,2])
        r<-cor(x$mest[,2],x$m[,2])
        xbar<-mean(x$mest[,2])
        ybar<-mean(x$m[,2])
        chk1=c(r*sdy/sdx,ybar-r*sdy/sdx*xbar)
        sdy<-sd(x$mest[,3])
        sdx<-sd(x$m[,3])
        r<-cor(x$mest[,3],x$m[,3])
        xbar<-mean(x$mest[,3])
        ybar<-mean(x$m[,3])
        chk2=c(r*sdy/sdx,ybar-r*sdy/sdx*xbar)
        return(cbind(chk1,chk2))
}

# rho.optimisation.homog <- function(PFS.PSA, OS.PSA,
#                                    lower.rho, upper.rho, no.rho.estimates,
#                                    method="euclidean_distance",
#                                    results.dir, file.identifier) {
#         iter <- nrow(PFS.PSA)
#         rhomatrix.optimised<-matrix(nrow=iter,ncol=1)
#         cat(paste0("\n*Starting rho optimisation - ", method, " . ",
#                    no.rho.estimates,
#                    " estimates of rho to be made for each PSA iteration. There are ",
#                    iter, " PSA iterations."))        
#         for(i in 1:iter){
#                 cat(paste0(i, "\n"))
#                 os<-OS.PSA[i,]
#                 pfs<-PFS.PSA[i,]
#                 rho.poss<-seq(lower.rho,upper.rho,length.out=no.rho.estimates)
#                 min_distance = Inf
#                 for (j in 1:no.rho.estimates) {
#                         if (method=="euclidean_distance") {
#                                 distance <- sum(check.rho.ed(os=os,pfs=pfs,rho.poss[j]))
#                                 if (distance < min_distance) {
#                                         min_distance <- distance
#                                         opt_rho <- rho.poss[j]
#                                 }
#                         } else if (method=="linear_model") {
#                                 distance <- sum((check.rho.lm(os=os,pfs=pfs,rho.poss[j])-c(1,0,1,0))^2)
#                                 if (distance < min_distance) {
#                                         min_distance <- distance
#                                         opt_rho <- rho.poss[j]
#                                 }                                
#                         }
#                 }
#                 rhomatrix.optimised[i,1] <- opt_rho
#         }
#         if (method=="euclidean_distance") {
#                 png(filename=paste(c(results.dir, "optimal_rho_plots/ED_", 
#                                      file.identifier, ".png"), collapse=""))
#                 plot(rhomatrix.optimised, xlab="iteration", ylab="optimal rho")
#                 title("Rho optimisation (E.D.)")
#         } else if (method=="linear_model") {
#                 png(filename=paste(c(results.dir, "optimal_rho_plots/LM_", 
#                                      file.identifier, ".png"), collapse=""))
#                 plot(rhomatrix.optimised, xlab="iteration", ylab="optimal rho")
#                 title("Rho optimisation (L.M.)")
#         }
#         dev.off()
#         return(rhomatrix.optimised)
# }

rho.optimisation.homog <- function(PFS.PSA, OS.PSA,
                                   lower.rho, upper.rho, no.rho.estimates,
                                   method="euclidean_distance",
                                   results.dir, file.identifier) {
        iter <- nrow(PFS.PSA)
        rhomatrix.optimised<-matrix(nrow=iter,ncol=1)
        cat(paste0("\n*Starting rho optimisation - ", method, " . ",
                   no.rho.estimates,
                   " estimates of rho to be made for each PSA iteration. There are ",
                   iter, " PSA iterations."))        
        for(i in 1:iter){
                cat(paste0(i, "\n"))
                os<-OS.PSA[i,]
                pfs<-PFS.PSA[i,]
                rho.poss<-seq(lower.rho,upper.rho,length.out=no.rho.estimates)
                min_distance <- Inf
                if (method=="euclidean_distance") {
                        for (j in 1:no.rho.estimates) {
                                distance <- sum(check.rho.ed(os=os,pfs=pfs,rho.poss[j]))
                                if (distance < min_distance) {
                                        min_distance <- distance
                                        opt_rho <- rho.poss[j]
                                }
                        }
                } else if (method=="linear_model") {
                        for (j in 1:no.rho.estimates) {
                                distance <- sum((check.rho.lm(os=os,pfs=pfs,rho.poss[j])-c(1,0,1,0))^2)
                                if (distance < min_distance) {
                                        min_distance <- distance
                                        opt_rho <- rho.poss[j]
                                }
                        }
                }
                rhomatrix.optimised[i,1] <- opt_rho
        }
        if (method=="euclidean_distance") {
                png(filename=paste(c(results.dir, "optimal_rho_plots/ED_", 
                                     file.identifier, ".png"), collapse=""))
                plot(rhomatrix.optimised, xlab="iteration", ylab="optimal rho")
                title("Rho optimisation (E.D.)")
        } else if (method=="linear_model") {
                png(filename=paste(c(results.dir, "optimal_rho_plots/LM_", 
                                     file.identifier, ".png"), collapse=""))
                plot(rhomatrix.optimised, xlab="iteration", ylab="optimal rho")
                title("Rho optimisation (L.M.)")
        }
        dev.off()
        return(rhomatrix.optimised)
}


# runs MM with optimal rho
recover.opt.transitions <- function(PFS.PSA, OS.PSA, optimal.rho) {
        PSA.simulations <- nrow(PFS.PSA)
        est.lambda = array(0,c(PSA.simulations, length(OS.PSA[1,])-1,3,3))        
        for (i in 1:PSA.simulations) {
                os <- OS.PSA[i,]
                pfs <- PFS.PSA[i,]
                est.props <- make.lambda(os=os, pfs=pfs, rho=optimal.rho[i])
                est.lambda[i,,,] <- est.props$lambda
        }
        return(list(est.props=est.props, est.lambda=est.lambda))
}



# # doesn't work well
# rho.optimisation.inhomog <- function(os,pfs,lower.rho,upper.rho,no.rho.estimates) {
#         m <- cbind(pfs,os-pfs,1-os)
#         t.length <- length(os)
#         opt.rho <- rep(1,t.length-1)
#         mest <- matrix(0,t.length,3)
#         mest[1,] <- c(1,0,0)
#         lambda <- array(0,c(t.length-1,3,3))
#         lambda[,3,3] <- 1
#         rho.poss <- seq(lower.rho,upper.rho,length.out=no.rho.estimates)
#         #time=1 - transition 1-2
#         min.distance <- Inf
#         for (j in 1:no.rho.estimates) {
#                 lambda[1,1,1] <- pfs[2]/pfs[1]
#                 lambda[1,1,2] <- (1-pfs[2]/pfs[1])*rho.poss[j]
#                 lambda[1,1,3] <- (1-pfs[2]/pfs[1])*(1-rho.poss[j])
#                 check <- which(lambda<0, arr.ind=TRUE)
#                 if(length(check)!=0){
#                         stop("Negative Probs!!!", print(i))
#                 }
#                 mest[2,] <- t(lambda[1,,])%*%mest[1,]
#                 distance <- sum((mest[2,]-m[2,])^2)
#                 if (distance < min.distance) {
#                         min.distance <- distance
#                         opt.rho[1] <- rho.poss[j]
#                 }
#         }
#         #time=2 - transition 2-3
#         #At time point 1 there is no one in post-progression
#         # t=2:T-1
#         for (i in 2:(t.length-1)) {
#                 min.distance <- Inf 
#                 for (j in 1:no.rho.estimates) {
#                         lambda[i,1,1] <- pfs[i+1]/pfs[i]
#                         lambda[i,1,2] <- (1-pfs[i+1]/pfs[i])*rho.poss[j]
#                         lambda[i,1,3] <- (1-pfs[i+1]/pfs[i])*(1-rho.poss[j])
#                         lambda[i,2,2] <- pmin(1,pmax(0, (os[i+1]-pfs[i+1]-pfs[i]*lambda[i-1,1,2]))/
#                                                       (os[i]-pfs[i]))
#                         lambda[i,2,3] <- 1-lambda[i,2,2]
#                         check <- which(lambda<0,arr.ind=TRUE)
#                         if(length(check)!=0) {
#                                 stop("Negative Probs!!!", print(i))
#                         }
#                         mest[i+1,] <- t(lambda[1,,])%*%mest[1,]
#                         distance <- sum(sum((mest[i+1,]-m[i+1,])^2)-c(1,0,1,0))^2
#                         if (distance < min.distance) {
#                                 min.distance <- distance
#                                 opt.rho[i] <- rho.poss[j]
#                         }
#                 }
#         }
#         return(opt.rho)
# }

panel_processing=function(sim_data, exact="both") {
        panel_states <- vector()
        panel_times <- vector()
        panel_ids <- vector()
        for (i in 1:nrow(sim_data)) {
                if (sim_data[i,]$tPFS < sim_data[i,]$tOS) {
                        if ((exact=="death_only")||(exact=="none")) {
                                panel_notprog <- rep(1,sim_data[i,]$tPFS)
                                panel_prog <- rep(2, ceiling(sim_data[i,]$tOS)-
                                                             sim_data[i,]$tPFS)
                                patient_ids <- rep(i, ceiling(sim_data[i,]$tOS))
                                patient_times = seq(0, ceiling(sim_data[i,]$tOS)-1)
                                panel_states <- c(panel_states, panel_notprog, panel_prog)
                        } else if (exact=="both") {
                                panel_notprog <- rep(1, ceiling(sim_data[i,]$tPFS))
                                panel_prog <- rep(2, ceiling(sim_data[i,]$tOS)-floor(sim_data[i,]$tPFS))
                                patient_ids <- rep(i, ceiling(sim_data[i,]$tOS)+ceiling(sim_data[i,]$tPFS)
                                                              -floor(sim_data[i,]$tPFS))
                                patient_times <- seq(0, floor(sim_data[i,])$tPFS)
                                if (sim_data[i,]$tPFS!=floor(sim_data[i,]$tPFS)) {
                                        patient_times <- c(patient_times, sim_data[i,]$tPFS)
                                }
                                if ((floor(sim_data[i,]$tPFS)+1)<=(ceiling(sim_data[i,]$tOS)-1)) {
                                        patient_times <- c(patient_times, seq(floor(sim_data[i,]$tPFS)+1,
                                                                              ceiling(sim_data[i,]$tOS)-1))
                                }
                                panel_states <- c(panel_states, panel_notprog, panel_prog)
                        }
                } else if (sim_data[i,]$tPFS == sim_data[i,]$tOS) {  
                        if ((exact=="death_only")||(exact=="both")||(exact=="none")) {
                                panel_notprog <- rep(1,ceiling(sim_data[i,]$tPFS))
                                patient_ids <- rep(i,ceiling(sim_data[i,]$tPFS))
                                patient_times <- seq(0, ceiling(sim_data[i,]$tPFS)-1)
                                panel_states <- c(panel_states, panel_notprog)
                        }
                }        
                if (sim_data[i,]$eOS==1) {
                        patient_ids <- c(patient_ids, i)
                        patient_times <- c(patient_times, sim_data[i,]$tOS)
                        panel_states <- c(panel_states, 3)
                }
                panel_ids <- c(panel_ids, patient_ids)
                panel_times <- c(panel_times, patient_times)
        }
        panel_data <- data.frame(cbind(panel_ids, panel_times, panel_states))
        colnames(panel_data) <- c("ID", "time", "state")
        return(panel_data) 
}

inhomog.true.probs <- function(survival.data, extrapolate.to, PSA.simulations) {
        inhomog.true.lambda.avg = array(0,c(extrapolate.to,3,3)) 
        inhomog.true.lambda.inf = array(0,c(extrapolate.to,3,3)) # lower CI
        inhomog.true.lambda.sup = array(0,c(extrapolate.to,3,3)) # upper CI
        survTP.obj <- with(survival.data, survTP(tPFS, ePFS, tOS, eOS))
        for (i in 1:extrapolate.to) {
                print(i)
                trans.probs <- transAJ(object=survTP.obj, s=as.numeric(i-1),
                                       t=as.numeric(i), conf=TRUE, 
                                       conf.level=0.95, n.boot=PSA.simulations)
                inhomog.true.lambda.avg[i,1,1] <- trans.probs$est[1]
                inhomog.true.lambda.avg[i,1,2] <- trans.probs$est[2]
                inhomog.true.lambda.avg[i,1,3] <- trans.probs$est[3]
                inhomog.true.lambda.avg[i,2,2] <- trans.probs$est[4]
                inhomog.true.lambda.avg[i,2,3] <- trans.probs$est[5]
                # lower CI (2.5%)
                inhomog.true.lambda.inf[i,1,1] <- trans.probs$inf[1]
                inhomog.true.lambda.inf[i,1,2] <- trans.probs$inf[2]
                inhomog.true.lambda.inf[i,1,3] <- trans.probs$inf[3]
                inhomog.true.lambda.inf[i,2,2] <- trans.probs$inf[4]
                inhomog.true.lambda.inf[i,2,3] <- trans.probs$inf[5]                
                # upper CI (97.5%)
                inhomog.true.lambda.sup[i,1,1] <- trans.probs$sup[1]
                inhomog.true.lambda.sup[i,1,2] <- trans.probs$sup[2]
                inhomog.true.lambda.sup[i,1,3] <- trans.probs$sup[3]
                inhomog.true.lambda.sup[i,2,2] <- trans.probs$sup[4]
                inhomog.true.lambda.sup[i,2,3] <- trans.probs$sup[5]
        }
        return(list(avg=inhomog.true.lambda.avg,
                    inf=inhomog.true.lambda.inf,
                    sup=inhomog.true.lambda.sup))
}

# runs MM for homogeneous empirical trans probs
run.MM.true.homog <- function(true.lambda, extrapolate.to) {
        est.props <- matrix(0, extrapolate.to, 3)
        est.props[1,] <- c(1,0,0)
        for (t in 2:nrow(est.props)) {
                est.props[t,] <- est.props[(t-1),]%*%true.lambda
        }
        return(est.props)
}

# runs MM for inhomogeneous empirical trans probs
run.MM.true.inhomog <- function(true.lambda, extrapolate.to) {
        est.props <- matrix(0, extrapolate.to, 3)
        est.props[1,]=c(1,0,0)
        for (t in 2:nrow(est.props)) {
                est.props[t,] <- est.props[(t-1),]%*%true.lambda[t-1,,]        
        }
        return(est.props)
}
  
plot_transitions=function(est.transitions.ed,
                          est.transitions.lm,
                          true.transitions.avg, 
                          true.transitions.inf,
                          true.transitions.sup,
                          y.lim, title.name, legend.position, extrapolate.to,
                          admin.cutoff) {
        plot(est.transitions.ed, col="red", type="s", pch=2, ylim=y.lim,
             xlab="t", ylab="probability", lwd=2)
        abline(h=true.transitions.avg, col="blue", lty="dashed",lwd=2)
        abline(h=true.transitions.inf, col="blue", lty="dotted", lwd=2)
        abline(h=true.transitions.sup, col="blue", lty="dotted", lwd=2)
        title(title.name)
        if (admin.cutoff < extrapolate.to) {
                abline(v=admin.cutoff, col="forestgreen", lty="dotted", lwd=2)
        }
        legend(x=legend.position, legend=c("empirical", "estimated (ED)"),
               col=c("blue", "red"), pch=c(NA,NA,NA), lty=c(2,1),
               lwd=c(2,2))
}

plot_proportions=function(est.props.ed,
                          est.props.lm,
                          props.partsa,
                          props.partsa.inf,
                          props.partsa.sup,
                          true.props.msm.avg, 
                          true.props.msm.inf,
                          true.props.msm.sup,
                          title.name, legend.position, extrapolate.to, admin.cutoff) {
        plot(est.props.lm, col="red", type="s", pch=2, ylim=c(0,1),
             xlab="t", ylab="probability", lwd=2)
        lines(props.partsa, col="orange", type="s", lty="dotdash", lwd=2)
        lines(true.props.msm.avg, col="blue", type="s", lty="dashed",lwd=2)
        lines(true.props.msm.inf, col="blue", type="s",lty="dotted", lwd=2)
        lines(true.props.msm.sup, col="blue", type="s", lty="dotted", lwd=2)
        title(title.name)
        if (admin.cutoff < extrapolate.to) {
                abline(v=admin.cutoff, col="forestgreen", lty="dotted", lwd=2)
        }
        legend(x=legend.position, legend=c("empirical (MM)", "fitted (PartSA)",
                                           "estimated (ED)"),
               col=c("blue", "orange", "red"), 
               lty=c(2, 4, 1), lwd=c(2,2,2))
}

