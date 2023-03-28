#------------------------------------------------------------------------------#
# the simulation for the three-arm equivalence trial                           #
# 2022 10 29                                                                   #
#------------------------------------------------------------------------------#
######part0: common setting ####################################################
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%`
path <- "E:\\BaiduNetdiskWorkspace\\research\\3armBE"
path_res <- path%+%"\\result"
source(path%+%"\\code\\function\\fun-BERR.R")
powersim <- function(pair,pait,paip,nr,nt,np,upper = 1.25,lower = 0.8,nsim = 10000,method = "FM",seed = 2){
  set.seed(seed)
  effnp <- rbinom(nsim,np,paip)
  effnt <- rbinom(nsim,nt,pait)
  effnr <- rbinom(nsim,nr,pair)
  
  env <- new.env()
  assign("nr",nr,envir = env)
  assign("nt",nt,envir = env)
  assign("np",np,envir = env)
  assign("lower",lower,envir = env)
  assign("upper",upper,envir = env)
  assign("method",method,envir = env)
  
  singletest <- function(x,envir){
    effnp <- x[1]
    effnt <- x[2]
    effnr <- x[3]
    nr <- get("nr",envir = envir)
    nt <- get("nt",envir = envir)
    np <- get("np",envir = envir)
    lower <- get("lower",envir = envir)
    upper <- get("upper",envir = envir)
    method <- get("method",envir = envir)
    
    x1 <- test_BERR(p1 = effnt/nt,p2 = effnr/nr,n1 = nt,n2 = nr,upper = upper,lower = lower,method = method)
    x2 <- test_superRR(p1 = effnt/nt,p2 = effnp/np,n1 = nt,n2 = np,margin = 1,method = method)
    x3 <- test_superRR(p1 = effnr/nr,p2 = effnp/np,n1 = nr,n2 = np,margin = 1,method = method)
    xx <- c(x1[1:2],x2[1],x3[1])
    any(xx>=c(0.05,0.05,0.025,0.025))
  }
  xdat <- cbind(effnp,effnt,effnr)
  res <- apply(xdat, 1, singletest,envir = env)
  1-mean(res)  
}
powersim_equ <- function(pair,pait,nr,nt,upper = 1.25,lower = 0.8,nsim = 10000,method = "FM",seed = 2){
  set.seed(seed)
  effnt <- rbinom(nsim,nt,pait)
  effnr <- rbinom(nsim,nr,pair)
  
  env <- new.env()
  assign("nr",nr,envir = env)
  assign("nt",nt,envir = env)
  assign("lower",lower,envir = env)
  assign("upper",upper,envir = env)
  assign("method",method,envir = env)
  
  singletest <- function(x,envir){
    effnt <- x[1]
    effnr <- x[2]
    nr <- get("nr",envir = envir)
    nt <- get("nt",envir = envir)
    lower <- get("lower",envir = envir)
    upper <- get("upper",envir = envir)
    method <- get("method",envir = envir)
    
    x1 <- test_BERR(p1 = effnt/nt,p2 = effnr/nr,n1 = nt,n2 = nr,upper = upper,lower = lower,method = method)
    xx <- x1[1:2]
    any(xx>c(0.05,0.05))
  }
  xdat <- cbind(effnt,effnr)
  res <- apply(xdat, 1, singletest,envir = env)
  1-mean(res)  
}
powersim_super <- function(paip,pait,np,nt,nsim = 10000,method = "FM",seed = 2){
  set.seed(seed)
  effnt <- rbinom(nsim,nt,pait)
  effnp <- rbinom(nsim,nr,paip)
  
  env <- new.env()
  assign("np",np,envir = env)
  assign("nt",nt,envir = env)
  assign("method",method,envir = env)
  
  singletest <- function(x,envir){
    effnt <- x[1]
    effnp <- x[2]
    np <- get("nr",envir = envir)
    nt <- get("nt",envir = envir)
    method <- get("method",envir = envir)
    
    x1 <- test_superBERR(p1 = effnt/nt,p2 = effnp/np,n1 = nt,n2 = np,method = method)
    xx <- x1[1]
    xx>0.05
  }
  xdat <- cbind(effnt,effnp)
  res <- apply(xdat, 1, singletest,envir = env)
  1-mean(res)  
}


######end part0#################################################################
######part1 fixratio############################################################
ndat <- xlsx::read.xlsx(path%+%"\\code\\r_fixratio-1-1-1.xlsx",sheetIndex = 1)
## FM ----
pai_R <- 0.5
method <- "FM"
parframe <- cbind(lower = 1/ndat$delta,
                  upper = ndat$delta,
                  pai_P = ndat$Placebo,
                  pai_R = pai_R,
                  pai_T = pai_R * 0.95,
                  N_R = ndat$n1_fm,
                  N_P = ndat$n1_fm,
                  N_T = ndat$n1_fm)
nsim <- 100000
powers <- NULL
for (i in 1:dim(parframe)[1]) {
  pari <- parframe[i,]
  lower <- pari[1]
  upper <- pari[2]
  pai_P <- pari[3]
  pai_R <- pari[4]
  pai_T <- pari[5]
  N_R <- pari[6]
  N_P <- pari[7]
  N_T <- pari[8]
  poweri <- powersim(pair = pai_R,pait = pai_T,paip = pai_P,nr = N_R,nt = N_T,np = N_P,
                     upper = upper,lower = lower,nsim = nsim,method = method,seed = 12)
  powers <- c(powers,poweri)
  cat(i,method,"\n")
}
power_FM_fixratio <- cbind(parframe,round(powers,4)) 
power_FM_fixratio
xlsx::write.xlsx(power_FM_fixratio,path_res%+%"\\power_sim.xlsx",sheetName = "FM-fixratio")

## log ----
method <- "log"
parframe <- cbind(lower = 1/ndat$delta,
                  upper = ndat$delta,
                  pai_P = ndat$Placebo,
                  pai_R = pai_R,
                  pai_T = pai_R * 0.95,
                  N_R = ndat$n1_log,
                  N_P = ndat$n1_log,
                  N_T = ndat$n1_log)
nsim <- 100000
powers <- NULL
for (i in 1:dim(parframe)[1]) {
  pari <- parframe[i,]
  lower <- pari[1]
  upper <- pari[2]
  pai_P <- pari[3]
  pai_R <- pari[4]
  pai_T <- pari[5]
  N_R <- pari[6]
  N_P <- pari[7]
  N_T <- pari[8]
  
  poweri <- powersim(pair = pai_R,pait = pai_T,paip = pai_P,nr = N_R,nt = N_T,np = N_P,
                   upper = upper,lower = lower,nsim = nsim,method = method,seed = 12)
 
  powers <- c(powers,poweri)
  cat(i,method,"\n")
}
power_log_fixratio <- cbind(parframe,powers) 
power_log_fixratio
xlsx::write.xlsx(power_log_fixratio,path_res%+%"\\power_sim.xlsx",sheetName = "log-fixratio",append = T)

###### end part1###############################################################
##### part2 noratio ###########################################################
ndat <- xlsx::read.xlsx(path%+%"\\code\\r_noratio.xlsx",sheetIndex = 1)
## FM ----
pai_R <- 0.5
method <- "FM"
parframe <- cbind(lower = 1/ndat$delta,
                  upper = ndat$delta,
                  pai_P = ndat$Placebo,
                  pai_R = pai_R,
                  pai_T = pai_R * 0.95,
                  N_R = ndat$N2_fm,
                  N_P = ndat$N3_fm,
                  N_T = ndat$N1_fm)
                                    
nsim <- 100000
powers <- NULL
for (i in 1:dim(parframe)[1]) {
  pari <- parframe[i,]
  lower <- pari[1]
  upper <- pari[2]
  pai_P <- pari[3]
  pai_R <- pari[4]
  pai_T <- pari[5]
  N_R <- pari[6]
  N_P <- pari[7]
  N_T <- pari[8]
  poweri <- powersim(pair = pai_R,pait = pai_T,paip = pai_P,nr = N_R,nt = N_T,np = N_P,
                     upper = upper,lower = lower,nsim = nsim,method = method,seed = 12)
  powers <- c(powers,poweri)
  cat(i,method,"\n")
}
power_FM_noratio <- cbind(parframe,round(powers,4)) 
power_FM_noratio
xlsx::write.xlsx(power_FM_noratio,path_res%+%"\\power_sim.xlsx",sheetName = "FM-noratio",append = T)
## log ----
method <- "log"
parframe <- cbind(lower = 1/ndat$delta,
                  upper = ndat$delta,
                  pai_P = ndat$Placebo,
                  pai_R = pai_R,
                  pai_T = pai_R * 0.95,
                  N_R = ndat$N1_log,
                  N_P = ndat$N3_log,
                  N_T = ndat$N2_log)

nsim <- 100000
powers <- NULL
for (i in 1:dim(parframe)[1]) {
  pari <- parframe[i,]
  lower <- pari[1]
  upper <- pari[2]
  pai_P <- pari[3]
  pai_R <- pari[4]
  pai_T <- pari[5]
  N_R <- pari[6]
  N_P <- pari[7]
  N_T <- pari[8]
  poweri <- powersim(pair = pai_R,pait = pai_T,paip = pai_P,nr = N_R,nt = N_T,np = N_P,
                     upper = upper,lower = lower,nsim = nsim,method = method,seed = 12)
  powers <- c(powers,poweri)
  
  cat(i,method,"\n")
}

power_log_noratio <- cbind(parframe,round(powers,4)) 
power_log_noratio
xlsx::write.xlsx(power_log_noratio,path_res%+%"\\power_sim.xlsx",sheetName = "log-noratio",append = T)
##### end part2 ###############################################################
#####part3 only equivalence####################################################
pai_R <- 0.5
pai_T <- pai_R*0.95
lower <- 0.7; upper <- 1/lower
method <- "log"
nsim <- 100000
N <- samplesize_BERR(trueratio = 0.95,p2 = pai_R,upper = upper,lower = lower,
                     alpha = 0.05,beta = 0.146,ratio = 1,method = method,step = 10)
N <- as.integer(N[1])
power_res <- powersim_equ(pair = pai_R,pait = pai_T,nr = N,nt = N,upper = upper,
                          lower = lower,nsim = 100000,method = method,seed = 7)
##### end part3 ###############################################################
samplesize_threearm_mvn(trueratio = 0.95,p2 = 0.5,p3 = 0.4,upper = 1/0.7,lower = 0.7,ratio = 1,verbose = T,method = method)
#####part4 only super####################################################
pai_R <- 0.3
pai_P <- 0.5
method <- "log"
nsim <- 100000
N <- samplesize_BERR(trueratio = 0.95,p2 = pai_R,upper = upper,lower = lower,
                     alpha = 0.05,beta = 0.146,ratio = 1,method = method,step = 10)
N <- as.integer(N[1])
power_res <- powersim_equ(pair = pai_R,pait = pai_T,nr = N,nt = N,upper = upper,
                          lower = lower,nsim = 100000,method = method,seed = 7)
##### end part3 ###############################################################
samplesize_threearm_mvn_noratio(grid = 0.01,maxgrid = 0.5,trueratio = 0.95,p2 = 0.6,
                                p3 = 0.4,upper = 1/0.7,lower = 0.7,verbose = F,method = method)
xx <- 0.08
test_superRR(0.5,xx,1860,1860,1,"FM") > test_superRR(0.5,xx,1860,1860,1,"log") 

