#------------------------------------------------------------------------------#
# The following functions is to calculate the sample sizes for three-arm       #
# equivalence trials with binary endpoints                                     #
# The last part is several examples.                                           #
# Date: 2022-10-11                                                              #
# Author: ------------------------------------                                 #
#------------------------------------------------------------------------------#
#############################function###########################################
#' @title test_BERR
#' @description This function,test_BERR, is used to implement the statistical test 
#' for the two-arm equivalence trial in term of the ratio of two effective rate, i.e.,RR.\n
#' H0: p1/p2 ≥ upper or p1/p2 ≤ lower;\n
#' H1: lower < p1/p2 < upper.
#' @param p1 The effective rate of T drug.
#' @param p2 The effective rate of R drug.
#' @param n1 The sample size for T drug group.
#' @param n2 The sample size for R drug group.
#' @param upper The upper boundary value of the equivalence hypothesis.
#' @param lower The lower boundary value of the equivalence hypothesis.
#' @param method A character indicating the statistic test method. "FM" is for the
#' #' Farrington Manning’s test. "log" is for the test in term of the logarithm of a ratio.
#' For "FM" and "log",the origin hypothesis can be transformed to the following, respectively.
#' For "FM", H0: p1 - upper · p2 ≥ 0 or p2 - p1/lower ≥ 0; \n
#'           H1: p1 - upper · p2 < 0 & p2 - p1/lower < 0.\n
#' For "log", H0: log(p1)-log(p2) ≥ log(upper) or log(p1)-log(p2) ≤ log(lower); \n
#'            H1: log(lower) < log(p1)-log(p2) < log(upper).
#' @return A vector of numeric. p_upper and p_lower are the test p-value for the 
#' upper null hypothesis,i.e.,p1/p2 ≥ upper, and the lower null hypothesis,i.e.,p1/p2 ≤ lower,
#' respectively. Only when both p_upper < alpha and p_lower < alpha are met. T and R drug can be 
#' considered equivalent. For method = "log" , the 90% CI of p1/p2 is given.   
#' @example      
#' test_BERR(p1 = 0.45,p2 = 0.5,n1 = 460,n2 = 460,upper = 1.25,lower = 0.8,method = "FM")
test_BERR <- function(p1,p2,n1,n2,upper,lower,method = "FM"){
  if (missing(lower)) lower <- 1/upper else lower <- lower
  ratio <- n2/n1
  
  a <- (1 + ratio)
  b <- -(upper * (1 + ratio * p2) + ratio + p1)
  c <- upper * (p1 + ratio * p2)
  p10upper <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p20upper <- p10upper/upper
  b <- -(lower * (1 + ratio * p2) + ratio + p1)
  c <- lower * (p1 + ratio * p2)
  p10lower <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p20lower <- p10lower/lower
  #browser()
  
  if(method == "FM"){
    sigmaupper <- sqrt((p10upper * (1 - p10upper) + 
                          upper^2 * p20upper * (1 - p20upper)/ratio)/n1)
    testupper <- (p1-p2*upper)/sigmaupper
    pupper <- pnorm(testupper)
    
    sigmalower <- sqrt((p10lower * (1 - p10lower) + 
                          lower^2 * p20lower * (1 - p20lower)/ratio)/n1)
    testlower <- (p1-p2*lower)/sigmalower
    plower <- 1 - pnorm(testlower)
    x <- c(pupper,plower)
    names(x) <- c("p_upper","p_lower")
  } else if(method == "log")  {
    sigmaupper <- sqrt((1-p10upper)/n1/p10upper + (1-p20upper)/n2/p20upper)
    testupper <- (log(p1)-log(p2)-log(upper))/sigmaupper
    pupper <- pnorm(testupper)
    sigmalower <- sqrt((1-p10lower)/n1/p10lower + (1-p20lower)/n2/p20lower)
    testlower <- (log(p1)-log(p2)-log(lower))/sigmalower
    plower <- 1 - pnorm(testlower)
    
    CIupper <- exp((log(p1)-log(p2)) + qnorm(0.9)*sigmaupper)
    CIlower <- exp((log(p1)-log(p2)) - qnorm(0.9)*sigmalower)
    x <- c(pupper,plower,CIupper,CIlower)
    names(x) <- c("p_upper","p_lower","90%CI-upper","90%CI-lower")
  } else stop("method must be log or FM")
  x
}

#' @title test_superRR
#' @description This function, test_superRR, is used to implement the statistical test 
#' for the two-arm superiority trial in term of the ratio of two effective rate.\n
#' H0: p1/p2 ≤ margin;\n
#' H1: p1/p2 > margin.
#' @param p1 The effective rate of T drug.
#' @param p2 The effective rate of placebo or control drug.
#' @param n1 The sample size for T drug group.
#' @param n2 The sample size for placebo or control drug.
#' @param margin The superior boundary value.
#' @param method A character indicating the statistic test method. "FM" is for 
#' Farrington Manning’s test. "log" is for the test in term of the logarithm of a ratio.
#' For "FM" and "log",the origin hypothesis can be transformed to the following, respectively.
#' For "FM", H0: p1 - p2 · margin ≤ 0 ; H1: p1 - p2 · margin > 0 .\n
#' For "log", H0: log(p1)-log(p2) ≤ log(margin);H1: log(p1)-log(p2) > log(margin).
#' @return a one side p-value
#' @example 
#' test_superRR(0.5,0.4,300,200)
test_superRR <- function(p1,p2,n1,n2,margin = 1,method = "FM"){
  lower <- margin
  ratio <- n2/n1
  
  d0 <- margin == 1
  a <- (1 + ratio)
  b <- -(lower * (1 + ratio * p2) + ratio + p1)
  c <- lower * (p1 + ratio * p2)
  p10lower <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p20lower <- p10lower/lower
  p10lower[d0] <- (p1[d0] + ratio[d0] * p2[d0])/(1 + ratio[d0])
  p20lower[d0] <- p10lower[d0]
  #browser()
  if(method == "FM"){
    sigmalower <- sqrt((p10lower * (1 - p10lower) + 
                          lower^2 * p20lower * (1 - p20lower)/ratio)/n1)
    testlower <- (p1-p2*lower)/sigmalower
    plower <- 1 - pnorm(testlower)
    x <- c(plower)
    names(x) <- c("p")
  } else if(method == "log")  {
    sigmalower <- sqrt((1-p10lower)/n1/p10lower + (1-p20lower)/n2/p20lower)
    testlower <- (log(p1)-log(p2)-log(lower))/sigmalower
    plower <- 1 - pnorm(testlower)
    
    CIupper <- exp((log(p1)-log(p2)) + qnorm(0.9)*sigmalower)
    CIlower <- exp((log(p1)-log(p2)) - qnorm(0.9)*sigmalower)
    x <- c(plower,CIupper,CIlower)
    names(x) <- c("p","90%CI-upper","90%CI-lower")
  } else stop("method must be log or FM")
  x
}


#' @title samplesize_BERR
#' @description This function is to calculate the sample size for two-arm equivalence trial 
#' in term of the ratio of two effective rates. The test statistic is seen in function, test_BERR.
#' @param trueratio The assumed true ratio of p1/p2. Commonly set as 0.95.
#' @param p2,upper,lower,method Same as that in function, test_BERR.
#' @param ratio The sample size ratio of T drug and R drug group,i.e.,n1/n2.
#' @param alpha The test critical p-value, i.e., type-I error.
#' @param beta The type-II error, i.e.,1-power.
#' @param step The length of step in iterative calculation. If missing, it will be automatically calculated.
#' @return A vector of numeric. n1 and n2 are the sample sizes of T and R drug, respectively. n is n1+n2.
#' power is the calculated power based on n1 and n2.
#' @example 
#' samplesize_BERR(0.95,0.5,1.25,0.8,alpha = 0.05,beta = 0.2,ratio = 1,method = "FM")
samplesize_BERR <- function(trueratio,p2,upper,lower,alpha = 0.05,beta = 0.2,
                            ratio = 1,method = "log",step){
  p1 <- p2 * trueratio
  if (missing(lower)) lower <- 1/upper else lower <- lower
  # give a initial sample n0 based on anticipate variance
  sigma001 <- (1-p1)/p1 + (1-p2)/p2/ratio
  zalpha <- qnorm(1-alpha)
  zbeta <- qnorm(1-beta/2)
  aa <- max((log(upper)-log(trueratio))^2,(-log(lower)+log(trueratio))^2)
  n0 <- floor((zalpha+zbeta)^2*sigma001/aa)
  if(n0 > 10000000) stop("The sample size is too large")
  if (missing(step)) step <- ceiling(n0 / 10) else step <- step
  if(step <2) stop("please set step again")
  a <- (1 + ratio)
  b <- -(upper * (1 + ratio * p2) + ratio + p1)
  c <- upper * (p1 + ratio * p2)
  p10upper <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p20upper <- p10upper/upper
  b <- -(lower * (1 + ratio * p2) + ratio + p1)
  c <- lower * (p1 + ratio * p2)
  p10lower <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p20lower <- p10lower/lower
  
  #browser()
  zalpha <- qnorm(1-alpha)
  if(method == "FM"){
    sigma1u <- sqrt(p1*(1-p1)+upper^2*(1-p2)*p2/ratio)
    sigma1l <- sqrt(p2*(1-p2)/ratio+p1*(1-p1)/lower^2)
    sigma0u <- sqrt(p10upper*(1-p10upper)+upper^2*(1-p20upper)*p20upper/ratio)
    sigma0l <- sqrt(p20lower*(1-p20lower)/ratio+p10lower*(1-p10lower)/lower^2)
    bbupper <- p1-upper*p2
    bblower <- p2-p1/lower
    limit1 <- -zalpha*sigma0u/sigma1u-bbupper/sigma1u*sqrt(n0)
    limit2 <- -zalpha*sigma0l/sigma1l-bblower/sigma1l*sqrt(n0)
    cor12 <- diag(1,nrow = 2)
    cor12[1,2] <- (-p1*(1-p1)/lower-p2*(1-p2)*upper/ratio)/sigma1u/sigma1l
    cor12[2,1] <- cor12[1,2]
    power0 <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf),upper = c(limit1,limit2),
                               mean = c(0,0),corr = cor12)
    fuhao <- power0 > 1-beta
    while (step > 0.25) {
      hh <- ifelse(fuhao,-1,1)
      n0 <- n0 + step*hh
      limit1 <- -zalpha*sigma0u/sigma1u-bbupper/sigma1u*sqrt(n0)
      limit2 <- -zalpha*sigma0l/sigma1l-bblower/sigma1l*sqrt(n0)
      power0 <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf),upper = c(limit1,limit2),
                                 mean = c(0,0),corr = cor12)
      if((power0 > 1-beta) != fuhao){
        step <- step/2
        fuhao <- !fuhao
      }
      #cat(n0,step,power0,"\n")
    }
    n0 <- round(n0)
    limit1 <- -zalpha*sigma0u/sigma1u-bbupper/sigma1u*sqrt(n0)
    limit2 <- -zalpha*sigma0l/sigma1l-bblower/sigma1l*sqrt(n0)
    power0 <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf),upper = c(limit1,limit2),
                               mean = c(0,0),corr = cor12)
    if(power0 < 0.8) {
      n0 <- n0 + 1
      limit1 <- -zalpha*sigma0u/sigma1u-bbupper/sigma1u*sqrt(n0)
      limit2 <- -zalpha*sigma0l/sigma1l-bblower/sigma1l*sqrt(n0)
      power0 <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf),upper = c(limit1,limit2),
                                 mean = c(0,0),corr = cor12)
    }
    n1 <- n0
    n2 <- round(n1 * ratio)
  } else if(method == "log")  {
    sigmaupper1 <- sqrt((1-p10upper)/p10upper + (1-p20upper)/p20upper/ratio)
    sigmalower1 <- sqrt((1-p10lower)/p10lower + (1-p20lower)/p20lower/ratio)
    aaupper <- sigmaupper1 * zalpha/sqrt(sigma001)
    aalower <- sigmalower1 * zalpha/sqrt(sigma001)
    bbupper <- (log(upper)-log(trueratio))/sqrt(sigma001)
    bblower <- (-log(lower)+log(trueratio))/sqrt(sigma001)
    power0 <- pnorm(bbupper*sqrt(n0)-aaupper) + pnorm(bblower*sqrt(n0)-aalower) - 1
    fuhao <- power0 > 1-beta
    while (step > 0.25) {
      hh <- ifelse(fuhao,-1,1)
      n0 <- n0 + step*hh
      power0 <- pnorm(bbupper*sqrt(n0)-aaupper) + pnorm(bblower*sqrt(n0)-aalower) - 1
      if((power0 > 1-beta) != fuhao){
        step <- step/2
        fuhao <- !fuhao
      }
    }
    n0 <- round(n0)
    power0 <- pnorm(bbupper*sqrt(n0)-aaupper) + pnorm(bblower*sqrt(n0)-aalower) - 1
    if(power0 < 1-beta) {
      n0 <- n0 + 1
      power0 <- pnorm(bbupper*sqrt(n0)-aaupper) + pnorm(bblower*sqrt(n0)-aalower) - 1
    }
    n1 <- n0
    n2 <- ceiling(n1 * ratio) 
  } else stop("method must be log or FM")
  data.frame(n1 = n1,n2 = n2,n = n1 + n2,power = power0)
}


#' @title samplesize_threearm_mvn
#' @description This function is to calculate the sample size for the three-arm equivalence trials
#' in term of one equivalence trial as in \code{\link{test_BERR}} and two superiority tests as 
#' in \code{\link{test_superRR}} when giving a fixed sample size allocation ratio for the three arms.
#' @param trueratio,p2,upper,lower,method,step,: same as \code{\link{samplesize_BERR}}.
#' @param alpha The overall type-I error. 
#' @param superside A character of either "one" or "two". "one" indicates the critical 
#' superiority test p-value is one side alpha. "two" indicates two side alpha.
#' @param p3 The assumed effective rate of Placebo drug.
#' @param ratio The sample size allocation ratio of Placebo drug. The sample sizes are set 
#' identical for T and R drug. So, the three sample size allocation ratio is n1:n2:n3 = 1:1:ratio.
#' @param verbose A logical value. If true, print the iterative calculation process.
#' @return n1, n2. n3 are the sample sizes of T drug, R drug and Placebo drug, respectively.
#' @example 
#' samplesize_threearm_mvn(0.95,0.5,0.4,1.25,0.8,1,margin = 1,alpha = 0.05,beta = 0.2,
#' method = "FM",verbose = F)
samplesize_threearm_mvn <- function(trueratio,p2,p3,upper,lower,ratio,margin = 1,
                                    beta = 0.2,alpha = 0.05,method = "log",
                                    verbose = F,superside = "two",step,...){
  #if(ratio > 1) stop("the ratio must not be larger than 1")
  p1 <- p2 * trueratio
  if (missing(lower)) lower <- 1/upper else lower <- lower
  r <- ratio
  ratio <- 1  # 等效性试验样本量设置相等
  a <- (1 + ratio)
  b <- -(upper * (1 + ratio * p2) + ratio + p1)
  c <- upper * (p1 + ratio * p2)
  p10upper <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p20upper <- p10upper/upper
  b <- -(lower * (1 + ratio * p2) + ratio + p1)
  c <- lower * (p1 + ratio * p2)
  p10lower <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p20lower <- p10lower/lower
    
  d0 <- margin == 1
  a <- (1 + r)
  b <- -(margin * (1 + r * p3) + r + p1)
  c <- margin * (p1 + r * p3)
  p102 <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p302 <- p102/margin
  p102[d0] <- (p1[d0] + r[d0] * p3[d0])/(1 + r[d0])
  p302[d0] <- p102[d0]
    
  a <- (1 + r)
  b <- -(margin * (1 + r * p3) + r + p2)
  c <- margin * (p2 + r * p3)
  p203 <- (-b - sqrt(b^2 - 4 * a * c))/2/a
  p303 <- p203/margin
  p203[d0] <- (p2[d0] + r[d0] * p3[d0])/(1 + r[d0])
  p303[d0] <- p203[d0]
  
  zalphabe <- qnorm(1-alpha)
  if(superside == "two") zalphasup <- qnorm(1-alpha/2)
    else if(superside == "one") zalphasup <- qnorm(1-alpha)
      else stop("superside must be 'one' or 'two'")
  
  n0 <- samplesize_BERR(trueratio,p2,upper,lower,alpha,beta,ratio,method = "log")
  n0 <- n0[1,1]
  if(missing(step)) step <- max(2,n0*0.1) else step <- step
  if(method == "FM"){
    #browser()
    sigma10upper <- sqrt(p10upper*(1-p10upper) + upper^2*p20upper*(1-p20upper)/ratio)
    sigma10lower <- sqrt(p20lower*(1-p20lower)/ratio + p10lower*(1-p10lower)/lower^2)
    sigma20 <- sqrt((1-p102)*p102 + (1-p302)*p302*margin^2/r)
    sigma30 <- sqrt((1-p203)*p203 + (1-p303)*p303*margin^2/r)
    sigma11upper <- sqrt(p1*(1-p1) + upper^2*p2*(1-p2)/ratio)
    sigma11lower <- sqrt(p2*(1-p2)/ratio + p1*(1-p1)/lower^2)
    sigma21 <- sqrt((1-p1)*p1 + (1-p3)*p3*margin^2/r)
    sigma31 <- sqrt((1-p2)*p2 + (1-p3)*p3*margin^2/r)
    r1u1l <- (-p1*(1-p1)/lower-p2*(1-p2)*upper/ratio)/sigma11upper/sigma11lower
    r1u2 <- p1*(1-p1)/sigma11upper/sigma21
    r1u3 <- -upper*p2*(1-p2)/sigma11upper/sigma31/ratio
    r1l2 <- -p1*(1-p1)/sigma11lower/sigma21/lower
    r1l3 <- p2*(1-p2)/sigma11lower/sigma31/lower/ratio
    r23 <- p3*(1-p3)*margin^2/sigma21/sigma31/r
    corx <- diag(1,nrow = 4)
    corx[upper.tri(corx)] <- c(r1u1l,r1u2,r1l2,r1u3,r1l3,r23)
    corx[lower.tri(corx)] <- c(r1u1l,r1u2,r1u3,r1l2,r1l3,r23)
    if(any(eigen(corx)[[1]]<1e-8)){
      warning("The non-positive-semidefinite  matrix has been transformed 
              as the nearest definite matrix ")
      corx <- as.matrix(Matrix::nearPD(corx,corr = T)[[1]]) 
    }
    eps1upper <- p1-upper*p2
    eps1lower <- p2-p1/lower
    eps2 <- p1-p3*margin
    eps3 <- p2-p3*margin
    
    limit1uperr <- -zalphabe*sigma10upper/sigma11upper-eps1upper/sigma11upper*sqrt(n0)
    limit1lower <- -zalphabe*sigma10lower/sigma11lower-eps1lower/sigma11lower*sqrt(n0)
    limit2 <- zalphasup*sigma20/sigma21-eps2/sigma21*sqrt(n0)
    limit3 <- zalphasup*sigma30/sigma31-eps3/sigma31*sqrt(n0)
    power0 <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf,limit2,limit3),
                               upper = c(limit1uperr,limit1lower,Inf,Inf),
                               mean = rep(0,4),corr = corx,...)
    fuhao <- power0 > 1-beta
    while (step > 0.25) {
      hh <- ifelse(fuhao,-1,1)
      n0 <- n0 + step*hh
      if(n0<1) stop("The initial step may be too large")
      limit1uperr <- -zalphabe*sigma10upper/sigma11upper-eps1upper/sigma11upper*sqrt(n0)
      limit1lower <- -zalphabe*sigma10lower/sigma11lower-eps1lower/sigma11lower*sqrt(n0)
      limit2 <- zalphasup*sigma20/sigma21-eps2/sigma21*sqrt(n0)
      limit3 <- zalphasup*sigma30/sigma31-eps3/sigma31*sqrt(n0)
      power0 <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf,limit2,limit3),
                                 upper = c(limit1uperr,limit1lower,Inf,Inf),
                                 mean = rep(0,4),corr = corx,...)
      if((power0 > 1-beta) != fuhao){
        step <- step/2
        fuhao <- !fuhao
      }
      if(verbose)  cat(n0,"  ",power0,"\n")
    }
    n0 <- round(n0)
    limit1uperr <- -zalphabe*sigma10upper/sigma11upper-eps1upper/sigma11upper*sqrt(n0)
    limit1lower <- -zalphabe*sigma10lower/sigma11lower-eps1lower/sigma11lower*sqrt(n0)
    limit2 <- zalphasup*sigma20/sigma21-eps2/sigma21*sqrt(n0)
    limit3 <- zalphasup*sigma30/sigma31-eps3/sigma31*sqrt(n0)
    power0 <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf,limit2,limit3),
                               upper = c(limit1uperr,limit1lower,Inf,Inf),
                               mean = rep(0,4),corr = corx,...)
    if(power0 < 1-beta) {
      n0 <- n0 + 1
      limit1uperr <- -zalphabe*sigma10upper/sigma11upper-eps1upper/sigma11upper*sqrt(n0)
      limit1lower <- -zalphabe*sigma10lower/sigma11lower-eps1lower/sigma11lower*sqrt(n0)
      limit2 <- zalphasup*sigma20/sigma21-eps2/sigma21*sqrt(n0)
      limit3 <- zalphasup*sigma30/sigma31-eps3/sigma31*sqrt(n0)
      power0 <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf,limit2,limit3),
                                 upper = c(limit1uperr,limit1lower,Inf,Inf),
                                 mean = rep(0,4),corr = corx,...)
    }
    
  } 
  else if(method == "log"){
    sigma10upper <- sqrt((1-p10upper)/p10upper + (1-p20upper)/p20upper/ratio)
    sigma10lower <- sqrt((1-p10lower)/p10lower + (1-p20lower)/p20lower/ratio)
    sigma20 <- sqrt((1-p102)/p102 + (1-p302)/p302/r)
    sigma30 <- sqrt((1-p203)/p203 + (1-p303)/p303/r)
    sigma11 <- sqrt((1-p1)/p1 + (1-p2)/p2/ratio)
    sigma21 <- sqrt((1-p1)/p1 + (1-p3)/p3/r)
    sigma31 <- sqrt((1-p2)/p2 + (1-p3)/p3/r)
    
    r12 <- (1-p1)/p1/sigma11/sigma21
    r13 <- -(1-p2)/p2/sigma11/sigma31/ratio
    r23 <- (1-p3)/p3/sigma21/sigma31/r
    corx123 <- diag(3)
    corx123[upper.tri(corx123)] <- c(r12,r13,r23)
    corx123[lower.tri(corx123)] <- c(r12,r13,r23)
    if(any(eigen(corx123)[[1]]<1e-8)){
      warning("The non-positive-semidefinite  matrix has been transformed 
              as the nearest definite matrix ")
      corx123 <- as.matrix(Matrix::nearPD(corx123,corr = T)[[1]]) 
    }
    logupper <-  log(upper)
    loglower <- log(lower)
    logmargin <- log(margin)
    logeps1 <- log(p1) - log(p2)
    logeps2 <- log(p1) - log(p3)
    logeps3 <- log(p2) - log(p3)
    
    limit1upper <- (logupper-logeps1)*sqrt(n0)/sigma11 - zalphabe*sigma10upper/sigma11
    limit1lower <- (loglower-logeps1)*sqrt(n0)/sigma11 + zalphabe*sigma10lower/sigma11
    limit2 <- (logmargin-logeps2)*sqrt(n0)/sigma21 + zalphasup*sigma20/sigma21
    limit3 <- (logmargin-logeps3)*sqrt(n0)/sigma31 + zalphasup*sigma30/sigma31
    power0 <- mvtnorm::pmvnorm(lower = c(limit1lower,limit2,limit3),
                               upper = c(limit1upper,Inf,Inf),
                               mean = rep(0,3),corr = corx123,...)
    fuhao <- power0 > 1-beta
    while (step > 0.25) {
      hh <- ifelse(fuhao,-1,1)
      n0 <- n0 + step*hh
      if(n0<1) stop("The initial step may be too large")
      limit1upper <- (logupper-logeps1)*sqrt(n0)/sigma11 - zalphabe*sigma10upper/sigma11
      limit1lower <- (loglower-logeps1)*sqrt(n0)/sigma11 + zalphabe*sigma10lower/sigma11
      limit2 <- (logmargin-logeps2)*sqrt(n0)/sigma21 + zalphasup*sigma20/sigma21
      limit3 <- (logmargin-logeps3)*sqrt(n0)/sigma31 + zalphasup*sigma30/sigma31
      if(limit1lower>=limit1upper) stop("The initial step may be too large")
      power0 <- mvtnorm::pmvnorm(lower = c(limit1lower,limit2,limit3),
                                 upper = c(limit1upper,Inf,Inf),
                                 mean = rep(0,3),corr = corx123,...)
      if((power0 > 1-beta) != fuhao){
        step <- step/2
        fuhao <- !fuhao
      }
      if(verbose)  cat(n0,"  ",power0,"\n")
    }
    n0 <- round(n0)
    limit1upper <- (logupper-logeps1)*sqrt(n0)/sigma11 - zalphabe*sigma10upper/sigma11
    limit1lower <- (loglower-logeps1)*sqrt(n0)/sigma11 + zalphabe*sigma10lower/sigma11
    limit2 <- (logmargin-logeps2)*sqrt(n0)/sigma21 + zalphasup*sigma20/sigma21
    limit3 <- (logmargin-logeps3)*sqrt(n0)/sigma31 + zalphasup*sigma30/sigma31
    power0 <- mvtnorm::pmvnorm(lower = c(limit1lower,limit2,limit3),
                               upper = c(limit1upper,Inf,Inf),
                               mean = rep(0,3),corr = corx123,...)
    if(power0 < 1-beta) {
      n0 <- n0 + 1
      limit1upper <- (logupper-logeps1)*sqrt(n0)/sigma11 - zalphabe*sigma10upper/sigma11
      limit1lower <- (loglower-logeps1)*sqrt(n0)/sigma11 + zalphabe*sigma10lower/sigma11
      limit2 <- (logmargin-logeps2)*sqrt(n0)/sigma21 + zalphasup*sigma20/sigma21
      limit3 <- (logmargin-logeps3)*sqrt(n0)/sigma31 + zalphasup*sigma30/sigma31
      power0 <- mvtnorm::pmvnorm(lower = c(limit1lower,limit2,limit3),
                                 upper = c(limit1upper,Inf,Inf),
                                 mean = rep(0,3),corr = corx123,...)
    }
  } else stop("method must be 'FM' or 'log'")
  
  data.frame(n1 = n0,n2 = n0,n3 = round(n0 * r),power = power0)
}


#' @title samplesize_threearm_mvn_noratio
#' @description This function adopt a grid search method to select an optimal sample size allocation ratio
#' for the three-arm equivalence trials. The following n1, n2. n3 indicate the sample sizes of T drug,
#'  R drug and Placebo drug, respectively.
#' @param grid The search length of grid. The search interval is detailed in details section.
#' @param maxgrid The maximal value of n3/(n1+n2+n3).
#' @param ... Seen in \code{\link{samplesize_threearm_mvn}}
#' @return 
#' @details The search interval of n3/(n1+n2+n3) is from grid to maxgrid with an increasement of grid.
#' @example 
#' xx <- samplesize_threearm_mvn_noratio(grid = 0.01,maxgrid = 0.5,trueratio = 0.95,
#' p2=0.5,p3=0.3,upper = 1.25,lower = 0.8,method = "FM",superside = "two",verbose = T)
#' xx$optimal_res
samplesize_threearm_mvn_noratio <- function(grid = 0.01,maxgrid = 0.5,verbose = T,...){
  res_all <- NULL
  grids <- seq(grid,maxgrid-grid,grid)
  for (i in 1:length(grids)) {
    ratio <- grids[i]/(1-grids[i])*2
    res0 <- samplesize_threearm_mvn(ratio = ratio,verbose = F,...)
    res_all <- rbind(res_all,unlist(c(ratio,sum(res0[1,1:3]),res0$power[1],res0[1,1:3])))
    if(verbose) cat("n1:n2:n3 =",round(ratio,3),"  ","total N =",sum(res0[1,1:3]),"\n")
  }
  res_all <- as.data.frame(res_all)
  colnames(res_all) <- c("ratio","N","Power","n1","n2","n3")
  rownames(res_all) <- paste0("1:1:",round(grids/(1-grids)*2,3))
  optimal_res <- res_all[which.min(res_all$N),]
  list(optimal_res = optimal_res,detail_res = res_all)
}



#############################example############################################
# calculate the sample size with p_R = 0.5, p_T = 0.95*0.5, p_P = 0.4. alpha = 0.05
# power = 0.8 with n1:n2:n3 = 1:1:x
samplesize_threearm_mvn(0.95,0.5,0.4,1.25,0.8,1,margin = 1,alpha = 0.05,beta = 0.2,
                        method = "FM",verbose = F)

# calculate the sample size with p_R = 0.5, p_T = 0.95*0.5, p_P = 0.3. alpha = 0.05
# power = 0.8 with n1:n2:n3 = 1:1:x
xx <- samplesize_threearm_mvn_noratio(grid = 0.01,maxgrid = 0.5,trueratio = 0.95,
      p2=0.5,p3=0.3,upper = 1.25,lower = 0.8,method = "FM",superside = "two",
      beta = 0.2,alpha = 0.05,verbose = T)
xx$optimal_res


