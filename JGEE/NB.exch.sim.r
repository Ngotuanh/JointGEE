library(MASS)
library(geepack)
NB.exch.sim <-function (){
  #Constant:
  N    <<- 260
  TT   <<- 20
  rho  <<- 0.9
  fixedPhi<<- 0
  NN = 1
  betaMat <- array(0, c(NN, 4))
  
  Beta <<- matrix(c(1, 1), c(2,1))
  Gama <<- matrix(c(0.5,0.1),c(2,1))
  for (k in 1:NN) {
    
    # Generating each subject, then combine by row:  
    X  <- array(0, c(1,2))
    W  <- array(0,c(1,2))
    Y  <- array(0, c(1,1))
    ID <- array(0, c(1,1))
    
    for (i in 1:N){# generate for each subject separately:
      Xi  <- genXi(i, N, TT)
      Mui <- calc_Mui(Xi, Beta)
      Wi  <- genWi(i, N, TT)
      Phi <- calc_Phi(Wi,Gama)
      Yi  <- genYi(Mui,Phi, TT, rho)
      
      IDi <- array(i, c(TT,1))
      X   <- rbind(X, Xi)
      W   <- rbind(W,Wi)
      Y   <- rbind(Y, Yi)
      ID  <- rbind(ID, IDi)
    }
    
    X   <- X[-1,] # delete first row
    Y   <- Y[-1,]
    W <- W[-1,]
    ID  <- ID[-1,]
    
    #Put data into frame and name the columns
    mydata <- data.frame(ID, X,W, Y)
    colnames(mydata) <- c('id', 'X0', 'X1','W0','W1', 'Y')
    
    fit <-  geese(  Y ~ 0 + X0 + X1 ,id = id,data = mydata,family = poisson,corstr = "exchangeable")  
    betaMat[k, 1] =  fit$beta[1]
    betaMat[k, 2] =  fit$beta[2]
    print(fit)
  }
  
  Res = SaveResult(betaMat)
  write.table(mydata, "mydata.csv",row.names=TRUE, 
              sep=',', col.names=TRUE,append=FALSE)
  
  return(mydata)
  
}


calc_Mui <- function(Xi, Beta){
 
  Mui <- (exp(Xi%*%Beta))
  return (Mui[1,1])
}

calc_Phi <- function(Wi, Gama){
  ones <- matrix(1,c(TT,1))
  Phi <- matrix(0,c(TT,1))
 #Phi <- mean(exp(Wi%*%Gama))
   Phi <- exp(Wi%*%Gama)
  Phi <- Phi 
  #print(dim(Phi))
  return (Phi[1,1])
}
# Generate Xi:
#=============
genXi <- function(i, N, TT){
  x0 <- gen_x0(i,N, TT)
  x1 <- gen_x1(TT)
  Xi <- cbind(x0, x1)
  return (Xi)
}


# Helpers for generating Xi:
gen_x0 <- function (i, N, TT){
  # fix x0 = -1, 0, 1
  if (i<=N/4) {x0 <- -1.0}
  else if (i<= (3*N)/4) {x0 <- 0.0}
  else x0 <- 1.0
  res = array(x0, c(TT,1))
  return (res)
}


gen_x1 <- function(TT){
  x <- uniform(TT)
  #x <- uniform2(TT)
  return(x)
}


uniform <- function(TT){
  # Generate TT values of the uniform distribution,
  # which are between min and max
  min <- 0.4
  max <- 0.9
  x <- runif(TT, min=min, max=max)
  x <- as.matrix(x)
  return (x)
}


uniform2 <- function(TT, mu=0.5, sd=0.01) { 
  # Generate TT values of the uniform distribution,
  # which have mean = mu, and standard deviation = sd
  x <- mu + sd * scale(rnorm(TT)) 
  return (x)
}

#Generate W given Gama
genWi <-function(i,N,TT){
  w0 <- gen_w0(i,N, TT)
  w1 <- gen_w1(TT)
  Wi <- cbind(w0, w1)
  return (Wi)
}
#helpers for generate Wi
gen_w0 <- function (i, N, TT){
  
  res = array(1, c(TT,1))
  return (res)
}


gen_w1 <- function(TT){
  w1 <- uniform(TT)
  #x <- uniform2(TT)
  return(w1)
}


# Generate Yi given Mui (scalar), TT, rho:
#================================
genYi <- function(Mui,Phi, TT, rho){
  #Yi <- poisson_exchangable(Mui, TT, rho)
  #Yi <- poisson_only(Mui, TT) 
  Yi <- NegBin_exchangable(Mui,Phi, TT, rho)
  Yi <- as.matrix(Yi)
  return (Yi)
}


poisson_exchangable <- function (Mui, TT, rho){
  lambda_eij <- Mui * (1 - sqrt(rho))
  eij <- rpois(TT, lambda_eij)
  
  lambda_yij <- Mui * sqrt(rho) 
  yij <- rpois(TT, lambda_yij)
  #Y0 = rpois(1,Mui)
  # yij <- binThin(sqrt(rho),Y0) 
  
  
  
  
  Yi <- yij + eij
  
  return (Yi)
}

NegBin_exchangable <- function(Mu_i,Phi_i,TT,rho){
  #Y_(t) = alpha * Y_0 + e_t
  Y_i = array(0,c(TT,1))
  # initial value Y_i0 = NB(mu_i,phi_i)
  # Y_i0 = rnbinom(1, size = 1/Phi_i, mu = Mu_i)
  # #get alpha = beta(shape1,shape2)
  # shape1 = sqrt(rho1) / Phi_i
  # shape2 = (1 - sqrt(rho1)) / Phi_i
  # alpha = rbeta(TT, shape1, shape2, ncp = 0)
  # alpha <- as.matrix(alpha)
  # d_t
  bij = rnbinom(TT, size= sqrt(rho)/Phi_i , mu= sqrt(rho)*Mu_i )
  bij <- as.matrix(bij)
  ei <- rnbinom(TT,size = (1-sqrt(rho))/Phi_i, mu = (1-sqrt(rho))*Mu_i)
  ei <- as.matrix(ei)
  # Get Y_1 first
  for (i in 1: TT)  {
    Y_i[i,1] = bij [i,1]+ ei[i,1]}
  return(Y_i)
  
}
binThin <- function(p,count){
  sum = 0
  for(i in 1: count){
    z = rbinom(1,1,p)
    sum = sum +z
  }
  return(sum)  
}

poisson_only <- function(Mui, TT){
  Yi <- rpois(TT, Mui)
}
##################### Result##############
SaveResult <- function(betaMat){
  Paramters <- data.frame(N,TT,rho,fixedPhi)
  
  Beta0.mean <- mean(betaMat[,1])
  Beta1.mean <- mean(betaMat[,2])
  Beta0.se = std(betaMat[,1])
  Beta1.se = std(betaMat[,2])
  
  Beta.res = data.frame(rbind(c(Beta0.mean,Beta1.mean),c(Beta0.se,Beta1.se)))
  Res = cbindPad(Paramters, Beta.res)
  
  # library(gridExtra)
  # pdf("data_output.pdf")#, height=8, width=8.5)
  # grid.table(Res)
  # dev.off()
  write.table(Res, "NB_sim.csv",row.names=FALSE, 
              sep=',', col.names=FALSE,append=TRUE)
  
  return(Res)
}

cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix("  ", mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}
std <- function(x) sd(x)/sqrt(length(x))

SaveResult <- function(betaMat){
  Paramters <- data.frame(N,TT,rho,fixedPhi)
  
  Beta0.mean <- mean(betaMat[,1])
  Beta1.mean <- mean(betaMat[,2])
  Beta0.se = std(betaMat[,1])
  Beta1.se = std(betaMat[,2])
  Beta0.compare = abs(Beta0.mean - Beta[1,1])
  Beta1.compare = abs(Beta1.mean - Beta[2,1])
  Beta.res = data.frame(rbind(c(Beta0.mean,Beta1.mean),c(Beta0.se,Beta1.se),c(Beta0.compare,Beta1.compare)))
  Res = cbindPad(Paramters, Beta.res)
  
  # library(gridExtra)
  # pdf("data_output.pdf")#, height=8, width=8.5)
  # grid.table(Res)
  # dev.off()
  write.table(Res, "ResultNegBin.csv",row.names=FALSE, 
              sep=',', col.names=FALSE,append=TRUE)
  
  return(Res)
}

