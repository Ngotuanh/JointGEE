
library(MASS)
library(corpcor)
JGEE1 <- function (data) {
  # TOTEST
  # get input params from data:
  input <- get_input(data)
  XMat <<- input$XMat
  WMat <<- input$WMat
  YMat <<- input$YMat
  p  <<- input$p
  q1 <<- input$q1
  N  <<- input$N
  TT <<- input$TT
  
  
  # init Beta, Gamma, epsilon, etc.:
  epsilon <- 0.0001
  tol <- 2
  Beta0 <- array(0, c(p , 1))
  Gama0 <- array(1, c(q1, 1))
  param_old <- array(0,c(p+q1,1))
  param_new <- array(0,c(p+q1,1))
  
  param_old = rbind( Beta0,Gama0)
  
  iter <- 1
  max_iter <- 300
  best_param <- param_old
  best_tol <- tol
  
  # iterating:
  while (tol > epsilon) {
    
    
    
    Mu  <- calc_mean(XMat, param_old[1:p,1])  # 2 columns: Mu$id, Mu$mean
    Phi <- calc_mean(WMat, param_old[(p+1) : (p+q1),1])  # same formula as Mu?
    r_res <- calc_r(YMat, Mu, Phi)
    s_res <- calc_s(YMat,Mu,Phi)
    #print(s_res)
    
    rho=calc_rho(r_res,s_res) # rho contains rho$rho1 and rho$rho2
    rho1 = rho$rho1
    rho2 = rho$rho2
    update_param = calculate_update_param(Mu, Phi,rho ) # sigma1 * sigma2
    
    param_new = param_old +update_param
    
    tol <- calc_tol(param_new, param_old)
    
    param_old = param_new # for next iteration if still tol > epsilon
    if (tol < best_tol){
      best_tol <- tol
      best_param <- param_new
    }
    
    if (iter > max_iter) {
      print ("Maximum number of iterations reached. Breaking the while loop...")
      print ("Best param so far is:")
      print (as.matrix(best_param))
      break
    }
    #print("tol")
    #print(tol)
    iter =iter +1
    res1 = t(as.matrix(param_old))
    res = cbind(iter,res1,tol)
    #print(res)
    #readline()
  }
  
  #result:
  if (tol <= epsilon){
  Beta <- as.matrix(best_param[1:p,1])
  Gama <- as.matrix(best_param[(p+1) : (p + q1),1])
  res = c(t(Beta),t(Gama))
  
  print(res)
  return (res)
  }
}


# Get input params from data: # ToDo & ToTest --> # TESTED if tested.
#============================
get_XMat <- function (data) {
  # ToTest
  
  keeps <- c("id", "X0", "X1")
  XMat <- data[, (names(data) %in% keeps)]
  XMat <- as.data.frame(XMat)
  return(XMat)
}

get_WMat <- function (data) {
  # ToTest
  
  keeps <- c("id", "W0", "W1")
  WMat <- data[, (names(data) %in% keeps)]
  WMat <- as.data.frame(WMat)
  return (WMat)
}

get_YMat <- function (data) {
  # ToTest
  
  keeps <- c("id", "Y")
  YMat <- data[, (names(data) %in% keeps)]
  YMat <- as.data.frame(YMat)
  
  return (YMat)
}


get_input <- function(data) {
  # TOTEST
  XMat <- get_XMat(data)
  WMat <- get_WMat(data)
  YMat <- get_YMat(data)
  ID <- data[, names(data) == "id"]
  ID <- as.matrix(ID)
  N  <- max(ID)
  TT <- nrow(ID) / N
  p  <- ncol(XMat) - 1 # exclude the id colum!
  q1 <- ncol(WMat) - 1 # exclude the id colum!
  
  output <-
    list(
      XMat = XMat,
      WMat = WMat,
      YMat = YMat,
      N = N,
      TT = TT,
      p = p,
      q1 = q1
    )
  #print (output)
  return (output)
} # TESTED?


# Caluclate tol given old and new param: # TOTEST
# ======================================
calc_tol <- function(param_new, param_old) {
  subtract = param_new - param_old  # same as update_param !
  subtract <- as.matrix(subtract)
  result = norm(subtract, type = "2") / norm(param_old, type = "2")
  return (result)
} # TESTED??


#Calculate Mu, Phi given XMat, WMat and Beta, Gama: # TOTEST
#==================================================
calc_mean <-  function(Mat, b) {
  X= Mat[,-1]
  MeanCol <- array(0,c(TT*N,1))
  MeanCol  = exp(as.matrix(X)%*%b)
  Mean = data.frame(Mat[,1], MeanCol)
  names(Mean) = c("id","Mean")
  return (Mean)
}



get_i <-
  function(XMat, i) {
    # Yi = get_i(YMat, i); Wi = get_i(WMat, i); Mui = get_i(Mu, i); Phii = get_i(Phi, i)
    Xi <- XMat[XMat$id == i,]
    Xi <- as.matrix(Xi[,-1]) # remove id column
    return (Xi)
  }



# Calculate rho$rho1, rho$rho2 given YMat, Mu, Phi:
# =================================================
# Note: XMat, YMat, WMat, Mu, Phi all have 'id' column
# To get values as matrix, use as example: Mui = get_i(Mu, id)

calc_rho <- function(r_res,s_res) {
  # TOTEST
  rho1 <- calc_rho1(r_res)
  rho2 <- calc_rho2(s_res)
  rho <- list(rho1 = rho1, rho2 = rho2)
  return (rho)
}

# Helper functions for calculating rho: # TODO & TOTEST
# -------------------------------------

calc_r <- function(Mat, Mu, Phi) {
  Y = as.matrix(Mat[,-1])
  MyMean = as.matrix(Mu[,-1])
  MyPhi = as.matrix(Phi[,-1])
  rcol = matrix(0,c(TT*N,1))
  rcol = Y - MyMean
  rcol = rcol/sqrt(MyMean + MyPhi * MyMean^2)
  rvec = data.frame(Mat[,1],rcol)
  names(rvec) =c("id","res")
  return(rvec)
}

rho_sum <- function(r_i) {
  # r_i is a vector
  sumsum = 0
  for (j in 2:TT)
    for (k in 1:(j - 1)) {
      sumsum = sumsum + r_i[j, 1] * r_i[k, 1]
    }
  return(sumsum)
}

calc_rho1 <- function(r_res) {
  s <- 0
  r_i <- array(0, c(TT, 1))
  
  #TODO, calc_r?
  for (i in 1:N) {
    
    r_i <- get_i(r_res,i)
    
    s = s + rho_sum(r_i)
  }
  
  # result:
  constant = ((N/2 * (TT - 1) * TT)) - p
  rho1 = s / constant
  return (rho1)
}

calc_s <- function(Mat, Mu,Phi) {
  Y <- as.matrix(Mat[,-1])
  Y2 = as.matrix(Y*Y)
  
  MyMu = as.matrix(Mu[,-1])
  MyPhi= as.matrix(Phi[,-1])
  
  s <- array(0,c(TT*N,1))
  varVec <-array(0,c(N*TT,1))
  m = as.matrix(mean_YtYt(MyMu,MyPhi))
  s = Y2 - m
  for (i in 1: TT*N){
    var= Var_YtYt(MyMu[i,1],MyPhi[i,1])
    varVec = sqrt(var)
  }
  s = s/(varVec)
  svec = data.frame(Mat[,1],s)
  names(svec) = c("id","sres")
  
  return(svec)
}

calc_rho2 <- function(s_res) {
  # Mu2 = mij
  s <- 0
  s_i <- array(0, c(TT, 1))
  #TODO, calc_s?
  for (i in 1:N) {
    s_i <- get_i(s_res,i)
    s = s + rho_sum(s_i)
  }
  # result:
  constant = ((N/2  * (TT - 1) * TT))- p - q1
  rho2 = s / constant
  return (rho2)
}




# Calculate the update_param value by calculating D,V,r... # TODO # TOTEST
# =========================================================

calculate_update_param <-
  function (Mu, Phi, rho) {
    # Gama = param[1], Beta = param[2] # TO TEST
    Sigma1 <- calc_sigma1(Mu, Phi, rho)
    Sigma2 <- calc_sigma2(Mu, Phi, rho)
    Sigma1_inv = pseudoinverse(Sigma1)
    update_param = Sigma1_inv %*% Sigma2
    return (update_param)
  }


# Helper functions for calculate update param: # TODO # TOTEST
# --------------------------------------------
calc_sigma1 <- function(Mu, Phi, rho) {
  sumMat <- array(0, c(p + q1, p + q1))
  
  for (i in 1:N) {
    Xi <- get_i(XMat, i)
    Wi <- get_i(WMat, i)
    Mui <- get_i(Mu, i)
    Phi_i <- get_i(Phi, i)
    
    Di <- calc_Di (Xi, Wi, Mui, Phi_i)
    Vi <- calc_Vi(Mui, Phi, rho)
    
    #Vi <- calc_Vi_2(Mui,Phi_i,rho)
    sumMat = sumMat + calc_sigma1_i(Di, Vi)
  }
  
  return(sumMat)
}

calc_sigma2 <- function(Mu, Phi, rho) {
  sumMat <- array(0, c(p + q1, 1))
  RMat <- calc_RMat(Mu,Phi)
  for (i in 1:N) {
    Xi <- get_i(XMat, i)
    Wi <- get_i(WMat, i)
    Yi <- get_i(YMat, i)
    Mui <- get_i(Mu, i)
    Phi_i <- get_i(Phi, i)
    R_i <- get_i(RMat,i)
    
    Di <- calc_Di (Xi, Wi, Mui, Phi_i)
    Vi <- calc_Vi(Mui, Phi_i, rho)
    sumMat = sumMat + calc_sigma2_i(Di, Vi, R_i)
  }
  return(sumMat)
  
}

calc_sigma2_i <- function(Di, Vi, Ri) {
  #sigma1 is matrix(p+q1 )x(p+q)
  mat <- array(0, c(p + q1, 1))
  V_inv = pseudoinverse(Vi)
  mat = t(Di) %*% V_inv %*% Ri
  
  return(mat)
}

calc_sigma1_i <- function(Di, Vi) {
  #sigma1 is matrix(p+q1 )x(p+q)
  mat <- array(0, c(p + q1, p + q1))
  V_inv = pseudoinverse(Vi)
  mat = t(Di) %*% V_inv %*% Di
  
  return(mat)
}


calc_Di <- function(Xi, Wi, Mui, Phi_i) {
  
  #TODO
  Di <- array(0, c(2 * TT, p + q1))
  Di1 <- calc_Di1(Mui, Xi)
  Di21 <- calc_Di21(Xi, Mui, Phi_i)
  Di2 <- calc_Di2(Wi, Mui, Phi_i)
  ZeroMat <- array(0, c(TT, q1))
  #Di21 = ZeroMat
  Di <- rbind(Di1, Di21)
  temp <- rbind(ZeroMat, Di2)
  
  Di <- cbind(Di, temp)
  
  return(Di)
}

calc_Di1 <- function(Mui, Xi) {
  mat <- array(0, c(TT, p))
  for (k in 1:TT) {
    for (j in 1:p) {
      mat[k, j] = Xi[k, j] * Mui[k, 1]
    }
  }
  
  return(mat)
}
calc_Di21 <- function(Xi, Mui, Phi_i) {
  mat <- array(0, c(TT, p))
  for (k in 1:TT) {
    for (j in 1:p) {
      mat[k, j] = Xi[k, j] * Mui[k, 1] * (1 + 2 * (Phi_i[k, 1] + 1) * Mui[k, 1])
    }
  }
  mat <- array(0, c(TT, p))
  
  return(mat)
}
calc_Di2 <- function(Wi, Mui, Phi_i) {
  mat <- array(0, c(TT, q1))
  for (k in 1:TT) {
    for (j in 1:q1) {
      mat[k, j] = Wi[k, j] * Mui[k, 1] ^ 2 * (Phi_i[k, 1])
    }
  }
  return(mat)
}

######



calc_Vi <- function (Mu, Phi, rho) {
  Vi1 <- calc_Vi1(Mu, Phi, rho$rho1)
  Vi2 <- calc_Vi2(Mu, Phi, rho$rho2)
  #Cov_i <- calc_Cov_i(Mu, Phi, rho$rho1)
  Cov_i <-array(0,c(TT,TT))
  
  Col1 <- rbind(Vi1, Cov_i)
  Col2 <- rbind(Cov_i, Vi2)
  mat <- cbind(Col1, Col2)
  
  return(mat)
  
}
calc_Vi1 <- function (Mu, Phi, rho1) {
  mat <- array(0, c(TT, TT))
  A <- array(0, c(TT, TT))
  
  for (t in 1:TT) {
    A[t, t] = Var_Yt(Mu[t, 1], Phi[t, 1])
  }
  R_rho1 <- array(rho1, c(TT, TT))
  diag(R_rho1) = 1
  mat = A^(1/2) %*% R_rho1 %*% (A^(1/2))
  return(mat)
  #TODO
}
calc_Vi2 <- function (Mu, Phi, rho2) {
  mat <- array(0, c(TT, TT))
  H <- array(0, c(TT, TT))
  
  for (t in 1:TT) {
    H[t, t] = Var_YtYt(Mu[t, 1], Phi[t, 1])
  }
  R_rho2 <- array(rho2, c(TT, TT))
  diag(R_rho2) = 1
  mat = H ^(1/2)%*%R_rho2 %*% H^(1/2)
  return(mat)
  #TODO
}
calc_Cov_i <- function (Mu, Phi, rho1) {
  mat <- array(0, c(TT, TT))
  for (j in 1:TT)
    for (k in j:TT)
    {
      if (j == k) {
        mat[j, k] = Cov_YtYtYt(Mu[j, 1], Phi[j, 1])
      }
      else{
        mat[j, k] = Cov_YtYsYs(Mu[j, 1], Phi[j, 1], Mu[k, 1], Phi[k, 1], rho1)
        mat[k, j] = mat[j, k]
      }
      
    }
  #TODO
  mat <- array(0, c(TT, TT))
  
  return(mat)
}

calc_Ri <- function (Mat, Mu,Phi) {
  # Ri = Yi - mui, Yi^2 - mij
  mat <- array(0, c(2 * TT, 1))
  for (i in 1:TT) {
    mat[i, 1] = Mat[i, 1] - Mu[i, 1]
    mij = mean_YtYt(Mu[i, 1], Phi[i, 1])
    mat[i + TT, 1] = Mat[i, 1] * Mat[i, 1] - mij
  }
  return(mat)
  
  #TODO
}

calc_RMat <- function(Mu,Phi){
  IDCol <- as.matrix(Mu[,1])
  Y <- as.matrix(YMat[,-1])
  Y2 = as.matrix(Y^2)
  
  MyMean <- as.matrix(Mu[,-1])
  MyPhi <- as.matrix(Phi[,-1])
  MeanY2 <- mean_YtYt(MyMean, MyPhi)
  
  FiPart <- array(0,array(N*TT,1))
  SePart <- array(0,array(N*TT,1))
  
  FiPart = Y - MyMean
  SePart = Y2 - MeanY2
  
  RPart = rbind(FiPart,SePart)
  RMat = data.frame(IDCol, RPart)
  names(RMat) = c("id","R")
  
  return(RMat)
}
####### Distribution : mean, variance, covariance
Var_Yt <- function(mean_t, phi_t) {
  output = mean_t + phi_t * mean_t
  return(output)
}
Var_YtYt <- function(mean_t, phi_t) {
  output =  mean_t +
    +mean_t ^ 2 * (6 + 7 * phi_t)+
    + mean_t ^ 3 * (4 + 16 * phi_t + 12 * phi_t ^ 2)+
    + mean_t ^ 4 * (4 + 10 * phi_t ^ 2  + 6 * phi_t ^ 3)
  
  return(output)
}
mean_YtYt <- function(mean_t, phi_t) {
  output= mean_t + (phi_t + 1) * mean_t^2
  return(output)
  
}
Cov_YtYtYt <- function(mean_t, phi_t) {
  #covariance cov(Yt,YtYt)
  output= mean_t * (1 + (2 + 3 * phi_t) * mean_t + 2 * phi_t * (1 + phi_t) * mean_t ^2)
  
  return(output)
  
}
Cov_Yt2Ys2 <- function(mean_t,phi_t,mean_s,phi_s,rho1){
  std2_t = Var_Yt(mean_t,phi_t)
  std2_s = Var_Yt(mean_s,phi_s)  
  m_t = mean_YtYt(mean_t,phi_t)
  m_s = mean_YtYt(mean_s,phi_s)
  output = 3*rho1^2 *std2_s*std2_t - 11* mean_t^2 *mean_s^2 + 
    + rho1 * (mean_s^2 *std2_t + 4*mean_t*mean_s*sqrt(std2_s*std2_t) +
                + mean_t^2 *std2_s) - m_t - m_s
  
  return(output)  
}
Cov_YtYsYs <- function(mean_t, phi_t, mean_s, phi_s, rho1) {
  # covariance cov(Yt, YsYs)
  std_Yt = sqrt(Var_Yt(mean_t, phi_t))
  std_Ys = sqrt(Var_Yt(mean_s, phi_s))
  
  output = 2 * rho1 * std_Yt * std_Ys * mean_s
  
  return(output)
}

