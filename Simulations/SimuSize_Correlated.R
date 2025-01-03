library(psych)
library(mvnfast)
library(Rcpp)
library(truncnorm)
library(glmnet)
library(scalreg)
library(doParallel)
library(foreach)

Rv = c(12, 16)
Cv = c(12, 16)

for (region in 1:length(Rv))
{
  set.seed(37)
  load('CMFD_Hebei.RData')
  #State = State[, -c(2, 5), , ]
  N = dim(PM25_Hebei)[1]; d1 = dim(PM25_Hebei)[2]; d2 = dim(PM25_Hebei)[3]
  
  n = 100
  selsample = sample(1:N, size = n, replace = F)
  
  TEMP = TEMP_Hebei[selsample, , ]
  WIND = WIND_Hebei[selsample, , ]
  SHUM = SHUM_Hebei[selsample, , ]
  PRES = PRES_Hebei[selsample, , ]
  
  R = Rv[region]; C = Cv[region]; p = 4
  
  totallen = p + R*C
  
  X = array(0, dim = c(n, R, C, p))
  
  s1 = floor(d1/2) - floor(R/2); e1 = s1 + R - 1
  s2 = floor(d2/2) - floor(C/2); e2 = s2 + C - 1
  
  X[ , , , 1] = 1
  X[ , , , 2] = WIND[, s1:e1, s2:e2]
  X[ , , , 3] = TEMP[, s1:e1, s2:e2] - 273.15
  X[ , , , 4] = PRES[, s1:e1, s2:e2]/1000
  
  CovNormal = matrix(0, R*C, R*C)
  for (j in 1:(R*C))
  {
    for (k in 1:(R*C))
    {
      CovNormal[j, k] = 0.5^abs(j - k)
    }
  }
  
  M = array(0, dim = c(n, R, C))
  for (i in 1:n)
  {
    tempNormal = matrix(rmvn(1, mu = rep(0, R*C), CovNormal), R, C)
    M[i, , ] = sign(tempNormal)
  }
  
  #M = array(sample(x = c(-1, 1), size = n*R*C, replace = T), dim = c(n, R, C))
  A = array(0, dim = c(n, R, C, R*C - 1))
  # 
  # Avector = t(apply(M, 1, as.vector))
  # for (r in 1:R)
  # {
  #   for (c in 1:C)
  #   {
  #     pos = (c - 1)*R + r
  #     
  #     A[, r, c, ] = Avector[, -pos]
  #   }
  # }
  # 
  # rm(Avector); gc()
  
  #####################################################################################
  
  Btrue = array(rnorm(R*C*p, 2, 1), dim = c(R, C, p))
  
  Ymid = array(0, dim = c(n, R, C))
  
  for (i in 1:n)
  {
    Ymid[i, , ] = apply(X[i, , , ]*Btrue, c(1, 2), sum)
  }
  ####################################################################################################
  lowp = 4; DEsignal = 0.01; #InFsignal = c(0, 0.002, 0.004, 0.006, 0.008, 0.01)
  InFsignal = c(0)
  
  DEMat = abs(apply(Ymid, c(2, 3), mean))
  DESVD = svd(DEMat)
  
  DEu = DESVD$u; DEv = DESVD$v; DEd = DESVD$d
  DEd[(lowp + 1):C] = 0
  
  Lbase = DEu %*% diag(DEd) %*% t(DEv)
  Ltrue = DEsignal*Lbase
  
  Sbase = array(0, dim = c(R, C, R*C - 1))
  startpoint = list()
  endpoint = list()
  tune = sqrt(log(1.8*totallen)/n)
  
  for (r in 1:R)
  {
    for (c in 1:C)
    {
      SrcTMat  = matrix(0, R, C)
      Matdis = matrix(0, R, C)
      for (j in 1:R)
      {
        for (k in 1:C)
        {
          Matdis[j, k] = (j - r)^2 + (k - c)^2
        }
      }
      
      disvec = sort(unique(as.vector(Matdis)), decreasing = F)[-1]
      clunum = length(disvec)
      
      Arc = NULL
      for (u in 1:clunum)
      {
        indexu = which(Matdis == disvec[u], arr.ind = T)
        
        for (v in 1:nrow(indexu))
        {
          Arc = cbind(Arc, M[, indexu[v, 1], indexu[v, 2]])
        }
      }
      
      A[, r, c, ] = Arc
      
      endp = rep(0, clunum)
      
      end = 0
      for (t in 1:clunum)
      {
        end = end + length(which(Matdis == disvec[t]))
        endp[t] = end
      }
      
      startp = c(1, endp[-length(endp)] + 1)
      
      regnum = sample(c(0, 1, 2), size = 1)
      
      if (regnum != 0)
      {
        regposind = sample(disvec[-c(clunum, clunum - 1)], regnum, replace = F)
        
        for (v in 1:regnum)
        {
          indexv = which(Matdis == regposind[v])
          numv = length(indexv)
          #meanv = abs(mean(Ymid[, r, c]))/numv/regnum
          meanv = abs(mean(Ymid[, r, c]))
          SrcTMat[indexv] = rtruncnorm(numv, a = 0, b = Inf, mean = meanv, sd = 1)
        }
      }
      
      Srctrue = NULL
      for (u in 1:clunum)
      {
        indexu = which(Matdis == disvec[u])
        Srctrue = c(Srctrue, SrcTMat[indexu])      
      }
      
      Sbase[r, c, ] = Srctrue
      startpoint = c(startpoint, list(startp))
      endpoint = c(endpoint, list(endp))
    }
  }
  
  rm(Arc); gc()
  
  StrueList = list()
  # LtrueList = list()
  for (num in 1:length(InFsignal))
  {
    Strue = InFsignal[num]*Sbase
    # Ltrue = DEsignal[num]*Lbase
    StrueList = c(StrueList, list(Strue))
    # LtrueList = c(LtrueList, list(Ltrue))
  }
  
  rm(list = ls()[-c(which(ls() == 'X'), which(ls() == 'A'), which(ls() == 'M'), 
                    which(ls() == 'Btrue'), which(ls() == 'Ltrue'), which(ls() == 'StrueList'), which(ls() == 'startpoint'), which(ls() == 'endpoint'), 
                    which(ls() == 'Rv'), which(ls() == 'Cv'), which(ls() == 'region'), which(ls() == 'tune'))])
  #############################################################################################################################
  #############################################################################################################################
  
  sellambdaL = rep(0, length(StrueList))
  source('BiRS_TEM.R')
  
  for (set in 1:length(StrueList))
  {
    Strue = StrueList[[set]]
    
    DIM = dim(X); n = DIM[1]; R = DIM[2]; C = DIM[3]; p = DIM[4]
    
    Delta = matrix(runif(R*C, 0.5, 1), R, C)
    Err = array(0, dim = c(n, R, C))
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        Err[, r, c] = rnorm(n, 0, Delta[r, c])
      }
    }
    
    Y = array(0, dim = c(n, R, C))
    
    for (i in 1:n)
    {
      Y[i, , ] = apply(X[i, , , ]*Btrue, c(1, 2), sum) + M[i, , ]*Ltrue + apply(A[i, , , ]*Strue, c(1, 2), sum) + Err[i, , ]
    }
    
    HSigma = matrix(0, R, C)
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        XMArc = cbind(X[, r, c, ], M[, r, c], A[, r, c, ])
        
        Estrc = scalreg(XMArc, Y[, r, c], tune)
        HSigma[r, c] = Estrc$hsigma
      }
    }
    LambdaMat = sqrt(2.1)*HSigma*sqrt(log(R*C - 1)/n)
    
    beta0 = array(0, dim = c(R, C, p))
    L0 = matrix(0, R, C)
    S0 = array(0, dim = c(R, C, R*C - 1))
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        Yrc = Y[ , r, c]
        XMArc = cbind(X[, r, c, -1], A[, r, c, ], M[, r, c])
        
        #bestlambda = cv.glmnet(x = XMArc, y = Yrc, alpha = 1)$lambda.min
        
        EstOLSrc = coef(glmnet(x = XMArc, y = Yrc, alpha = 1, lambda = LambdaMat[r, c]))
        
        beta0[r, c, ] = EstOLSrc[1:p]
        S0[r, c, ] = EstOLSrc[(p + 1):(p + R*C - 1)]
        L0[r, c] = EstOLSrc[p + R*C]
      }
    }
    
    tol = 1e-04; iter_max = 100; deltaL = 2
    
    lambdaLv = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
    numcore = min(length(lambdaLv), 50)
    #ratio = c(0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5)
    #lambdaCandidate = cbind(rep(ratio, rep(length(lambdaLv), length(ratio))), rep(lambdaLv, length(ratio)))
    
    foldlen = floor(n/5)
    permu = sample(1:n, n, replace = F)
    
    CV = function(r)
    {
      set.seed(37)
      lambdaLr = lambdaLv[r]
      
      ResdFold = rep(0, 5)
      ResdSpar = rep(0, 5)
      for (fold in 1:5)
      {
        if (fold != 5)
        {
          indf = permu[((fold - 1)*foldlen + 1):(fold*foldlen)]
        }
        if (fold == 5)
        {
          indf = permu[((fold - 1)*foldlen + 1):n]
        }
        
        Ytrain = Y[-indf, , ]; Ytest = Y[indf, , ]
        Xtrain = X[-indf, , , ]; Xtest = X[indf, , , ]
        Mtrain = M[-indf, , ]; Mtest = M[indf, , ]
        Atrain = A[-indf, , , ]; Atest = A[indf, , , ]
        
        TrainEst = LAFISTA(Ytrain, Xtrain, Mtrain, Atrain, beta0, L0, S0, deltaL, lambdaLr, LambdaMat, tol, iter_max, method = 'glmnet')
        Strain = TrainEst$Shat; Ltrain = TrainEst$Lhat; betatrain = TrainEst$betahat
        
        ErrMat = array(0, dim = c(length(indf), R, C))
        for (j in 1:length(indf))
        {
          ErrMat[j, , ] = Ytest[j, , ] - apply(Xtest[j, , , ]*betatrain, c(1, 2), sum) - Mtest[j, , ]*Ltrain - apply(Atest[j, , , ]*Strain, c(1, 2), sum)
        }
        
        ResdFold[fold] = sum(ErrMat^2)
        
        print(paste0('Fold', fold, '-----------complete'))
      }
      
      Resd = mean(ResdFold)
      
      return(Resd)
    }
    
    cl = makeCluster(numcore)
    registerDoParallel(cl)
    
    RESD = foreach(r = 1:length(lambdaLv), .packages = c("mvnfast", "glmnet", "scalreg"), .combine = 'c') %dopar% CV(r)
    
    stopImplicitCluster()
    stopCluster(cl)
    
    sellambdaL[set] = lambdaLv[which.min(RESD)]
  }
  
  rm(list = ls()[-c(which(ls() == 'X'), which(ls() == 'A'), which(ls() == 'M'),
                    which(ls() == 'Btrue'), which(ls() == 'Ltrue'), which(ls() == 'StrueList'), which(ls() == 'startpoint'), which(ls() == 'endpoint'), which(ls() == 'sellambdaL'),
                    which(ls() == 'Rv'), which(ls() == 'Cv'), which(ls() == 'region'), which(ls() == 'tune'))])
  #########################################################################################################################################################################################
  
  #selratio = matrix(1, 5, 2)
  for (set in 1:length(StrueList))
  {
    Strue = StrueList[[set]]
    R = dim(Strue)[1]; C = dim(Strue)[2]
    
    source('BiRS_TEM.R')
    source('Knockoff_TEM.R')
    
    Simu = function(s)
    {
      DIM = dim(X); n = DIM[1]; R = DIM[2]; C = DIM[3]; p = DIM[4]
      
      set.seed(s)
      
      Delta = matrix(runif(R*C, 0.5, 1), R, C)
      Err = array(0, dim = c(n, R, C))
      for (r in 1:R)
      {
        for (c in 1:C)
        {
          Err[, r, c] = rnorm(n, 0, Delta[r, c])
        }
      }
      
      Y = array(0, dim = c(n, R, C))
      
      for (i in 1:n)
      {
        Y[i, , ] = apply(X[i, , , ]*Btrue, c(1, 2), sum) + M[i, , ]*Ltrue + apply(A[i, , , ]*Strue, c(1, 2), sum) + Err[i, , ]
      }
      
      LambdaMat = matrix(0, R, C)
      HSigma = matrix(0, R, C)
      for (r in 1:R)
      {
        for (c in 1:C)
        {
          XMArc = cbind(X[, r, c, ], M[, r, c], A[, r, c, ])
          
          Estrc = scalreg(XMArc, Y[, r, c], tune)
          HSigma[r, c] = Estrc$hsigma
        }
      }
      
      LambdaMat = sqrt(2.1)*HSigma*sqrt(log(R*C - 1)/n)
      
      # effectratio1 = matrix(0, R, C)
      # for (r in 1:R)
      # {
      #   for (c in 1:C)
      #   {
      #     effectratio1[r, c] = abs((sum(Strue[r, c, ])/sum(Strue[r, c, ]!=0))/mean(Y[, r, c]))
      #   }
      # }
      # 
      # effectratio2 = matrix(0, R, C)
      # for (r in 1:R)
      # {
      #   for (c in 1:C)
      #   {
      #     effectratio2[r, c] = abs(Ltrue[r, c]/mean(Y[, r, c]))
      #   }
      # }
      
      #LambdaMat1 = matrix(0, R, C)
      beta0 = array(0, dim = c(R, C, p))
      L0 = matrix(0, R, C)
      S0 = array(0, dim = c(R, C, R*C - 1))
      for (r in 1:R)
      {
        for (c in 1:C)
        {
          Yrc = Y[ , r, c]
          XMArc = cbind(X[, r, c, -1], A[, r, c, ], M[, r, c])
          
          #bestlambda = cv.glmnet(x = XMArc, y = Yrc, alpha = 1)$lambda.min
          
          EstOLSrc = coef(glmnet(x = XMArc, y = Yrc, alpha = 1, lambda = LambdaMat[r, c]))
          
          #LambdaMat1[r, c] = bestlambda
          
          beta0[r, c, ] = EstOLSrc[1:p]
          S0[r, c, ] = EstOLSrc[(p + 1):(p + R*C - 1)]
          L0[r, c] = EstOLSrc[p + R*C]
        }
      }
      
      MB = 1000; alpha = 0.05; tol = 1e-04; iter_max = 100; deltaL = 2
      lambdaL = sellambdaL[set]
      
      #LambdaMat = LambdaMat*selratio[set, 1]
      resBoot = BiRS_IF(Y, X, M, A, startpoint, endpoint, beta0, L0, S0, deltaL, lambdaL, LambdaMat, HSigma, MB, alpha, tol, iter_max, ReMax = 10, method = 'Boot')
      resDebias = BiRS_IF(Y, X, M, A, startpoint, endpoint, beta0, L0, S0, deltaL, lambdaL, LambdaMat, HSigma, MB, alpha, tol, iter_max, ReMax = 10, method = 'debiased')
      #resLasso = BiRS_IF(Y, X, M, A, startpoint, endpoint, beta0 = NULL, L0 = NULL, S0 = NULL, deltaL = NULL, lambdaL = NULL, LambdaMat, HSigma, MB, alpha, tol = NULL, iter_max = NULL, ReMax = 10, method = 'Lasso')
      #resKnockoff = KnockoffSelect(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL, LambdaMat, tol, iter_max)
      #res1 = LAFISTA(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL = 1, LambdaMat)
      
      return(list(resBoot = resBoot, resDebias = resDebias))
    }
    
    nsimu = 1000
    
    cl = makeCluster(100)
    registerDoParallel(cl)
    
    RES = foreach(s = 1:nsimu, .packages = c("mvnfast", "glmnet", "scalreg", "MASS", "selectiveInference", "knockoff")) %dopar% Simu(s)
    
    stopImplicitCluster()
    stopCluster(cl)
    
    save(list = c('RES', 'Strue', 'sellambdaL'), file = paste0('Results_Correlated/R_', R, 'C_', C, 'Signal', set, '.RData'))
    rm(RES); gc()
  }
  
  rm(list = ls()[-c(which(ls() == 'Rv'), which(ls() == 'Cv'), which(ls() == 'region'))])
}