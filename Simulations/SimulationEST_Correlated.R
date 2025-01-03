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
  tune = sqrt(log(totallen)/n)
  
  # L1 = function(x)
  # {
  #   return((1/pnorm(1 - x/totallen))^4 + 2*(1/pnorm(1 - x/totallen))^2 - x)
  # }
  # 
  # tunnpara = uniroot(L1, interval = c(1, 10))$root
  # tune = sqrt(2)*n^(-1/2)/pnorm(1 - tunnpara/totallen)
  
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
  lowp = 4; DEsignal = 0.01; InFsignal = c(0.002, 0.004, 0.006, 0.008, 0.01)
  
  DEMat = abs(apply(Ymid, c(2, 3), mean))
  DESVD = svd(DEMat)
  
  DEu = DESVD$u; DEv = DESVD$v; DEd = DESVD$d
  DEd[(lowp + 1):C] = 0
  
  Lbase = DEu %*% diag(DEd) %*% t(DEv)
  Ltrue = DEsignal*Lbase
  
  Sbase = array(0, dim = c(R, C, R*C - 1))
  startpoint = list()
  endpoint = list()
  
  num_neighbor = matrix(0, R, C)
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
      
      num_neighbor[r, c] = sum(Matdis <= 2) - 1
      
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
  
  rm(list = ls()[-c(which(ls() == 'X'), which(ls() == 'A'), which(ls() == 'M'), which(ls() == 'num_neighbor'),
                    which(ls() == 'Btrue'), which(ls() == 'Ltrue'), which(ls() == 'StrueList'), which(ls() == 'Rv'), which(ls() == 'Cv'), which(ls() == 'region'))])
  #############################################################################################################################
  #############################################################################################################################
  
  postEst = function(Yrc, Xrc, Mrc, Arc, indrc)
  {
    Zrc = cbind(Xrc, Mrc, Arc[, indrc])
    
    d = ncol(Xrc); p = ncol(Zrc)
    Estrc = lm(Yrc ~ Zrc - 1)$coefficients[(d+1):p]
    
    return(Estrc)
  }
  
  for (set in 1:length(StrueList))
  {
    Strue = StrueList[[set]]
    R = dim(Strue)[1]; C = dim(Strue)[2]
    
    filename = paste0('Results_Correlated/R_', R, 'C_', C, 'Signal', set + 1, '.RData')
    load(filename)
    
    ESTsimu = function(s)
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
      
      TE_Boot = matrix(0, R, C); TE_Neighbor = matrix(0, R, C); TE_True = matrix(0, R, C)
      for (r in 1:R)
      {
        for (c in 1:C)
        {
          indrc_est = NULL
          
          if (!is.character(RES[[s]]$resBoot))
          {
            indrc_est = which(RES[[s]]$resBoot$Signal[r, c, ] != 0)
          }
          
          Yrc = Y[, r, c]; Xrc = X[, r, c, ]; Mrc = M[, r, c]; Arc = A[, r, c, ]
          indrc_neighbor = 1:num_neighbor[r, c]
          
          Bootrc = postEst(Yrc, Xrc, Mrc, Arc, indrc_est)
          Neighborrc = postEst(Yrc, Xrc, Mrc, Arc, indrc_neighbor)
          
          TE_Boot[r, c] = sum(Bootrc); TE_Neighbor[r, c] = sum(Neighborrc); TE_True[r, c] = Ltrue[r, c] + sum(Strue[r, c, ])
        }
      }
      
      ATE_Boot = mean(TE_Boot); ATE_Neighbor = mean(TE_Neighbor); ATE_True = mean(TE_True)
      
      return(list(ATE_Boot = ATE_Boot, ATE_Neighbor = ATE_Neighbor, ATE_True = ATE_True))
    }
    
    nsimu = 500
    
    cl = makeCluster(50)
    registerDoParallel(cl)
    
    RES = foreach(s = 1:nsimu, .packages = c("mvnfast", "glmnet", "MASS")) %dopar% ESTsimu(s)
    
    stopImplicitCluster()
    stopCluster(cl)
    
    save(list = c('RES'), file = paste0('Results_Correlated/Post_Estimation_R_', R, 'C_', C, 'Signal', set + 1, '.RData'))
    rm(RES); gc()
  }
  
  rm(list = ls()[-c(which(ls() == 'Rv'), which(ls() == 'Cv'), which(ls() == 'region'))])
}

RMSE_Boot12 = rep(0, 5); RMSE_Neighbor12 = rep(0, 5)
RMSE_Boot16 = rep(0, 5); RMSE_Neighbor16 = rep(0, 5)

nsimu = 500
for (set in 1:5)
{
  filename = paste0('Results_Correlated/Post_Estimation_R_12C_12Signal', set + 1, '.RData')
  load(filename)
  
  SE_Boot = rep(0, nsimu); SE_Neighbor = rep(0, nsimu)
  
  for (i in 1:nsimu)
  {
    SE_Boot[i] = (RES[[i]]$ATE_Boot - RES[[i]]$ATE_True)^2
    SE_Neighbor[i] = (RES[[i]]$ATE_Neighbor - RES[[i]]$ATE_True)^2
  }
  
  RMSE_Boot12[set] = sqrt(mean(SE_Boot))/RES[[i]]$ATE_True; RMSE_Neighbor12[set] = sqrt(mean(SE_Neighbor))/RES[[i]]$ATE_True
  #######################################################################################################
  
  filename = paste0('Results_Correlated/Post_Estimation_R_16C_16Signal', set + 1, '.RData')
  load(filename)
  
  SE_Boot = rep(0, nsimu); SE_Neighbor = rep(0, nsimu)
  
  for (i in 1:nsimu)
  {
    SE_Boot[i] = (RES[[i]]$ATE_Boot - RES[[i]]$ATE_True)^2
    SE_Neighbor[i] = (RES[[i]]$ATE_Neighbor - RES[[i]]$ATE_True)^2
  }
  
  RMSE_Boot16[set] = sqrt(mean(SE_Boot))/RES[[i]]$ATE_True; RMSE_Neighbor16[set] = sqrt(mean(SE_Neighbor))/RES[[i]]$ATE_True
}

delta = c(0.002, 0.004, 0.006, 0.008, 0.010)
jpeg(filename = "Figures/Correlated_MSE.png", width = 1200, height = 900)
plot(delta, RMSE_Neighbor16, main = 'Relative RMSE, Correlated Design', ylab = '', cex.main = 3.5, cex.lab = 3, cex.axis = 3, lty = 2, col = 'red', type = 'l', lwd = 3, ylim = c(0, 1))
lines(delta, RMSE_Neighbor12, type = 'l', lty = 1, lwd = 3, col = 'red')
lines(delta, RMSE_Boot12, type = 'l', lty = 1, lwd = 3, col = 'black')
lines(delta, RMSE_Boot16, type = 'l', lty = 2, lwd = 3, col = 'black')
legend('topleft', legend = c('ATE-Post, R=12, C=12', 'ATE-Post, R=16, C=16', 'ATE-MeanField, R=12, C=12', 'ATE-MeanField, R=16, C=16'), 
       lty = c(1, 2, 1, 2), col = c('black', 'black', 'red', 'red'), lwd = 3, cex = 3)
dev.off()
