library(knockoff)
library(glmnet)
library(scalreg)

Grad = function(Y, M, L)
{
  n = dim(Y)[1]; R = dim(Y)[2]; C = dim(Y)[3]
  
  Err = array(0, dim = c(n, R, C))
  
  for (i in 1:n)
  {
    Err[i, , ] = Y[i, , ] - M[i, , ]*L 
  }
  
  Grad = -2*apply(M*Err, c(2, 3), mean)
  
  return(Grad)
}

Soft = function(L, tau)
{
  SVDMat = svd(L)
  Um = SVDMat$u; Vm = SVDMat$v; Dm = SVDMat$d
  
  correct = rep(0, length(Dm))
  
  indp = which((Dm - tau) > 0)
  correct[indp] = Dm[indp] - tau
  
  SoftMat = Um %*% diag(correct) %*% t(Vm)
  
  return(SoftMat)
}

LAFISTA = function(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL, LambdaMat, tol = 1e-04, iter_max = 200, method = "glmnet")
{
  L1 = L0; S1 = S0; beta1 = beta0; t0 = 1; t1 = t0
  
  n = dim(Y)[1]; R = dim(Y)[2]; C = dim(Y)[3]
  
  diff = 1000000
  index = 0
  
  while(diff > tol)
  {
    T0 = L1 + (t0 - 1)/t1*(L1 - L0)
    
    L0 = L1; S0 = S1; beta0 = beta1; t0 = t1
    
    Ytilde = array(0, dim = c(n, R, C))
    
    for (i in 1:n)
    {
      Ytilde[i, , ] = Y[i, , ] - apply(X[i, , , ]*beta0, c(1, 2), sum) - apply(A[i, , , ]*S0, c(1, 2), sum)
    }
    
    GradT0 = Grad(Ytilde, M, T0)
    
    L1 = Soft(T0 - GradT0/deltaL, lambdaL/deltaL)
    ############################################################################
    
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        Ydotrc = Y[ , r, c] - A[, r, c, ] %*% S0[r, c, ] - M[, r, c]*L1[r, c]
        beta1[r, c, ] = lm(Ydotrc ~ X[, r, c, ] - 1)$coefficients
        
        Ybarrc = Y[ , r, c] - X[ , r, c, ] %*% beta1[r, c, ] - M[, r, c]*L1[r, c]
        
        if (method == "glmnet")
        {
          S1[r, c, ] = coef(glmnet(A[, r, c, ], Ybarrc, alpha = 1, lambda = LambdaMat[r, c], intercept = F))[-1]
        }
        
        if (method == "scale")
        {
          S1[r, c, ] = scalreg(A[, r, c, ], Ybarrc, lam0 = LambdaMat[r, c])$coefficients
        }
        
      }
    }
    
    diff = sqrt(sum((L1 - L0)^2))/sqrt(sum(L0^2)) 
    + sqrt(sum((S1 - S0)^2))/sqrt(sum(S0^2)) + sqrt(sum((beta1 - beta0)^2))/sqrt(sum(beta0^2))
    
    t1 = (1 + sqrt(1 + 4*t0^2))/2
    
    index = index + 1
    
    print(paste0('Iter', index, '----difference:', diff))
    
    if (index > iter_max) break
  }
  
  print(paste0('Converges---------------------difference:', diff))
  
  return(list(betahat = beta1, Lhat = L1, Shat = S1, diff = diff, index = index))
}

KnockoffSelect = function(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL, LambdaMat, tol = 1e-04, iter_max = 200)
{
  n = dim(Y)[1]; R = dim(Y)[2]; C = dim(Y)[3]; d = dim(S0)[3]
  
  Est = LAFISTA(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL, LambdaMat, tol, iter_max, method = 'glmnet')
  Shat = Est$Shat; betahat = Est$betahat; Lhat = Est$Lhat
  
  Signal = array(0, dim = c(R, C, R*C - 1))
  index = 0
  for (r in 1:R)
  {
    for (c in 1:C)
    {
      index = index + 1
      
      Ybarrc = Y[ , r, c] - X[ , r, c, ] %*% betahat[r, c, ] - M[, r, c]*Lhat[r, c]
      
      detectrc = knockoff.filter(X = A[, r, c, ], y = Ybarrc, fdr = 0.05)$selected
      Signal[r, c, detectrc] = 1
      
      print(index)
    }
  }
  
  return(Signal)
}