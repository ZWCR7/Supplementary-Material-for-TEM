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


LFISTA = function(Y, X, M, beta0, L0, deltaL, lambdaL, tol = 1e-04, iter_max = 200)
{
  L1 = L0; beta1 = beta0; t0 = 1; t1 = t0
  
  n = dim(Y)[1]; R = dim(Y)[2]; C = dim(Y)[3]
  
  diff = 1000000
  index = 0
  
  while(diff > tol)
  {
    T0 = L1 + (t0 - 1)/t1*(L1 - L0)
    
    L0 = L1; beta0 = beta1; t0 = t1
    
    Ytilde = array(0, dim = c(n, R, C))
    
    for (i in 1:n)
    {
      Ytilde[i, , ] = Y[i, , ] - apply(X[i, , , ]*beta0, c(1, 2), sum)
    }
    
    GradT0 = Grad(Ytilde, M, T0)
    
    L1 = Soft(T0 - GradT0/deltaL, lambdaL/deltaL)
    ############################################################################
    
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        Ydotrc = Y[ , r, c] - M[, r, c]*L1[r, c]
        beta1[r, c, ] = lm(Ydotrc ~ X[, r, c, ] - 1)$coefficients
      }
    }
    
    diff = sqrt(sum((L1 - L0)^2))/sqrt(sum(L0^2)) + sqrt(sum((beta1 - beta0)^2))/sqrt(sum(beta0^2))
    
    t1 = (1 + sqrt(1 + 4*t0^2))/2
    
    index = index + 1
    
    print(paste0('Iter', index, '----difference:', diff))
    
    if (index > iter_max) break
  }
  
  print(paste0('Converges---------------------difference:', diff))
  
  return(list(betahat = beta1, Lhat = L1, diff = diff, index = index))
}

SFISTA = function(Y, X, M, A, beta0, L0, S0, LambdaMat, tol = 1e-04, iter_max = 200)
{
  n = dim(Y)[1]; R = dim(Y)[2]; C = dim(Y)[3]; p = dim(beta0)[3]
  
  betaL0 = array(0, dim = c(R, C, p + 1))
  betaL0[, , 1:p] = beta0; betaL0[, , p+1] = L0
  
  betaL1 = array(0, dim = c(R, C, p + 1))
  S1 = array(0, dim = c(R, C, R*C - 1))
  
  XM = array(0, dim = c(n, R, C, p + 1))
  XM[, , , 1:p] = X; XM[, , , p+1] = M
  
  diff = 1000000
  index = 0
  
  while(diff > tol)
  {
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        Ydotrc = Y[ , r, c] - A[, r, c, ] %*% S0[r, c, ]
        betaL1[r, c, ] = lm(Ydotrc ~ XM[, r, c, ] - 1)$coefficients
        
        Ybarrc = Y[ , r, c] - XM[ , r, c, ] %*% betaL1[r, c, ]
        
        S1[r, c, ] = coef(glmnet(A[, r, c, ], Ybarrc, alpha = 1, lambda = LambdaMat[r, c], intercept = F))[-1]
         
      }
    }
    
    diff = sqrt(sum((S1 - S0)^2))/sqrt(sum(S0^2)) + sqrt(sum((betaL1 - betaL0)^2))/sqrt(sum(betaL0^2))
    
    index = index + 1
    
    betaL0 = betaL1
    S0 = S1
    
    print(paste0('Iter', index, '----difference:', diff))
    
    if (index > iter_max) break
  }
  
  print(paste0('Converges---------------------difference:', diff))
  
  beta1 = betaL1[, , 1:p]; L1 = betaL1[, , p+1]
  
  return(list(betahat = beta1, Lhat = L1, Shat = S1, diff = diff, index = index))
  
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


BootRTest = function(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL, LambdaMat, HSigma, MB = 1000, alpha = 0.05, tol = 1e-04, iter_max = 200)
{
  n = dim(Y)[1]; R = dim(Y)[2]; C = dim(Y)[3]; d = dim(S0)[3]
  
  Est = LAFISTA(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL, LambdaMat, tol, iter_max, method = 'glmnet')
  Shat = Est$Shat; betahat = Est$betahat; Lhat = Est$Lhat
  
  # HSigma = matrix(0, R, C)
  # for (r in 1:R)
  # {
  #   for (c in 1:C)
  #   {
  #     XMArc = cbind(X[, r, c, ], M[, r, c], A[, r, c, ])
  #     
  #     Estrc = scalreg(XMArc, Y[, r, c])
  #     HSigma[r, c] = Estrc$hsigma
  #   }
  # }
  
  Un = sqrt(n)*abs(Shat); Tn = max(Un)
  
  UnbMat = array(0, dim = c(R, C, d, MB))
  Tnb = rep(0, MB)
  
  #LambdaMatb = matrix(0, R, C)
  LambdaMatb = LambdaMat
  I = diag(rep(1, n))
  for (b in 1:MB)
  {
    Yb = array(0, dim = c(n, R, C))
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        erc = rnorm(n, 0, HSigma[r, c])
        PXrc = X[, r, c, ] %*% solve(t(X[, r, c, ]) %*% X[, r, c, ]) %*% t(X[, r, c, ])
        Yb[, r, c] = (I - PXrc) %*% erc
      }
    }
    
    Shatb = array(0, dim = c(R, C, R*C - 1))
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        # if (b == 1)
        # {
        #   lambdarc = cv.glmnet(x = A[, r, c, ], y = Yb[, r, c], alpha = 1, intercept = F)
        #   LambdaMatb[r, c] = lambdarc$lambda.min
        #   Shatb[r, c, ] = coef(lambdarc, s = "lambda.min")[-1]
        # }
        # else
        # {
        #   Shatb[r, c, ] = coef(glmnet(x = A[, r, c, ], y = Yb[, r, c], alpha = 1, lambda = LambdaMatb[r, c], intercept = F))[-1]
        # }
        
        Shatb[r, c, ] = coef(glmnet(x = A[, r, c, ], y = Yb[, r, c], alpha = 1, lambda = LambdaMatb[r, c], intercept = F))[-1]
      }
    }
    
    UnbMat[ , , , b] = sqrt(n)*abs(Shatb)
    Tnb[b] = max(sqrt(n)*abs(Shatb))
    
    print(b)
  }
  
  CriV = quantile(Tnb, 1 - alpha)
  reject = (Tn > CriV)
  
  return(list(Tn = Tn, CriV = CriV, Un = Un, UnbMat = UnbMat, reject = reject, Shat = Shat, Lhat = Lhat, betahat = betahat, HSigma = HSigma))
}

LassoTest = function(Y, X, M, A, LambdaMat, HSigma, MB = 1000, alpha = 0.05)
{
  n = dim(Y)[1]; R = dim(Y)[2]; C = dim(Y)[3]; d = dim(A)[4]; p = dim(X)[4]
  
  betahat = array(0, dim = c(R, C, p))
  Lhat = matrix(0, R, C)
  Shat = array(0, dim = c(R, C, R*C - 1))
  for (r in 1:R)
  {
    for (c in 1:C)
    {
      Yrc = Y[ , r, c]
      XMArc = cbind(X[, r, c, -1], A[, r, c, ], M[, r, c])
      
      EstOLSrc = coef(glmnet(x = XMArc, y = Yrc, alpha = 1, lambda = LambdaMat[r, c]))
      
      betahat[r, c, ] = EstOLSrc[1:p]
      Shat[r, c, ] = EstOLSrc[(p + 1):(p + R*C - 1)]
      Lhat[r, c] = EstOLSrc[p + R*C]
    }
  }
  
  Un = sqrt(n)*abs(Shat); Tn = max(Un)
  
  UnbMat = array(0, dim = c(R, C, d, MB))
  Tnb = rep(0, MB)
  
  #LambdaMatb = matrix(0, R, C)
  LambdaMatb = LambdaMat
  I = diag(rep(1, n))
  for (b in 1:MB)
  {
    Yb = array(0, dim = c(n, R, C))
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        erc = rnorm(n, 0, HSigma[r, c])
        PXrc = X[, r, c, ] %*% solve(t(X[, r, c, ]) %*% X[, r, c, ]) %*% t(X[, r, c, ])
        Yb[, r, c] = (I - PXrc) %*% erc
      }
    }
    
    Shatb = array(0, dim = c(R, C, R*C - 1))
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        # if (b == 1)
        # {
        #   lambdarc = cv.glmnet(x = A[, r, c, ], y = Yb[, r, c], alpha = 1, intercept = F)
        #   LambdaMatb[r, c] = lambdarc$lambda.min
        #   Shatb[r, c, ] = coef(lambdarc, s = "lambda.min")[-1]
        # }
        # else
        # {
        #   Shatb[r, c, ] = coef(glmnet(x = A[, r, c, ], y = Yb[, r, c], alpha = 1, lambda = LambdaMatb[r, c], intercept = F))[-1]
        # }
        
        Shatb[r, c, ] = coef(glmnet(x = A[, r, c, ], y = Yb[, r, c], alpha = 1, lambda = LambdaMatb[r, c], intercept = F))[-1]
      }
    }
    
    UnbMat[ , , , b] = sqrt(n)*abs(Shatb)
    Tnb[b] = max(sqrt(n)*abs(Shatb))
    
    print(b)
  }
  
  CriV = quantile(Tnb, 1 - alpha)
  reject = (Tn > CriV)
  
  return(list(Tn = Tn, CriV = CriV, Un = Un, UnbMat = UnbMat, reject = reject, Shat = Shat, Lhat = Lhat, betahat = betahat, HSigma = HSigma))
}