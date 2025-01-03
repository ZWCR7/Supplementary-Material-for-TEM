source('BootTest.R')
source('DebiasedTest.R')

BiRS_IF = function(Y, X, M, A, startpoint, endpoint, beta0, L0, S0, deltaL, lambdaL, LambdaMat, HSigma, MB = 1000, alpha = 0.05, tol = 1e-04, iter_max = 200, ReMax = 10, method = 'Boot')
{
  n = dim(Y)[1]; R = dim(Y)[2]; C = dim(Y)[3]
  
  if (method == 'Boot')
  {
    TestGlo = BootRTest(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL, LambdaMat, HSigma, MB, alpha, tol, iter_max)
    Un = TestGlo$Un; UnbMat = TestGlo$UnbMat; reject = TestGlo$reject
  }
  
  if (method == 'debiased')
  {
    TestGlo = DeLassoTest(Y, X, M, A, beta0, L0, S0, deltaL, lambdaL, LambdaMat, HSigma, MB, alpha, tol, iter_max)
    Un = TestGlo$TnArr; UnbMat = TestGlo$BootArr; reject = TestGlo$reject
  }
  
  if (method == 'Lasso')
  {
    TestGlo = LassoTest(Y, X, M, A, LambdaMat, HSigma, MB, alpha)
    Un = TestGlo$Un; UnbMat = TestGlo$UnbMat; reject = TestGlo$reject
  }
  
  rm(TestGlo); gc()
  
  Signal = array(0, dim = c(R, C, R*C - 1))
  if (reject == F) return('there is no signal regions')
  
  Signal_step = stepdown(Un, UnbMat)
  
  Slist = list()
  
  for (r in 1:R)
  {
    for (c in 1:C)
    {
      l = (r - 1)*C + c
      
      num_regionrc = length(startpoint[[l]])
      mrc = ceiling(num_regionrc/2)
      
      Slistrc = list(1:mrc, (mrc + 1):num_regionrc)
      Slist[[l]] = Slistrc
    }
  }

  Signal = BiSearch(Slist, startpoint, endpoint, Signal, Un, UnbMat, alpha)
  
  indest = which(Signal != 0, arr.ind = T)
  for (j in 1:nrow(indest))
  {
    indj = indest[j, ]
    UnbMat[indj[1], indj[2], indj[3], ] = 0
    Un[indj[1], indj[2], indj[3]] = 0
  }
  
  thres_re = quantile(apply(UnbMat, 4, max), 1 - alpha)
  ind_re = (max(Un) > thres_re)
  
  index = 0
  
  while(ind_re && (index < ReMax))
  {
    Signal = BiSearch(Slist, startpoint, endpoint, Signal, Un, UnbMat, alpha)
    
    indest = which(Signal != 0, arr.ind = T)
    for (j in 1:nrow(indest))
    {
      indj = indest[j, ]
      UnbMat[indj[1], indj[2], indj[3], ] = 0
      Un[indj[1], indj[2], indj[3]] = 0
    }
    
    thres_re = quantile(apply(UnbMat, 4, max), 1 - alpha)
    ind_re = (max(Un) > thres_re)
    
    index = index + 1
  }
  
  return(list(Signal = Signal, Signal_step = Signal_step, index = index))
}

BiSearch = function(Slist, startpoint, endpoint, Signal, Un, UnbMat, alpha)
{
  R = dim(Signal)[1]; C = dim(Signal)[2]
  
  loop.ind = 1
  while(loop.ind == 1)
  {
    Slist1 = list()
    
    thres = quantile(apply(UnbMat, 4, max), 1 - alpha)
    
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        l = (r - 1)*C + c
        regIndrc = Slist[[l]]
        
        if (length(regIndrc) > 0)
        {
          Slistrc = list()
          for (i in 1:length(regIndrc))
          {
            regIndrci = regIndrc[[i]]
            pi = length(regIndrci); mi = ceiling(pi/2)
            regionrci = startpoint[[l]][regIndrci[1]]:endpoint[[l]][regIndrci[pi]]
            
            Tnrci = max(Un[r, c, regionrci]); initrci = (Tnrci > thres)
            
            if (initrci == 0)
            {
              UnbMat[r, c, regionrci, ] = 0
            }
            
            if (initrci == 1 && pi > 1) 
            {
              Slistrc = c(Slistrc, list(regIndrci[1:mi]))
              Slistrc = c(Slistrc, list(regIndrci[(mi + 1):pi]))
            }
            
            if (initrci == 1 && pi == 1)
            {
              Signal[r, c, regionrci] = 1
              UnbMat[r, c, regionrci, ] = 0
            }
          }
          
          Slist1[[l]] = Slistrc
        }
        else
        {
          Slist1[[l]] = list()
        }
      }
    }
    
    Slist = Slist1
    if (sum(lengths(Slist)) == 0) loop.ind = 0
  }
  
  return(Signal)
}

stepdown = function(Un, UnbMat)
{
  R = dim(Un)[1]; C = dim(Un)[2]; d = dim(Un)[3]
  
  Signal = array(0, dim = c(R, C, d))
  
  thres = quantile(apply(UnbMat, 4, max), 0.95)
  rej = (max(Un) > thres)
  
  while(rej)
  {
    for (r in 1:R)
    {
      for (c in 1:C)
      {
        Unrc = Un[r, c, ]
        indrc = which(Unrc > thres)
        
        if (length(indrc) > 0)
        {
          Signal[r, c, indrc] = 1
          
          Un[r, c, indrc] = 0
          UnbMat[r, c, indrc, ] = 0
        }
      }
    }
    
    thres = quantile(apply(UnbMat, 4, max), 0.95)
    rej = (max(Un) > thres)
  }
  
  return(Signal)
}



