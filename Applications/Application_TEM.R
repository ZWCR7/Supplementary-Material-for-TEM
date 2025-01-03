set.seed(37)
load('CMFD_PM10_Huabei.RData')
#State = State[, -c(2, 5), , ]
R = dim(PM10_Huabei)[2]; C = dim(PM10_Huabei)[3]

# WIND_Huabei = WIND_Huabei[, , 6:21]
# TEMP_Huabei = TEMP_Huabei[, , 6:21]
# PRES_Huabei = PRES_Huabei[, , 6:21]
# PM10_Huabei = PM10_Huabei[, , 6:21]
# PREC_Huabei = PREC_Huabei[, , 6:21]

#indx = which(apply(PM10_Huabei, 1, mean) > 75)
indx = 1:100
indx1 = 2:101 

n = length(indx); p = 3
X = array(0, dim = c(n, R, C, p))

X[ , , , 1] = 1
X[ , , , 2] = PREC_Huabei[indx1 , , ]
X[ , , , 3] = TEMP_Huabei[indx1 , , ] - 273.15

WIND_Huabei = sign(WIND_Huabei - 3.5)
 
M = WIND_Huabei[indx, , ]
A = array(0, dim = c(n, R, C, R*C - 1))

Avector = t(apply(M, 1, as.vector))
startpoint = list()
endpoint = list()

for (c in 1:C)
{
  for (r in 1:R)
  {
    pos = (c - 1)*R + r

    A[, r, c, ] = Avector[, -pos]
    startpoint = c(startpoint, list(1:(R*C - 1)))
    endpoint = c(endpoint, list(1:(R*C - 1)))
  }
}

rm(Avector); gc()

Y = PM10_Huabei[indx1, , ]

rm(list = c('PM10_Huabei', 'PRES_Huabei', 'SHUM_Huabei', 'TEMP_Huabei', 'WIND_Huabei', 'PREC_Huabei'))
gc()

source('BiRS_TEM.R')

totallen = p + R*C
tune = sqrt(log(1.8*totallen)/n)

# L1 = function(x)
# {
#   return((1/pnorm(1 - x/totallen))^4 + 2*(1/pnorm(1 - x/totallen))^2 - x)
# }
# 
# tunnpara = uniroot(L1, interval = c(1, 10))$root
# tune = sqrt(2)*n^(-1/2)/pnorm(1 - tunnpara/totallen)

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

#load('HSigma.RData')
LambdaMat = sqrt(3)*HSigma*sqrt(log(R*C - 1)/n)

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
    #LambdaMat[r, c] = bestlambda
    
    beta0[r, c, ] = EstOLSrc[1:p]
    S0[r, c, ] = EstOLSrc[(p + 1):(p + R*C - 1)]
    L0[r, c] = EstOLSrc[p + R*C]
  }
}

tol = 1e-04; iter_max = 100; deltaL = 2
lambdaLv = seq(10, 100, 10)
numcore = min(length(lambdaLv), 50)

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

sellambdaL = lambdaLv[which.min(RESD)]
#########################################################################################################################################################################################
#sellambdaL = 1.5
MB = 1000; alpha = 0.05; tol = 1e-04; iter_max = 200; lambdaL = sellambdaL; deltaL = 2

res = BiRS_IF(Y, X, M, A, startpoint, endpoint, beta0, L0, S0, deltaL, lambdaL, LambdaMat, HSigma, MB, alpha, tol, iter_max, ReMax = 10, method = 'Boot')



   