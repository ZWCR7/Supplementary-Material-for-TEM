
TDR = rep(0, 6); FDR = rep(0, 6)
TDR_step = rep(0, 6); FDR_step = rep(0, 6)
TDR_lasso = rep(0, 6); FDR_lasso = rep(0, 6)
TDR_debias = rep(0, 6); FDR_debias = rep(0, 6)
TDR_db_step = rep(0, 6); FDR_db_step = rep(0, 6)
TDR_kf = rep(0, 6); FDR_kf = rep(0, 6)

Power_LS = rep(0, 6); Power_DS = rep(0, 6)

for (set in 1:6)
{
  load(paste0('Results_Independent/R_12C_12Signal', set, '.RData'))
  
  nsimu = 500
  TP = rep(0, nsimu); FP = rep(0, nsimu)
  TP_step = rep(0, nsimu); FP_step = rep(0, nsimu)
  TP_lasso = rep(0, nsimu); FP_lasso = rep(0, nsimu)
  TP_debias = rep(0, nsimu); FP_debias = rep(0, nsimu)
  TP_db_step = rep(0, nsimu); FP_db_step = rep(0, nsimu)
  TP_kf = rep(0, nsimu); FP_kf = rep(0, nsimu)
  
  indtrue = which(Strue != 0)
  for (i in 1:nsimu)
  {
    if (!is.character(RES[[i]]$resBoot))
    {
      Power_LS[set] = Power_LS[set] + 1
      
      signal = RES[[i]]$resBoot$Signal
      signal_step = RES[[i]]$resBoot$Signal_step
      
      indest = which(signal != 0)
      indest_step = which(signal_step != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP[i] = TP[i] + 1
          }
          else
          {
            FP[i] = FP[i] + 1
          }
        }
      }
      
      if (length(indest_step) > 0)
      {
        for (j in 1:length(indest_step))
        {
          if (sum(indtrue == indest_step[j]) > 0)
          {
            TP_step[i] = TP_step[i] + 1
          }
          else
          {
            FP_step[i] = FP_step[i] + 1
          }
        }
      }
    }
    
    # if (!is.character(RES[[i]]$resLasso))
    # {
    #   signal = RES[[i]]$resLasso$Signal
    #   indest = which(signal != 0)
    # 
    #   if (length(indest) > 0)
    #   {
    #     for (j in 1:length(indest))
    #     {
    #       if (sum(indtrue == indest[j]) > 0)
    #       {
    #         TP_lasso[i] = TP_lasso[i] + 1
    #       }
    #       else
    #       {
    #         FP_lasso[i] = FP_lasso[i] + 1
    #       }
    #     }
    #   }
    # }

    if (!is.character(RES[[i]]$resDebias))
    {
      Power_DS[set] = Power_DS[set] + 1

      signal = RES[[i]]$resDebias$Signal
      signal_step = RES[[i]]$resDebias$Signal_step

      indest = which(signal != 0)
      indest_step = which(signal_step != 0)

      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP_debias[i] = TP_debias[i] + 1
          }
          else
          {
            FP_debias[i] = FP_debias[i] + 1
          }
        }
      }

      if (length(indest_step) > 0)
      {
        for (j in 1:length(indest_step))
        {
          if (sum(indtrue == indest_step[j]) > 0)
          {
            TP_db_step[i] = TP_db_step[i] + 1
          }
          else
          {
            FP_db_step[i] = FP_db_step[i] + 1
          }
        }
      }
    }

    if (!is.character(RES[[i]]$resKnockoff))
    {
      signal = RES[[i]]$resKnockoff
      indest = which(signal != 0)

      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP_kf[i] = TP_kf[i] + 1
          }
          else
          {
            FP_kf[i] = FP_kf[i] + 1
          }
        }
      }
    }
    
    
  }
  
  TDR[set] = mean(TP)/length(indtrue)
  FDR[set] = mean(FP)/mean(FP + TP)
  
  TDR_step[set] = mean(TP_step)/length(indtrue)
  FDR_step[set] = mean(FP_step)/mean(TP_step + FP_step)
  
  # TDR_lasso[set] = mean(TP_lasso)/length(indtrue)
  # FDR_lasso[set] = mean(FP_lasso)/mean(TP_lasso + FP_lasso)
  # 
  TDR_debias[set] = mean(TP_debias)/length(indtrue)
  FDR_debias[set] = mean(FP_debias)/mean(TP_debias + FP_debias)

  TDR_db_step[set] = mean(TP_db_step)/length(indtrue)
  FDR_db_step[set] = mean(FP_db_step)/mean(TP_db_step + FP_db_step)

  TDR_kf[set] = mean(TP_kf)/length(indtrue)
  FDR_kf[set] = mean(FP_kf)/mean(TP_kf + FP_kf)
  
}

delta = c(0.002, 0.004, 0.006, 0.008, 0.010)

jpeg(filename = "Figures/Independent_R12_C12_TPR.png", width = 1200, height = 900)
#par(mai=c(0.6,0.6,0.3,0.3),omi=c(0.1,0.1,0,0))
plot(delta, TDR[2:6], type = 'l', lty = 1, col = 1, lwd = 3, main = 'TPRs, R = 12, C = 12', ylim = c(0, 1), ylab = '', cex.main = 3.5, cex.lab = 2.5, cex.axis = 2.5)
lines(delta, TDR_step[2:6], type = 'l', lty = 1, col = 2, lwd = 3)
lines(delta, TDR_debias[2:6], type = 'l', lty = 1, col = 3, lwd = 3)
lines(delta, TDR_db_step[2:6], type = 'l', lty = 1, col = 4, lwd = 3)
lines(delta, TDR_kf[2:6], type = 'l', lty = 1, col = 5, lwd = 3)
#legend('topleft', legend = c('Proposed-BiRS', 'Proposed-Step', 'Debias-BiRS', 'Debias-Step', 'Knockoff'), lwd = 2, col = c(1, 2, 3, 4, 5), cex = 0.8)
dev.off()

jpeg(filename = "Figures/Independent_R12_C12_FDR.png", width = 1200, height = 900)
plot(delta, FDR[2:6], type = 'l', lty = 1, col = 1, lwd = 3, main = 'FDRs, R = 12, C = 12', ylim = c(0, 1), ylab = '', cex.main = 3.5, cex.lab = 2.5, cex.axis = 2.5)
lines(delta, FDR_step[2:6], type = 'l', lty = 1, col = 2, lwd = 3)
lines(delta, FDR_debias[2:6], type = 'l', lty = 1, col = 3, lwd = 3)
lines(delta, FDR_db_step[2:6], type = 'l', lty = 1, col = 4, lwd = 3)
lines(delta, FDR_kf[2:6], type = 'l', lty = 1, col = 5, lwd = 3)
#legend('topleft', legend = c('Proposed-BiRS', 'Proposed-Step', 'Debias-BiRS', 'Debias-Step', 'Knockoff'), lwd = 3, col = c(1, 2, 3, 4, 5), cex = 4)
dev.off()

rm(list = ls())
#############################################################################################################################################################################

TDR = rep(0, 6); FDR = rep(0, 6)
TDR_step = rep(0, 6); FDR_step = rep(0, 6)
TDR_lasso = rep(0, 6); FDR_lasso = rep(0, 6)
TDR_debias = rep(0, 6); FDR_debias = rep(0, 6)
TDR_db_step = rep(0, 6); FDR_db_step = rep(0, 6)
TDR_kf = rep(0, 6); FDR_kf = rep(0, 6)

Power_LS = rep(0, 6); Power_DS = rep(0, 6)

for (set in 1:6)
{
  load(paste0('Results_Independent/R_16C_16Signal', set, '.RData'))
  
  nsimu = 500
  TP = rep(0, nsimu); FP = rep(0, nsimu)
  TP_step = rep(0, nsimu); FP_step = rep(0, nsimu)
  TP_lasso = rep(0, nsimu); FP_lasso = rep(0, nsimu)
  TP_debias = rep(0, nsimu); FP_debias = rep(0, nsimu)
  TP_db_step = rep(0, nsimu); FP_db_step = rep(0, nsimu)
  TP_kf = rep(0, nsimu); FP_kf = rep(0, nsimu)
  
  indtrue = which(Strue != 0)
  for (i in 1:nsimu)
  {
    if (!is.character(RES[[i]]$resBoot))
    {
      Power_LS[set] = Power_LS[set] + 1
      
      signal = RES[[i]]$resBoot$Signal
      signal_step = RES[[i]]$resBoot$Signal_step
      
      indest = which(signal != 0)
      indest_step = which(signal_step != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP[i] = TP[i] + 1
          }
          else
          {
            FP[i] = FP[i] + 1
          }
        }
      }
      
      if (length(indest_step) > 0)
      {
        for (j in 1:length(indest_step))
        {
          if (sum(indtrue == indest_step[j]) > 0)
          {
            TP_step[i] = TP_step[i] + 1
          }
          else
          {
            FP_step[i] = FP_step[i] + 1
          }
        }
      }
    }
    
    if (!is.character(RES[[i]]$resDebias))
    {
      Power_DS[set] = Power_DS[set] + 1
      
      signal = RES[[i]]$resDebias$Signal
      signal_step = RES[[i]]$resDebias$Signal_step
      
      indest = which(signal != 0)
      indest_step = which(signal_step != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP_debias[i] = TP_debias[i] + 1
          }
          else
          {
            FP_debias[i] = FP_debias[i] + 1
          }
        }
      }
      
      if (length(indest_step) > 0)
      {
        for (j in 1:length(indest_step))
        {
          if (sum(indtrue == indest_step[j]) > 0)
          {
            TP_db_step[i] = TP_db_step[i] + 1
          }
          else
          {
            FP_db_step[i] = FP_db_step[i] + 1
          }
        }
      }
    }
    
    if (!is.character(RES[[i]]$resKnockoff))
    {
      signal = RES[[i]]$resKnockoff
      indest = which(signal != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP_kf[i] = TP_kf[i] + 1
          }
          else
          {
            FP_kf[i] = FP_kf[i] + 1
          }
        }
      }
    }
    
    
  }
  
  TDR[set] = mean(TP)/length(indtrue)
  FDR[set] = mean(FP)/mean(FP + TP)
  
  TDR_step[set] = mean(TP_step)/length(indtrue)
  FDR_step[set] = mean(FP_step)/mean(TP_step + FP_step)
  
  TDR_debias[set] = mean(TP_debias)/length(indtrue)
  FDR_debias[set] = mean(FP_debias)/mean(TP_debias + FP_debias)
  
  TDR_db_step[set] = mean(TP_db_step)/length(indtrue)
  FDR_db_step[set] = mean(FP_db_step)/mean(TP_db_step + FP_db_step)
  
  TDR_kf[set] = mean(TP_kf)/length(indtrue)
  FDR_kf[set] = mean(FP_kf)/mean(TP_kf + FP_kf)
  
}

delta = c(0.002, 0.004, 0.006, 0.008, 0.010)

jpeg(filename = "Figures/Independent_R16_C16_TPR.png", width = 1200, height = 900)
#par(mai=c(0.6,0.6,0.3,0.3),omi=c(0.1,0.1,0,0))
plot(delta, TDR[2:6], type = 'l', lty = 1, col = 1, lwd = 3, main = 'TPRs, R = 16, C = 16', ylim = c(0, 1), ylab = '', cex.main = 3.5, cex.lab = 2.5, cex.axis = 2.5)
lines(delta, TDR_step[2:6], type = 'l', lty = 1, col = 2, lwd = 3)
lines(delta, TDR_debias[2:6], type = 'l', lty = 1, col = 3, lwd = 3)
lines(delta, TDR_db_step[2:6], type = 'l', lty = 1, col = 4, lwd = 3)
lines(delta, TDR_kf[2:6], type = 'l', lty = 1, col = 5, lwd = 3)
#legend('topleft', legend = c('Proposed-BiRS', 'Proposed-Step', 'Debias-BiRS', 'Debias-Step', 'Knockoff'), lwd = 2, col = c(1, 2, 3, 4, 5), cex = 0.8)
dev.off()

jpeg(filename = "Figures/Independent_R16_C16_FDR.png", width = 1200, height = 900)
plot(delta, FDR[2:6], type = 'l', lty = 1, col = 1, lwd = 3, main = 'FDRs, R = 16, C = 16', ylim = c(0, 1), ylab = '', cex.main = 3.5, cex.lab = 2.5, cex.axis = 2.5)
lines(delta, FDR_step[2:6], type = 'l', lty = 1, col = 2, lwd = 3)
lines(delta, FDR_debias[2:6], type = 'l', lty = 1, col = 3, lwd = 3)
lines(delta, FDR_db_step[2:6], type = 'l', lty = 1, col = 4, lwd = 3)
lines(delta, FDR_kf[2:6], type = 'l', lty = 1, col = 5, lwd = 3)
legend('topleft', legend = c('Proposed-BiRS', 'Proposed-Step', 'Debias-BiRS', 'Debias-Step', 'Knockoff'), lwd = 3, col = c(1, 2, 3, 4, 5), cex = 4)
dev.off()

rm(list = ls())
################################################################################################################################################################################
################################################################################################################################################################################


TDR = rep(0, 6); FDR = rep(0, 6)
TDR_step = rep(0, 6); FDR_step = rep(0, 6)
TDR_lasso = rep(0, 6); FDR_lasso = rep(0, 6)
TDR_debias = rep(0, 6); FDR_debias = rep(0, 6)
TDR_db_step = rep(0, 6); FDR_db_step = rep(0, 6)
TDR_kf = rep(0, 6); FDR_kf = rep(0, 6)

Power_LS = rep(0, 6); Power_DS = rep(0, 6)

for (set in 1:6)
{
  load(paste0('Results_Correlated/R_12C_12Signal', set, '.RData'))
  
  nsimu = 500
  TP = rep(0, nsimu); FP = rep(0, nsimu)
  TP_step = rep(0, nsimu); FP_step = rep(0, nsimu)
  TP_lasso = rep(0, nsimu); FP_lasso = rep(0, nsimu)
  TP_debias = rep(0, nsimu); FP_debias = rep(0, nsimu)
  TP_db_step = rep(0, nsimu); FP_db_step = rep(0, nsimu)
  TP_kf = rep(0, nsimu); FP_kf = rep(0, nsimu)
  
  indtrue = which(Strue != 0)
  for (i in 1:nsimu)
  {
    if (!is.character(RES[[i]]$resBoot))
    {
      Power_LS[set] = Power_LS[set] + 1
      
      signal = RES[[i]]$resBoot$Signal
      signal_step = RES[[i]]$resBoot$Signal_step
      
      indest = which(signal != 0)
      indest_step = which(signal_step != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP[i] = TP[i] + 1
          }
          else
          {
            FP[i] = FP[i] + 1
          }
        }
      }
      
      if (length(indest_step) > 0)
      {
        for (j in 1:length(indest_step))
        {
          if (sum(indtrue == indest_step[j]) > 0)
          {
            TP_step[i] = TP_step[i] + 1
          }
          else
          {
            FP_step[i] = FP_step[i] + 1
          }
        }
      }
    }
    
    # if (!is.character(RES[[i]]$resLasso))
    # {
    #   signal = RES[[i]]$resLasso$Signal
    #   indest = which(signal != 0)
    # 
    #   if (length(indest) > 0)
    #   {
    #     for (j in 1:length(indest))
    #     {
    #       if (sum(indtrue == indest[j]) > 0)
    #       {
    #         TP_lasso[i] = TP_lasso[i] + 1
    #       }
    #       else
    #       {
    #         FP_lasso[i] = FP_lasso[i] + 1
    #       }
    #     }
    #   }
    # }
    
    if (!is.character(RES[[i]]$resDebias))
    {
      Power_DS[set] = Power_DS[set] + 1
      
      signal = RES[[i]]$resDebias$Signal
      signal_step = RES[[i]]$resDebias$Signal_step
      
      indest = which(signal != 0)
      indest_step = which(signal_step != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP_debias[i] = TP_debias[i] + 1
          }
          else
          {
            FP_debias[i] = FP_debias[i] + 1
          }
        }
      }
      
      if (length(indest_step) > 0)
      {
        for (j in 1:length(indest_step))
        {
          if (sum(indtrue == indest_step[j]) > 0)
          {
            TP_db_step[i] = TP_db_step[i] + 1
          }
          else
          {
            FP_db_step[i] = FP_db_step[i] + 1
          }
        }
      }
    }
    
    if (!is.character(RES[[i]]$resKnockoff))
    {
      signal = RES[[i]]$resKnockoff
      indest = which(signal != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP_kf[i] = TP_kf[i] + 1
          }
          else
          {
            FP_kf[i] = FP_kf[i] + 1
          }
        }
      }
    }
    
    
  }
  
  TDR[set] = mean(TP)/length(indtrue)
  FDR[set] = mean(FP)/mean(FP + TP)
  
  TDR_step[set] = mean(TP_step)/length(indtrue)
  FDR_step[set] = mean(FP_step)/mean(TP_step + FP_step)
  
  # TDR_lasso[set] = mean(TP_lasso)/length(indtrue)
  # FDR_lasso[set] = mean(FP_lasso)/mean(TP_lasso + FP_lasso)
  # 
  TDR_debias[set] = mean(TP_debias)/length(indtrue)
  FDR_debias[set] = mean(FP_debias)/mean(TP_debias + FP_debias)
  
  TDR_db_step[set] = mean(TP_db_step)/length(indtrue)
  FDR_db_step[set] = mean(FP_db_step)/mean(TP_db_step + FP_db_step)
  
  TDR_kf[set] = mean(TP_kf)/length(indtrue)
  FDR_kf[set] = mean(FP_kf)/mean(TP_kf + FP_kf)
  
}

delta = c(0.002, 0.004, 0.006, 0.008, 0.010)

jpeg(filename = "Figures/Correlated_R12_C12_TPR.png", width = 1200, height = 900)
#par(mai=c(0.6,0.6,0.3,0.3),omi=c(0.1,0.1,0,0))
plot(delta, TDR[2:6], type = 'l', lty = 1, col = 1, lwd = 3, main = 'TPRs, R = 12, C = 12', ylim = c(0, 1), ylab = '', cex.main = 3.5, cex.lab = 2.5, cex.axis = 2.5)
lines(delta, TDR_step[2:6], type = 'l', lty = 1, col = 2, lwd = 3)
lines(delta, TDR_debias[2:6], type = 'l', lty = 1, col = 3, lwd = 3)
lines(delta, TDR_db_step[2:6], type = 'l', lty = 1, col = 4, lwd = 3)
lines(delta, TDR_kf[2:6], type = 'l', lty = 1, col = 5, lwd = 3)
#legend('topleft', legend = c('Proposed-BiRS', 'Proposed-Step', 'Debias-BiRS', 'Debias-Step', 'Knockoff'), lwd = 2, col = c(1, 2, 3, 4, 5), cex = 0.8)
dev.off()

jpeg(filename = "Figures/Correlated_R12_C12_FDR.png", width = 1200, height = 900)
plot(delta, FDR[2:6], type = 'l', lty = 1, col = 1, lwd = 3, main = 'FDRs, R = 12, C = 12', ylim = c(0, 1), ylab = '', cex.main = 3.5, cex.lab = 2.5, cex.axis = 2.5)
lines(delta, FDR_step[2:6], type = 'l', lty = 1, col = 2, lwd = 3)
lines(delta, FDR_debias[2:6], type = 'l', lty = 1, col = 3, lwd = 3)
lines(delta, FDR_db_step[2:6], type = 'l', lty = 1, col = 4, lwd = 3)
lines(delta, FDR_kf[2:6], type = 'l', lty = 1, col = 5, lwd = 3)
#legend('topleft', legend = c('Proposed-BiRS', 'Proposed-Step', 'Debias-BiRS', 'Debias-Step', 'Knockoff'), lwd = 3, col = c(1, 2, 3, 4, 5), cex = 4)
dev.off()

rm(list = ls())
#############################################################################################################################################################################

TDR = rep(0, 6); FDR = rep(0, 6)
TDR_step = rep(0, 6); FDR_step = rep(0, 6)
TDR_lasso = rep(0, 6); FDR_lasso = rep(0, 6)
TDR_debias = rep(0, 6); FDR_debias = rep(0, 6)
TDR_db_step = rep(0, 6); FDR_db_step = rep(0, 6)
TDR_kf = rep(0, 6); FDR_kf = rep(0, 6)

Power_LS = rep(0, 6); Power_DS = rep(0, 6)

for (set in 1:6)
{
  load(paste0('Results_Correlated/R_16C_16Signal', set, '.RData'))
  
  nsimu = 500
  TP = rep(0, nsimu); FP = rep(0, nsimu)
  TP_step = rep(0, nsimu); FP_step = rep(0, nsimu)
  TP_lasso = rep(0, nsimu); FP_lasso = rep(0, nsimu)
  TP_debias = rep(0, nsimu); FP_debias = rep(0, nsimu)
  TP_db_step = rep(0, nsimu); FP_db_step = rep(0, nsimu)
  TP_kf = rep(0, nsimu); FP_kf = rep(0, nsimu)
  
  indtrue = which(Strue != 0)
  for (i in 1:nsimu)
  {
    if (!is.character(RES[[i]]$resBoot))
    {
      Power_LS[set] = Power_LS[set] + 1
      
      signal = RES[[i]]$resBoot$Signal
      signal_step = RES[[i]]$resBoot$Signal_step
      
      indest = which(signal != 0)
      indest_step = which(signal_step != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP[i] = TP[i] + 1
          }
          else
          {
            FP[i] = FP[i] + 1
          }
        }
      }
      
      if (length(indest_step) > 0)
      {
        for (j in 1:length(indest_step))
        {
          if (sum(indtrue == indest_step[j]) > 0)
          {
            TP_step[i] = TP_step[i] + 1
          }
          else
          {
            FP_step[i] = FP_step[i] + 1
          }
        }
      }
    }
    
    if (!is.character(RES[[i]]$resDebias))
    {
      Power_DS[set] = Power_DS[set] + 1
      
      signal = RES[[i]]$resDebias$Signal
      signal_step = RES[[i]]$resDebias$Signal_step
      
      indest = which(signal != 0)
      indest_step = which(signal_step != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP_debias[i] = TP_debias[i] + 1
          }
          else
          {
            FP_debias[i] = FP_debias[i] + 1
          }
        }
      }
      
      if (length(indest_step) > 0)
      {
        for (j in 1:length(indest_step))
        {
          if (sum(indtrue == indest_step[j]) > 0)
          {
            TP_db_step[i] = TP_db_step[i] + 1
          }
          else
          {
            FP_db_step[i] = FP_db_step[i] + 1
          }
        }
      }
    }
    
    if (!is.character(RES[[i]]$resKnockoff))
    {
      signal = RES[[i]]$resKnockoff
      indest = which(signal != 0)
      
      if (length(indest) > 0)
      {
        for (j in 1:length(indest))
        {
          if (sum(indtrue == indest[j]) > 0)
          {
            TP_kf[i] = TP_kf[i] + 1
          }
          else
          {
            FP_kf[i] = FP_kf[i] + 1
          }
        }
      }
    }
    
    
  }
  
  TDR[set] = mean(TP)/length(indtrue)
  FDR[set] = mean(FP)/mean(FP + TP)
  
  TDR_step[set] = mean(TP_step)/length(indtrue)
  FDR_step[set] = mean(FP_step)/mean(TP_step + FP_step)
  
  TDR_debias[set] = mean(TP_debias)/length(indtrue)
  FDR_debias[set] = mean(FP_debias)/mean(TP_debias + FP_debias)
  
  TDR_db_step[set] = mean(TP_db_step)/length(indtrue)
  FDR_db_step[set] = mean(FP_db_step)/mean(TP_db_step + FP_db_step)
  
  TDR_kf[set] = mean(TP_kf)/length(indtrue)
  FDR_kf[set] = mean(FP_kf)/mean(TP_kf + FP_kf)
  
}

delta = c(0.002, 0.004, 0.006, 0.008, 0.010)

jpeg(filename = "Figures/Correlated_R16_C16_TPR.png", width = 1200, height = 900)
#par(mai=c(0.6,0.6,0.3,0.3),omi=c(0.1,0.1,0,0))
plot(delta, TDR[2:6], type = 'l', lty = 1, col = 1, lwd = 3, main = 'TPRs, R = 16, C = 16', ylim = c(0, 1), ylab = '', cex.main = 3.5, cex.lab = 2.5, cex.axis = 2.5)
lines(delta, TDR_step[2:6], type = 'l', lty = 1, col = 2, lwd = 3)
lines(delta, TDR_debias[2:6], type = 'l', lty = 1, col = 3, lwd = 3)
lines(delta, TDR_db_step[2:6], type = 'l', lty = 1, col = 4, lwd = 3)
lines(delta, TDR_kf[2:6], type = 'l', lty = 1, col = 5, lwd = 3)
#legend('topleft', legend = c('Proposed-BiRS', 'Proposed-Step', 'Debias-BiRS', 'Debias-Step', 'Knockoff'), lwd = 2, col = c(1, 2, 3, 4, 5), cex = 0.8)
dev.off()

jpeg(filename = "Figures/Correlated_R16_C16_FDR.png", width = 1200, height = 900)
plot(delta, FDR[2:6], type = 'l', lty = 1, col = 1, lwd = 3, main = 'FDRs, R = 16, C = 16', ylim = c(0, 1), ylab = '', cex.main = 3.5, cex.lab = 2.5, cex.axis = 2.5)
lines(delta, FDR_step[2:6], type = 'l', lty = 1, col = 2, lwd = 3)
lines(delta, FDR_debias[2:6], type = 'l', lty = 1, col = 3, lwd = 3)
lines(delta, FDR_db_step[2:6], type = 'l', lty = 1, col = 4, lwd = 3)
lines(delta, FDR_kf[2:6], type = 'l', lty = 1, col = 5, lwd = 3)
legend('topleft', legend = c('Proposed-BiRS', 'Proposed-Step', 'Debias-BiRS', 'Debias-Step', 'Knockoff'), lwd = 3, col = c(1, 2, 3, 4, 5), cex = 4)
dev.off()

rm(list = ls())
################################################################################################################################################################################