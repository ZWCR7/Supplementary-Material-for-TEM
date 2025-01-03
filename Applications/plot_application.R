library(ggplot2)
library(ggmap)

load('Application_PM10.RData')
L = R*C - 1

latisel = seq(35.05, 40.05, length = R)
longisel = seq(107.05, 117.05, length = C)

Map = NULL
for (c in 1:C)
{
  for (r in 1:R)
  {
    l = (c - 1)*R + r
    
    
    if (l == 1)
    {
      INTF_Vec = c(0, res$Signal[r, c, ])
    }
    
    if (l == R*C)
    {
      INTF_Vec = c(res$Signal[r, c, ], 0)
    }
    
    if ((l != 1) & (l != R*C))
    {
      INTF_Vec = c(res$Signal[r, c, 1:(l - 1)], 0, res$Signal[r, c, (l:L)])
    }
    
    INTF_Map = matrix(INTF_Vec, R, C)
    
    if (sum(INTF_Map) != 0)
    {
      INFT_Point = which(INTF_Map != 0, arr.ind = T)
      lattitude = c(r, INFT_Point[, 1])
      longitude = c(c, INFT_Point[, 2])
      
      labels = c('origin', rep('effect', nrow(INFT_Point)))
      
      Maprc = data.frame(lattitude, longitude, labels)
      # plot_Map = ggplot(Map, aes(x = longitude, y = lattitude, colour = labels)) +
      #   geom_point(size = 10) +
      #   scale_y_continuous(limits = c(1, R)) +
      #   scale_x_continuous(limits = c(1, C))
      # 
      # png(filename = paste0('Figure_CMFD/Interference_Map_longi = ', c, '_latti = ', r, '.png'), width = 1200, height = 600)
      # print(plot_Map)
      # dev.off()
      
      Map = rbind(Map, Maprc)
    }
  }
}

# map_Huabei = get_stamenmap(bbox = c(left = longisel[1], bottom = latisel[1], right = longisel[C], top = latisel[R]), zoom = 8)
# plot_Map = ggmap(map_Huabei) + 
#   geom_point(data = Map, aes(x = longitude, y = lattitude, colour = labels), size = 10)
# 
# png(filename = paste0('Figure_CMFD/Interference_Map.png'), width = 2400, height = 1200)
# print(plot_Map)
# dev.off()

plot_Map = ggplot() + geom_point(data = Map, aes(x = longitude, y = lattitude, colour = labels), size = 10)
png(filename = paste0('Figure_CMFD/Interference_Map.png'), width = 2400, height = 1200)
print(plot_Map)
dev.off()