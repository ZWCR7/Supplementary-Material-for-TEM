.rs.restartR()

library(raster)
library(ncdf4)
library(rasterVis)
library(lattice)

input_wind = "Data_forcing_01dy_010deg_2018/wind_ITPCAS-CMFD_V0106_B-01_01dy_010deg_201801-201812.nc"
windfile = nc_open(input_wind)
wind_name = names(windfile$var)
ras_wind = stack(input_wind, varname = wind_name)
longitude = ncvar_get(windfile,'lon')
latitude = ncvar_get(windfile,'lat')
nc_close(windfile)

n = dim(ras_wind)[3]; p1 = dim(ras_wind)[1]; p2 = dim(ras_wind)[2]
WIND = array(0, dim = c(n, p1, p2))
for (i in 1:n)
{
  WIND[i, , ] = as.matrix(ras_wind[[i]])
}

#WIND = as.data.frame(WIND)

rm(list = ls()[-c(which(ls() == 'WIND'), which(ls() == 'longitude'), which(ls() == 'latitude'))])
########################################################################################

input_prec = "Data_forcing_01dy_010deg_2018/prec_ITPCAS-CMFD_V0106_B-01_01dy_010deg_201801-201812.nc"
precfile = nc_open(input_prec)
prec_name = names(precfile$var)
ras_prec = stack(input_prec, varname = prec_name)
nc_close(precfile)

n = dim(ras_prec)[3]; p1 = dim(ras_prec)[1]; p2 = dim(ras_prec)[2]
PREC = array(0, dim = c(n, p1, p2))
for (i in 1:n)
{
  PREC[i, , ] = as.matrix(ras_prec[[i]])
}

rm(list = ls()[-c(which(ls() == 'WIND'), which(ls() == 'PREC'), which(ls() == 'longitude'), which(ls() == 'latitude'))])
#####################################################################################

input_temp = "Data_forcing_01dy_010deg_2018/temp_ITPCAS-CMFD_V0106_B-01_01dy_010deg_201801-201812.nc"
tempfile = nc_open(input_temp)
temp_name = names(tempfile$var)
ras_temp = stack(input_temp, varname = temp_name)
nc_close(tempfile)

n = dim(ras_temp)[3]; p1 = dim(ras_temp)[1]; p2 = dim(ras_temp)[2]
TEMP = array(0, dim = c(n, p1, p2))
for (i in 1:n)
{
  TEMP[i, , ] = as.matrix(ras_temp[[i]])
}

rm(list = ls()[-c(which(ls() == 'WIND'), which(ls() == 'PREC'), which(ls() == 'TEMP'), which(ls() == 'longitude'), which(ls() == 'latitude'))])
##########################################################################################

input_shum = "Data_forcing_01dy_010deg_2018/shum_ITPCAS-CMFD_V0106_B-01_01dy_010deg_201801-201812.nc"
shumfile = nc_open(input_shum)
shum_name = names(shumfile$var)
ras_shum = stack(input_shum, varname = shum_name)
nc_close(shumfile)

n = dim(ras_shum)[3]; p1 = dim(ras_shum)[1]; p2 = dim(ras_shum)[2]
SHUM = array(0, dim = c(n, p1, p2))
for (i in 1:n)
{
  SHUM[i, , ] = as.matrix(ras_shum[[i]])
}

rm(list = ls()[-c(which(ls() == 'WIND'), which(ls() == 'PREC'), which(ls() == 'TEMP'), 
                  which(ls() == 'SHUM'), which(ls() == 'longitude'), which(ls() == 'latitude'))])
######################################################################################

input_pres = "Data_forcing_01dy_010deg_2018/pres_ITPCAS-CMFD_V0106_B-01_01dy_010deg_201801-201812.nc"
presfile = nc_open(input_pres)
pres_name = names(presfile$var)
ras_pres = stack(input_pres, varname = pres_name)
nc_close(presfile)

n = dim(ras_pres)[3]; p1 = dim(ras_pres)[1]; p2 = dim(ras_pres)[2]
PRES = array(0, dim = c(n, p1, p2))
for (i in 1:n)
{
  PRES[i, , ] = as.matrix(ras_pres[[i]])
}

rm(list = ls()[-c(which(ls() == 'WIND'), which(ls() == 'PREC'), which(ls() == 'TEMP'), 
                  which(ls() == 'SHUM'), which(ls() == 'PRES'), which(ls() == 'longitude'), which(ls() == 'latitude'))])

# n = dim(WIND)[1]; p1 = dim(WIND)[2]; p2 = dim(WIND)[3]
# State = array(0, dim = c(n, 5, p1, p2))
# State[ , 1, , ] = WIND
# State[ , 2, , ] = PRES
# State[ , 3, , ] = SHUM
# State[ , 4, , ] = TEMP
# State[ , 5, , ] = PREC

#rm("PREC", "WIND", "PRES", "SHUM", "TEMP", "n", "p1", "p2")

save.image(file = "CMFD_day.RData")

#################################################################################################
load('CMFD_day.RData')

latisel = round(seq(35.05, 40.05, 0.5), 2)
longisel = round(seq(107.05, 117.05, 0.5), 2)

latitude = round(latitude, 2)
longitude = round(longitude, 2)

indlati = rep(0, length(latisel))
indlongi = rep(0, length(longisel))

for (j in 1:length(latisel))
{
  indlati[j] = which(latitude == latisel[j])
}

for (k in 1:length(longisel))
{
  indlongi[k] = which(longitude == longisel[k])
}

WIND_Huabei = WIND[, indlati, indlongi]
PREC_Huabei = PREC[, indlati, indlongi]
TEMP_Huabei = TEMP[, indlati, indlongi]
SHUM_Huabei = SHUM[, indlati, indlongi]
PRES_Huabei = PRES[, indlati, indlongi]

dayname = NULL

M1 = c(1, 3, 5, 7, 8, 10, 12); M2 = c(4, 6, 9, 11)

for (m in 1:9)
{
  if (sum(M1 == m) == 1)
  {
    for (d in 1:9)
    {
      dayname = c(dayname, paste0('20180', m, '0', d))
    }
    
    for (d in 10:31)
    {
      dayname = c(dayname, paste0('20180', m, d))
    }
  }
  
  if (m == 2)
  {
    for (d in 1:9)
    {
      dayname = c(dayname, paste0('20180', m, '0', d))
    }
    
    for (d in 10:28)
    {
      dayname = c(dayname, paste0('20180', m, d))
    }
  }
  
  if (sum(M2 == m) == 1)
  {
    for (d in 1:9)
    {
      dayname = c(dayname, paste0('20180', m, '0', d))
    }
    
    for (d in 10:30)
    {
      dayname = c(dayname, paste0('20180', m, d))
    }
  } 
   
}

for (m in 10:12)
{
  if (sum(M1 == m) == 1)
  {
    for (d in 1:9)
    {
      dayname = c(dayname, paste0('2018', m, '0', d))
    }
    
    for (d in 10:31)
    {
      dayname = c(dayname, paste0('2018', m, d))
    }
  }
  
  if (sum(M2 == m) == 1)
  {
    for (d in 1:9)
    {
      dayname = c(dayname, paste0('2018', m, '0', d))
    }
    
    for (d in 10:30)
    {
      dayname = c(dayname, paste0('2018', m, d))
    }
  } 
  
}

PM10_Huabei_upp = array(0, dim = c(365, length(latisel), length(longisel)))
PM10_Huabei_low = array(0, dim = c(365, length(latisel), length(longisel)))
for (i in 1:365)
{
  filename = paste0("PM10DATA2018/CHAP_PM10_D1K_", dayname[i], '_V4.nc')
  ncdata = nc_open(filename)
  
  # extract date from filename
  # date = substr(filename,regexpr('\\d{8}',filename),regexpr('\\d{8}',filename)+7)
  
  # read longitude and latitude from NetCDF file
  lon = ncvar_get(ncdata, 'lon')
  lat = ncvar_get(ncdata, 'lat')
  
  # read PM10 data
  PM10 = t(ncvar_get(ncdata, 'PM10'))
  latlow = rep(0, length(latisel))
  latupp = rep(0, length(latisel))
  for (j in 1:length(latisel))
  {
    latlow[j] = min(which(lat <= latisel[j]))
    latupp[j] = max(which(lat >= latisel[j]))
  }
  
  lonlow = rep(0, length(longisel))
  lonupp = rep(0, length(longisel))
  for (k in 1:length(longisel))
  {
    lonlow[k] = max(which(lon <= longisel[k]))
    lonupp[k] = min(which(lon >= longisel[k]))
  }
  
  PM10_Huabei_low[i, , ] = PM10[latlow, lonlow]
  PM10_Huabei_upp[i, , ] = PM10[latupp, lonupp]
}

PM10_Huabei = (PM10_Huabei_low + PM10_Huabei_upp)/2
save(list = c('WIND_Huabei', 'PREC_Huabei', 'TEMP_Huabei', 'SHUM_Huabei', 'PRES_Huabei', 'PM10_Huabei', 'latisel', 'longisel'), file = 'CMFD_PM10_Huabei_Test.RData')





