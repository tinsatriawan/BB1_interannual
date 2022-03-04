# T. Satriawan 2022-01-31
# This script is for data wrangling, including: 
# (1) Filtering data according to the dates needed for the analysis, 
# (2) calculate daily fluxes (gapfilled and non-gapfilled) to get gC/m2/day
# (3) Calculate daily met variables
# to run the codes in another script, use "source('data_wrangling.R')"



library(plyr)
library(ggplot2)
library(lubridate)
library(naniar)
library(tidyr)
library(plotly)
library(knitr)
library(data.table)
library(kableExtra)
library(ggsci)
library(effects)
library(tidyverse)
library(zoo)
library(caret)
library(ggpubr)
library(gridExtra)
library(lognorm)
library(reshape2)
library(matrixStats)
library(gtools)
library(dplyr)



#####------------- FUNCTIONS -------------####

# function for setting year period
set_period <- function(data) {
  data$period <- NA 
  
  y_2015_2016 <- data$date >= "2015-10-01" & data$date < "2016-10-01"
  y_2016_2017 <- data$date >= "2016-10-01" & data$date < "2017-10-01"
  y_2017_2018 <- data$date >= "2017-10-01" & data$date < "2018-10-01"
  y_2018_2019 <- data$date >= "2018-10-01" & data$date < "2019-10-01"
  y_2019_2020 <- data$date >= "2019-10-01" & data$date < "2020-10-01"
  y_2020_2021 <- data$date >= "2020-10-01" & data$date < "2021-10-01"
  
  data$period[y_2015_2016] <- "2015-2016"
  data$period[y_2016_2017] <- "2016-2017"
  data$period[y_2017_2018] <- "2017-2018"
  data$period[y_2018_2019] <- "2018-2019"
  data$period[y_2019_2020] <- "2019-2020"
  data$period[y_2020_2021] <- "2020-2021"
  
  data
}

# set season 
set_season <- function(input, date_vector){
  input$month_local <- month(date_vector)
  GS <- input$month_local %in% c(4, 5, 6, 7, 8, 9)
  NGS <- input$month_local %in% c(10, 11, 12, 1, 2, 3)
  input$season <- NA
  input$season[GS] <- "growing season"
  input$season[NGS] <- "non-growing season"
  input$season <- factor(input$season, levels = c("non-growing season", "growing season"))
  input
}



#####------------------- READ DATA ---------------------#####

# define directory

dir_rf <- "G:\\.shortcut-targets-by-id\\1txCh9lZ7VGCujXGvBaJCMVnuxT-65q4K\\Micromet Lab\\Projects\\2014-BB1 Burns Bog\\Flux-tower (1)\\flux_data\\for_uncertainty\\RF_results"
dir_unc <- "G:\\.shortcut-targets-by-id\\1txCh9lZ7VGCujXGvBaJCMVnuxT-65q4K\\Micromet Lab\\Projects\\2014-BB1 Burns Bog\\Flux-tower (1)\\flux_data\\for_uncertainty"

# dir_rf <- "G:\\.shortcut-targets-by-id\\1txCh9lZ7VGCujXGvBaJCMVnuxT-65q4K\\Micromet Lab\\People\\2020-Tin Satriawan\\for_uncertainty\\RF_results"
# dir_unc <- "G:\\.shortcut-targets-by-id\\1txCh9lZ7VGCujXGvBaJCMVnuxT-65q4K\\Micromet Lab\\People\\2020-Tin Satriawan\\for_uncertainty"



## read NEE RF gapfilling with 20x iteration
result_RF_NEE <- readRDS(paste0(dir_rf,"/BB_NEE_rf_result"))
# average over 20 iteration
RF_NEE <- result_RF_NEE %>% pivot_wider(id_cols = DateTime, names_from = iteration, values_from = NEE_RF_filled)
RF_NEE <- data.frame(DateTime = RF_NEE$DateTime, 
                     NEE_RF_filled = rowMeans(RF_NEE[, 2:ncol(RF_NEE)]), 
                     NEE_RF_fsd = rowSds(as.matrix(RF_NEE[, 2:ncol(RF_NEE)])))

## read CH4 RF gapfilling with 20x iteration
result_RF_FCH4 <- read_csv(paste0(dir_rf,"/BB_FCH4_rf_result.csv"))

# read flux & met data
BB1 <- fread(paste0(dir_unc,"/BB_L3.csv"))
BB1[BB1 == -9999] <- NA
BB1$DATE <- as.POSIXct(BB1$DATE, tz = "UTC")
BB1$date <- as.Date(BB1$DATE, format="%Y-%m-%d %H:%M:%S")

BB1 <- BB1 %>%
  mutate(period = NA, 
         season = NA)


# fill long gaps after MDS gapfilling with values from random forest gapfilling
long_gap2 <- BB1$DATE >= as.POSIXct("2015-12-01 00:30:00", tz = "UTC") &
  BB1$DATE < as.POSIXct("2016-02-21 00:30:00", tz = "UTC")
long_gap <- BB1$DATE >= as.POSIXct("2016-10-22 00:30:00", tz = "UTC") &
  BB1$DATE < as.POSIXct("2017-01-21 00:30:00", tz = "UTC")


BB1$NEE_RF <- RF_NEE$NEE_RF_filled  # random forest gapfill
BB1$NEE_f_RF <- BB1$NEE_f
BB1$NEE_RF_fsd <- NA
BB1$NEE_RF_fsd <- RF_NEE$NEE_RF_filled_SD

# fill all NA values with values from random forest
BB1$NEE_f_RF[long_gap] <- BB1$NEE_RF[long_gap]
BB1$NEE_f_RF[long_gap2] <- BB1$NEE_RF[long_gap2]


# check if there are still gaps
which(is.na(BB1$NEE_f_RF))


BB1$GPP_f_RF <- BB1$Reco - BB1$NEE_f_RF


# filter date
BB1 <- BB1 %>% 
  filter(DATE >= as.POSIXct("2015-10-01 00:30:00", tz = "UTC") & DATE < as.POSIXct("2021-10-01 00:00:00", tz = "UTC")) 
result_RF_FCH4 <- result_RF_FCH4 %>% 
  filter(DateTime >= as.POSIXct("2015-10-01 00:30:00", tz = "UTC") & DateTime < as.POSIXct("2021-10-01 00:00:00", tz = "UTC")) 
result_RF_NEE <- result_RF_NEE %>% 
  filter(DateTime >= as.POSIXct("2015-10-01 00:30:00", tz = "UTC") & DateTime < as.POSIXct("2021-10-01 00:00:00", tz = "UTC"))





#####-------------------- CALCULATION for FLUX DATA -----------------#####

#Conversion factors
ch4_conv <- 12.01/(10^6)
co2_conv <- 12.01/(10^6) 
d.avg <- 1800 * 48  # 60s * 30min * 2 hours * 24 hours

BB_flux <- BB1 %>%
  select (c(DATE, date, jday, month_local, Year_local, time, period, season,
            Tau, ET, RH, u., jday, co2_flux, ch4_flux, NEE_f,
            Reco, GPP_f, GPP_DT, Reco_DT, FCH4_f, FCH4_gf_RF,
            NEE_f_RF, GPP_f_RF)) %>%
  dplyr::rename(year_local = Year_local)


#### DAILY NON-GAPFILLED FLUX DATA ####

# filter non gap-filled data to get >30% of half-hourly measurements
daily.nGF <- BB_flux %>%
  dplyr::select(c(year_local, jday, month_local, co2_flux, ch4_flux)) %>%
  drop_na()

obs_count <- daily.nGF %>% 
  group_by(year_local, jday) %>% 
  tally() %>% 
  ungroup()

daily.nGF <- daily.nGF %>%
  group_by(year_local, jday) %>% 
  summarize(co2 = mean(co2_flux),
            ch4 = mean(ch4_flux)) %>%
  mutate(co2gC = co2 * d.avg * co2_conv) %>%
  mutate(ch4gC = ch4* d.avg * co2_conv) %>%
  mutate(ch4mgC = ch4gC * 1000) %>% 
  ungroup() %>% 
  mutate(n = obs_count$n) %>% 
  # filter(n > 48 * 0.6) %>%
  mutate(date = as.POSIXct(paste(year_local, jday, sep = "-"), format = "%Y-%j"))

# filter non gap-filled data to get >30% of half-hourly measurements
daily.nGF[which(daily.nGF$n < 48 * 0.3), 2:7] <- NA



#### DAILY GAPFILLED FLUX DATA ####

daily <- BB_flux %>%
  ddply("date", summarize,
        NEE_f = mean(NEE_f),
        NEE_f_RF = mean(NEE_f_RF),
        GPP_DT = mean(GPP_DT),
        GPP_f = mean(GPP_f),
        GPP_f_RF = mean(GPP_f_RF),
        Reco = mean(Reco),
        Reco_DT = mean(Reco_DT),
        CH4_f_RF = mean(FCH4_gf_RF),
        CH4_f = mean(FCH4_f)) %>%
  mutate(NEEgC_f = NEE_f * d.avg * co2_conv, 
         NEEgC_f_RF = NEE_f_RF * d.avg * co2_conv,
         GPPgC_DT = GPP_DT * d.avg * co2_conv, 
         GPPgC_f = GPP_f *d.avg * co2_conv,
         GPPgC_f_RF = GPP_f_RF * d.avg * co2_conv,
         RecogC = Reco * d.avg * co2_conv, 
         RecogC_DT = Reco_DT * d.avg * co2_conv,
         CH4gC_f = CH4_f * d.avg * ch4_conv, 
         CH4mgC_f = CH4gC_f * 1000,
         CH4gC_f_RF = CH4_f_RF * d.avg * ch4_conv,
         CH4mgC_f_RF = CH4gC_f_RF * 1000,
         year_local = year(date), 
         month_local = month(date), 
         floor_date = floor_date(date, "month"))



# set year period
daily <- daily %>% set_period() %>% set_season(date = daily$date)



#####-------------------- CALCULATION FOR MET DATA ------------------#####


# DAILY MET DATA
BB_met <- BB1 %>%  
  rename("TS.5" = "SOIL_TEMP_5CM", 
         "TS.10" = "SOIL_TEMP_10CM",
         "TS.50" = "SOIL_TEMP_50CM") %>% 
  group_by(date) %>% 
  summarize(
    TS.5 = mean(TS.5, na.rm = T),  # soil temp
    TS.10 = mean(TS.10, na.rm = T),
    TS.50 = mean(TS.50, na.rm = T),
    RH = mean(RH_2M, na.rm = T),  # relative humidity
    PARin = mean(INCOMING_PAR, na.rm = T) * d.avg * 1e-6,  # mol m2 day -1 
    SWin = mean(SHORTWAVE_IN, na.rm = T),
    VPD.y = mean(VPD.y, na.rm = T),
    VPD.x = mean(VPD.x, na.rm = T),
    Precip = sum(PRECIP, na.rm = T),
    PA_EC = mean(PA_EC2_AIR, na.rm = T),  # barrometric pressure from EC
    PA_2M = mean(PA_2M, na.rm = T),  # barrometric pressure from met
    Ta = mean(AIR_TEMP_2M, na.rm = T),  # air temp at 2 m
    # WTH2 = mean(4.606 -(-0.89*WTH2+82)/100),  # water table height before 2018
    wtd = mean(WTH_absolute, na.rm = T),  # water Marion
    SWC = mean(SVWC, na.rm = T), # soil water content
    LE = mean(LE, na.rm = T), 
    H = mean(H, na.rm = T), 
    br = mean(bowen_ratio, na.rm = T), 
    ustar = mean(u., na.rm = T)
  ) %>%  
  mutate(month_local = month(date), 
         year_local = year(date), 
         floor_date = floor_date(date, "month"))%>%
  ungroup()


# set period
BB_met <- BB_met %>% set_period() %>% set_season(date = BB_met$date)




#####------------------------- WEEKLY DATA --------------------------####

##### daily non gapfilled flux data (averaged per week) #####

weekly.nGF <- BB_flux %>%
  dplyr::select(c(year_local, jday, month_local, co2_flux, ch4_flux, Reco, GPP_f_RF, date))
# 
# weekly.nGF[is.na(weekly.nGF$co2_flux) , "Reco"] <- NA
# weekly.nGF[is.na(weekly.nGF$co2_flux) , "GPP_f_RF"] <- NA


weekly.nGF <- weekly.nGF  %>% 
  mutate(floor_date = floor_date(date, "week", week_start =  getOption("lubridate.week.start", 4)))  ## week

obs_count <- weekly.nGF %>%
  group_by(floor_date) %>% 
  summarise(n = sum(!is.na(co2_flux)))


weekly.nGF <- weekly.nGF %>%
  group_by(floor_date) %>% 
  summarize(co2 = mean(co2_flux, na.rm = T),
            ch4 = mean(ch4_flux, na.rm = T),
            GPP_f_RF = mean(GPP_f_RF, na.rm = T), 
            Reco = mean(Reco, na.rm = T)) %>%
  mutate(co2gC = co2 * d.avg * co2_conv, 
         ch4gC = ch4 * d.avg  * ch4_conv,
         GPP_f_RF = GPP_f_RF * d.avg * co2_conv, 
         Reco = Reco * d.avg * co2_conv) %>%  # not multiplied by 7 because we want the unit to be /day
  ungroup()%>% 
  mutate(n = obs_count$n)

# filter non gap-filled data to get >30% of half-hourly measurements
weekly.nGF[which(weekly.nGF$n < 48 * 7 * 0.3), 2:ncol(weekly.nGF)] <- NA

percent_week <- (weekly.nGF %>% drop_na() %>% tally()) / nrow(weekly.nGF) * 100

print(paste("percent of weekly data points available in 6 years =", percent_week, "%"))


#### daily gapfilled flux data (averaged per week) #####

weekly.GF <- BB_flux %>%
  dplyr::select(c(NEE_f_RF, Reco, FCH4_gf_RF, GPP_f_RF, date)) %>% 
  mutate(floor_date = floor_date(date, "week", week_start =  getOption("lubridate.week.start", 4)))%>%
  group_by(floor_date) %>% 
  summarize(GPP_f_RF = mean(GPP_f_RF, na.rm = T), 
            NEE_f_RF = mean(NEE_f_RF, na.rm = T), 
            Reco = mean(Reco, na.rm = T), 
            CH4_f_RF = mean(FCH4_gf_RF, na.rm = T))%>%
  mutate(GPP_f_RF = GPP_f_RF * d.avg * co2_conv, 
         NEE_f_RF = NEE_f_RF * d.avg * co2_conv, 
         Reco = Reco * d.avg * co2_conv,
         CH4_f_RF = CH4_f_RF * d.avg * ch4_conv) %>%  # not multiplied by 7 because we want the unit to be /day
  ungroup()


##### weekly met data #####

# DAILY MET DATA (averaged per week)
weekly.met <- BB1 %>%  
  rename("TS.5" = "SOIL_TEMP_5CM", 
         "TS.10" = "SOIL_TEMP_10CM",
         "TS.50" = "SOIL_TEMP_50CM") %>% 
  mutate(floor_date = floor_date(date, "week", week_start =  getOption("lubridate.week.start", 4))) %>%
  group_by(floor_date) %>% 
  summarize(
    TS.5 = mean(TS.5, na.rm = T),  # soil temp
    TS.10 = mean(TS.10, na.rm = T),
    TS.50 = mean(TS.50, na.rm = T),
    RH = mean(RH_2M, na.rm = T),  # relative humidity
    PARin = mean(INCOMING_PAR, na.rm = T) * d.avg * 1e-6,  # mol m2 day -1 
    SWin = mean(SHORTWAVE_IN, na.rm = T),
    VPD.y = mean(VPD.y, na.rm = T),
    Precip = sum(PRECIP, na.rm = T),
    PA_2M = mean(PA_2M, na.rm = T),  # barrometric pressure from met
    Ta = mean(AIR_TEMP_2M, na.rm = T),  # air temp at 2 m
    wtd = mean(WTH_absolute, na.rm = T),  # water Marion
    SWC = mean(SVWC, na.rm = T), # soil water content
    LE = mean(LE, na.rm = T), 
    H = mean(H, na.rm = T), 
    ustar = mean(u., na.rm = T)
  ) %>%
  ungroup()
