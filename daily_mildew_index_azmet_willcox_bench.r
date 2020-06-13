

# This code generates an exploratory analysis of daily mildew index (dmi) values 
# based on historical temperature observations at the AZMET Willcox Bench 
# station

# The dmi_ref algorithm, 'Model 2 of 2', is described at:
# http://ipm.ucanr.edu/DISEASE/DATABASE/grapepowderymildew.html

# AZMET data are at: https://cals.arizona.edu/azmet/

# Author:
# Jeremy Weiss, Climate and Geospatial Extension Scientist
# School of Natural Resources and the Environment
# University of Arizona
# 520-626-8063, jlweiss@email.arizona.edu


# SETUP --------------------


# Load needed libraries
library("reshape2")
library("dplyr")
library("ggplot2")


library("tidyverse")
library("extrafont")
library("lubridate")


# Set the AZMET station name and years of interest (if different from actual
# station data coverage)
stn_name <- "Willcox Bench"
yr_start <- 2017
yr_end <- 2020

# Load AZMET station list
stn_list <- read.csv("azmet-station-list.csv", sep = ",")

# Load grape budbreak dates
budbreak_dates <- read.csv("budbreak-dates.csv")
budbreak_dates["date"] <- as.Date.character(paste(
  budbreak_dates$year, budbreak_dates$month, budbreak_dates$day), 
  format = "%Y %m %d")
budbreak_dates["doy"] <- as.numeric(format(budbreak_dates$date, "%j"))

# Load and transform dmi_ref values
dmi_ref <- read.csv("daily-mildew-index.csv")
colnames(dmi_ref) <- c("Tmin_range", 
                   # Tmax ranges
                   "55-60", "60-65", "65-70", "70-75", "75-80", "80-85",
                   "85-90", "90-95", "95-100", "100-105", "105-110")
dmi_ref <- melt(data = dmi_ref, na.rm = FALSE,value.name = "dmi_value")
colnames(dmi_ref)[2] <- "Tmax_range"

dmi_ref["Tmin_low"] <- as.character(dmi_ref$Tmin_range) %>%
  strsplit("-") %>%
  sapply("[", 1) %>%
  as.numeric()
dmi_ref["Tmin_high"] <- as.character(dmi_ref$Tmin_range) %>%
  strsplit("-") %>%
  sapply("[", 2) %>%
  as.numeric()
dmi_ref["Tmax_low"] <- as.character(dmi_ref$Tmax_range) %>%
  strsplit("-") %>%
  sapply("[", 1) %>%
  as.numeric()
dmi_ref["Tmax_high"] <- as.character(dmi_ref$Tmax_range) %>%
  strsplit("-") %>%
  sapply("[", 2) %>%
  as.numeric()

# Load function to download and transform daily AZMET data
source("azmet.daily.data.download.R")


# DOWNLOAD AND TRANSFORM DAILY AZMET DATA --------------------


stn_data <- azmet.daily.data.download(stn_name)

# Retain necessary variables and data
stn_data <- select(stn_data, date, year, month, day, doy, Tmax, Tmin, PRCtot)
stn_data <- filter(stn_data, year >= yr_start & year <= yr_end)

# Convert temperatures from Celsius to Fahrenheit, and precipitation from 
# millimeters to inches
stn_data$Tmax <- (1.8 * stn_data$Tmax) + 32
stn_data$Tmin <- (1.8 * stn_data$Tmin) + 32
stn_data$PRCtot <- stn_data$PRCtot / 25.4


# CALCULATE DMI --------------------


stn_data["dmi"] <- NA
stn_data["dmi_notes"] <- NA

# For nodata case
stn_data$dmi[
  which(is.na(stn_data$Tmin) == TRUE | is.na(stn_data$Tmax) == TRUE)
  ] <- NA

# For outside lower bounds of dmi reference table
stn_data$dmi[
  which(stn_data$Tmin <= min(dmi_ref$Tmin_low) | 
          stn_data$Tmax <= min(dmi_ref$Tmax_low))
  ] <- 0

# For within bounds of dmi reference table
breaks_tmin <- unique(c(dmi_ref$Tmin_low, dmi_ref$Tmin_high))
breaks_tmax <- unique(c(dmi_ref$Tmax_low, dmi_ref$Tmax_high))
for (btmin in 1:(length(breaks_tmin) - 1)) {
  for (btmax in 1:(length(breaks_tmax) - 1)) {
    
    stn_data$dmi[which(
      stn_data$Tmin >= breaks_tmin[btmin] & 
        stn_data$Tmin < breaks_tmin[btmin + 1] &
        stn_data$Tmax >= breaks_tmax[btmax] &
        stn_data$Tmax < breaks_tmax[btmax + 1]
        )] <- 
      dmi_ref$dmi_value[which(
        dmi_ref$Tmin_low == breaks_tmin[btmin] & 
          dmi_ref$Tmin_high == breaks_tmin[btmin + 1] &
          dmi_ref$Tmax_low == breaks_tmax[btmax] &
          dmi_ref$Tmax_high == breaks_tmax[btmax + 1]
      )]
    
  }
}
rm(btmin, btmax)

# For excessive heat and leaf burn issue
stn_data$dmi_notes[which(
  stn_data$Tmin >= 75 & stn_data$Tmin < 80 & stn_data$Tmax >= 105
  )] <- "excessive heat"
stn_data$dmi_notes[which(
  stn_data$Tmin >= 80 & stn_data$Tmax >= 95
)] <- "excessive heat"


# CALCULATE PMI --------------------


# PMI (Powdery Mildex Index) is based on accumulation of DMI values from the 
# start of the growing season, or budbreak.

# According to the model, the first dusting is twelve days after initial leaf 
# appearance or 6-inch shoot growth, whichever comes first. Subsequent dustings 
# should occur when the difference between the current PMI and the PMI on the 
# last dusting date equals or exceeds 1.0. When precipitation exceeds 0.10 inch, 
# the vineyard should be re-dusted.

stn_data["pmi"] <- NA
stn_data["dusting"] <- NA

for (yr in yr_start:yr_end) {
  iyr_entries <- which(stn_data$year == yr)
  gs_start <- min(budbreak_dates$doy[which(budbreak_dates$year == yr)])
  
  # PMI
  for(entry in iyr_entries) {
    if (stn_data$doy[entry] <= gs_start + 12) {
      stn_data$pmi[entry] <- 0
    } else {
      stn_data$pmi[entry] <- stn_data$dmi[entry] + stn_data$pmi[entry - 1]
    }
  }
  rm(entry)
  
  # Dusting
  for(entry in iyr_entries) {
    if (stn_data$doy[entry] < gs_start + 12) {
      stn_data$dusting[entry] <- 0
    } else if (stn_data$doy[entry] == gs_start + 12) {
      stn_data$dusting[entry] <- 1
    } else {
      ilast_dust <- max(which(stn_data$dusting == 1))
      if (stn_data$pmi[entry] - stn_data$pmi[ilast_dust] >= 1 |
          stn_data$PRCtot[entry] > 0.10) {
            stn_data$dusting[entry] <- 1
      } else {
        stn_data$dusting[entry] <- 0
      }
    
    }
  }
  rm(entry)
  
  rm(iyr_entries, gs_start)
}
rm(yr)


for (yr in yr_start:yr_end) {
  iyr_entries <- which(stn_data$year == yr)
  gs_start <- min(budbreak_dates$doy[which(budbreak_dates$year == yr)])
  
  for(entry in iyr_entries) {
    if (stn_data$doy[entry] <= gs_start + 12) {
      stn_data$pmi[entry] <- 0
    } else {
      stn_data$pmi[entry] <- stn_data$dmi[entry] + stn_data$pmi[entry - 1]
    }
  }
  rm(entry)
  
  rm(iyr_entries, gs_start)
}
rm(yr)

















# Load the .csv file that contains bud break dates observed at Buhl Memorial 
# Vineyard
bmv <- read.csv(
  "buhl-memorial-vineyard-budbreak-dates.csv", header = TRUE, sep = ","
 )

# Generate a new column in the bud break data dataframe that contains the day of 
# year for each row. This requires that the year, month, and day are first 
# combined into a new variable for 'date'. The 'yday' function is from the 
# library 'lubridate'.
bmv["date"] <- as.Date(paste(bmv$year, bmv$month, bmv$day, sep = "-"))
bmv["doy"] <- as.numeric(yday(bmv$date))


# LOAD AND TRANSFORM AZMET TEMPERATURE DATA  --------------------


# Data from: https://cals.arizona.edu/azmet/az-data.htm. The following is set 
# for data from the Willcox Bench station.
data_dir <- 
  "C:/Users/jlweiss/OneDrive - University of Arizona/Documents/R/datastore/azmet/"

# Data collection at the Willcox Bench station started in middle 2016 and 
# continues at present. 'obs_yrs' format matches data filename convention.
obs_yrs <- seq(15, 20)

# Data column names are from: https://cals.arizona.edu/azmet/raw2003.htm
col_names <- 
  c("Year", "JDay", "Hour", "Temp", "rh", "vpd", "sr", "prcp", "4sm", "20sm",
    "wavg", "wvm", "wvd", "wdsd", "mws", "reto", "avp", "dp")

# Load data, which is in files for individual years
for (iyr in 1:length(obs_yrs)) {
  obs <- read.table(paste0(data_dir, "09", obs_yrs[iyr], "rh.txt"),
                     header = FALSE,
                     sep = ',',
                     fill = TRUE)
  colnames(obs) <- col_names
  obs <- select(obs, Year, JDay, Hour, Temp)
  obs_date <- as.Date(obs$JDay-1, origin = paste0("20", obs_yrs[iyr], "-01-01"))
  obs["Date"] <- as.character(obs_date)
  obs["Month"] <- as.numeric(format(obs_date, "%m"))
  obs["Day"]  <- as.numeric(format(obs_date, "%d"))
  rm(obs_date)
  
  # Concatenate station data from different years
  if (iyr == 1) {
    stn_data <- obs
  }
  else {
    stn_data <- rbind(stn_data, obs)
  }
  
  rm(obs)
}

rm(col_names, data_dir, iyr)


# ANALYSIS WITH AZMET TEMPERATURE DATA  --------------------


# Set the start and end days-of-year for the period of interest for analysis by 
# non-leap-year values. These values should match the overall period of interest 
# (pertaining to both chill and heat units), in order to facilitate plotting 
# below. Heat units will start accumulating on December 1 for a given 
# winter year.
doy_start <- 305  # November 1
doy_end <- 181  # June 30

# Aggregate hourly data to daily data
stn_data_daily <- stn_data %>%
  group_by(Year, JDay, Date, Month, Day) %>%
  summarize(avgTempC = mean(Temp, na.rm = TRUE))


# Missing value, quick fix
x <- c(
  stn_data_daily$avgTempC[(which(stn_data_daily$avgTempC == 999)) - 1],
  stn_data_daily$avgTempC[(which(stn_data_daily$avgTempC == 999)) + 1]
)
stn_data_daily$avgTempC[which(stn_data_daily$avgTempC == 999)] <- mean(x)


# Convert daily average temperature data from degrees C to degrees F, as we will
# follow NPN Viz Tool approach to calculate GDDs: 
# https://www.usanpn.org/data/agdd_maps
stn_data_daily["avgTempF"] <- (1.8 * stn_data_daily$avgTempC) + 32

# Calculate growing degree days (GDDs), based on 'avgTempF'. If the Tmean value 
# for a given day is present and less than or equal to 'gdd_base_temp', GDD 
# equals 0. Otherwise, GDD equals the difference between the daily average 
# temperature and 'gdd_base_temp'.
stn_data_daily["GDDs"] <- NA

for (d in 1:nrow(stn_data_daily)) {
  if (stn_data_daily$avgTempF[d] < gdd_base_temp) {
    stn_data_daily$GDDs[d] <- 0
  } 
  else {
    stn_data_daily$GDDs[d] <- stn_data_daily$avgTempF[d] - gdd_base_temp
  }
}
rm(d)

# Calculate accumulated growing degree days (AGDDs) and put into new dataframe, 
# accounting for leap years. Start date for calculations is December 1 of a 
# given winter year.
years <- seq(
  as.numeric(paste0(20, obs_yrs[1])),
  as.numeric(paste0(20, obs_yrs[length(obs_yrs)])) 
 )

for (iyr in 1:(length(years) - 1)) {
  
  winter_yr1 <- years[iyr]
  winter_yr2 <- years[iyr+1]
  
  ##### CONDITION 1
  
  # Condition 1: when first calendar year of a winter year is a leap year
  if (leap_year(winter_yr1) == TRUE & leap_year(winter_yr2) == FALSE) {
    
    # Filter data by given winter year and setup to store AGDD data
    winter_data <- rbind(stn_data_daily %>%
                            filter(Year == winter_yr1 & JDay >= doy_start+1),
                         stn_data_daily %>%
                            filter(Year == winter_yr2 & JDay <= doy_end))
    
    winter_data["Winter"] <- paste(years[iyr], years[iyr+1], sep = "-")
    winter_data["iDay"] <- seq(1, nrow(winter_data))
    winter_data["AGDDs"] <- NA
    
    # Calculate AGDD values for the given winter year
    for (d in 1:nrow(winter_data)) {
      # November values
      if (
        winter_data$JDay[d] >= doy_start+1 & winter_data$JDay[d] <= doy_start+30
        ) {
        winter_data$AGDDs[d] <- NA
      } 
      # December 1 value
      else if (winter_data$JDay[d] == doy_start+30+1) {
        winter_data$AGDDs[d] <- winter_data$GDDs[d]
      } 
      # all other day values
      else {
        winter_data$AGDDs[d] <- winter_data$GDDs[d] + winter_data$AGDDs[d-1]
      }
    }
    
  }
  
  ##### CONDITION 2
  
  # Condition 2: when neither calendar years of a winter year are leap years
  else if (leap_year(winter_yr1) == FALSE & leap_year(winter_yr2) == FALSE) {
    
    # Filter data by given winter year and setup to store AGDD data
    winter_data <- rbind(stn_data_daily %>%
                            filter(Year == winter_yr1 & JDay >= doy_start),
                         stn_data_daily %>%
                            filter(Year == winter_yr2 & JDay <= doy_end))
    
    winter_data["Winter"] <- paste(years[iyr], years[iyr+1], sep = "-")
    winter_data["iDay"] <- seq(1, nrow(winter_data))
    winter_data["AGDDs"] <- NA
    
    # Calculate AGDD values for the given winter year
    for (d in 1:nrow(winter_data)) {
      # November values
      if (
        winter_data$JDay[d] >= doy_start & winter_data$JDay[d] <= doy_start+30-1
        ) {
        winter_data$AGDDs[d] = NA
      } 
      # December 1 value
      else if (winter_data$JDay[d] == doy_start+30) {
        winter_data$AGDDs[d] <- winter_data$GDDs[d]  
      } 
      # all other day values
      else {
        winter_data$AGDDs[d] <- winter_data$GDDs[d] + winter_data$AGDDs[d-1]
      }
    }
    
  }
  
  ##### CONDITION 3
  
  # Condition 3: when second calendar year of a winter year is a leap year
  else if (leap_year(winter_yr1) == FALSE & leap_year(winter_yr2) == TRUE) {
    
    # Filter data by given winter year and setup to store AGDD data
    winter_data <- rbind(stn_data_daily %>%
                            filter(Year == winter_yr1 & JDay >= doy_start),
                         stn_data_daily %>%
                            filter(Year == winter_yr2 & JDay <= doy_end+1))
    
    winter_data["Winter"] <- paste(years[iyr], years[iyr+1], sep = "-")
    winter_data["iDay"] <- seq(1, nrow(winter_data))
    winter_data["AGDDs"] <- NA
    
    # Calculate AGDD values for the given winter year
    for (d in 1:nrow(winter_data)) {
      # November values
      if (
        winter_data$JDay[d] >= doy_start & winter_data$JDay[d] <= doy_start+30-1
        ) {
        winter_data$AGDDs[d] = NA
      } 
      # December 1 value
      else if (winter_data$JDay[d] == doy_start+30) {
        winter_data$AGDDs[d] <- winter_data$GDDs[d]  
      } 
      # all other day values
      else {
        winter_data$AGDDs[d] <- winter_data$GDDs[d] + winter_data$AGDDs[d-1]
      }
    }
    
  }
  
  rm(d)
  
  # Concatenate AGDD data from different winters
  if (iyr == 1) {
    agdd <- winter_data
  } 
  else {
    agdd <- rbind(agdd, winter_data)
  }
  
  rm(winter_data, winter_yr1, winter_yr2)
  
}
rm(iyr)



# NEEDS CLEANING UP


#####

#  calculate daily chill portions for winter 2015-2016
stn_data_x <- filter( stn_data, Year == 2015 | Year == 2016 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 305,  # non-leap year value for November 1
               End_JDay = doy_end )  #  ***UPDATE***

# initialize and label columns of dataframe
chill1516 <- data.frame( matrix( data = NA,
                                 nrow = x$Data_days[ which( x$Season == "2015/2016" ) ],
                                 ncol = 7 ) )
colnames( chill1516 ) <- c( "Winter",
                            "iDay",
                            "JDay",
                            "Chilling_Hours",
                            "Utah_Model",
                            "Chill_portions",
                            "GDH" )
chill1516$Winter <- "2015-2016"
rm( x )

#  calculate daily chill portions for 2015 and write to 
#  dataframe
d <- 1
for ( jday in 305:365 ) {
  chill1516$iDay[ d ] <- d
  chill1516$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1516$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2015 ) ]
  chill1516$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2015 ) ]
  chill1516$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2015 ) ]
  chill1516$GDH[ d ] <- x$GDH[ which( x$End_year == 2015 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2016 and write to 
#  dataframe
for ( jday in 1:doy_end ) {  #  ***UPDATE***
  chill1516$iDay[ d ] <- d
  chill1516$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1516$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2016 ) ]
  chill1516$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2016 ) ]
  chill1516$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2016 ) ]
  chill1516$GDH[ d ] <- x$GDH[ which( x$End_year == 2016 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )







##### 

#  calculate daily chill portions for winter 2016-2017
stn_data_x <- filter( stn_data, Year == 2016 | Year == 2017 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 306,  # leap year value for November 1
               End_JDay = doy_end )

# initialize and label columns of dataframe
chill1617 <- data.frame( matrix( data = NA,
                                 nrow = x$Data_days[ which( x$Season == "2016/2017" ) ],
                                 ncol = 7 ) )
colnames( chill1617 ) <- c( "Winter",
                            "iDay",
                            "JDay",
                            "Chilling_Hours",
                            "Utah_Model",
                            "Chill_portions",
                            "GDH" )
chill1617$Winter <- "2016-2017"
rm( x )

#  calculate daily chill portions for 2016 and write to 
#  dataframe
d <- 1
for ( jday in 306:366 ) {
  chill1617$iDay[ d ] <- d
  chill1617$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 305,  # leap year value for October 31
                 End_JDay = jday )
  chill1617$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2016 ) ]
  chill1617$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2016 ) ]
  chill1617$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2016 ) ]
  chill1617$GDH[ d ] <- x$GDH[ which( x$End_year == 2016 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2017 and write to 
#  dataframe
for ( jday in 1:doy_end ) {
  chill1617$iDay[ d ] <- d
  chill1617$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 305,  # leap year value for October 31
                 End_JDay = jday )
  chill1617$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2017 ) ]
  chill1617$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2017 ) ]
  chill1617$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2017 ) ]
  chill1617$GDH[ d ] <- x$GDH[ which( x$End_year == 2017 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )



#####

#  calculate daily chill portions for winter 2017-2018
stn_data_x <- filter( stn_data, Year == 2017 | Year == 2018 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 305,  # non-leap year value for November 1
               End_JDay = doy_end )

# initialize and label columns of dataframe
chill1718 <- data.frame( matrix( data = NA,
                                 nrow = x$Data_days[ which( x$Season == "2017/2018" ) ],
                                 ncol = 7 ) )
colnames( chill1718 ) <- c( "Winter",
                            "iDay",
                            "JDay",
                            "Chilling_Hours",
                            "Utah_Model",
                            "Chill_portions",
                            "GDH" )
chill1718$Winter <- "2017-2018"
rm( x )

#  calculate daily chill portions for 2017 and write to 
#  dataframe
d <- 1
for ( jday in 305:365 ) {
  chill1718$iDay[ d ] <- d
  chill1718$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1718$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2017 ) ]
  chill1718$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2017 ) ]
  chill1718$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2017 ) ]
  chill1718$GDH[ d ] <- x$GDH[ which( x$End_year == 2017 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2018 and write to 
#  dataframe
for ( jday in 1:doy_end ) {
  chill1718$iDay[ d ] <- d
  chill1718$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1718$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2018 ) ]
  chill1718$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2018 ) ]
  chill1718$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2018 ) ]
  chill1718$GDH[ d ] <- x$GDH[ which( x$End_year == 2018 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )


#####

#  calculate daily chill portions for winter 2018-2019
stn_data_x <- filter( stn_data, Year == 2018 | Year == 2019 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 305,  # non-leap year value for November 1
               End_JDay = doy_end )

# initialize and label columns of dataframe
chill1819 <- data.frame( matrix( data = NA,
                                 nrow = x$Data_days[ which( x$Season == "2018/2019" ) ],
                                 ncol = 7 ) )
colnames( chill1819 ) <- c( "Winter",
                            "iDay",
                            "JDay",
                            "Chilling_Hours",
                            "Utah_Model",
                            "Chill_portions",
                            "GDH" )
chill1819$Winter <- "2018-2019"
rm( x )

#  calculate daily chill portions for 2018 and write to 
#  dataframe
d <- 1
for ( jday in 305:365 ) {
  chill1819$iDay[ d ] <- d
  chill1819$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1819$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2018 ) ]
  chill1819$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2018 ) ]
  chill1819$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2018 ) ]
  chill1819$GDH[ d ] <- x$GDH[ which( x$End_year == 2018 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2019 and write to 
#  dataframe
for ( jday in 1:doy_end ) {
  chill1819$iDay[ d ] <- d
  chill1819$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1819$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2019 ) ]
  chill1819$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2019 ) ]
  chill1819$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2019 ) ]
  chill1819$GDH[ d ] <- x$GDH[ which( x$End_year == 2019 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )


#####

#  calculate daily chill portions for winter 2019-2020
stn_data_x <- filter( stn_data, Year == 2019 | Year == 2020 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 305,  # non-leap year value for November 1
               End_JDay = doy_end )  #  ***UPDATE***

# initialize and label columns of dataframe
chill1920 <- data.frame( matrix( data = NA,
                                 nrow = x$Data_days[ which( x$Season == "2019/2020" ) ],
                                 ncol = 7 ) )
colnames( chill1920 ) <- c( "Winter",
                            "iDay",
                            "JDay",
                            "Chilling_Hours",
                            "Utah_Model",
                            "Chill_portions",
                            "GDH" )
chill1920$Winter <- "2019-2020"
rm( x )

#  calculate daily chill portions for 2019 and write to 
#  dataframe
d <- 1
for ( jday in 305:365 ) {
  chill1920$iDay[ d ] <- d
  chill1920$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1920$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2019 ) ]
  chill1920$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2019 ) ]
  chill1920$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2019 ) ]
  chill1920$GDH[ d ] <- x$GDH[ which( x$End_year == 2019 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2020 and write to 
#  dataframe
for ( jday in 1:141 ) {  #  ***UPDATE***
  chill1920$iDay[ d ] <- d
  chill1920$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1920$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2020 ) ]
  chill1920$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2020 ) ]
  chill1920$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2020 ) ]
  chill1920$GDH[ d ] <- x$GDH[ which( x$End_year == 2020 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )


#####


#  Concatenate chill portions data from different winters
cp <- rbind(chill1516, chill1617, chill1718, chill1819, chill1920)
rm(chill1516, chill1617, chill1718, chill1819, chill1920)


# PLOT CHILL AND HEAT ACCUMULATIONS AT BUD BREAK  --------------------


# Prepare 'bmv' dataframe
colnames(bmv) <- 
  c("Section", "Block", "Rows", "Year", "Month", "Day", "Date", "JDay")

bmv["Winter"] <- NA
bmv["AGDDs"] <- NA
bmv["CPs"] <- NA

# Construct 'bmv$Winter' entries
for (d in 1:nrow(bmv)) {
  # First calendar year in winter year case
  if (bmv$JDay[d] >= doy_start) {
    bmv$Winter[d] <- 
      paste(as.character(bmv$Year[d]), as.character(bmv$Year[d] + 1), sep = "-")
  }
  # Second calendar year in winter year case
  else {
    bmv$Winter[d] <- 
      paste(as.character(bmv$Year[d] - 1), as.character(bmv$Year[d]), sep = "-")
  }
}
rm(d)

# Match the date between the individual bud break observations and corresponding 
# AGDDs and chill portion values
for (obs in 1:nrow(bmv)) {
  
  # AGDDs
  if (length(which(agdd$Date == bmv$Date[obs], arr.ind = TRUE)) > 0) {
    bmv$AGDDs[obs] <- agdd$AGDDs[
      which(agdd$Date == bmv$Date[obs], arr.ind = TRUE)
      ]
  }
  
  # Chill portions
  if (length(which(cp$Winter == bmv$Winter[obs] & cp$JDay == bmv$JDay[obs],
                     arr.ind = TRUE)) > 0) {
    bmv$CPs[obs] <- cp$Chill_portions[
      which(cp$Winter == bmv$Winter[obs] & cp$JDay == bmv$JDay[obs], 
            arr.ind = TRUE)
      ]
  }
  
}
rm(obs)  

# Load additional font options
font_import()
y
loadfonts(device = "win")

# Create a 'ggplot' object for the data
p <- ggplot(
  data = filter(bmv, Winter != "2015-2016"), aes(x = CPs, y = AGDDs)
  ) +
  
  # Add data for individual years as differently colored circles for data
  geom_point(data = filter(bmv, Winter == "2015-2016"),  # circle fill
             alpha = 0.5, color = "#400639", size = 5, stroke = 1) +
  geom_point(data = filter(bmv, Winter == "2015-2016"),  # circle outline
             alpha = 1.0, pch = 21, color = "#400639", size = 5, stroke = 0.5) +
  
  geom_point(data = filter(bmv, Winter == "2016-2017"),  # circle fill
             alpha = 0.5, color = "#E38C27", size = 5, stroke = 1) +
  geom_point(data = filter(bmv, Winter == "2016-2017"),  # circle outline
             alpha = 1.0, pch = 21, color = "#E38C27", size = 5, stroke = 0.5) +
  
  geom_point(data = filter(bmv, Winter == "2017-2018"),  # circle fill
             alpha = 0.5, color = "#66584A", size = 5, stroke = 1) +
  geom_point(data = filter(bmv, Winter == "2017-2018"),  # circle outline
             alpha = 1.0, pch = 21, color = "#66584A", size = 5, stroke = 0.5) +
  
  geom_point(data = filter(bmv, Winter == "2018-2019"),  # circle fill
             alpha = 0.5, color = "#8F1915", size = 5, stroke = 1) +
  geom_point(data = filter(bmv, Winter == "2018-2019"),  # circle outline
             alpha = 1.0, pch = 21, color = "#8F1915", size = 5, stroke = 0.5) +
  
  geom_point(data = filter(bmv, Winter == "2019-2020"),  # circle fill
             alpha = 0.5, color = "#4BA692", size = 5, stroke = 1) +
  geom_point(data = filter(bmv, Winter == "2019-2020"),  # circle outline
             alpha = 1.0, pch = 21, color = "#4BA692", size = 5, stroke = 0.5) +
  
  # Add the graph title, subtitle, and axis labels
  ggtitle("Chill and Heat Accumulations at Bud Break, \n2016-2020") +
  labs(subtitle = "Buhl Memorial Vineyard, Cochise County, Arizona",
       x = "\nChill Portions",
       y = "Cumulative  Growing  Degree  Days\n") +
  
  # Specify axis breaks, gridlines, and limits
  scale_x_continuous(
    breaks = seq(0, max(bmv$CPs, na.rm = TRUE), 5),
    limits = c(min(bmv$CPs, na.rm = TRUE), max(bmv$CPs, na.rm = TRUE))
    ) +
  scale_y_continuous(
    breaks = seq(0, max(bmv$AGDDs, na.rm = TRUE), 100),
    limits = c(min(bmv$AGDDs, na.rm = TRUE), max(bmv$AGDDs, na.rm = TRUE))
  ) +
  
  # Further customize the figure appearance
  theme_light(base_family = "Source Sans Pro") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(color = "gray40", size = 10),
        axis.text.y = element_text(color = "gray40", size = 10),
        axis.ticks.length = unit(0.0, "mm"),
        axis.title.x = element_text(color = "gray40", size = 10),
        axis.title.y = element_text(color = "gray40", size = 10),
        legend.title = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_line(color = "gray80", size = 0.25),
        panel.grid.major.y = element_line(color = "gray80", size = 0.25),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color = "gray80", size = 0.25),
        plot.margin = unit(c(1, 1 ,1, 1), "mm"),
        plot.subtitle = (element_text(family = "Source Serif Pro", size = 12)), 
        plot.title = (element_text(face = "bold", family = "Source Serif Pro", size = 16))
  )
  
p

# Save the figure as a .png file in a current sub-directory
ggsave(
  "./figures/buhl-memorial-vineyard-bud-break-cgdds-cps-2016-2020-202005-sub-azmet-bonita.png",
  plot = p, device = "png", path = NULL, scale = 1,
  width = 6, height = 4, units = "in", dpi = 300
  ) 


# EXTRA CODE  --------------------


# Determine length of bud break period by year
for (yr in years){
  x <- filter(bmv, Year == yr)
  len <- max(x$JDay,na.rm = TRUE) - min(x$JDay, na.rm = TRUE) 
  print(paste0("Length of bud break period in ",yr,": ",len))
  flush.console()
}
rm(yr)

#
agdd["Chill_portions"] <- NA

# Match the date between corresponding AGDDs and chill portion values
for (obs in 1:nrow(agdd)) {
  
  if (length(which(cp$Winter == agdd$Winter[obs] & cp$JDay == agdd$JDay[obs],
                   arr.ind = TRUE)) > 0) {
    agdd$Chill_portions[obs] <- cp$Chill_portions[
      which(cp$Winter == agdd$Winter[obs] & cp$JDay == agdd$JDay[obs], 
            arr.ind = TRUE)
      ]
  }
  
}
rm(obs)
  
# Create a 'ggplot' object for the data
ggplot(data = filter(agdd, Winter != "2015-2016"),
       aes(x = Chill_portions, y = AGDDs, color = Winter)) +
  geom_point()

