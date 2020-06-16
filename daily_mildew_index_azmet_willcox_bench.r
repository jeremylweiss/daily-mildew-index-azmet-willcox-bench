

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
library("lubridate")
library("extrafont")
library("ggplot2")

# Load additional font options for plotting
font_import()
y
loadfonts(device = "postscript")

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


# CALCULATE PMI & DUSTING --------------------


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


# PLOT --------------------


# Make PMI values earlier than 12 days prior to the first dusting NAs to more 
# accurately depict accumulation
for (yr in yr_start:yr_end) {
  iyr_entries <- which(stn_data$year == yr)
  gs_start <- min(budbreak_dates$doy[which(budbreak_dates$year == yr)])
  
  for(entry in iyr_entries) {
    if (stn_data$doy[entry] < gs_start) {
      stn_data$pmi[entry] <- NA
    }
  }
  rm(entry)
  
  rm(iyr_entries, gs_start)
}
rm(yr)

# For leap years, subtract 1 from the doy values in order to line up data as a 
# function of calendar day, which works here as plot will start after February
stn_data["doy_plot"] <- stn_data$doy
for (yr in yr_start:yr_end) {
  if (leap_year(yr) == TRUE) {
    iyr_entries <- which(stn_data$year == yr)
    for(entry in iyr_entries) {
      stn_data$doy_plot[entry] <- stn_data$doy[entry] - 1 
    }
    rm(entry)
  }
}
rm(yr)

# Make new variables to help plot dusting dates
stn_data["dusting_dates"] <- NA
stn_data["dusting_code"] <- NA
for (i in 1:length(stn_data$dusting)) {
  if (stn_data$dusting[i] == 1) {
    stn_data$dusting_dates[i] <- stn_data$pmi[i]
  }
  if (stn_data$dusting[i] == 1 & stn_data$PRCtot[i] > 0.10) {
    # precipitation forcing
    stn_data$dusting_code[i] <- 2
  } else if (stn_data$dusting[i] == 1) {
    # temperature forcing
    stn_data$dusting_code[i] <- 1
  }
}
rm(i)

# Create a 'ggplot' object
p <- ggplot() +
  
  # Add PMI data as time series lines with circles for dusting dates coded by 
  # temperature or precipitation forcing, by individual years
  geom_line(data = stn_data, 
            aes(x = doy_plot, y = pmi), 
            color = "gray60", size = 1) +
  geom_point(data = filter(stn_data, dusting_code == 1),
             aes(x = doy_plot, y = dusting_dates, size = 2),
             alpha = 0.5, color = "#f46d43") +
  geom_point(data = filter(stn_data, dusting_code == 2),
             aes(x = doy_plot, y = dusting_dates, size = 2),
             alpha = 0.5, color = "#74add1") +
  facet_wrap(~ year, ncol = 1) +
  
  # Specify axis breaks, gridlines, and limits
  scale_x_continuous(
    #breaks = c(1, 15, 32, 46, 60, 74, 
    #           91, 105, 121, 135, 152, 166,
    #           182, 196, 213, 227, 244),
    #labels = c("Jan 1", "Jan 15", "Feb 1", "Feb 15", "Mar 1", "Mar 15",
    #           "Apr 1", "Apr 15", "May 1", "May 15", "Jun 1", "Jun 15",
    #           "Jul 1", "Jul 15", "Aug 1", "Aug 15", "Sep 1"),
    breaks = c(1, 8, 15, 22, 32, 39, 46, 53, 
               60, 67, 74, 81, 91, 98, 105, 112,
               121, 128, 135, 142, 152, 159, 166, 173,
               182, 189, 196, 203, 213, 220, 227, 234, 244),
    labels = c("1/1", "", "1/15", "", "2/1", "", "2/15", "", 
               "3/1", "", "3/15", "", "4/1", "", "4/15", "",
               "5/1", "", "5/15", "", "6/1", "", "6/15", "",
               "7/1", "", "7/15", "", "8/1", "", "8/15", "", "9/1"),
    limits = c((min(filter(budbreak_dates, year >= yr_start)$doy) - 1), 245),
    #minor_breaks = c(1, 8, 15, 22, 32, 39, 46, 53, 60, 67, 74, 81, 
    #                 91, 98, 105, 112, 121, 128, 135, 142, 152, 159, 166, 173,
    #                 182, 189, 196, 203, 213, 220, 227, 234, 244),
    expand = c(0.0, 0.0)
    ) +
  
  scale_y_continuous(
    breaks = seq(from = 0, to = max(stn_data$pmi, na.rm = TRUE), by = 4),
    limits = c(0, max(filter(stn_data, doy_plot <= 245)$pmi, na.rm = TRUE)),
    minor_breaks = seq(from = 0, to = max(stn_data$pmi, na.rm = TRUE), by = 1),
    expand = c(0.06, 0.0)
  ) +
  
  # Add the graph title, subtitle, and axis labels
  ggtitle("Model-based Timing of Sulfur Dustings \nfor Powdery Mildew, 2017-2020") +
  labs(subtitle = "Southcentral Willcox AVA, Arizona",
       x = "\nDate",
       y = "Powdery Mildew Index\n",
       caption = paste0("\nmodel: UC IPM model 2 (ipm.ucanr.edu/DISEASE/DATABASE/grapepowderymildew.html)",
                        "\nmeteorological data: AZMET Willcox Bench station (cals.arizona.edu/azmet)",
                        "\nbudbreak data: Jesse Noble, Merkin Vineyards")) +
  
  # Further customize the figure appearance
  theme_light(base_family = "Source Sans Pro") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(color = "gray40", size = 10),
        axis.text.y = element_text(color = "gray40", size = 10),
        axis.ticks.x.bottom = element_line(color = "gray80", size = 0.25),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0.0, "mm"),
        axis.ticks.length.y = unit(0.0, "mm"),
        axis.title.x = element_text(color = "gray40", size = 10),
        axis.title.y = element_text(color = "gray40", size = 10),
        legend.title = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_line(color = "gray80", size = 0.25),
        panel.grid.major.y = element_line(color = "gray80", size = 0.25),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color = "gray80", size = 0.25),
        plot.caption = element_text(color = "gray40", hjust = 0.0, size = 7),
        plot.caption.position = "plot",
        plot.margin = unit(c(1, 1 ,1, 1), "mm"),
        plot.subtitle = (element_text(family = "Source Serif Pro", size = 12)), 
        plot.title = (element_text(face = "bold", family = "Source Serif Pro", size = 16)),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(color = "gray40", size = 12, face = "bold")
  )

#  Save the figure
ggsave(file = "daily_mildew_index_azmet_willcox_bench-20200615.eps",
       plot = p, device = cairo_pdf, path = NULL, scale = 1,
       width = 6, height = 9, units = "in", dpi = 300) 

