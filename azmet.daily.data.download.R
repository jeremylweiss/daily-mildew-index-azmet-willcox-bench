

# FUNCTION FOR DOWNLOADING AND FORMATTING DAILY AZMET DATA 

# Authors:
# Jeremy Weiss, Climate and Geospatial Extension Scientist
# School of Natural Resources and the Environment
# University of Arizona
# 520-626-8063, jlweiss@email.arizona.edu
#
# Michael Crimmins, Climate Science Extension Specialist
# Department of Soil, Water, and Environmental Science
# University of Arizona
# 520-626-4244, crimmins@email.arizona.edu

# This function downloads daily AZMET data for an individual station, formats
# data into a dataframe, checks for missing or duplicate dates or other
# oddities, and writes the station data dataframe to the current environment


azmet.daily.data.download <- function(stn_name) {
  
  
  # SETUP  -------------------- 
  
  
  # AZMET data format changes between the periods 1987-2002 and 2003-present, as 
  # the number of variables measured / reported and their order in the data file 
  # are slightly different. We will set up a column name string that matches the 
  # variables and variable order of the latter period.
  
  # Set column name string for the 2003-present period. This list can be found 
  # at http://ag.arizona.edu/azmet/raw2003.htm. Note that the soil temperature 
  # depths change between the 1987-2002 and 2003-present periods. We use the 
  # depths from the latter to name these columns instead of generating new 
  # columns for the different depths between the two periods. As we do not 
  # anticipate using soil temperature data, this is of no consequence. However, 
  # this code will need to be changed in order to address this issue if soil 
  # temperature data becomes of interest. 
  col_names <- c("year", "doy", "stn_no", "Tmax", "Tmin", "Tmean", "RHmax", 
                 "RHmin", "RHmean", "VPDmean", "SORADtot", "PRCtot", "4STmax",
                 "4STmin", "4STmean", "20STmax", "20STmin", "20STmean", 
                 "WSmean", "WVmag", "WVdir", "Wdirstd", "WSmax", "HU8555",
                 "ETref", "ETrefPM", "AVPmean", "DPTmean")
  
  # Set the string elements that together will build the full URL where 
  # individual AZMET station data reside. Note that daily station data are 
  # available by individual years.
  
  # Extract the row of information (station name, station number, start year, 
  # and end year) tied to the selected AZMET station
  stn_info <- subset(x = stn_list, subset = stn == stn_name )
  
  # Set the station number based on the information extracted from 'stn_list' in 
  # the previous command. The station number will need to be converted to a 
  # character string in order to be put together with the other full URL string 
  # elements. Also, if the station number is less than 10, the station number 
  # character string will need to have a '0' preceeding it, in order to match
  # the AZMET daily data file name format.
  stn_no <- as.character(select(stn_info, stn_no))
  if (as.integer(select(stn_info, stn_no)) < 10) { 
    stn_no <- paste0("0", stn_no)
  }
  
  # Set the range of years for which to download data for the selected station
  stn_yrs <- as.integer(select(stn_info, start_yr)):
    as.integer(select(stn_info, end_yr))
  
  # Set the base URL of the AZMET data
  baseurl <- "http://ag.arizona.edu/azmet/data/"
  
  # Set the suffix of the data file to be downloaded
  suffix <- "rd.txt"
  
  
  # DOWNLOAD DATA  -------------------- 
  
  
  # Recall that AZMET data are provided year-by-year. We will need to 
  # iteratively download the annual files.
  
  # Loop through the 'stn_yrs' integer vector in order to build the full URL 
  # where the AZMET daily data for individual years reside. We will treat the 
  # 1987-2002 and 2003-present periods differently within the loop.
  for (i in 1:length(stn_yrs)) {
    
    # Set the data URL
    url <- paste0(baseurl,
                  stn_no,
                  substr(as.character(stn_yrs[i]), 3, 4 ),
                  suffix)
    
    # Test for the condition of a year falling in the 1987-2002 period. If true, 
    # switch the last two columns and add three new empty columns after the 
    # existing columns that will contain new variables that start in 2003. These 
    # changes are described at http://ag.arizona.edu/azmet/raw2003.htm.
    if (stn_yrs[i] <= 2002) {
      ann_data <- read.table(url, header = FALSE, sep = ',', fill = TRUE)
      ann_data <- ann_data[c(1:23, 25, 24)]
      ann_data[, 26:28] <- NA
    } else {
      # If the year falls in the 2003-present period instead, simply read in the 
      # data as is
      ann_data <- read.table(url, header = FALSE, sep = ',', fill = TRUE )
    }
    
    # Years prior to 2000 are to be two-digit values instead of four-digit 
    # values. Overwrite the first column for all years with four-digit values.
    ann_data[,1] <- rep(stn_yrs[i], nrow(ann_data))
    
    # Concatenate the data in the row dimension as it is downloaded year-by-year
    if (i == 1) {
      stn_data <- ann_data
    } else {
      stn_data <- rbind(stn_data, ann_data)
    }
    
  }
  rm(i)
  
  
  # FORMAT DATA  --------------------
  
  
  # Set the column names for the downloaded data
  colnames(stn_data) <- col_names
  
  # Populate new 'date', 'month', and 'day' columns
  stn_data["date"] <- as.Date.character(paste(stn_data$year, stn_data$doy),
                                        format = "%Y %j")
  stn_data["month"] <- as.numeric(format(stn_data$date, "%m"))
  stn_data["day"] <- as.numeric(format(stn_data$date, "%d"))
  
  
  
  
  
  
  # Based on previous work with AZMET data, there are several known
  # formatting bugs in the original / downloaded data files. We will
  # address these individually.
  
  # An odd character (".") appears at the end of some data files for 
  # some years and some stations. In the R dataframe, this results
  # in a row of NAs. Find and remove these rows.
  stn_data <- stn_data[ rowSums( is.na( stn_data ) ) != ncol( stn_data ), ]
  
  # Replace 'nodata' values in the downloaded AZMET data with 'NA'.
  # Values for 'nodata' in AZMET data are designated as '999'.
  # However, other similar values also appear (e.g., 999.9).
  stn_data[ stn_data == 999 ] <- NA
  stn_data[ stn_data == 999.9 ] <- NA
  stn_data[ stn_data == 9999 ] <- NA
  
  # Some AZMET data files have duplicate row entries. Find and
  # remove these rows.
  stn_data <- distinct( stn_data )
  
  
  ##### ADDRESS MISSING VALUES
  
  
  # Create a dataframe that mimics the actual station data, but
  # has a full YYYYMMDD list. AZMET data can have missing date
  # entries.
  
  # Create a new object that fully expands dates between the 
  # first and last dates in the station data. 
  date_full <- seq( first( stn_data$date ),last( stn_data$date ),by="days" )
  
  # Convert this new object to a dataframe that will allow us to
  # join with the station data. Dataframe column names must match
  # the column names in the station data in order to join.
  stn_data_full <- data.frame( matrix( NA,
                                       nrow=length( date_full ) ) )
  colnames( stn_data_full ) <- "date"
  stn_data_full$date <- date_full
  
  # Join the complete dates dataframe with the station data by
  # using 'date' as the key.
  stn_data_full <- left_join( stn_data_full,
                              stn_data,
                              by="date" )
  
  stn_data <- as_tibble( stn_data_full )
  rm( stn_data_full )
  
  # Fill in values for year, month, day, and day-of-year for any
  # date entries that may be missing in the downloaded original 
  # data.
  stn_data$year <- as.numeric( format( stn_data$date,"%Y" ) )
  stn_data$month <- as.numeric( format( stn_data$date,"%m" ) )
  stn_data$day <- as.numeric( format( stn_data$date,"%d" ) )
  stn_data$doy <- as.numeric( format( stn_data$date,"%j" ) )
  
  # Do similarly with the station number value.
  stn_data$stn_no <- stn_info$stn_no
  
  
  ##### RETURN THE DATA AND CLOSE THE FUNCTION
  
  
  return( stn_data )
}



