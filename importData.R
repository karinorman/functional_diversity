library(tidyverse)
library(stringr)

#### Get scripts from Weecology bbs import ####

source_https <- function(u, unlink.tmp.certs = FALSE) {
  # load package
  require(RCurl)
  
  # read script lines from website using a security certificate
  if(!file.exists("cacert.pem")) download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.pem")
  script <- getURL(u, followlocation = TRUE, cainfo = "cacert.pem")
  if(unlink.tmp.certs) unlink("cacert.pem")
  
  # parase lines and evealuate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}

#create own get_species_data() function to use DBI/table method for connecting to database
get_species_data = function() {
  data_path <- paste('./data/', 'bbs', '_species.csv', sep = "")
  if (file.exists(data_path)) {
    return(read.csv(data_path))
  }else{
    write.csv(species, file = data_path, row.names = FALSE, quote = FALSE)
    return(species)
  }
}

source_https("https://raw.githubusercontent.com/weecology/bbs-forecasting/master/R/forecast-bbs-core.R")
source_https("https://raw.githubusercontent.com/weecology/bbs-forecasting/master/R/save_provenance.R")

#Adapted from get_bbs_data() from sourced scripts
get_bbs <- function(){
  data_path <- paste('./data/', 'bbs', '_data.csv', sep="")
  if (file.exists(data_path)){
    return(read_csv(data_path))
  }
  else{
    
    if (!db_engine(action='check', db = "~/Dropbox/Data/functional-diversity/bbsforecasting.sqlite",
                   table_to_check = 'breed_bird_survey_counts')){
      install_dataset('breed-bird-survey')
    }
    
    birds <- DBI::dbConnect(RSQLite::SQLite(), "~/Dropbox/Data/functional-diversity/bbsforecasting_old.sqlite")
    
    #save database tables as table in R to use with tidyverse commands
    counts <- tbl(birds, "breed_bird_survey_counts")
    weather <- tbl(birds, "breed_bird_survey_weather")
    routes <- tbl(birds, "breed_bird_survey_routes")
    species <- tbl(birds, "breed_bird_survey_species") %>%
      select(-species_id) #drop the column that BBS is calling species_id (not the same as our species ID which is the AOU code)
    
    #join all data into one table
    bbs <- left_join(weather, counts, by = c("year", "statenum", "route", "rpid", "year")) %>%
      left_join(routes, by = c("statenum", "route")) %>%
      left_join(species, by = "aou") %>%
      ###
      dplyr::filter(runtype == 1 & rpid == 101) %>%
      mutate(site_id = (statenum*1000) + route) %>%
      select(site_id, latitude, longitude, aou, year, speciestotal) %>%
      rename(species_id = aou, abundance = speciestotal, lat = latitude, long = longitude) %>%
      collect() 
    
    #create own get_species_data() function to use DBI/table method for connecting to database
    get_species_data = function() {
      data_path <- paste('./data/', 'bbs', '_species.csv', sep = "")
      if (file.exists(data_path)) {
        return(read.csv(data_path))
      }else{
        write.csv(species, file = data_path, row.names = FALSE, quote = FALSE)
        return(species)
      }
    }
    
    #clean up specie(i.e. combine subspecies, exclude poorly sampled species), see source script for details - probably doesn't work
    bbs_clean <- bbs %>% 
      filter_species() %>%
      group_by(site_id) %>%
      combine_subspecies() %>%
      #add taxonomy
      left_join(collect(species), by = c("species_id" = "aou")) %>%
      select (site_id, year, species_id, lat, long, abundance, genus, species, english_common_name) %>%
      rename (common_name = english_common_name)
    
    write.csv(bbs_clean, file = data_path, row.names = FALSE, quote = FALSE)
    return(bbs_clean)
    
  }
}

bbs <- get_bbs()

###################
####Trait Data#####
###################
get_trait_data <- function(){
  data_path <- paste('./data/elton_traits/', 'elton_traits', '_BirdFuncDat.csv', sep = "")
  if (file.exists(data_path)){
    return(read_csv(data_path))
  }else{
    # dir.create("data/elton_traits")
    # rdataretriever::install("elton-traits", 'csv', data_dir = "data/elton_traits")
  }
}

trait_data <- get_trait_data()
