library(tidyverse)

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

source_https("https://raw.githubusercontent.com/weecology/bbs-forecasting/master/R/forecast-bbs-core.R")

##Not currently working, installs the database but doesn't produce the csv. May work on pre 3.5.0 R verison
bbs <- get_bbs_data()



# ##testing
# data_path <- paste('./data/', 'bbs', '_data.csv', sep="")
# if (file.exists(data_path)){
#   return(read_csv(data_path))
# }
# 
# if (!db_engine(action='check', table_to_check = 'breed_bird_survey_counts')){
#   print('yes')
# }
# 
# bbs_query ="SELECT
#                   (counts.statenum*1000) + counts.Route AS site_id,
# Latitude AS lat,
# Longitude AS long,
# aou AS species_id,
# counts.Year AS year,
# speciestotal AS abundance
# FROM
# breed_bird_survey_counts AS counts
# JOIN breed_bird_survey_weather
# ON counts.statenum=breed_bird_survey_weather.statenum
# AND counts.route=breed_bird_survey_weather.route
# AND counts.rpid=breed_bird_survey_weather.rpid
# AND counts.year=breed_bird_survey_weather.year
# JOIN breed_bird_survey_routes
# ON counts.statenum=breed_bird_survey_routes.statenum
# AND counts.route=breed_bird_survey_routes.route
# WHERE breed_bird_survey_weather.runtype=1 AND breed_bird_survey_weather.rpid=101"
# 
# bbs_data=db_engine(action='read', sql_query = bbs_query) %>%
#   filter_species() %>%
#   group_by(site_id) %>%
#   combine_subspecies()
# save_provenance(bbs_data)
# write.csv(bbs_data, file = data_path, row.names = FALSE, quote = FALSE)
# return(bbs_data)


###################
####Trait Data#####
###################
dir.create("data/elton_traits")
rdataretriever::install("elton-traits", 'csv', data_dir = "data/elton_traits")


