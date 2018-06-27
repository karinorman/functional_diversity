### Copy of bbs_forecasting_functions.R from Weecology's bbs-forecasting (6/1/18) found here:
### (https://github.com/weecology/bbs-forecasting/blob/master/R/forecast-bbs-core.R)

#' Install a particular dataset via rdataretreiver
#'
#' @param dataset name
install_dataset <- function(dataset){
  # Install a dataset using the rdataretriever
  
  rdataretriever::install(dataset, 'sqlite', db_file='./data/bbsforecasting.sqlite')
}

#' Single wrapper for all database actions
#'
#' We require only a few simple sql methods. They are 1. Writing an entire dataframe
#' directly to a database as it's own table, 2. Reading the same tables as dataframes,
#' possibly with some modification using SQL statements, 3. Checking to see if a
#' table exists. If a particular table does it exists it's assumed it has all data
#' required.
#'
#' read returns a dataframe
#' write returns nothing
#' check returns boolean
#'
#' @param action Action to perform in db call. Either read, write, or check
#' @param db name of database. A file if using sqlite
#' @param sql_query SQL statement if action is read
#' @param df Dataframe of data if action is write. Will copy the dataframe verbatim to it's own table with name new_table_name
#' @param new_table_name Table name for new data being written
#' @param table_to_check Table name to check if it exists for when action is check
#' @importFrom dplyr copy_to src_sqlite src_tbls collect tbl

db_engine=function(action, db='./data/bbsforecasting.sqlite', sql_query=NULL,
                   df=NULL, new_table_name=NULL, table_to_check=NULL){
  
  if(!dir.exists("data")){dir.create("data")}
  
  con <- src_sqlite(db, create=TRUE)
  
  if(action=='read'){
    to_return=collect(tbl(con, sql(sql_query)), n=Inf)
    
  } else if(action=='write') {
    copy_to(con, df, name=new_table_name, temporary = FALSE)
    to_return=NA
    
  } else if(action=='check') {
    #Only works with sqlite for now.
    to_return=tolower(table_to_check) %in% tolower(src_tbls(con))
    
  } else {
    stop(paste0('DB action: ',action,' not found'))
  }
  
  #Close the connection before returning results.
  rm(con)
  return(to_return)
}

#' Filter poorly sampled BBS species
#'
#' Removes waterbirds, shorebirds, owls, kingfishers, knightjars,
#' dippers. These species are poorly sampled due to their aquatic or
#' noctural nature. Also removes taxa that were either partially unidentified
#' (e.g. "sp.") or were considered hybrids (e.g. "A x B") or were listed as more
#' than one species (e.g. "A / B")
#'
#' @param df dataframe containing an species_id column
#'
#' @return dataframe, filtered version of initial dataframe
#' @importFrom dplyr "%>%" inner_join do rowwise select filter group_by ungroup full_join n_distinct semi_join left_join
filter_species <- function(df){
  species_table = get_species_data()
  
  is_unidentified = function(names) {
    #Before filtering, account for this one hybrid of 2 subspecies so it's kept
    names[names=='auratus auratus x auratus cafer']='auratus auratus'
    grepl('sp\\.| x |\\/', names)
  }
  
  valid_taxa = species_table %>%
    filter(!is_unidentified(species)) %>%
    filter(aou > 2880) %>%
    filter(aou < 3650 | aou > 3810) %>%
    filter(aou < 3900 | aou > 3910) %>%
    filter(aou < 4160 | aou > 4210) %>%
    filter(aou != 7010)
  
  filter(df, species_id %in% valid_taxa$aou)
}

#' Combine subspecies into their common species
#'
#' @importFrom dplyr "%>%" filter slice group_by summarise ungroup pull
#' @importFrom stringr word
combine_subspecies = function(df){
  
  species_table = get_species_data()
  
  # Subspecies have two spaces separated by non-spaces
  subspecies_names = species_table %>%
    filter(aou %in% unique(df$species_id)) %>%
    pull(spanish_common_name) %>%
    grep(" [^ ]+ ", ., value = TRUE)
  
  subspecies_ids = species_table %>%
    filter(spanish_common_name %in% subspecies_names) %>%
    pull(aou)
  
  # Drop all but the first two words to get the root species name,
  # then find the AOU code
  new_subspecies_ids = species_table %>%
    slice(match(word(subspecies_names, 1,2),
                species_table$spanish_common_name)) %>%
    pull(aou)
  
  # replace the full subspecies names with species-level names
  for (i in seq_along(subspecies_ids)) {
    df$species_id[df$species_id == subspecies_ids[i]] = new_subspecies_ids[i]
  }
  
  df %>%
    group_by(site_id, year, species_id, lat, long) %>%
    summarise(abundance = sum(abundance)) %>%
    ungroup()
}

##My version
get_species_data = function() {
  data_path <- paste('./data/', 'bbs', '_species.csv', sep = "")
  if (file.exists(data_path)) {
    return(read.csv(data_path))
  }else{
    write.csv(species, file = data_path, row.names = FALSE, quote = FALSE)
    return(species)
  }
}

##Ethan Version
# get_species_data = function() {
#   data_path <- paste('./data/', 'bbs', '_species.csv', sep = "")
#   if (file.exists(data_path)) {
#     return(read.csv(data_path))
#   }else{
#     species_table=db_engine(action = 'read', sql_query = 'SELECT * FROM breed_bird_survey_species')
#     write.csv(species_table, file = data_path, row.names = FALSE, quote = FALSE)
#     save_provenance(species_table)
#     return(species_table)
#   }
# }

#' Get the primary bbs data file which compiles the counts, route info, and weather
#' data. Install it via rdataretriever if needed.
#'
#' @export
#' @importFrom dplyr "%>%" group_by
#' @importFrom readr read_csv
get_bbs_data <- function(){
  
  data_path <- paste('./data/', 'bbs', '_data.csv', sep="")
  if (file.exists(data_path)){
    return(read_csv(data_path))
  }
  else{
    
    if (!db_engine(action='check', table_to_check = 'breed_bird_survey_counts')){
      install_dataset('breed-bird-survey')
    }
    
    #Primary BBS dataframe
    bbs_query ="SELECT
    (counts.statenum*1000) + counts.Route AS site_id,
    Latitude AS lat,
    Longitude AS long,
    aou AS species_id,
    counts.Year AS year,
    speciestotal AS abundance
    FROM
    breed_bird_survey_counts AS counts
    JOIN breed_bird_survey_weather
    ON counts.statenum=breed_bird_survey_weather.statenum
    AND counts.route=breed_bird_survey_weather.route
    AND counts.rpid=breed_bird_survey_weather.rpid
    AND counts.year=breed_bird_survey_weather.year
    JOIN breed_bird_survey_routes
    ON counts.statenum=breed_bird_survey_routes.statenum
    AND counts.route=breed_bird_survey_routes.route
    WHERE breed_bird_survey_weather.runtype=1 AND breed_bird_survey_weather.rpid=101"
    
    bbs_data=db_engine(action='read', sql_query = bbs_query) %>%
      filter_species() %>%
      group_by(site_id) %>%
      combine_subspecies()
    save_provenance(bbs_data)
    write.csv(bbs_data, file = data_path, row.names = FALSE, quote = FALSE)
    return(bbs_data)
  }
}

#' Get route locations in a SpatialPointsDataFrame
#'
#' @param projection string projection for route data
#'
#' @return a spatial data frame including site_id, long, and lat
#' @importFrom sp SpatialPointsDataFrame CRS
#' @importFrom dplyr collect copy_to src_sqlite src_tbls tbl %>%
get_route_data <- function(){
  p=CRS('+proj=longlat +datum=WGS84')
  bbs_data <- get_bbs_data()
  route_locations <- unique(dplyr::select(bbs_data, site_id, long, lat))
  spatial_routes <- route_locations %>%
    dplyr::select(long, lat) %>%
    SpatialPointsDataFrame(data=route_locations, proj4string=p)
}
