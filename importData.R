library(tidyverse)
library(stringr)

###################
###### BBS ########
###################

source("bbs_forecasting_functions.R")

#Adapted from get_bbs_data() from sourced scripts
get_bbs <- function(){
  data_path <- paste('./data/', 'bbs', '_data.csv', sep="")
  if (file.exists(data_path)){
    return(read_csv(data_path))
    print("yes csv")
  }
  else{
    if (!db_engine(action='check', db = "~/Dropbox/Data/functional-diversity/bbsforecasting.sqlite",
                   table_to_check = 'breed_bird_survey_counts')){
      print("no database")
      install_dataset('breed-bird-survey')
    }
    
    birds <- DBI::dbConnect(RSQLite::SQLite(), "~/Dropbox/Data/functional-diversity/bbsforecasting_old.sqlite")
    
    #save database tables as table in R to use with tidyverse commands
    counts <- tbl(birds, "breed_bird_survey_counts")
    weather <- tbl(birds, "breed_bird_survey_weather")
    routes <- tbl(birds, "breed_bird_survey_routes")
    species <- tbl(birds, "breed_bird_survey_species") %>%
      select(-species_id) #drop the column that BBS is calling species_id (not the same as our species ID which is the AOU code)
    
    print('tables ran')
    
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

    #clean up specie(i.e. combine subspecies, exclude poorly sampled species), see source script for details - probably doesn't work
    bbs_clean <- bbs %>% 
      filter_species() %>%
      group_by(site_id) %>%
      combine_subspecies() %>%
      #add taxonomy
      left_join(collect(species), by = c("species_id" = "aou")) %>%
      select (site_id, year, species_id, lat, long, abundance, genus, species, english_common_name) %>%
      rename (common_name = english_common_name) %>%
      unite(scientific, genus, species, sep = " ")
    
    write.csv(bbs_clean, file = data_path, row.names = FALSE, quote = FALSE)
    return(bbs_clean)
    
  }
}

bbs <- get_bbs()

###################
####Trait Data#####
###################
get_elton_trait <- function(){
  data_path <- paste('./data/elton_traits/', 'elton_traits', '_BirdFuncDat.csv', sep = "")
  if (file.exists(data_path)){
    return(read_csv(data_path))
  }else{
    dir.create("data/elton_traits")
    rdataretriever::install("elton-traits", 'csv', data_dir = "data/elton_traits")
  }
}

elton_trait <- get_elton_trait()

####################
### BBS & Trait ####
### Master CSV  ####
####################

get_bbs_w_traits <- function(){
  
  data_path <- paste('./data/', 'bbsTraits_master.csv', sep = "")
  if (file.exists(data_path)){
    return(read_csv(data_path))
  }else{
    #join bbs and traits on scientific name to find taxonomic mismatches, 
    #and get a dataframe of taxonomic equivalancies based on common names
    sci_equivalent <- bbs %>% 
      select(scientific, common_name) %>%
      unique() %>%
      left_join(select(elton_trait, scientific, english), by = "scientific") %>% #join on scientific name
      subset(is.na(english)) %>% #get rows where the scientic names didn't match
      select(scientific, common_name) %>% #select bbs sci name and common name
      left_join(select(elton_trait, scientific, english), 
                by = c("common_name" = "english")) %>% #join bbs and trait data on common name to see taxanomic equivalents
      rename(bbs_sci = scientific.x, trait_sci = scientific.y) 
    
    
    ## Yellow-rumped warbler = yellow-rumped warbler w/o commentary
    ## Black-throated Gray Warbler = Black-throated Grey Warbler
    ## Gray Hawk = grey hawk (technically different subspecies?)
    
    ###No trait data for: Pacific Wren, Eastern Yellow Wagtail, Sagebrush sparrow, Woodhouse's scrub jay, Bell's sparrow
    
    
    get_compatible_sci_names <- function(data, sci_equiv){
      #Create a new column called compat_sci that replaces BBS sci names with their trait data equivalent#
      
      #case when formula - one for each row of equivalencies 
      replace_form <- lapply(1:dim(sci_equiv)[1],function(var) 
        formula(paste0('data$scientific == as.character(sci_equiv$bbs_sci[',
                       var, ']) ~ as.character(sci_equiv$trait_sci[', var,'])')))
      
      #add special cases where neither common or scientific names match, but the species are the same
      ##still don't work, not sure why
      replace_form <- append(replace_form, 
                             c(formula('data$scientific == \'Setophaga coronata\' ~ \'Dendroica coronata\''),
                               formula('data$scientific == \'Setophaga nigrescens\' ~ \'Dendroica nigrescens\''),
                               formula('data$scientific == as.character(\'Buteo plagiatus\') ~ \'Buteo nitidus\'')
                             )
      )
      
      #add case for when there is no equivalence and we keep the original name
      replace_form <- append(replace_form, formula(paste0("TRUE ~ as.character(data$scientific)")))
      
      #add column
      data %>%
        mutate(compat_sci = case_when(!!!replace_form))
      
    }
    
    #get equivalence column for BBS data
    bbs_compat <- get_compatible_sci_names(bbs, sci_equivalent)
    
    
    #Merge bbs
    bbs_trait <- elton_trait %>%
      select(-specid, -passnonpass, -iocorder, -blfamilylatin, -blfamilyenglish, -blfamsequid, -taxo, -bodymass_speclevel, -english,
             -ends_with("source"), -ends_with("comment"), -ends_with("enteredby")) %>%
      right_join(select(bbs_compat, -scientific), by = c("scientific" = "compat_sci")) %>%
      filter(!is.na(scientific))
    
    write.csv(bbs_trait, file = data_path, row.names = FALSE, quote = FALSE)
    return(bbs_trait)
  }
}

bbs_trait <- get_bbs_w_traits()