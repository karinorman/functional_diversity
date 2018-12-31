########################
### Dayment Climate#####
######   Data  #########
########################

library(daymetr)
library(here)
library(tidyverse)
library(lubridate)

data_path <- here::here("data", "modis", "daymet_climate.rdata")
if (file.exists(data_path)){
  load(data_path)
  print("yes climate")
}else{
  climate <- download_daymet_batch(file_location = here("data", 'bbs_sites.csv'),
                      start = 2007,
                      end = 2016,
                      internal = TRUE)
  save(climate, data_path)
}

error_list <- which(map(climate, class) == 'try-error')

#climate_df <- as.data.frame(climate[-error_list])

extract_daymet <- function(daymet_object) {
  site_id <- daymet_object$site
  
  df <- as.data.frame(daymet_object$data) %>%
    mutate(site = site_id)
  
  
}

climate_df <- map_df(climate[-error_list], extract_daymet) %>%
  mutate(site = as.factor(site))

samp_period_avg <- climate_df %>%
  filter(yday > 121, yday < 181) %>%
  select(-yday) %>%
  group_by(year, site) %>%
  summarise_at(.vars = vars(avg_cols), .funs = mean) %>%
  rename_all(paste0, "samp") %>%
  rename(year = yearsamp, site = sitesamp)
  
year_avg <- climate_df %>%
  select(-yday) %>%
  group_by(year, site) %>%
  summarise_all(mean) %>%
  rename_all(paste0, "yearly") %>%
  rename(year = yearyearly, site = siteyearly)

climate_means <- left_join(samp_period_avg, year_avg, by = c("year", "site"))

#######################
###### MODIS ##########
#######################

library(MODISTools)

sites <- read_csv(here("data", "bbs_sites.csv")) %>%
  mutate(site_id = as.factor(site_id)) %>%
  rename(lon = long) %>%
  as.data.frame()

extract_MODIS <- function(site_name, modis_object){
  as.data.frame(modis_object$data) %>% mutate(site_id = site_name)
}

#### YEARLY DATA  #####

clean_yearly_df <- function(data, data_name){
  
  #get all dates in the same format
  data$date <- gsub("[^[:alnum:][:blank:]?&/\\-]", "", data$calendar_date)
  data$month_first <- mdy(data$date)
  
  data_clean <- data %>% subset(is.na(month_first)) %>%
    select(-month_first) %>%
    mutate(month_first = ydm(date)) %>% 
    bind_rows(subset(data, !is.na(month_first))) %>%
    subset(!is.na(band)) %>%
    mutate(year = year(month_first), site_id = as.factor(site_id)) %>%
    select(site_id, year, !!quo_name(data_name) := data)
  
  return(data_clean)
    
}

#FAO land cover type
lc <- mt_batch_subset(df = sites,
                      product = "MCD12Q1",
                      band = "LC_Type4",
                      internal = TRUE,
                      start = "2007-01-01",
                      end = "2016-12-31")

lc_df <- map2(sites$site_id, lc, extract_MODIS) %>%
  bind_rows()

write_csv(lc_df, "data/modis/land_cover.csv")
  
#NPP
npp <- mt_batch_subset(df = sites,
                       product = "MYD17A3H",
                       band = "Npp_500m",
                       internal = TRUE,
                       start = "2007-01-01",
                       end = "2016-12-31")

npp_df <- map2(sites$site_id, npp, extract_MODIS) %>%
  bind_rows()

write_csv(npp_df, "data/modis/npp.csv")


### MONTHLY DATA

get_month_subset <- function(year, site_df = sites, product_name, band_name, folder){
   modis_data <- mt_batch_subset(df = site_df, 
                                product = product_name,
                                band = band_name,
                                internal = TRUE,
                                start = paste0(year, "-05-01"),
                                end = paste0(year, "-06-30"))
  
  modis_df <- map2(site_df$site_id, modis_data, extract_MODIS) %>% bind_rows()
  write_csv(modis_df, paste0("data/modis/", folder, "/", folder, "_", year, ".csv"))
  
  return(modis_df)
}


#Leaf area index

dir.create("data/modis/lai")
system.time(lai <- map_df(2007:2016, get_month_subset,
                                    site_df = site_sub,
                                    product_name = "MCD15A2H",
                                    band_name = "Lai_500m", 
                                    folder = "lai"))
            

#NDVI
dir.create("data/modis/ndvi")
system.time(ndvi <- map_df(2012:2016, get_month_subset,
                          product_name = "MOD13Q1",
                          band_name = "250m_16_days_NDVI", 
                          folder = "ndvi"))

ndvi_df <- map_df(paste0("data/modis/ndvi/", dir(path = "data/modis/ndvi", pattern = "*.csv")), read_csv) %>%
  select(site_id, calendar_date, data) %>%
  mutate(year = year(calendar_date)) %>%
  group_by(site_id, year) %>%
  summarise(ndvi = mean(data)) %>%
  ungroup() %>%
  mutate(site_id = as.factor(site_id))


site_data <- climate_means %>% 
  mutate(site_id = site) %>%
  left_join(clean_yearly_df(lc_df, "land_cover"), by = c("site_id", "year")) %>%
  left_join(clean_yearly_df(npp_df, "npp"), by = c("site_id", "year")) %>%
  left_join(ndvi_df, by = c("site_id", "year"))

write_tsv(site_data, "data/site_variables.tsv.bz2")
