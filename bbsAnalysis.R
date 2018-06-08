library(tidyverse)
library(dplyr)
library(FD)
library(spData)
library(sf)
library(tmap)

trait <- read_csv("data/elton_traits/elton_traits_BirdFuncDat.csv")
bbs <- read_csv("data/bbs_data.csv") #422 species
bbs_trait <- read_csv("data/bbsTraits_master.csv") #414 species, 8 species not merged with trait data, see importData.R for which species and why

min_year = 2006 #define the minimum year of sampling to include
p <-  4326 # +proj=longlat +datum=WGS84: Projection for spatial analyses

#Species Matrix  
get_species_matrix <- function(){
  species <- bbs_trait %>% #401 species, 4176 sites, some not seen in the time period
    filter(year > min_year) %>%
    dplyr::select(scientific, site_id, abundance) %>%
    group_by(scientific, site_id) %>%
    summarize(m = mean(abundance)) %>% #think about something other than mean?
    spread(scientific, m) %>%
    column_to_rownames(var = "site_id")
}

#Trait Matrix
get_trait_matrix <- function(species_list = colnames(species)){ #species_list should be scientific names derived from trait or bbs_trait, bbs isn't neccesarily compatible
  traits <- trait %>%
    filter(scientific %in% species_list) %>%
    dplyr::select(-specid, -passnonpass, -iocorder, -blfamilylatin, -blfamilyenglish, -blfamsequid, -taxo, -bodymass_speclevel, -english, -diet_certainty,
           -ends_with("source"), -ends_with("comment"), -ends_with("enteredby")) %>%
    arrange(scientific) %>%
    column_to_rownames(var = "scientific")
}

#Get Functional Diversity Metrics  
get_site_FD <- function(){  
  data_path <- paste('./data/', 'FD_stats.RData', sep="")
  if (file.exists(data_path)){
    print("FD present")
    FD_file <- load(data_path)
    return(as.data.frame(get(FD_file)))
  }else{
    species <- get_species_matrix()
    traits <- get_trait_matrix()

    FD <- dbFD(traits, species, w.abun = TRUE)
    save(FD, file = "FD_stats.RData")
  }
}

#Get Ecoregions for each site

get_ecoreg_shp <- function(){
  bcr <- st_read("data/bcr_shp/BCR_Terrestrial_master.shp") %>%
    st_transform(crs = p) %>%
    filter(REGION %in% c("CANADA", "USA"))
}

#get bbs sites for our time period as a spatial dataframe
get_sites_sf <- function(){
  spatial_routes <- bbs_trait %>%
    filter(year > min_year) %>%
    dplyr::select(site_id, long, lat) %>%
    unique() %>%
    st_as_sf(coords = c("long", "lat"), crs = p)
}


get_sites_w_region_FD <- function(){
  bcr <- get_ecoreg_shp()
  bbs_sites <- get_sites_sf()
  bcr_names <- unique(bcr$BCRNAME)
  
  region_sites <- matrix(ncol = length(bcr_names), nrow = dim(bbs_sites)[1])
  for (i in 1:length(bcr_names)){
    reg_name <- bcr %>%
      filter(BCRNAME == bcr_names[i])
    int_mat <- st_intersects(reg_name, bbs_sites, sparse = FALSE)
    
    # case_when(
    #   dim(int_mat)[1] > 0 & sum(int_mat) > 0 ~ as.logical(colSums(int_mat)),
    #   sum(int_mat) == 0 ~ rep(FALSE, dim(bbs_sites)[1]),
    #   TRUE ~ int_mat
    # )
    
    #case when there is more than one polygon, and some intersections are found
    if (dim(int_mat)[1] > 0 & sum(int_mat > 0)){
      int_mat <- as.logical(colSums(int_mat))
    }
    
    #case when no intersections are found
    if (sum(int_mat) == 0){
      int_mat <- rep(FALSE, dim(bbs_sites)[1])
    }
    region_sites[,i] <- (int_mat)
  }
  
  region_sites <- as.data.frame(region_sites)
  colnames(region_sites) <- bcr_names
  region_sites <- cbind(site_id = bbs_sites$site_id, region_sites)
  
  #Check if any sites were classified in two 
  bad_sites <- apply(dplyr::select(region_sites, -site_id), 1, function(x) length(which(x)))
  
  if(max(bad_sites) > 1){
    warning("One or more sites has been classified to multiple regions")
  }
  
  FD <- get_site_FD() %>%
    rownames_to_column() %>%
    mutate(site_id = as.integer(rowname)) %>%
    dplyr::select(-rowname)
  
  region_sites[region_sites == FALSE ] <- NA
  site_labels <- region_sites %>% 
    gather(region, value, -site_id) %>% 
    na.omit() %>% 
    dplyr::select(-value) %>%
    left_join(., dplyr::select(region_sites, site_id)) %>%
    left_join(., bbs_sites) %>%
    left_join(., FD, by = "site_id") %>%
    arrange(site_id)
  ## ^^^ above joins result in multiple regions for each site - not sure what's going on 
}

bbs_site_FD <- get_sites_w_region_FD()

sites_in_region <- bbs_site_FD %>% 
  group_by(region) %>%
  summarise(n = n()) %>%
  arrange(n)


#Map dropped sites
bbs_sites <- get_sites_sf()

dropped_site_ids <- dplyr::setdiff(bbs_sites$site_id, bbs_site_FD$site_id)
dropped_sites <- bbs_sites %>% filter(site_id %in% dropped_site_ids)

bcr <- get_ecoreg_shp()
map_dropped = tm_shape(bcr) + tm_borders()
map_dropped + tm_shape(dropped_sites) + tm_dots(col = "red")


#Simulate null model for one region
n_rockies_FD <- filter(bbs_site_FD, region == "NORTHERN ROCKIES")
n_rockies <- bbs_trait %>% 
  filter(year > min_year & site_id %in% unique(n_rockies_FD$site_id))

species_pool <- unique(n_rockies$scientific)

#for loop option, still need to add database for iterative storage  
# for(i in 1:length(species_pool)){
#   samp_trait_mat <- get_trait_matrix(sample(species_pool, i))
#   sample_FD <- dbFD(samp_trait_mat)
# }

#try to do it the R way with an apply
get_sample_fd <- function(x){
  samp_trait_mat <- get_trait_matrix(sample(species_pool, x))
  sample_FD <- dbFD(samp_trait_mat)
  return(c(richness = x, head(sample_FD, -1))) #remove last element, which is the CWM for each trait - maybe add back in later?
}

test_sim <- plyr::ldply(1:length(species_pool), get_sample_fd()) #should work 
