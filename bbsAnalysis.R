source("installPackages.R")
library(tidyverse)
library(dplyr)
library(FD)
library(spData)
library(sf)
library(tmap)

trait <- read_csv("~/Dropbox/functional-diversity/elton_traits/elton_traits_BirdFuncDat.csv")
bbs <- read_csv("~/Dropbox/functional-diversity/bbs_data_compatible.csv") #422 species

min_year = 2006 #define the minimum year of sampling to include
p <-  102003 # USA Contiguous Albers Equal Area Conic - planar projection for st_ functions

#Species Matrix  
get_species_matrix <- function(){
  species <- bbs %>% #401 species, 4176 sites, some not seen in the time period
    filter(year > min_year) %>%
    dplyr::select(scientific, site_id, abundance) %>%
    group_by(scientific, site_id) %>%
    summarize(m = mean(abundance)) %>% #think about something other than mean?
    spread(scientific, m) %>%
    column_to_rownames(var = "site_id")
}

#Trait Matrix
get_trait_matrix <- function(species_list = colnames(species)){ 
  traits <- trait %>%
    filter(scientific %in% species_list) %>%
    dplyr::select(-specid, -passnonpass, -iocorder, -blfamilylatin, -blfamilyenglish, -blfamsequid, -taxo, -bodymass_speclevel, -english, -diet_certainty,
           -ends_with("source"), -ends_with("comment"), -ends_with("enteredby")) %>%
    arrange(scientific) %>%
    column_to_rownames(var = "scientific")
}

#Get Functional Diversity Metrics  
get_site_FD <- function(){  
  data_path <- paste('~/Dropbox/functional-diversity/', 'FD_stats.RData', sep="")
  if (file.exists(data_path)){
    print("FD present")
    FD_file <- load(data_path)
    return(as.data.frame(get(FD_file)))
  }else{
    print("No FD")
    species <- get_species_matrix()
    traits <- get_trait_matrix()

    FD <- as.data.frame(dbFD(traits, species, w.abun = TRUE))
    save(FD, file = data_path)
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
  spatial_routes <- bbs %>%
    filter(year > min_year) %>%
    dplyr::select(site_id, long, lat) %>%
    unique() %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = p)
}


get_sites_w_region <- function(sites = FALSE, method = c("intersect", "dist")){
  
  if(!hasArg(sites)) sites <- get_sites_sf()
  
  bcr <- get_ecoreg_shp()
  bcr_names <- unique(bcr$BCRNAME)

  if (method == "intersect"){
    region_sites <- matrix(ncol = length(bcr_names), nrow = dim(sites)[1])
    for (i in 1:length(bcr_names)){
      reg_name <- bcr %>%
        filter(BCRNAME == bcr_names[i])
      
      int_mat <- st_intersects(reg_name, sites, sparse = FALSE)
      
      #case when there is more than one polygon, and some intersections are found
      if (dim(int_mat)[1] > 0 & sum(int_mat > 0)){
        int_mat <- as.logical(colSums(int_mat))
      }
      
      #case when no intersections are found
      if (sum(int_mat) == 0){
        int_mat <- rep(FALSE, dim(sites)[1])
      }
      region_sites[,i] <- (int_mat)
    }
    
    region_sites <- as.data.frame(region_sites)
    colnames(region_sites) <- bcr_names
    region_sites <- cbind(site_id = sites$site_id, region_sites)
    
    #Check if any sites were classified in two 
    bad_sites <- apply(dplyr::select(region_sites, -site_id), 1, function(x) length(which(x)))
    
    if(max(bad_sites) > 1){
      warning("One or more sites has been classified to multiple regions")
    }
    
    region_sites[region_sites == FALSE ] <- NA
    site_labels <- region_sites %>% 
      gather(region, value, -site_id) %>% 
      na.omit() %>% 
      dplyr::select(-value) %>%
      left_join(., dplyr::select(region_sites, site_id), by = "site_id") %>%
      left_join(., sites, by = "site_id") %>%
      arrange(site_id) %>%
      st_sf()
  }
  
  if (method == "dist"){
    dist <- st_distance(sites, bcr)
    poly_index <- apply(dist, 1, which.min) #get index of the minimum distance polygon for each site
    region <- bcr$BCRNAME[poly_index] #get the name of the nearest region by the index
    site_labels <- cbind(sites, region)
  }
  return(site_labels)
}

get_complete_site_data <- function(){
  continent <- get_sites_w_region(method = "intersect") #get sites that intersect
  
  bbs_sites <- get_sites_sf()
  dropped_site_ids <- dplyr::setdiff(bbs_sites$site_id, continent$site_id) #find dropped
  dropped_sites <- bbs_sites %>% filter(site_id %in% dropped_site_ids)
  coast <- get_sites_w_region(sites = dropped_sites, method = "dist") #get sites that were dropped (on the coast)
  
  all <- rbind(continent, coast)
  
  FD <- get_site_FD() %>%
    rownames_to_column() %>%
    mutate(site_id = as.integer(rowname)) %>%
    dplyr::select(-rowname) %>%
    left_join(., all, by = "site_id")
}

bbs_site_FD <- get_complete_site_data()


#Simulate null model for one region
n_rockies_FD <- filter(bbs_site_FD, region == "NORTHERN ROCKIES")
n_rockies <- bbs %>% 
  filter(year > min_year & site_id %in% unique(n_rockies_FD$site_id))

species_pool <- unique(n_rockies$scientific)

get_sample_fd <- function(x, ...){
  samp_trait_mat <- get_trait_matrix(sample(species_pool, x))
  samp_species <- rownames(samp_trait_mat)
  sample_FD <- dbFD(x = samp_trait_mat, ...)
  #return(c(richness = x, head(sample_FD, -1))) #remove last element, which is the CWM for each trait - maybe add back in later?
  return(list("species" = samp_species, "FD" = head(sample_FD, -1)))
}

#test_sim <- plyr::ldply(100:length(species_pool), get_sample_fd(calc.FRic = FALSE)$species) #would work if dbFD didn't error out

FDdf <- data.frame()
rich_vals <- c()
#for loop option, still need to add database for iterative storage
for(i in 106:length(species_pool)){
  possibleError <- tryCatch(
    samp_fd <- get_sample_fd(i),
    error=function(e)e
  )
  if(inherits(possibleError, "error")) next
  
  #rich_vals <- c(rich_vals, i)
  FDdf <- rbind(FDdf, samp_fd$FD)
}

#recalculate broken case from for loop ^^^
test_trait_mat <- get_trait_matrix(samp_fd$species)
test_samp_fd <- dbFD(x = test_trait_mat)


#preliminary plot of null curve
FDdf %>% ggplot(aes(x = nbsp, y = FDiv)) + geom_smooth() +
  theme_classic()
