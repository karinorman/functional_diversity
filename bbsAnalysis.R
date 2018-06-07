library(tidyverse)
library(FD)
library(spData)
library(sf)

trait <- read_csv("data/elton_traits/elton_traits_BirdFuncDat.csv")
bbs <- read_csv("data/bbs_data.csv") #422 species
bbs_trait <- read_csv("data/bbsTraits_master.csv") #414 species
  
#Species Matrix  
species <- bbs_trait %>% #401 species - where are they going??
  filter(year > 2006) %>%
  select(scientific, site_id, abundance) %>%
  group_by(scientific, site_id) %>%
  summarize(m = mean(abundance)) %>%
  spread(scientific, m) %>%
  column_to_rownames(var = "site_id")

#Trait Matrix  
traits <- trait %>%
  filter(scientific %in% colnames(species)) %>%
  select(-specid, -passnonpass, -iocorder, -blfamilylatin, -blfamilyenglish, -blfamsequid, -taxo, -bodymass_speclevel, -english, -diet_certainty,
         -ends_with("source"), -ends_with("comment"), -ends_with("enteredby")) %>%
  arrange(scientific) %>%
  column_to_rownames(var = "scientific")

#Get Functional Diversity Metrics  
testFD <- dbFD(traits, species, w.abun = TRUE)
save(testFD, file = "FD_stats.RData")


#Get Ecoregion for each site
p <-  4326 # +proj=longlat +datum=WGS84
bcr <- st_read("data/bcr_shp/BCR_Terrestrial_master.shp") %>%
  st_transform(crs = p) %>%
  filter(REGION %in% c("CANADA", "USA"))

get_route_data <- function(){
  route_locations <- unique(dplyr::select(bbs, site_id, long, lat))
  spatial_routes <- route_locations %>%
    #dplyr::select(long, lat) %>%
    st_as_sf(coords = c("long", "lat"), crs = p)
}

bbs_routes <- get_route_data()

bcr_names <- unique(bcr$BCRNAME)
region_sites <- matrix(ncol = length(bcr_names), nrow = dim(bbs_routes)[1])

for (i in 1:length(bcr_names)){
  reg_name <- bcr %>%
    filter(BCRNAME == bcr_names[i])
  int_mat <- st_intersects(reg_name, bbs_routes, sparse = FALSE)
  
  if (dim(int_mat)[1] > 0 & sum(int_mat > 0)){
    print(1)
    int_mat <- as.logical(colSums(int_mat))
  }else{
    print(c(2, i))
    int_mat <- rep(FALSE, dim(bbs_routes)[1])
  }
  print(dim(int_mat))
  region_sites[,i] <- (int_mat)
  region_sites <- as.data.frame(region_sites)
  colnames(region_sites) <- bcr_names
}

site_region_num <- apply(region_sites, 1, function(x) min(which(x == TRUE)))

site_region_num_min <- apply(region_sites, 1, function(x) max(which(x == TRUE)))


