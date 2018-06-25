
#' Get a species by site matrix from a dataframe of observations
#' @export
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

#' Get the Functional Diversity metrics for all bbs sites, aggregating over a defined time period.
#' Will use previously calculated FD if it is in the appropriate path
#' @export
get_site_FD <- function(){  
  data_path <- paste('./data/', 'FD_stats.RData', sep="")
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

#' Read in shapefile for bird conservation regions
#' @export
get_ecoreg_shp <- function(){
  bcr <- st_read("data/bcr_shp/BCR_Terrestrial_master.shp") %>%
    st_transform(crs = p) %>%
    filter(REGION %in% c("CANADA", "USA"))
}

#' Get bbs sites for our time period as a spatial dataframe
#' @export
get_sites_sf <- function(){
  spatial_routes <- bbs %>%
    filter(year > min_year) %>%
    dplyr::select(site_id, long, lat) %>%
    unique() %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = p)
}


#' Label sites with their conservation region
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

#' Get site-level FD metrics with their region label
#' @export
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
