region <- "APPALACHIAN MOUNTAINS"

#' Check if the richness levels are the same across original BBS, the site-level FD output, and the simulated regional richness
#' @export
test_richness <- function(original, FD, region_name){
  
  fd_rich <- FD %>%
    filter(region == region_name) %>%
    select(nbsp, site_id) %>%
    unique()
  
  orig_rich <- original %>%
    filter(year > min_year, site_id %in% fd_rich$site_id) %>%
    select(site_id, species_id) %>%
    group_by(site_id) %>%
    summarise(nbsp = n_distinct(species_id))
  
  org_fd_test <- sum(!orig_rich$nbsp %in% fd_rich$nbsp)
  
  if (org_fd_test != 0){
    warning("original richness count does not match output from site-level FD calculations")
  }
  
  sim_rich <- read_tsv(here("data", "stat", paste0(get_safe_name(region_name), ".tsv.bz2"))) %>%
    select(richness) %>%
    unique()
  
  sim_fd_test <- sum(!sim_rich$richness %in% fd_rich$nbsp)
  
  if(sim_fd_test != 0){
    warning("simulated richness does not match output from site-level FD calculations")
  }
}

#test_richness(original = bbs, FD = bbs_site_FD, region_name = region)
