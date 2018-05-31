library(tidyverse)

bird_trait <- read_csv("~/Documents/Projects/functional_diversity/data/elton_traits/elton_traits_BirdFuncDat.csv")
bbs <- read_csv("data/bbs_data.csv")

#join bbs and traits on scientific name to find taxonomic mismatches, and get a dataframe of taxonomic equivalancies based on common names
sci_equivalent <- bbs %>% 
  select(scientific, common_name) %>%
  unique() %>%
  left_join(select(bird_trait, scientific, english), by = "scientific") %>% #join on scientific name
  subset(is.na(english)) %>% #get rows where the scientic names didn't match
  select(scientific, common_name) %>% #select bbs sci name and common name
  left_join(select(bird_trait, scientific, english), 
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

#test case
test_bbs <- data.frame(scientific = c("a", "b", "c", "d"), abundance = c(1,2,3,4))
test_sci <- data.frame(bbs_sci = c("a", "b"), trait_sci = c("a1", "b1"))
get_compatible_sci_names(data = test_bbs, sci_equiv = test_sci)

#get equivalence column for BBS data
bbs_compat <- get_compatible_sci_names(bbs, sci_equivalent)


#Merge bbs
bbs_trait <- bird_trait %>%
  select(-specid, -passnonpass, -iocorder, -blfamilylatin, -blfamilyenglish, -blfamsequid, -taxo, -bodymass_speclevel, -english,
         -ends_with("source"), -ends_with("comment"), -ends_with("enteredby")) %>%
  right_join(select(bbs_compat, -scientific), by = c("scientific" = "compat_sci")) %>%
  filter(!is.na(scientific))
  

