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
                                                                                        

## Yellow-warbler = yellow-warbler w/o commentary
## Black-throated Gray Warbler = Black-throated Grey Warbler
## Gray Hawk = grey hawk (technically different subspecies?)

###No trait data for: Pacific Wren, Eastern Yellow Wagtail, Sagebrush sparrow, Woodhouse's scrub jay, Bell's sparrow


### Replace old sci names with new sci names - test case
test_bbs <- data.frame(scientific = c("a", "b", "c", "d"), abundance = c(1,2,3,4))
test_sci <- data.frame(bbs_sci = c("a", "b"), trait_sci = c("a1", "b1"))

#case when formula
replace_form <- lapply(1:dim(test_sci)[1],function(var) 
  formula(paste0('test_bbs$scientific == as.character(test_sci$bbs_sci[',var, ']) ~ as.character(test_sci$trait_sci[', var,'])'),
          env=globalenv()))
replace_form <- append(replace_form, formula(paste0("TRUE ~ as.character(test_bbs$scientific)")))

test_bbs %>%
  mutate(compat_sci = case_when(!!!replace_form))




bbs_trait <- bird_trait %>%
  select(-specid, -passnonpass, -iocorder, -blfamilylatin, -blfamilyenglish, -blfamsequid, -taxo, -bodymass_speclevel,
         -ends_with("source"), -ends_with("comment")) %>%
  right_join(bbs, by = "scientific")


species_joined_com <- bbs %>% 
  select(scientific, common_name) %>%
  unique() %>%
  full_join(select(bird_trait, scientific, english), by = c("common_name" = "english"))
