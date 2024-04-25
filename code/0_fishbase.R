
library(rfishbase)
library(dplyr)

#===============================================================================
#Main
#===============================================================================

# search habitat info from FishBase

# import and search
data <- read.csv("./input/dataset.csv", row.names = 1)
species_list <- rownames(data)

res <- species(species_list) %>%
  select(Species, SpecCode, Fresh, Brack, Saltwater)

# 88 species were not matched.
species_list_2 <- species_list[!species_list %in% res$Species]

# search edited species
edited <- read.csv("./input/fishbase/edited_species_list.csv")

res2 <- edited %>%
  pull(To) %>%
  species() %>%
  select(Species, SpecCode, Fresh, Brack, Saltwater) %>%
  rename(To = Species) %>%
  left_join(edited, ., by = "To") %>%
  filter(!is.na(SpecCode)) %>%
  bind_rows(res)

# 15 more speices
species_list_3 <- species_list[!species_list %in% res2$Species]

# for synonyms of original
res_syn <- synonyms(species_list_3) %>%
  select(synonym, Species, SpecCode) %>%
  filter(!is.na(SpecCode))

res3 <- species(res_syn$Species) %>%
  select(Species, SpecCode, Fresh, Brack, Saltwater) %>%
  left_join(res_syn[c("synonym", "SpecCode")], by = "SpecCode") %>%
  rename(To = Species, Species = synonym) %>%
  bind_rows(res2, .)

# 5 more speices
species_list_4 <- species_list[!species_list %in% res3$Species]

# for synonyms of edited
tmp <- edited %>%
  filter(Species %in% species_list_4)

res_syn <- tmp %>%
  pull(To) %>%
  synonyms() %>%
  select(synonym, Species, SpecCode) %>%
  filter(!is.na(SpecCode))

res4 <- species(res_syn$Species) %>%
  select(Species, SpecCode, Fresh, Brack, Saltwater) %>%
  left_join(res_syn[c("synonym", "SpecCode")], by = "SpecCode") %>%
  rename(To = synonym, syn = Species) %>%
  select(-syn) %>%
  left_join(., tmp, by = "To") %>%
  bind_rows(res3, .) %>%
  arrange(Species)
  
# ================================================================
# assign flags by habitat info.
# ================================================================

res_final <- res4 %>%
  mutate(habitat_class = case_when(
    Fresh == 1 & Brack == 0 & Saltwater == 0 ~ 1, # F
    Fresh == 1 & Brack == 1 & Saltwater == 0 ~ 2, # F-B
    Fresh == 1 & Brack == 1 & Saltwater == 1 ~ 3, # F-B-M
    Fresh == 0 & Brack == 1 & Saltwater == 0 ~ 4, # B
    Fresh == 0 & Brack == 1 & Saltwater == 1 ~ 5, # B-M
    Fresh == 0 & Brack == 0 & Saltwater == 1 ~ 6, # M
    TRUE ~ NA_integer_))  %>% # NA when not matched
  select(Species, habitat_class)

write.csv(res_final, "./output/habitat.csv")
