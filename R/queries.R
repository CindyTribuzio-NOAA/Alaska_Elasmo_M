# age data from fishery and survey as available from AFSC

# Set up ----
libs <- c("tidyverse", "RODBC", "lubridate", "janitor")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

dbname <- "akfin"
db <- read_csv('database.csv')
database_akfin=db %>% filter(database == dbname) %>% select(database) #need to add filter for AKFIN user/pass only
username_akfin=db %>% filter(database == dbname) %>% select(username)
password_akfin=db %>% filter(database == dbname) %>% select(password)

channel_akfin <- odbcConnect(dbname, uid = username_akfin, pwd = password_akfin, believeNRows=FALSE)


# Create a year subdirectory to store annual data
dat_path <- paste0("data")
dir.create(dat_path)
raw_path <- paste0(dat_path, "/raw") # raw data
dir.create(raw_path) 

# race species_codes ----
# not really necessary, but a handy look-up table for quick reference
spp <- sqlQuery(channel_akfin, query = ("
            select * 
            from gap_products.akfin_taxonomic_classification
            where species_code between 230 and 495
            and ID_rank = 'species' 
            and database_id is not null")) %>% 
   clean_names() 

# in case the query needs to be run in SQL developer (e.g., if RODBC is being a pill)
#spp <- read_csv(paste0(getwd(), '/', raw_path, '/race_spp.csv.')) %>% 
#  clean_names()

write_csv(spp, paste0(dat_path, '/race_spp.csv'))

# survey age data ----
srv_ages <- sqlQuery(channel_akfin, query = ("
                select b.cruisejoin, a.hauljoin, c.year, c.survey_name, a.specimen_id, a.species_code, a.length_mm, a.sex, a.age
                from gap_products.akfin_specimen a
                left join gap_products.akfin_haul b
                on a.hauljoin = b.hauljoin
                left join gap_products.akfin_cruise c
                on b.cruisejoin = c.cruisejoin
                where species_code between 230 and 495 and age is not null")) %>% 
  clean_names() 

# in case the query needs to be run in SQL developer (e.g., if RODBC is being a pill)
#srv_ages <- read_csv(paste0(getwd(), '/', raw_path, '/survey_ages.csv.')) %>% 
#  clean_names()

# NOTE: no AI samples, so not added in this mutate, may need to fix for future use
srv_ages <- srv_ages %>% 
  mutate(area = ifelse(survey_name == 'Gulf of Alaska Bottom Trawl Survey', 'GOA', 'BS'))

write_csv(srv_ages, paste0(dat_path, "/survey_ages_akfin.csv"))

# fishery data ----

# akr species code translation for observer data
akrspp <- read_csv("data/species_translation_akr_obs.csv")

akrspp <- akrspp %>%
  filter(grepl('Skate|Shark', akr_name)) %>% 
  distinct(akr_name, akr_code)

akrspp <- akrspp %>% 
  mutate(common_name = tolower(akr_name)) %>% 
  left_join(spp)

myakrspp <- akrspp %>% pull(akr_code)
myakrspp_string <- toString(sprintf("'%s'", myakrspp)) # allows you to pass vector into sql query
  
fshbio <- sqlQuery(channel_akfin, query = ("
                select   *
                from     norpac.debriefed_age_flat_mv
                where    akr_species_codes between 689 and 705 and
                        age IS NOT NULL")) %>% 
  clean_names() 

# in case the query needs to be run in SQL developer (e.g., if RODBC is being a pill)
#fshbio <- read_csv(paste0(getwd(), '/', raw_path, '/fishery_ages.csv.')) %>% 
#  clean_names()

#there is only spiny dogfish on the west coast, but this step is just in case other data show up for future use
fshbio <- fshbio %>% 
  mutate(area = ifelse(fmp_subarea %in% c('AI', 'BS'), fmp_subarea,
                       ifelse(fmp_subarea == 'WOC', fmp_subarea, 'GOA')))

fshbio <- fshbio %>% 
  left_join(akrspp %>% 
              select(akr_species_codes = akr_code, common_name))


write_csv(fshbio, paste0(dat_path, "/fishery_ages.csv"))
