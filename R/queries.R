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

write_csv(srv_ages, paste0(dat_path, "/survey_ages_akfin.csv"))

# fishery data ----

# akr species code translation for observer data
akrspp <- read_csv("data/species_translation_akr_obs.csv")

akrspp <- akrspp %>%
  filter(grepl('Skate|Shark', akr_name)) %>% 
  distinct(akr_name, akr_code)


#########################
###########################
#start here once RODBC is fixed
akrspp <- akrspp %>% 
  mutate(common_name = tolower(akr_name)) %>% 
  left_join(spp)

myakrspp <- akrspp %>% pull(akr_code)
myakrspp_string <- toString(sprintf("'%s'", myakrspp)) # allows you to pass vector into sql query
  
fshbio <- fshbio %>% 
  mutate(area = ifelse(fmp_subarea %in% c('AI', 'BS'), fmp_subarea, 'GOA'))

fshbio <- fshbio %>% 
  left_join(akrspp %>% 
              select(akr_species_codes = akr_code, common_name))

fish_age <- sqlQuery(channel_akfin, query = ("
                select   *
                from     norpac.debriefed_age_flat_mv
                where    akr_species_codes between 689 and 705 and
                        age IS NOT NULL")) %>% 
  clean_names() 

write_csv(fshbio, paste0(dat_path, "/fishery_ages.csv"))
