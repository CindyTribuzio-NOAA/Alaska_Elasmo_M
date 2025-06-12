# get tmax using data available from Alaska-based surveys and fisheries
# contact cindy.tribuzio@noaa.gov
# last updated 2025-6-12

# get three options for tmax 
# 1) max age in sample - use
# 2) mean of top 5 ages in sample - use
# 3) 99th percentile in sample - don't use

# setup ----
libs <- c("tidyverse")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# use queries.R to get updated data

out_path <- paste0("results")
dir.create(out_path)
dat_path <- paste0("data")

# data ----

srvbio <- read_csv(paste0(getwd(), '/', dat_path, "/survey_ages_akfin.csv")) # use survey_ages_akfin.csv to reproduce what's in results/YEAR/using_akfin_survey_specimen_data. RACE only delivers data once to AKFIN per year
fshbio <- read_csv(paste0(getwd(), '/', dat_path, "/fishery_ages.csv"))
spp <- read_csv(paste0(dat_path, '/race_spp.csv'))

# survey ----

srvbio <- srvbio %>% 
  left_join(spp %>% select(-species_name)) %>% 
  select(area, source = survey_name, year, common_name, sex, age, length = length_mm) 
  
# just a test for errors
srvbio %>% 
  group_by(source, year, area) %>% 
  summarize(n = n()) %>% 
  arrange(n) %>% 
  print(n=Inf)

# below step from rockfish work, not relevant to elasmos at this time
# get rid of the oddball observations that lack relevant sampling plan
# documentation, including specimens collected (1) on the west coast in a
# cooperative hake survey in 1988 and 1990, and (2) in the BS in a US-Japanese
# cooperative survey in 1979. Keep pre-survey data (e.g. GRFSH & RFSH KODK and AI)
#srvbio <- srvbio %>% filter(!area %in% c('WC', 'HWC') & !source %in% c('USJPN COOP'))

tmpsrv <- srvbio %>% 
  mutate(length = length / 10,
         sex = case_when(sex == 1 ~ 'M',
                         sex == 2 ~ 'F',
                         sex == 3 ~ 'U'),
         gear = 'AFSC bottom trawl survey',
         source = 'RACEBASE',
         common_name = ifelse(common_name == 'Alaska skate', 'alaska skate', common_name))

names(tmpsrv)

# fishery ----

tmpfsh <- fshbio %>%
  mutate(common_name = ifelse(common_name == 'shark, spiny dogfish', 'spiny dogfish', common_name),
         gear = case_when(gear_description == 'NON PELAGIC' ~ 'Bottom trawl fishery',
                            gear_description == 'PELAGIC' ~ 'Pelagic trawl fishery',
                            gear_description == 'LONGLINER' ~ 'Hook-and-line fishery'), 
         source = 'NORPAC') %>% 
  select(area, source, year, common_name, sex, age, length, gear) 

# combined ----

comb <- bind_rows(tmpsrv, tmpfsh) %>% 
  mutate(common_name = factor(common_name,
                              labels = c('alaska skate', 'longnose skate',
                                         'big skate','spiny dogfish'),
                              levels = c('alaska skate', 'longnose skate',
                                         'big skate','spiny dogfish'),
                              ordered = TRUE)) %>% 
  arrange(common_name, source, year)

write_csv(comb, paste0(dat_path, '/all_age_data_used.csv'))

combsum <- comb %>% 
  group_by(common_name, area) %>% 
  summarise(n_total = n(),
            n_survey = length(which(source == 'RACEBASE')),
            n_fishery = length(which(source == 'NORPAC')),
            q99 = quantile(age, 0.99),
            max_age = max(age)) %>% 
  left_join(comb %>% 
              group_by(common_name, area) %>% 
              slice_max(order_by = age, n = 5) %>% 
              summarize(mean_top5 = mean(age))) %>% 
  left_join(comb %>%
              group_by(common_name, area) %>% 
              summarise(years_survey = paste(sort(unique(year[which(source == 'RACEBASE')])), collapse = ', '),
                        years_fishery = paste(sort(unique(year[which(source == 'NORPAC')])), collapse = ', '))) %>% 
  arrange(common_name, area) %>% 
  select(`Species` = common_name, Region = area, `Total N` = n_total,
         `Survey N` = n_survey, `Fishery N` = n_fishery,
         `Maximum age observed` = max_age,
         `Mean top 5 ages` = mean_top5,
         `Survey years` = years_survey, `Fishery years` = years_fishery) %>% 
  mutate(Species2 = ifelse(Species == "Alaska skate", 'alaska skate', Species))
  

combsum %>% print(n=Inf)

write_csv(combsum, paste0(out_path, '/tmax_summary.csv'))

top5 <- comb %>% 
  group_by(common_name, area) %>% 
  slice_max(order_by = age, n = 5) %>% 
  arrange(common_name, area, -age) %>%
  select(`Species` = common_name, Region = area, `Year sampled` = year, Sex = sex,
         `Age (yr)` = age, `Fork length (cm)` = length, Gear = gear,  Source = source) 

top5 %>% print(n=Inf)

write_csv(top5, paste0(out_path, '/top5_tmax_detailed.csv'))

# graphic Cindy requested 
comb <- read_csv(paste0(dat_path, '/all_age_data_used.csv'))
df <- comb %>% 
  droplevels() %>% 
  group_by(common_name, area, source, year) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  tidyr::complete(common_name, area, source, year, fill = list(n = 0)) 

# https://github.com/ben-williams/funcr/blob/master/R/tickr.R
tickr <- function(data, var, to = 5, start = NULL, end = NULL, min = NULL){
  
  data %>%
    dplyr::summarise(min = min({{var}}, na.rm = T),
                     max = max({{var}}, na.rm = T)) -> out
  
  if(is.null(start) & is.null(end)){
    
    data.frame(breaks = out$min:out$max) %>%
      dplyr::mutate(labels = ifelse(breaks %in%
                                      seq(to * min(breaks) / to,
                                          max(breaks), by = to), breaks, ""))
    
  } else if(!is.null(start) & is.null(end) & is.null(min)){
    
    data.frame(breaks = start:out$max) %>%
      dplyr::mutate(labels = ifelse(breaks %in%
                                      seq(to * start / to, max(breaks),
                                          by = to), breaks, ""))
    
  } else if(!is.null(start) & is.null(end) & !is.null(min)){
    data.frame(breaks = start:out$max) %>%
      dplyr::mutate(labels = ifelse(breaks %in%
                                      seq(to * start / to, max(breaks),
                                          by = to), breaks, "")) %>%
      dplyr::filter(breaks >= min) -> lb
    lb$labels[1] <- lb$breaks[1]
    lb
    
  } else if(is.null(start) & !is.null(end)){
    
    data.frame(breaks = out$min:end) %>%
      dplyr::mutate(labels = ifelse(breaks %in%
                                      seq(to * min(breaks) / to, end, by = to),
                                    breaks, ""))
  } else {
    
    data.frame(breaks = start:end) %>%
      dplyr::mutate(labels = ifelse(breaks %in%
                                      seq(to * start / to, end, by = to),
                                    breaks, ""))
  }
}

internet_theme <- theme_set(axis.text.x = theme_text(angle = 90,
                                                     hjust = 1), 
                            panel.grid.major = theme_line(colour = "grey90"),
                            panel.grid.minor = theme_blank(), 
                            panel.background = theme_blank(),
                            axis.ticks = theme_blank(), 
                            legend.position = "none")

tmp <- df %>%  filter(area == 'GOA') 

xaxis <- tickr(tmp, year, 10, min = 1978, start = 1980)
tmp %>% 
  ggplot(aes(x = year, y = n, fill = source)) +
  geom_col(width = 0.7) +
  facet_wrap(~ common_name, scales = 'free') +
  scale_x_continuous(labels = xaxis$labels, breaks = xaxis$breaks) +
  labs(x = NULL, y = 'Number of observations', 
       title = 'GOA specimen data summary')

