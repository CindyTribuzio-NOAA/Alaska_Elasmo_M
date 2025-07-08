# BSAI skate specific Updates ----
# Species specific updates to M values
# copied from OROX>Outside analytics>2022_M_approach_updates>M_updates.R

# Setup ----
libs <- c("tidyverse", "janitor", "readr", "rstan", "fishmethods", "truncnorm", "shiny", "googlesheets4")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

datadir<-paste(getwd(),"/Data/",sep="")
outdir<-paste(getwd(),"/Results/",sep="")

priorcomb<-function(medianvec,sdvec,log,interval){
  #medianvec = vector of medians in real space = means in real space and gives means in log space when transformed
  #sd vec = standard deviations in either real or log space
  #log indicates if in real (0) or log space (anything else - I put in 1)
  #interval = confidence interval coverage - typcially 0.9 or 0.95
  n <- length(medianvec)
  if(length(sdvec)==1){
    sdvec=rep(sdvec,n)}
  invar <- 1/(sdvec^2)
  wt<- invar/sum(invar)
  if(log==0){
    mean <- sum(wt*medianvec)
    sd <- sqrt(1/sum(invar))
    upper<-max(medianvec+3*sdvec)
    lower<-min(medianvec-3*sdvec)
    range<-upper-lower
    x <- lower+c(0:1000)*range/1000
    y <- dnorm(x,mean,sd)
    ui <- qnorm((1+interval)/2,mean,sd)
    li <- qnorm((1-interval)/2,mean,sd)
    plot(x,y,type = "l",lty=4,xlab ="M", ylab="density")
    for(i in 1:n){
      m<-medianvec[i]
      sdm<-sdvec[i]
      y<-dnorm(x,m,sdm)
      lines(x,y)}
    col1<-c("mean","sd","CIlower","CIupper")
    return(cbind(col1,c(mean,sd,li,ui)))
  }
  else
  {
    lmedianvec<-log(medianvec)
    lmin <- min(lmedianvec)
    lmax <- max(lmedianvec)
    meanl <- sum(wt*lmedianvec)
    sd <- sqrt(1/sum(invar))
    upper<-max(lmedianvec+3*sdvec)
    lower<-min(lmedianvec-3*sdvec)
    range<-upper-lower
    x <- lower+c(0:1000)*range/1000
    y <- dnorm(x,meanl,sd)
    upperl<-exp(max(lmedianvec+sdvec))
    xx<-c(0:1000)*upperl/1000
    yy <- dlnorm(xx,meanl,sd)
    ui <- qnorm((1+interval)/2,meanl,sd)
    li <- qnorm((1-interval)/2,meanl,sd)
    par(mfrow=c(2,1))
    plot(x,y,type = "l",lty=4,xlab ="ln M", ylab="density")
    for(i in 1:n){
      m<-lmedianvec[i]
      sdm<-sdvec[i]
      y<-dnorm(x,m,sdm)
      lines(x,y)}
    plot(xx,yy,type = "l",lty=4,xlab ="M", ylab="density")
    lines(xx,yy,type="l",lty=4)
    for(i in 1:n){
      m<-lmedianvec[i]
      sdm<-sdvec[i]
      yy<-dlnorm(xx,m,sdm)
      lines(xx,yy)}
    col1<-c("logmean","logsd","min","max","Mean","Median","CIlower","CIupper")
    return(data.frame(Label=col1,Value=signif(c(meanl,sd,exp(lmin), exp(lmax), exp(meanl+(sd^2)/2),exp(meanl),exp(li),exp(ui)),3)))
  }
}

#Species list w/codes
#spec_list <- read_csv("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/AFSC_GOA_OROX_Assessment/Data/Static/GOA_OROX_codes.csv") %>% 
#  filter(RACE_name != "dusky rockfish")
#
#spec_cmn <- c(spec_list$RACE_name)

# Get updated M values ----
#from github.com/CindyTribuzio-NOAA

urlfile_raw <- "https://raw.githubusercontent.com/CindyTribuzio-NOAA/Alaska_Elasmo_M/refs/heads/main/results/skates/M_estimates_skates.csv"
urlfile_formatted <- "https://raw.githubusercontent.com/CindyTribuzio-NOAA/Alaska_Elasmo_M/refs/heads/main/results/skates/formatted_M_results_skates.csv"

Mup_raw <- read_csv(url(urlfile_raw)) %>% 
  pivot_longer(cols = !c("species", "version", "M_method"), 
               names_to =  "area",
               values_to = "M_estimate") %>% 
  filter(!is.na(M_estimate))
M_method <- unique(Mup_raw$M_method) %>% 
  sort()

Mup_form <- read_csv(url(urlfile_formatted)) %>% 
  select(!M_estimate) %>% 
  mutate(version = as.numeric(unlist(str_extract_all(M_method, '(?<=.v)\\d+'))),
         M_method = gsub('\\..*', '', M_method))

M_meth2 <- Mup_form %>% 
  mutate(M_method = gsub('\\..*', '', M_method)) %>% 
  group_by(M_method) %>% 
  summarise(outtest = unique(data_input)) %>% 
  rename(data_input = outtest)

Mup_raw <- Mup_raw %>% 
  left_join(M_meth2, by = 'M_method') %>% 
  left_join(Mup_form, by = c('species', 'area', 'version', 'data_input', 'M_method'))

# Current assessment M values ----
Mspecs <- unique(Mup_raw$species)
#NOTE: all skates in Tier 5 M = 0.1, the AK skate Tier 3 model M = 0.13


#M_assess <- read_csv(paste(datadir, "OROX_natM.csv", sep = "")) %>%
#  clean_names() %>% 
#  mutate(spec_name = tolower(spec_name)) %>% 
#  filter(spec_name %in% Mspecs) %>% 
#  rename(species = spec_name,
#         M_assess = m) %>% 
#  select(!species_code)

# Weighting decisions ----
# general decision rules (these are in the elasmo_M_updates google sheet: 
# 1) when multiple inputs for a given method, weights should sum to 1
# 2) if multiple inputs, and method needs to be downweighted, apply that after step 1
# 3) data downweight by 0.5 if assumed to be biased
# 2) data downweight by 0.5 if borrowed from another region
# 3) data downweight by 0.25 if input data are highly uncertain (e.g., dry weight/temp)
# 4) method downweight by 0.5 if reliant on M/K ratio (holdover from rockfish decision, open for discussion)

data_wts <- read_sheet("https://docs.google.com/spreadsheets/d/1Aw0wCgkyBacbN87FJBGN8P6Iim_nDyn9DyjCPsC5pqk/edit?gid=0#gid=0") %>% 
  filter(species %nin% c('spiny_dogfish', 'blue_shark', 'pacific_sleeper_shark', 'salmon_shark')) %>% 
  mutate(version = as.numeric(unlist(str_extract_all(M_method, '(?<=.v)\\d+'))),
         M_method = gsub('\\..*', '', M_method)) %>% 
  select(species, M_method, version, area, data_weight, data_weight_reason)

Mwt <- Mup_raw %>% 
  select(species, M_method, version, area, data_input_values, references) %>% 
  left_join(data_wts)

# Alaska skate ----

########
#harlequin
# amax = 47, downweight by 0.25 because this is quite out of range from other amax values
# amax = 79, full weight because best value for all areas
# amax = 41, downweight by 0.25  because this is quite out of range from other amax values
# amax = 63, downweight by 0.5 because borrowed from neighboring area
# gsi = 0.0270, full weight because recent, updated research from same region
# linf = 30.9 and k = 0.167, downweight by 0.5 because k based method
# temp = 5.5 and drywt = 226, downweight by 0.25 because need to re-evaluate input data here

# rescale amax weights
h_norm <- 1/sum(c(0.25, 1, 0.25, 0.5))
h_tmax <- c(0.25, 1, 0.25, 0.5) * h_norm

# data weights
h_dwt <- c(h_tmax, 1, 1, 0.25)
# method weights
h_methwt <- c(rep(1,length(h_tmax)), 1, 0.5, 1)
# final weights
h_wt <- h_dwt * h_methwt

hwts <- OROX_Mup %>% 
  filter(species == "harlequin rockfish") %>% 
  mutate(dwt = h_dwt, 
         methwt = h_methwt,
         f_wt = h_wt)

###########
# redbanded
# amax = 106, full weight because local value
# linf = 54.8, k = 0.05, downweight by 0.5 because neighboring region
# temp = 5.5, drywt = 1960, downweight by 0.5 because need to re-evaluate input data

# data weights
rb_dwt <- c(1, 0.5, 0.5)
# method weights
rb_methwt <- c(1, 0.5, 1)
# final weights
rb_wt <- rb_dwt * rb_methwt

rbwts <- OROX_Mup %>% 
  filter(species == "redbanded rockfish") %>% 
  mutate(dwt = rb_dwt, 
         methwt = rb_methwt,
         f_wt = rb_wt)

#####
# redstripe
# amax = 46, full weight because local value
# amax = 55, downweight by 0.5 because neighboring area
# amax = 39, full weight because local value

# rescale amax weights
rs_norm <- 1/sum(c(1, 0.5, 1))
rs_tmax <- c(1, 0.5, 1) * rs_norm

# data weights
rs_dwt <- rs_tmax
# method weights
rs_methwt <- rep(1,length(rs_tmax))
# final weights
rs_wt <- rs_dwt * rs_methwt

rswts <- OROX_Mup %>% 
  filter(species == "redstripe rockfish") %>% 
  mutate(dwt = rs_dwt, 
         methwt = rs_methwt,
         f_wt = rs_wt)

#####
# sharpchin
# amax = 48, full weight because local value
# amax = 58, keep this as full weight, even though from WC because they actually age this species regularly
# amax = 43, full weight because local value
# Linf = 32.6, k = 0.131, downweight by 0.5 because M/K ratio assumptions
# Linf = 34.9, k = 0.095, downweight by 0.25 because M/k ratio AND neighboring area
# Linf = 33.2, k = 0.17, downweight by 0.25 because M/k ratio and two areas away
# temp = 6, drywt = 533, downweight by 0.25 because need to re-evaluate input data

# rescale amax weights
sc_tnorm <- 1/sum(c(1, 1, 1))
sc_tmax <- c(1, 1, 1) * sc_tnorm

# rescale vb weights
sc_vbnorm <- 1/sum(c(1, 0.5, 0.5))
sc_vb <- c(1, 0.5, 0.5)* sc_vbnorm

# data weights
sc_dwt <- c(sc_tmax, sc_vb, 0.25)
# method weights
sc_methwt <- c(rep(1,length(sc_tmax)), rep(0.5, length(sc_vb)), 1)
# final weights
sc_wt <- sc_dwt * sc_methwt

scwts <- OROX_Mup %>% 
  filter(species == "sharpchin rockfish") %>% 
  mutate(dwt = sc_dwt, 
         methwt = sc_methwt,
         f_wt = sc_wt)

#####
# silvergray
# amax = 79, full weight because local value
# amax = 81, downweight by 0.5 because neighboring value
# amax = 71, full weight because local value
# amax = 75, full weight because local value
# temp = 7, drywt = 1960, downweight by 0.25 because input values need to be re-evaluated

# rescale amax weights
sg_norm <- 1/sum(c(1, 0.5, 1, 1))
sg_tmax <- c(1, 0.5, 1, 1) * sg_norm

# data weights
sg_dwt <- c(sg_tmax, 0.25)
# method weights
sg_methwt <- c(rep(1,length(sg_tmax)), 1)
# final weights
sg_wt <- sg_dwt * sg_methwt

sgwts <- OROX_Mup %>% 
  filter(species == "silvergray rockfish") %>% 
  mutate(dwt = sg_dwt, 
         methwt = sg_methwt,
         f_wt = sg_wt)

#####
# yelloweye
# amax = 122, downweight by 0.5 because of K Munk has history of being biased high
# amax = 115, downweight by 0.5 because different region
# amax = 114, full weight because local values, even though old
# amax = 118, full weight because local values, even though old
# gsi = 0.0285, full weight because based on recent local research
# temp = 6 and drywt = 4200, downweight by 0.25 because need to re-evaluate input data here
yewt <- c(0.25, 0.5, 0.5, 0.5, 1, 0.25)

# rescale amax weights
ye_norm <- 1/sum(c(0.5, 0.5, 1, 1))
ye_tmax <- c(0.5, 0.5, 1, 1) * ye_norm

# data weights
ye_dwt <- c(ye_tmax, 1, 0.25)
# method weights
ye_methwt <- c(rep(1,length(ye_tmax)), 1, 1)
# final weights
ye_wt <- ye_dwt * ye_methwt

yewts <- OROX_Mup %>% 
  filter(species == "yelloweye rockfish") %>% 
  mutate(dwt = ye_dwt, 
         methwt = ye_methwt,
         f_wt = ye_wt)

# add wts to OROX_Mup
OROX_Mup <- bind_rows(hwts, rbwts, rswts, scwts, sgwts, yewts)

# Additional M values ----
# created after the Sullivan et al. analysis

new_vals <- data.frame(matrix(ncol = ncol(OROX_Mup), nrow = 2))
colnames(new_vals) <- names(OROX_Mup)

# harlequin
# add 4 years to the largest max age to account for demonstrated underageing (Kastelle et al. 2020), amax = 85
new_vals$species <- "harlequin rockfish"
new_vals$version <-4
new_vals$M_method <- "amax"
new_vals$area <- c("AI", "GOA")
new_vals$M_estimate <- c(Then_M(83)[3], Then_M(51)[3])
new_vals$data_input <- "Max age (y)"
new_vals$data_input_values <- c(83, 51)
new_vals$references <- "TenBrink with Kastelle correction"
new_vals$dwt <- 0
new_vals$methwt <- 0
new_vals$f_wt <- 0

OROX_Mup <- rbind(OROX_Mup, new_vals)

# rescale with this value downweighted by 0.5
# rescale amax weights
h_norm <- 1/sum(c(0.25, 1, 0.25, 0.5, 0.75, 0.25))
h_tmax <- c(0.25, 1, 0.25, 0.5, 0.75, 0.25) * h_norm

# data weights
h_dwt <- c(h_tmax, 1, 1, 0.25)
# method weights
h_methwt <- c(rep(1,length(h_tmax)), 1, 0.5, 1)
# final weights
h_wt <- h_dwt * h_methwt

hwts <- OROX_Mup %>% 
  filter(species == "harlequin rockfish") %>% 
  arrange(M_method) %>% 
  mutate(dwt = h_dwt, 
         methwt = h_methwt,
         f_wt = h_wt)

OROX_Mup <- OROX_Mup %>% 
  filter(species != "harlequin rockfish") %>% 
  bind_rows(hwts) %>% 
  filter(f_wt > 0)

# Composite Estimates by species ----
Comp_M <- OROX_Mup %>%
  group_by(species) %>% 
  summarise(M_Low = min(M_estimate),
            M_High = max(M_estimate),
            M_Med = median(M_estimate),
            M_Mean = mean(M_estimate),
            M_comp = weighted.mean(M_estimate, f_wt), #note the weighted.mean command doesn't handle cases of n=1
            M_rec = M_comp) %>% 
  mutate(Method = "Wt_M") 

# Hamel and Cope Prior distributions ----
# H and C report that a CV of 0.31 is appropriate
# for CVs use following rules:
# amax CV = 0.31 as per Hamel and Cope, noting a direct link between age and M
# lvb CV = 0.85 as per Hamel and Cope, noting this is a secondary relationship between age and M
# gsi CV = 0.31 because method is a direct link between reproductive energy and M, using amax as a proxy "good" value
# temp CV = 0.85 using lvb as a proxy "bad" value, not going to re-analyze all of the McG data

OROX_Mup <- OROX_Mup %>% 
  mutate(CV = if_else(M_method == "amax", 0.31,
                      if_else(M_method == "lvb", 0.85,
                              if_else(M_method == "gsi", 0.31, 0.85))))

HC_M_prior<-OROX_Mup %>% 
  group_by(species) %>% 
  summarize(M_Low = priorcomb(M_estimate, CV, "lognormal",0.95)[3,2],
            M_High = priorcomb(M_estimate, CV, "lognormal",0.95)[4,2],
            M_Mean = priorcomb(M_estimate, CV, "lognormal",0.95)[5,2],
            M_Med = priorcomb(M_estimate, CV, "lognormal",0.95)[6,2],
            LL = priorcomb(M_estimate, CV, "lognormal",0.95)[7,2],
            UL = priorcomb(M_estimate, CV, "lognormal",0.95)[8,2]) %>% 
  mutate(Method = "HCCV")

HC_M_prior_CV31<-OROX_Mup %>% 
  group_by(species) %>% 
  summarize(M_Low = priorcomb(M_estimate, rep(0.31, length(M_estimate)), "lognormal",0.95)[3,2],
            M_High = priorcomb(M_estimate, rep(0.31, length(M_estimate)), "lognormal",0.95)[4,2],
            M_Mean = priorcomb(M_estimate, rep(0.31, length(M_estimate)), "lognormal",0.95)[5,2],
            M_Med = priorcomb(M_estimate, rep(0.31, length(M_estimate)), "lognormal",0.95)[6,2],
            LL = priorcomb(M_estimate, rep(0.31, length(M_estimate)), "lognormal",0.95)[7,2],
            UL = priorcomb(M_estimate, rep(0.31, length(M_estimate)), "lognormal",0.95)[8,2]) %>% 
  mutate(Method = "HC")

# To avoid psuedo replication, use the weighted values as input 
HC_input <- OROX_Mup %>% 
  mutate(wtm = M_estimate*f_wt) %>% 
  group_by(species, M_method) %>% 
  summarise(M_comp = sum(wtm)) %>% 
  mutate(CV = if_else(M_method == "amax", 0.31,
                      if_else(M_method == "lvb", 0.85,
                              if_else(M_method == "gsi", 0.31, 0.85))))

# but some species only have one method for those species, use all inputs, filter those out and use values from above
#spec_1 <- HC_input %>% 
#  group_by(species) %>% 
#  summarise(tot_N = length(M_method)) %>% 
#  filter(tot_N == 1)

#spec_1_HC <- HC_M_prior %>% 
#  filter(species %in% spec_1$species) %>% 
#  transform(Method = "HCCV_Wt_M")

HC_wtM_prior <- HC_input %>% 
  #filter(species %nin% spec_1$species) %>% 
  group_by(species) %>% 
  summarize(M_Low = priorcomb(M_comp, CV, "lognormal",0.95)[3,2],
            M_High = priorcomb(M_comp, CV, "lognormal",0.95)[4,2],
            M_Mean = priorcomb(M_comp, CV, "lognormal",0.95)[5,2],
            M_Med = priorcomb(M_comp, CV, "lognormal",0.95)[6,2],
            LL = priorcomb(M_comp, CV, "lognormal",0.95)[7,2],
            UL = priorcomb(M_comp, CV, "lognormal",0.95)[8,2]) %>% 
  mutate(Method = "HCCV_Wt_M") 
  
#spec_1_HCCV31 <- HC_M_prior_CV31 %>% 
#  filter(species %in% spec_1$species) %>% 
#  transform(Method = "HC_Wt_M")

HC_wtM_priorCV31 <- HC_input %>% 
  #filter(species %nin% spec_1$species) %>% 
  group_by(species) %>% 
  summarize(M_Low = priorcomb(M_comp, rep(0.31, length(M_comp)), "lognormal",0.95)[3,2],
            M_High = priorcomb(M_comp, rep(0.31, length(M_comp)), "lognormal",0.95)[4,2],
            M_Mean = priorcomb(M_comp, rep(0.31, length(M_comp)), "lognormal",0.95)[5,2],
            M_Med = priorcomb(M_comp, rep(0.31, length(M_comp)), "lognormal",0.95)[6,2],
            LL = priorcomb(M_comp, rep(0.31, length(M_comp)), "lognormal",0.95)[7,2],
            UL = priorcomb(M_comp, rep(0.31, length(M_comp)), "lognormal",0.95)[8,2]) %>% 
  mutate(Method = "HC_Wt_M") 

HC_M_out <- HC_M_prior %>% 
  bind_rows(HC_M_prior_CV31, HC_wtM_prior, HC_wtM_priorCV31) %>% 
  mutate(M_rec = M_Med)

# Summary ----
OROX_M_Summary <- Comp_M %>% 
  bind_rows(HC_M_out) %>% 
  left_join(M_assess) %>%
  mutate(Prop_change = round((M_rec - M_assess)/M_assess*100, 0))%>% 
  transform(species = gsub(" rockfish", "", species),
            M_Low = round(M_Low, 3),
            M_High = round(M_High, 3),
            M_Med = round(M_Med, 3),
            M_Mean = round(M_Mean, 3),
            LL = round(LL, 3),
            UL = round(UL, 3),
            M_comp = round(M_comp, 3),
            M_rec = (round(M_rec, 3))) %>% 
  mutate(rank = if_else(Method == "Wt_M", 1,
                      if_else(Method == "HC", 2,
                              if_else(Method == "HC_Wt_M", 3,
                                      if_else(Method == "HCCV", 4, 5))))) %>% 
  arrange(species, rank) %>% 
  select(-c(rank, M_comp)) %>% 
  relocate(Method, .after = species)

write.csv(OROX_M_Summary, paste(outdir, "OROX_updatedM.csv", sep = ""), row.names = F)

plot_dat <- OROX_M_Summary %>% 
  mutate(rank = if_else(Method == "Wt_M", 1,
                        if_else(Method == "HC", 2,
                                if_else(Method == "HC_Wt_M", 3,
                                        if_else(Method == "HCCV", 4, 5))))) %>% 
  arrange(rank) %>% 
  #mutate(Method = factor(Method, levels = Method)) %>% 
  #filter(Method %in% c("Wt_M", "HCCV_Wt_M")) %>% 
  transform(species = gsub(" rockfish", "", species))

OROX_Mup <- OROX_Mup %>% 
  transform(species = gsub(" rockfish", "", species))
write.csv(OROX_Mup, paste(outdir, "OROX_Minputs.csv", sep = ""), row.names = F)

M_plot <- ggplot(plot_dat, aes(x = as.factor(Method), y = M_rec))+
  #geom_point(data = OROX_Mup, aes(x= "Composite", y = M_estimate, size = f_wt), color = "grey75")+
  geom_jitter(data = OROX_Mup, aes(x= "Wt_M", y = M_estimate, size = f_wt, shape = M_method), color = "grey50", alpha = 0.5)+
  geom_segment(aes(x = Method, y = M_assess, xend = Method, yend = M_rec), size = 1, color = "red")+
  geom_hline(aes(yintercept = M_assess))+
  geom_errorbar(aes(ymin = LL, ymax = UL), width = 0.5)+
  geom_point(size = 4, aes(color = Prop_change))+
  #scale_x_discrete(limits = rev)+
  xlab("Composite Method")+
  ylab("Composite M Value")+
  labs(shape = "M Type", size = "M Wt", color = "% Change")+
  facet_grid(species~., scales = "free")+
  theme_bw()  

ggsave("GOAOROX_updatedM.png",path = outdir, plot= M_plot,dpi=600, bg="transparent",width=7,height=8)


