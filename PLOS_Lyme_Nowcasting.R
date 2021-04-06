###############################################
# Lyme Nowcasting code for PLOS One submission
# Created on 1/12/2021
###############################################

#### DATA PREPERATION ####
setwd() #Allows user to set their own working directory

#### LOAD NEEDED PACKAGES ####
install.packages(c("tidyverse", "gtrendsR", "caret", "glmnet", "cowplot"))
library(tidyverse)
library(gtrendsR)
library(caret)
library(glmnet)
library(cowplot)

#### LOAD LYME REGIONAL LYME COUNT AND STATE POPULATION DATA ####
#Here we load in the Lyme disease incidence as counted from the CDC survillance data
#Population data is from the 2010 US Census
regional_sum <- read.csv("Regional_Sums.csv")

Regional_pop <- readxl::read_xlsx("Sate_populations.xlsx") %>% select(stname, Pop_2010)

#Creating Regional population counts to join with Regional Lyme disease incidence data. 
Region_agg <- Regional_pop %>% 
  mutate(Region = if_else(stname %in% c("Arkansas", "Louisiana", "Mississippi", "Tennessee",
                                        "Kentucky", "West Virginia", "Maryland", "Delaware",
                                        "District of Columbia", "Virginia", "North Carolina",
                                        "South Carolina", "Georgia", "Alabama", "Florida"), "Southeast", #Southeast
                          
                          if_else(stname %in% c("Ohio", "Indiana", "Michigan", "Illinois", "Wisconsin",
                                                "Minnesota", "Iowa", "Missouri", "North Dakota", "South Dakota",
                                                "Nebraska", "Kansas"), "Midwest", #Midwest
                                  
                                  if_else(stname %in% c("Washington", "Oregon", "California", "Nevada", "Idaho",
                                                        "Montana", "Wyoming", "Utah", "Colorado", "Alaska", 
                                                        "Hawaii"), "West", #West
                                          
                                          if_else(stname %in% c("Arizona", "New Mexico", "Texas", "Oklahoma"), "Southwest", #Southwest
                                                  
                                                  if_else(stname %in% c("Pennsylvania", "New Jersey", "New York", "Connecticut",
                                                                        "Rhode Island", "Massachusetts", "Vermont", "New Hampshire",
                                                                        "Maine"), "Northeast", "Other")))))) %>%  #Other is used here to round out the ifelse function
  group_by(Region) %>% 
  summarize(region_pop = sum(Pop_2010))
table(Region_agg) #To double check that there is no "Other" in the data
#Joining 
Region_rates <- regional_sum %>% 
  inner_join(Region_agg) %>% 
  mutate(Lyme_Rate = 100000*(total_lyme / region_pop))
#write_csv(Region_rates, path = "Region_Rates.csv") These regional rates will be used for modeling and collecting of search term data

#### SUBSETING DATA ####
#This part of code subsets the data region and training time period.
#These data subsets were used in Google Correlate to identify search terms to collect data on.
#Filtering data to only training period to identify candidate serach terms. If use the whole data, could lead to lack of accuracy of hold-out data
Region_rates <- read_csv("Region_Rates.csv")
Northeast_rate <- Region_rates %>% 
  group_by(Region) %>% 
  filter(Region == "Northeast")
Northeast_rate_train <- Northeast_rate %>% filter(date <= "2012-12-01") #This is the end of the training period for the models. 

Midwest_rate <- Region_rates %>% 
  group_by(Region) %>% 
  filter(Region == "Midwest")
Midwest_rate_train <- Midwest_rate %>% filter(date <= "2012-12-01")

Southeast_rate <- Region_rates %>% 
  group_by(Region) %>% 
  filter(Region == "Southeast")
Southeast_rate_train <- Southeast_rate %>% filter(date <= "2012-12-01")

Southwest_rate <- Region_rates %>% 
  group_by(Region) %>% 
  filter(Region == "Southwest")
Southwest_rate_train <- Southwest_rate %>% filter(date <= "2012-12-01")

West_rate <- Region_rates %>% 
  group_by(Region) %>% 
  filter(Region == "West")
West_rate_train <- West_rate %>% filter(date <= "2012-12-01")

# write_csv(Northeast_rate_train, path = "Correlate Data/Northeast_rates_train.csv")
# write_csv(Midwest_rate_train, path = "Correlate Data/Midwest_rates_train.csv")
# write_csv(Southeast_rate_train, path = "Correlate Data/Southeast_rate_train.csv")
# write_csv(Southwest_rate_train, path = "Correlate Data/Southwest_rate_train.csv")
# write_csv(West_rate_train, path = "Correlate Data/West_rate_train.csv")

#### COLLECTING SEARCH TERM HIT DATA FROM GOOGLE TRENDS ####
#Here we pull the Search terms identified via Google Correlate and add them to the list of Vector and disease specific Terms
# Once terms are combined, for each region we pull search traffic for each term.
# We then aggregate to the average search traffic for each term in each region.
# Aggregated search traffic is combined with Lyme Rate data for each region. 
# Saving Final and Model Data sets.
# Final data sets have all the search term data up to the date that Google data was collected (Sept. 2019)
# Model data ends December 2017 because that is when CDC surveillance data ends. Can not model with NAs in the outcome therefore we cut the data off there for modeling purposes. 

# Function used to collect serach term data for each state
State_SearchTimeseries <- function(term, geo){
  print(paste(term, geo))
  Sys.sleep(2)
  gt_data_1 <- gtrends(term, geo = geo, time = "all")
  if(is.null(gt_data_1$interest_over_time)){
    out <- tibble(date=character(),
                  hits=integer(),
                  geo = character())
    return(out)
  } else{
    out <- gt_data_1$interest_over_time %>% 
      as_data_frame() %>% 
      select(date,hits)
    return(out)
  }
}
#Load in Google Correlate data
NorthEast_gt_lyme <- read_csv("Correlate Data/correlate-Northeast_Rates_Corr.csv", comment = "#")

NorthEast_search_terms <- c(names(NorthEast_gt_lyme)[3:102],"tick", "black tick", "lyme", "lyme disease", "rash", "bullseye rash",
                            "bell's palsy", "facial paralysis", "side of face paralyzed", "knee pain", "swollen knee", "swollen knees",
                            "swollen joint", "swollen joints", "joint pain", "fever", "tired", "deer tick", "black-legged tick", "black legged tick",
                            "black leg tick")
NorthEast_Table <- tibble(geo = rep(c("US-PA", "US-NY", "US-VT", "US-NH", "US-ME", "US-MA",
                                      "US-RI", "US-CT", "US-NJ"), #Using state aberviations here becuase that is what Google Trends implements to pull data from the state level
                                    length(NorthEast_search_terms))) %>% 
  arrange(geo) %>% mutate(term = rep(NorthEast_search_terms, 9))
table(NorthEast_Table$geo)#making sure each state has the same number of search terms used. 

NorthEast_timeseries <- NorthEast_Table %>% 
  mutate(data=map2(term, geo, ~State_SearchTimeseries(.x,.y)))

NorthEast_timeseries_wide <- NorthEast_timeseries  %>% 
  mutate(data=map(data,~mutate(.,hits=as.integer(hits)))) %>% 
  unnest() %>% 
  spread(key = term, value = hits)
#write_csv(NorthEast_timeseries_wide, path = "Google Hit Data/NorthEast Search Hits.csv")
#NorthEast_timeseries_wide <- read_csv("Google Hit Data/NorthEast Search Hits.csv", comment = "#") %>% select(-geo1)

NorthEast_timeseries_wide2 <- NorthEast_timeseries_wide %>% 
  group_by(date) %>% 
  summarize_at(vars("akron racers":"wild raspberries"), mean, na.rm = T)

NorthEast_final <- Northeast_rate %>% ungroup() %>% 
  select(date, Lyme_Rate) %>% 
  right_join(., NorthEast_timeseries_wide2, by = "date") %>% 
  mutate_at(vars(-date, -Lyme_Rate),funs(ifelse(is.na(.),0L,.)))

NorthEast_ModelData <- NorthEast_final %>% filter(date <= "2017-12-01")
#Double checking missing data 
sum(is.na(NorthEast_final$Lyme_Rate)) #21 missing outcome points, but this is due to lack of surveillance data between Dec. 2017 and Sept. 2019. do not need to worry about this
sum(is.na(NorthEast_ModelData)) #0 Missing data
# write_csv(NorthEast_final, path = "Final Datasets/NorthEast_Final.csv")
# write_csv(NorthEast_ModelData, path = "Model Datasets/NorthEast_Model.csv")

MidWest_gt_lyme <- read_csv("../../Correlate Data/correlate-Midwest_Rates_Corr.csv", comment = "#")

Midwest_search_terms <- c(names(MidWest_gt_lyme)[3:102],"tick", "black tick", "lyme", "lyme disease", "rash", "bullseye rash",
                          "bell's palsy", "facial paralysis", "side of face paralyzed", "knee pain", "swollen knee", "swollen knees",
                          "swollen joint", "swollen joints", "joint pain", "fever", "tired", "deer tick", "black-legged tick", "black legged tick",
                          "black leg tick")
MidWest_Table <- tibble(geo = rep(c("US-ND", "US-SD", "US-NE", "US-KS", "US-MN", "US-IA",
                                    "US-MO", "US-WI", "US-IL", "US-MI", "US-IN", "US-OH"), 
                                  length(Midwest_search_terms))) %>% 
  arrange(geo) %>% mutate(term = rep(Midwest_search_terms, 12))
table(MidWest_Table$geo)

MidWest_timeseries <- MidWest_Table %>% 
  mutate(data=map2(term, geo, ~State_SearchTimeseries(.x,.y)))

MidWest_timeseries_wide <- MidWest_timeseries  %>% 
  mutate(data=map(data,~mutate(.,hits=as.integer(hits)))) %>% 
  unnest() %>% 
  spread(key = term, value = hits)
#write_csv(MidWest_timeseries_wide, path = "MidWest Search Hits.csv")
#MidWest_timeseries_wide <- read_csv("Google Hit Data/MidWest Search Hits.csv", comment = "#") %>% select(-geo1)

MidWest_timeseries_wide2 <- MidWest_timeseries_wide %>% 
  group_by(date) %>% 
  summarize_at(vars("49 drive in":"worlds of fun discount"), mean, na.rm = T)

MidWest_final <- Midwest_rate %>% ungroup() %>% 
  select(date, Lyme_Rate) %>% 
  right_join(., MidWest_timeseries_wide2, by = "date") %>% 
  mutate_at(vars(-date, -Lyme_Rate),funs(ifelse(is.na(.),0L,.)))

MidWest_ModelData <- MidWest_final %>% filter(date <= "2017-12-01")

#Double checking missing data 
sum(is.na(MidWest_final$Lyme_Rate)) #21 missing outcome points, but this is due to lack of surveillance data between Dec. 2017 and Sept. 2019. do not need to worry about this
sum(is.na(MidWest_ModelData))
# write_csv(MidWest_final, path = "Final Datasets/MidWest_Final.csv")
# write_csv(MidWest_ModelData, path = "Model Datasets/MidWest_Model.csv")

SouthEast_gt_lyme <- read_csv("../../Correlate Data/correlate-Southeast_Rates_Corr.csv", comment = "#")

SouthEast_search_terms <- c(names(SouthEast_gt_lyme)[3:102],"tick", "black tick", "lyme", "lyme disease", "rash", "bullseye rash",
                            "bell's palsy", "facial paralysis", "side of face paralyzed", "knee pain", "swollen knee", "swollen knees",
                            "swollen joint", "swollen joints", "joint pain", "fever", "tired", "deer tick", "black-legged tick", "black legged tick",
                            "black leg tick")
SouthEast_Table1 <- tibble(geo = rep(c("US-AR", "US-LA", "US-MS", "US-TN", "US-KY", "US-WV",
                                       "US-VA", "US-NC", "US-SC", "US-GA"), 
                                     length(SouthEast_search_terms))) %>% 
  arrange(geo) %>% mutate(term = rep(SouthEast_search_terms, 10))
table(SouthEast_Table1$geo)

SouthEast_Table2 <- tibble(geo = rep(c("US-FL", "US-AL","US-MD", "US-DE", "US-DC"), 
                                     length(SouthEast_search_terms))) %>% 
  arrange(geo) %>% mutate(term = rep(SouthEast_search_terms, 5))
table(SouthEast_Table2$geo)

SouthEast_timeseries1 <- SouthEast_Table1 %>% 
  mutate(data=map2(term, geo, ~State_SearchTimeseries(.x,.y)))

SouthEast_timeseries2 <- SouthEast_Table2 %>% 
  mutate(data=map2(term, geo, ~State_SearchTimeseries(.x,.y)))

SouthEast_timeseries <- rbind(SouthEast_timeseries1, SouthEast_timeseries2)

SouthEast_timeseries_wide <- SouthEast_timeseries %>% 
  mutate(data=map(data,~mutate(.,hits=as.integer(hits),
                               date = as.Date(date)))) %>% 
  unnest() %>% 
  spread(key = term, value = hits) %>% select(-geo1)
#write_csv(SouthEast_timeseries_wide, path = "SouthEast Search Hits.csv")
#SouthEast_timeseries_wide <- read_csv("Google Hit Data/SouthEast Search Hits.csv", comment = "#")

SouthEast_timeseries_wide2 <- SouthEast_timeseries_wide %>%
  group_by(date) %>%
  summarize_at(vars("alabama water park":"wood bat tournament"), mean, na.rm = T)

SouthEast_final <- Southeast_rate %>% ungroup() %>%
  select(date, Lyme_Rate) %>%
  right_join(., SouthEast_timeseries_wide2, by = "date") %>%
  mutate_at(vars(-date, -Lyme_Rate),funs(ifelse(is.na(.),0L,.)))

SouthEast_ModelData <- SouthEast_final %>% filter(date <= "2017-12-01")
#Double checking missing data 
sum(is.na(SouthEast_final$Lyme_Rate)) #21 missing outcome points, but this is due to lack of surveillance data between Dec. 2017 and Sept. 2019. do not need to worry about this
sum(is.na(SouthEast_ModelData))
# write_csv(SouthEast_final, path = "Final Datasets/SouthEast_Final.csv")
# write_csv(SouthEast_ModelData, path = "Model Datasets/SouthEast_Model.csv")

SouthWest_gt_lyme <- read_csv("../../Correlate Data/correlate-Southwest_Rates_Corr.csv", comment = "#")

SouthWest_search_terms <- c(names(SouthWest_gt_lyme)[3:102],"tick", "black tick", "lyme", "lyme disease", "rash", "bullseye rash",
                            "bell's palsy", "facial paralysis", "side of face paralyzed", "knee pain", "swollen knee", "swollen knees",
                            "swollen joint", "swollen joints", "joint pain", "fever", "tired", "deer tick", "black-legged tick", "black legged tick",
                            "black leg tick")
SouthWest_Table <- tibble(geo = rep(c("US-AZ", "US-NM", "US-TX", "US-OK"), 
                                    length(SouthWest_search_terms))) %>% 
  arrange(geo) %>% mutate(term = rep(SouthWest_search_terms, 4))
table(SouthWest_Table$geo)

SouthWest_timeseries <- SouthWest_Table %>% 
  mutate(data=map2(term, geo, ~State_SearchTimeseries(.x,.y)))

SouthWest_timeseries_wide <- SouthWest_timeseries  %>% 
  mutate(data=map(data,~mutate(.,hits=as.integer(hits),
                               date = as.Date(date)))) %>% 
  unnest() %>% 
  spread(key = term, value = hits) %>% select(-geo1)

#write_csv(SouthWest_timeseries_wide, path = "SouthWest Search Hits.csv")
#SouthWest_timeseries_wide <- read_csv("Google Hit Data/SouthWest Search Hits.csv", comment = "#")

SouthWest_timeseries_wide2 <- SouthWest_timeseries_wide %>% 
  group_by(date) %>% 
  summarize_at(vars("bell's palsy":"world rv"), mean, na.rm = T)

SouthWest_final <- Southwest_rate %>% ungroup() %>% 
  select(date, Lyme_Rate) %>% 
  right_join(., SouthWest_timeseries_wide2, by = "date") %>% 
  mutate_at(vars(-date, -Lyme_Rate),funs(ifelse(is.na(.),0L,.)))

SouthWest_ModelData <- SouthWest_final %>% filter(date <= "2017-12-01")

#Checking of missing data in outcome
sum(is.na(SouthWest_final$Lyme_Rate)) #39 missing data points. this indicates that there is more missing than just data after surveillance data ends
sum(is.na(SouthWest_ModelData$Lyme_Rate)) #18 missing outcome data points. 
#The missing here is due to no incident cases reported for:
# Feb, March, Sept, and Dec of 2009; Nov and Dec 2005; Jan, March, April, June, and Dec 2007; Nov 2010; Jan, Feb, and Oct 2014; Jan, March, and Dec  2015
#To double check this we run the code below
SouthWest_ModelData$date[is.na(SouthWest_ModelData$Lyme_Rate)]
#Want to make sure that there is no missing data in our outcome variables.
#Filling in these missing values of no reported cases in a month with 0
SouthWest_ModelData$Lyme_Rate[is.na(SouthWest_ModelData$Lyme_Rate)] <- 0
#Double checking no missing data.
sum(is.na(SouthWest_ModelData$Lyme_Rate)) #0 missing data no
sum(is.na(SouthWest_ModelData))#0 missing data overall

na.index <- is.na(SouthWest_final$Lyme_Rate)
na.index_surv <- c(na.index[1:168],rep(as.logical(FALSE), 21)) #Need to add these 21 Falses at the end to no override NAs where we do not have surveillance data
SouthWest_final$Lyme_Rate[na.index_surv]<-0
sum(is.na(SouthWest_final$Lyme_Rate))  #21 missing outcome points, but this is due to lack of surveillance data between Dec. 2017 and Sept. 2019. do not need to worry about this

# write_csv(SouthWest_final, path = "Final Datasets/SouthWest_Final.csv")
# write_csv(SouthWest_ModelData, path = "Model Datasets/SouthWest_Model.csv")

West_gt_lyme <- read_csv("../../Correlate Data/correlate-West_Rates_Corr.csv", comment = "#")

West_search_terms <- c(names(West_gt_lyme)[3:102],"tick", "black tick", "lyme", "lyme disease", "rash", "bullseye rash",
                       "bell's palsy", "facial paralysis", "side of face paralyzed", "knee pain", "swollen knee", "swollen knees",
                       "swollen joint", "swollen joints", "joint pain", "fever", "tired", "deer tick", "black-legged tick", "black legged tick",
                       "black leg tick")
West_Table <- tibble(geo = rep(c("US-WA", "US-OR", "US-CA", "US-NV", "US-ID", "US-UT",
                                 "US-MT", "US-WY", "US-CO"), 
                               length(West_search_terms))) %>% 
  arrange(geo) %>% mutate(term = rep(West_search_terms, 9))
table(West_Table$geo)

West_timeseries <- West_Table %>% 
  mutate(data=map2(term, geo, ~State_SearchTimeseries(.x,.y)))

West_timeseries_wide <- West_timeseries  %>% 
  mutate(data=map(data,~mutate(.,hits=as.integer(hits),
                               date = as.Date(date)))) %>% 
  unnest() %>% 
  spread(key = term, value = hits)

#write_csv(West_timeseries_wide, path = "West Search Hits.csv")
#West_timeseries_wide <- read_csv("Google Hit Data/West Search Hits.csv", comment = "#") %>% select(-geo1)

West_timeseries_wide2 <- West_timeseries_wide %>% 
  group_by(date) %>% 
  summarize_at(vars( "bandshell":"zucchini flower" ), mean, na.rm = T)

West_final <- West_rate %>% ungroup() %>% 
  select(date, Lyme_Rate) %>% 
  right_join(., West_timeseries_wide2, by = "date") %>% 
  mutate_at(vars(-date, -Lyme_Rate),funs(ifelse(is.na(.),0L,.)))

West_ModelData <- West_final %>% filter(date <= "2017-12-01")

#Checking for missing values in outcome
sum(is.na(West_final$Lyme_Rate)) #26 missing data points
sum(is.na(West_ModelData$Lyme_Rate)) # 5 missing outcome data points. 
#The missing here is due to no incident cases reported for:
# Feb and Oct 2005; Dec 2006; March and Dec 2011 
#To double check this we run the code below
West_ModelData$date[is.na(West_ModelData$Lyme_Rate)]
#Filling in these missing values of no reported cases in a month with 0
West_ModelData$Lyme_Rate[is.na(West_ModelData$Lyme_Rate)] <- 0
#Double checking no missing data.
sum(is.na(West_ModelData$Lyme_Rate)) #0 missing data no
sum(is.na(West_ModelData))#0 missing data overall

na.index <- is.na(West_final$Lyme_Rate)
na.index_surv <- c(na.index[1:168],rep(as.logical(FALSE), 21)) #Need to add these 21 Falses at the end to no override NAs where we do not have surveillance data
West_final$Lyme_Rate[na.index_surv]<-0
sum(is.na(West_final$Lyme_Rate))  #21 missing outcome points, but this is due to lack of surveillance data between Dec. 2017 and Sept. 2019. do not need to worry about this

# write_csv(West_final, path = "Final Datasets/West_Final.csv")
# write_csv(West_ModelData, path = "Model Datasets/West_Model.csv")

#### LOAD IN  DATA SETS ####
#Northeast
Northeast_final <- read_csv(file = "Final Datasets/NorthEast_Final.csv")
Northeast_model <- read_csv(file = "Model Datasets/NorthEast_Model.csv")

#Midwest
Midwest_final <- read_csv(file = "Final Datasets/MidWest_Final.csv")
Midwest_model <- read_csv(file = "Model Datasets/MidWest_Model.csv")

#Southeast
Southeast_final <- read_csv(file = "Final Datasets/SouthEast_Final.csv")
Southeast_model <- read_csv(file = "Model Datasets/SouthEast_Model.csv")

#Southwest
Southwest_final <- read_csv(file = "Final Datasets/SouthWest_Final.csv")
Southwest_model <- read_csv(file = "Model Datasets/SouthWest_Model.csv")

#West
West_final <- read_csv(file = "Final Datasets/West_Final.csv")
West_model <- read_csv(file = "Model Datasets/West_Model.csv")

#### TESTING FOR ZERO VARIANCE VARIABLES ####
#only using the training part of the data set as this is what will be used to train the models
Northeast_train <- Northeast_final %>% 
  filter(as.Date(date) < "2015-01-01")
Midwest_train <- Midwest_final %>% 
  filter(as.Date(date) < "2015-01-01")
Southeast_train <- Southeast_final %>% 
  filter(as.Date(date) < "2015-01-01")
Southwest_train <- Southwest_final %>% 
  filter(as.Date(date) < "2015-01-01")
West_train <- West_final %>% 
  filter(as.Date(date) < "2015-01-01")

Northeast_test <- anti_join(Northeast_model, Northeast_train, by = "date")
Midwest_test <- anti_join(Midwest_model, Midwest_train, by = "date")
Southeast_test <- anti_join(Southeast_model, Southeast_train, by = "date")
Southwest_test <- anti_join(Southwest_model, Southwest_train, by = "date")
West_test <- anti_join(West_model, West_train, by = "date")

nearZeroVar(Northeast_train, saveMetrics = T)
nearZeroVar(Midwest_train, saveMetrics = T)
nearZeroVar(Southeast_train, saveMetrics = T)
nearZeroVar(Southwest_train, saveMetrics = T)
nearZeroVar(West_train, saveMetrics = T)



#### MODELING WITH ENVIRONMENT TERMS ####
eval_ctrl <- trainControl(method = "timeslice",
                          initialWindow = 12,   # 12 months (1 year) of training
                          horizon = 1,         # 1 month of data to validate. Use one year worth of data to predict the current month.        
                          fixedWindow = TRUE,  # Set this to TRUE for a rolling window         
                          savePredictions = TRUE,
                          allowParallel = TRUE)

# Decided to do a rolling window approach with a shorter training window to help prevent over fitting
# this approach will allow for multiple validation periods and thus may help prevent over fitting
#NORTHEAST ELASTIC NET MODEL 1 
set.seed(1234)
Northeast_els1_full <- train(Lyme_Rate ~ .,                             
                        data = Northeast_train %>% select(-date, -`bible school songs`, -`vbs games`),    
                        method = "glmnet",                        
                        preProc = c("center", "scale"),
                        tuneLength = 250,
                        trControl = eval_ctrl)
save(Northeast_els1_full, file = "Model RDAs/Full Search term Models/Northeast ElsNet 1 Full.rda")


# NORTHEAST ELASTIC NET MODEL 2
set.seed(1234)
Northeast_els2_full <- train(Lyme_Rate ~ .,                             
                        data = Northeast_train %>% select(-date, -`bible school songs`, -`vbs games`),    
                        method = "glmnet",                        
                        preProc = c("center", "scale"),
                        tuneLength = 150,
                        trControl = eval_ctrl)
save(Northeast_els2_full, file = "Model RDAs/Full Search term Models/Northeast ElsNet 2 Full.rda")

# MIDWEST ELASTIC NET MODEL 1
set.seed(1234)
Midwest_els1_full <- train(Lyme_Rate ~ .,                             
                      data = Midwest_train %>% select(-date, -`river float trips`),    
                      method = "glmnet",                        
                      preProc = c("center", "scale"),
                      tuneLength = 50,
                      trControl = eval_ctrl)
save(Midwest_els1_full, file = "Model RDAs/Full Search term Models/Midwest ElsNet 1 Full.rda")

# MIDWEST ELASTIC NET MODEL 2
set.seed(1234)
Midwest_els2_full <- train(Lyme_Rate ~ .,                             
                      data = Midwest_train %>% select(-date, -`river float trips`),    
                      method = "glmnet",                        
                      preProc = c("center", "scale"),
                      tuneLength = 150,
                      trControl = eval_ctrl)
save(Midwest_els2_full, file = "Model RDAs/Full Search term Models/Midwest ElsNet 2 Full.rda")

# SOUTHEAST ELASTIC NET MODEL 1
set.seed(1234)
Southeast_els1_full <- train(Lyme_Rate ~ .,                             
                        data = Southeast_train %>% select(-date, -`alive at five`, -`asheboro copperheads`, -`danville dans`, 
                                                              -`diamond devils`, -`fayetteville swampdogs`,-`free summer kids movies`, 
                                                              -`herndon braves`, -`palmetto falls`, -`wilson tobs`, -`wood bat tournament`),    
                        method = "glmnet",                        
                        preProc = c("center", "scale"),
                        tuneLength = 50,
                        trControl = eval_ctrl)
save(Southeast_els1_full, file = "Model RDAs/Full Search term Models/SouthEast ElsNet 1 Full.rda")

# SOUTHEAST ELASTIC NET MODEL 2
set.seed(1234)
Southeast_els2_full <- train(Lyme_Rate ~ .,                             
                        data = Southeast_train %>% select(-date, -`alive at five`, -`asheboro copperheads`, -`danville dans`, 
                                                              -`diamond devils`, -`fayetteville swampdogs`,-`free summer kids movies`, 
                                                              -`herndon braves`, -`palmetto falls`, -`wilson tobs`, -`wood bat tournament`),    
                        method = "glmnet",                        
                        preProc = c("center", "scale"),
                        tuneLength = 150,
                        trControl = eval_ctrl)
save(Southeast_els2_full, file = "Model RDAs/Full Search term Models/SouthEast ElsNet 2 Full.rda") 


# SOUTHWEST ELASTIC NET MODEL 1
set.seed(1234)
Southwest_els1_full <- train(Lyme_Rate ~ .,
                        data = Southwest_train %>% select(-date),
                        method = "glmnet",
                        preProc = c("center", "scale"),
                        tuneLength = 50,
                        trControl = eval_ctrl)
save(Southwest_els1_full, file = "Model RDAs/Full Search term Models/SouthWest ElsNet 1 Full.rda")

# SOUTHWEST ELASTIC NET MODEL 2
set.seed(1234)
Southwest_els2_full <- train(Lyme_Rate ~ .,
                        data = Southwest_train %>% select(-date),
                        method = "glmnet",
                        preProc = c("center", "scale"),
                        tuneLength = 150,
                        trControl = eval_ctrl)
save(Southwest_els2_full, file = "Model RDAs/Full Search term Models/SouthWest ElsNet 2 Full.rda")

# WEST ELASTIC NET MODEL 1
set.seed(1234)
West_els1_full <- train(Lyme_Rate ~ .,
                   data = West_train %>% select(-date),
                   method = "glmnet",
                   preProc = c("center", "scale"),
                   tuneLength = 50,
                   trControl = eval_ctrl)
save(West_els1_full, file = "Model RDAs/Full Search term Models/West ElsNet 1 Full.rda")

# WEST ELASTIC NET MODEL 2
set.seed(1234)
West_els2_full <- train(Lyme_Rate ~ .,
                   data = West_train %>% select(-date),
                   method = "glmnet",
                   preProc = c("center", "scale"),
                   tuneLength = 150,
                   trControl = eval_ctrl)
save(West_els2_full, file = "Model RDAs/Full Search term Models/West ElsNet 2 Full.rda")


#### SELECTING VECTOR AND DISEASE TERMS ONLY ####
VSD_Northeast_model <- Northeast_model %>% 
                          select(date, Lyme_Rate, tick, `black tick`, lyme, `lyme disease`, rash, `bullseye rash`,
                                 `bell's palsy`, `facial paralysis`, `knee pain`, `swollen knee`,
                                 `swollen knees`, `swollen joint`, `swollen joints`, `joint pain`, fever, tired, `deer tick`,
                                 `black legged tick`)

VSD_Midwest_model<- Midwest_model %>% 
                       select(date, Lyme_Rate, tick, `black tick`, lyme, `lyme disease`, rash, `bullseye rash`,
                              `bell's palsy`, `facial paralysis`, `knee pain`, `swollen knee`,
                              `swollen knees`, `swollen joint`, `swollen joints`, `joint pain`, fever, tired, `deer tick`,
                              `black legged tick`)

VSD_Southeast_model<- Southeast_model %>% 
                       select(date, Lyme_Rate, tick, `black tick`, lyme, `lyme disease`, rash, `bullseye rash`,
                              `bell's palsy`, `facial paralysis`, `knee pain`, `swollen knee`,
                              `swollen knees`, `swollen joint`, `swollen joints`, `joint pain`, fever, tired, `deer tick`,
                              `black legged tick`)

VSD_Southwest_model<- Southwest_model %>% 
                        select(date, Lyme_Rate, tick, `black tick`, lyme, `lyme disease`, rash, `bullseye rash`,
                               `bell's palsy`, `facial paralysis`, `knee pain`, `swollen knee`,
                               `swollen knees`, `swollen joint`, `swollen joints`, `joint pain`, fever, tired, `deer tick`,
                               `black legged tick`)

VSD_West_model<- West_model %>% 
                   select(date, Lyme_Rate, tick, `black tick`, lyme, `lyme disease`, rash, `bullseye rash`,
                          `bell's palsy`, `facial paralysis`, `knee pain`, `swollen knee`,
                          `swollen knees`, `swollen joint`, `swollen joints`, `joint pain`, fever, tired, `deer tick`,
                          `black legged tick`)

Northeast_train_VSD <- VSD_Northeast_model %>% 
  filter(as.Date(date) < "2015-01-01")
Midwest_train_VSD <- VSD_Midwest_model %>% 
  filter(as.Date(date) < "2015-01-01")
Southeast_train_VSD <- VSD_Southeast_model %>% 
  filter(as.Date(date) < "2015-01-01")
Southwest_train_VSD <- VSD_Southwest_model %>% 
  filter(as.Date(date) < "2015-01-01")
West_train_VSD <- VSD_West_model %>% 
  filter(as.Date(date) < "2015-01-01")

Northeast_test_VSD <- anti_join(VSD_Northeast_model, Northeast_train_VSD, by = "date")
Midwest_test_VSD <- anti_join(VSD_Midwest_model, Midwest_train_VSD, by = "date")
Southeast_test_VSD <- anti_join(VSD_Southeast_model, Southeast_train_VSD, by = "date")
Southwest_test_VSD <- anti_join(VSD_Southwest_model, Southwest_train_VSD, by = "date")
West_test_VSD <- anti_join(VSD_West_model, West_train_VSD, by = "date")

#### MODELING WITH VECTOR AND DISEASE TERMS ONLY  ####

#NORTHEAST ELASTIC NET MODEL 1 
eval_ctrl <- trainControl(method = "timeslice",
                          initialWindow = 12,   # 12 months (1 year) of training
                          horizon = 1,         # 1 month of data to validate. Use one year worth of data to predict the current month.        
                          fixedWindow = TRUE,  # Set this to TRUE for a rolling window         
                          savePredictions = TRUE,
                          allowParallel = TRUE)

# This is a change in approached compared to the BMC submission
# Decided to do a rolling window approach with a shorter training window to help prevent over fitting
# this approach will allow for multiple validation periods and thus may help prevent over fitting
#NORTHEAST ELASTIC NET MODEL 1 
set.seed(1234)
Northeast_els1_VSD <- train(Lyme_Rate ~ .,                             
                             data = Northeast_train_VSD %>% select(-date),    
                             method = "glmnet",                        
                             preProc = c("center", "scale"),
                             tuneLength = 50,
                             trControl = eval_ctrl)
save(Northeast_els1_VSD, file = "Model RDAs/Disease And Vector term only Models/Northeast ElsNet 1 VSD.rda")


# NORTHEAST ELASTIC NET MODEL 2
set.seed(1234)
Northeast_els2_VSD <- train(Lyme_Rate ~ .,                             
                             data = Northeast_train_VSD %>% select(-date),    
                             method = "glmnet",                        
                             preProc = c("center", "scale"),
                             tuneLength = 150,
                             trControl = eval_ctrl)
save(Northeast_els2_VSD, file = "Model RDAs/Disease And Vector term only Models/Northeast ElsNet 2 VSD.rda")

# MIDWEST ELASTIC NET MODEL 1
set.seed(1234)
Midwest_els1_VSD <- train(Lyme_Rate ~ .,                             
                           data = Midwest_train_VSD %>% select(-date),    
                           method = "glmnet",                        
                           preProc = c("center", "scale"),
                           tuneLength = 50,
                           trControl = eval_ctrl)
save(Midwest_els1_VSD, file = "Model RDAs/Disease And Vector term only Models/Midwest ElsNet 1 VSD.rda")

# MIDWEST ELASTIC NET MODEL 2
set.seed(1234)
Midwest_els2_VSD <- train(Lyme_Rate ~ .,                             
                           data = Midwest_train_VSD %>% select(-date),    
                           method = "glmnet",                        
                           preProc = c("center", "scale"),
                           tuneLength = 150,
                           trControl = eval_ctrl)
save(Midwest_els2_VSD, file = "Model RDAs/Disease And Vector term only Models/Midwest ElsNet 2 VSD.rda")

# SOUTHEAST ELASTIC NET MODEL 1
set.seed(1234)
Southeast_els1_VSD <- train(Lyme_Rate ~ .,                             
                             data = Southeast_train_VSD %>% select(-date),    
                             method = "glmnet",                        
                             preProc = c("center", "scale"),
                             tuneLength = 50,
                             trControl = eval_ctrl)
save(Southeast_els1_VSD, file = "Model RDAs/Disease And Vector term only Models/SouthEast ElsNet 1 VSD.rda")

# SOUTHEAST ELASTIC NET MODEL 2
set.seed(1234)
Southeast_els2_VSD <- train(Lyme_Rate ~ .,                             
                             data = Southeast_train_VSD %>% select(-date),    
                             method = "glmnet",                        
                             preProc = c("center", "scale"),
                             tuneLength = 150,
                             trControl = eval_ctrl)
save(Southeast_els2_VSD, file = "Model RDAs/Disease And Vector term only Models/SouthEast ElsNet 2 VSD.rda") 


# SOUTHWEST ELASTIC NET MODEL 1
set.seed(1234)
Southwest_els1_VSD <- train(Lyme_Rate ~ .,
                             data = Southwest_train_VSD %>% select(-date),
                             method = "glmnet",
                             preProc = c("center", "scale"),
                             tuneLength = 50,
                             trControl = eval_ctrl)
save(Southwest_els1_VSD, file = "Model RDAs/Disease And Vector term only Models/SouthWest ElsNet 1 VSD.rda")

# SOUTHWEST ELASTIC NET MODEL 2
set.seed(1234)
Southwest_els2_VSD <- train(Lyme_Rate ~ .,
                             data = Southwest_train_VSD %>% select(-date),
                             method = "glmnet",
                             preProc = c("center", "scale"),
                             tuneLength = 150,
                             trControl = eval_ctrl)
save(Southwest_els2_VSD, file = "Model RDAs/Disease And Vector term only Models/SouthWest ElsNet 2 VSD.rda")

# WEST ELASTIC NET MODEL 1
set.seed(1234)
West_els1_VSD <- train(Lyme_Rate ~ .,
                        data = West_train_VSD %>% select(-date),
                        method = "glmnet",
                        preProc = c("center", "scale"),
                        tuneLength = 50,
                        trControl = eval_ctrl)
save(West_els1_VSD, file = "Model RDAs/Disease And Vector term only Models/West ElsNet 1 VSD.rda")

# WEST ELASTIC NET MODEL 2
set.seed(1234)
West_els2_VSD <- train(Lyme_Rate ~ .,
                        data = West_train_VSD %>% select(-date),
                        method = "glmnet",
                        preProc = c("center", "scale"),
                        tuneLength = 150,
                        trControl = eval_ctrl)
save(West_els2_VSD, file = "Model RDAs/Disease And Vector term only Models/West ElsNet 2 VSD.rda")

#### ENV MODEL EVALUATION ####
load(file = "Model RDAs/Full Search term Models/Northeast ElsNet 1 Full.rda")
load(file = "Model RDAs/Full Search term Models/Northeast ElsNet 2 Full.rda")
load(file = "Model RDAs/Full Search term Models/Midwest ElsNet 1 Full.rda")
load(file = "Model RDAs/Full Search term Models/Midwest ElsNet 2 Full.rda")
load(file = "Model RDAs/Full Search term Models/Southeast ElsNet 1 Full.rda")
load(file = "Model RDAs/Full Search term Models/Southeast ElsNet 2 Full.rda")
load(file = "Model RDAs/Full Search term Models/Southwest ElsNet 1 Full.rda")
load(file = "Model RDAs/Full Search term Models/Southwest ElsNet 2 Full.rda")
load(file = "Model RDAs/Full Search term Models/West ElsNet 1 Full.rda")
load(file = "Model RDAs/Full Search term Models/West ElsNet 2 Full.rda")
##Northeast
#in-sample data prediction error and accuracy
postResample(pred = predict(Northeast_els1_full, newdata = Northeast_els1_full$trainingData),
             obs = Northeast_els1_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Northeast_els1_full$pred %>% 
               filter(lambda == Northeast_els1_full$bestTune$lambda,
                      alpha == Northeast_els1_full$bestTune$alpha) %>% select(pred), 
             obs = Northeast_els1_full$pred %>% 
               filter(lambda == Northeast_els1_full$bestTune$lambda,
                      alpha == Northeast_els1_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Northeast_els1_full, newdata = Northeast_test),
             obs = Northeast_test$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(Northeast_els2_full, newdata = Northeast_els2_full$trainingData),
             obs = Northeast_els2_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Northeast_els2_full$pred %>% 
               filter(lambda == Northeast_els2_full$bestTune$lambda,
                      alpha == Northeast_els2_full$bestTune$alpha) %>% select(pred), 
             obs = Northeast_els2_full$pred %>% 
               filter(lambda == Northeast_els2_full$bestTune$lambda,
                      alpha == Northeast_els2_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Northeast_els2_full, newdata = Northeast_test),
             obs = Northeast_test$Lyme_Rate)

##Midwest
#in-sample data prediction error and accuracy
postResample(pred = predict(Midwest_els1_full, newdata = Midwest_els1_full$trainingData),
             obs = Midwest_els1_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Midwest_els1_full$pred %>% 
               filter(lambda == Midwest_els1_full$bestTune$lambda,
                      alpha == Midwest_els1_full$bestTune$alpha) %>% select(pred), 
             obs = Midwest_els1_full$pred %>% 
               filter(lambda == Midwest_els1_full$bestTune$lambda,
                      alpha == Midwest_els1_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Midwest_els1_full, newdata = Midwest_test),
             obs = Midwest_test$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(Midwest_els2_full, newdata = Midwest_els2_full$trainingData),
             obs = Midwest_els2_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Midwest_els2_full$pred %>% 
               filter(lambda == Midwest_els2_full$bestTune$lambda,
                      alpha == Midwest_els2_full$bestTune$alpha) %>% select(pred), 
             obs = Midwest_els2_full$pred %>% 
               filter(lambda == Midwest_els2_full$bestTune$lambda,
                      alpha == Midwest_els2_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Midwest_els2_full, newdata = Midwest_test),
             obs = Midwest_test$Lyme_Rate)

##Southeast
#in-sample data prediction error and accuracy
postResample(pred = predict(Southeast_els1_full, newdata = Southeast_els1_full$trainingData),
             obs = Southeast_els1_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Southeast_els1_full$pred %>% 
               filter(lambda == Southeast_els1_full$bestTune$lambda,
                      alpha == Southeast_els1_full$bestTune$alpha) %>% select(pred), 
             obs = Southeast_els1_full$pred %>% 
               filter(lambda == Southeast_els1_full$bestTune$lambda,
                      alpha == Southeast_els1_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Southeast_els1_full, newdata = Southeast_test),
             obs = Southeast_test$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(Southeast_els2_full, newdata = Southeast_els2_full$trainingData),
             obs = Southeast_els2_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Southeast_els2_full$pred %>% 
               filter(lambda == Southeast_els2_full$bestTune$lambda,
                      alpha == Southeast_els2_full$bestTune$alpha) %>% select(pred), 
             obs = Southeast_els2_full$pred %>% 
               filter(lambda == Southeast_els2_full$bestTune$lambda,
                      alpha == Southeast_els2_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Southeast_els2_full, newdata = Southeast_test),
             obs = Southeast_test$Lyme_Rate)

##Southwest
#in-sample data prediction error and accuracy
postResample(pred = predict(Southwest_els1_full, newdata = Southwest_els1_full$trainingData),
             obs = Southwest_els1_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Southwest_els1_full$pred %>% 
               filter(lambda == Southwest_els1_full$bestTune$lambda,
                      alpha == Southwest_els1_full$bestTune$alpha) %>% select(pred), 
             obs = Southwest_els1_full$pred %>% 
               filter(lambda == Southwest_els1_full$bestTune$lambda,
                      alpha == Southwest_els1_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Southwest_els1_full, newdata = Southwest_test),
             obs = Southwest_test$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(Southwest_els2_full, newdata = Southwest_els2_full$trainingData),
             obs = Southwest_els2_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Southwest_els2_full$pred %>% 
               filter(lambda == Southwest_els2_full$bestTune$lambda,
                      alpha == Southwest_els2_full$bestTune$alpha) %>% select(pred), 
             obs = Southwest_els2_full$pred %>% 
               filter(lambda == Southwest_els2_full$bestTune$lambda,
                      alpha == Southwest_els2_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Southwest_els2_full, newdata = Southwest_test),
             obs = Southwest_test$Lyme_Rate)

##West
#in-sample data prediction error and accuracy
postResample(pred = predict(West_els1_full, newdata = West_els1_full$trainingData),
             obs = West_els1_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =West_els1_full$pred %>% 
               filter(lambda == West_els1_full$bestTune$lambda,
                      alpha == West_els1_full$bestTune$alpha) %>% select(pred), 
             obs = West_els1_full$pred %>% 
               filter(lambda == West_els1_full$bestTune$lambda,
                      alpha == West_els1_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(West_els1_full, newdata = West_test),
             obs = West_test$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(West_els2_full, newdata = West_els2_full$trainingData),
             obs = West_els2_full$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =West_els2_full$pred %>% 
               filter(lambda == West_els2_full$bestTune$lambda,
                      alpha == West_els2_full$bestTune$alpha) %>% select(pred), 
             obs = West_els2_full$pred %>% 
               filter(lambda == West_els2_full$bestTune$lambda,
                      alpha == West_els2_full$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(West_els2_full, newdata = West_test),
             obs = West_test$Lyme_Rate)

#### VECTOR, SYMPTOM, and DISEASE MODEL EVALUATION ####
load(file = "Model RDAs/Disease And Vector term only Models/Northeast ElsNet 1 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/Northeast ElsNet 2 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/Midwest ElsNet 1 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/Midwest ElsNet 2 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/Southeast ElsNet 1 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/Southeast ElsNet 2 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/Southwest ElsNet 1 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/Southwest ElsNet 2 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/West ElsNet 1 VSD.rda")
load(file = "Model RDAs/Disease And Vector term only Models/West ElsNet 2 VSD.rda")

##Northeast
#in-sample data prediction error and accuracy
postResample(pred = predict(Northeast_els1_VSD, newdata = Northeast_els1_VSD$trainingData),
             obs = Northeast_els1_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Northeast_els1_VSD$pred %>% 
               filter(lambda == Northeast_els1_VSD$bestTune$lambda,
                      alpha == Northeast_els1_VSD$bestTune$alpha) %>% select(pred), 
             obs = Northeast_els1_VSD$pred %>% 
               filter(lambda == Northeast_els1_VSD$bestTune$lambda,
                      alpha == Northeast_els1_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Northeast_els1_VSD, newdata = Northeast_test_VSD),
             obs = Northeast_test_VSD$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(Northeast_els2_VSD, newdata = Northeast_els2_VSD$trainingData),
             obs = Northeast_els2_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Northeast_els2_VSD$pred %>% 
               filter(lambda == Northeast_els2_VSD$bestTune$lambda,
                      alpha == Northeast_els2_VSD$bestTune$alpha) %>% select(pred), 
             obs = Northeast_els2_VSD$pred %>% 
               filter(lambda == Northeast_els2_VSD$bestTune$lambda,
                      alpha == Northeast_els2_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Northeast_els2_VSD, newdata = Northeast_test_VSD),
             obs = Northeast_test_VSD$Lyme_Rate)

##Midwest
#in-sample data prediction error and accuracy
postResample(pred = predict(Midwest_els1_VSD, newdata = Midwest_els1_VSD$trainingData),
             obs = Midwest_els1_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Midwest_els1_VSD$pred %>% 
               filter(lambda == Midwest_els1_VSD$bestTune$lambda,
                      alpha == Midwest_els1_VSD$bestTune$alpha) %>% select(pred), 
             obs = Midwest_els1_VSD$pred %>% 
               filter(lambda == Midwest_els1_VSD$bestTune$lambda,
                      alpha == Midwest_els1_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Midwest_els1_VSD, newdata = Midwest_test_VSD),
             obs = Midwest_test_VSD$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(Midwest_els2_VSD, newdata = Midwest_els2_VSD$trainingData),
             obs = Midwest_els2_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Midwest_els2_VSD$pred %>% 
               filter(lambda == Midwest_els2_VSD$bestTune$lambda,
                      alpha == Midwest_els2_VSD$bestTune$alpha) %>% select(pred), 
             obs = Midwest_els2_VSD$pred %>% 
               filter(lambda == Midwest_els2_VSD$bestTune$lambda,
                      alpha == Midwest_els2_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Midwest_els2_VSD, newdata = Midwest_test_VSD),
             obs = Midwest_test_VSD$Lyme_Rate)

##Southeast
#in-sample data prediction error and accuracy
postResample(pred = predict(Southeast_els1_VSD, newdata = Southeast_els1_VSD$trainingData),
             obs = Southeast_els1_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Southeast_els1_VSD$pred %>% 
               filter(lambda == Southeast_els1_VSD$bestTune$lambda,
                      alpha == Southeast_els1_VSD$bestTune$alpha) %>% select(pred), 
             obs = Southeast_els1_VSD$pred %>% 
               filter(lambda == Southeast_els1_VSD$bestTune$lambda,
                      alpha == Southeast_els1_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Southeast_els1_VSD, newdata = Southeast_test_VSD),
             obs = Southeast_test_VSD$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(Southeast_els2_VSD, newdata = Southeast_els2_VSD$trainingData),
             obs = Southeast_els2_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Southeast_els2_VSD$pred %>% 
               filter(lambda == Southeast_els2_VSD$bestTune$lambda,
                      alpha == Southeast_els2_VSD$bestTune$alpha) %>% select(pred), 
             obs = Southeast_els2_VSD$pred %>% 
               filter(lambda == Southeast_els2_VSD$bestTune$lambda,
                      alpha == Southeast_els2_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Southeast_els2_VSD, newdata = Southeast_test_VSD),
             obs = Southeast_test_VSD$Lyme_Rate)

##Southwest
#in-sample data prediction error and accuracy
postResample(pred = predict(Southwest_els1_VSD, newdata = Southwest_els1_VSD$trainingData),
             obs = Southwest_els1_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Southwest_els1_VSD$pred %>% 
               filter(lambda == Southwest_els1_VSD$bestTune$lambda,
                      alpha == Southwest_els1_VSD$bestTune$alpha) %>% select(pred), 
             obs = Southwest_els1_VSD$pred %>% 
               filter(lambda == Southwest_els1_VSD$bestTune$lambda,
                      alpha == Southwest_els1_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Southwest_els1_VSD, newdata = Southwest_test_VSD),
             obs = Southwest_test_VSD$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(Southwest_els2_VSD, newdata = Southwest_els2_VSD$trainingData),
             obs = Southwest_els2_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =Southwest_els2_VSD$pred %>% 
               filter(lambda == Southwest_els2_VSD$bestTune$lambda,
                      alpha == Southwest_els2_VSD$bestTune$alpha) %>% select(pred), 
             obs = Southwest_els2_VSD$pred %>% 
               filter(lambda == Southwest_els2_VSD$bestTune$lambda,
                      alpha == Southwest_els2_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(Southwest_els2_VSD, newdata = Southwest_test_VSD),
             obs = Southwest_test_VSD$Lyme_Rate)

##West
#in-sample data prediction error and accuracy
postResample(pred = predict(West_els1_VSD, newdata = West_els1_VSD$trainingData),
             obs = West_els1_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =West_els1_VSD$pred %>% 
               filter(lambda == West_els1_VSD$bestTune$lambda,
                      alpha == West_els1_VSD$bestTune$alpha) %>% select(pred), 
             obs = West_els1_VSD$pred %>% 
               filter(lambda == West_els1_VSD$bestTune$lambda,
                      alpha == West_els1_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(West_els1_VSD, newdata = West_test_VSD),
             obs = West_test_VSD$Lyme_Rate)

#in-sample data prediction error and accuracy
postResample(pred = predict(West_els2_VSD, newdata = West_els2_VSD$trainingData),
             obs = West_els2_VSD$trainingData$.outcome)

#Validation 1 month prediction error and accuracy
postResample(pred =West_els2_VSD$pred %>% 
               filter(lambda == West_els2_VSD$bestTune$lambda,
                      alpha == West_els2_VSD$bestTune$alpha) %>% select(pred), 
             obs = West_els2_VSD$pred %>% 
               filter(lambda == West_els2_VSD$bestTune$lambda,
                      alpha == West_els2_VSD$bestTune$alpha) %>% select(obs)) 
#True out of sample prediction error and accuracy
postResample(pred = predict(West_els2_VSD, newdata = West_test_VSD),
             obs = West_test_VSD$Lyme_Rate)

#### SEARCH TERM CORRELATION ####
NE__cor <- Northeast_model %>% 
  filter(as.Date(date) < "2015-01-01") %>%  
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(NE__cor$value, na.rm = T)
mean(NE__cor$value, na.rm = T)
median(NE__cor$value, na.rm = T)

MW__cor <- Midwest_model %>% 
  filter(as.Date(date) < "2015-01-01") %>%  
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(MW__cor$value, na.rm = T)
mean(MW__cor$value, na.rm = T)
median(MW__cor$value, na.rm = T)

SE__cor <- Southeast_model %>% 
  filter(as.Date(date) < "2015-01-01") %>%  
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(SE__cor$value, na.rm = T)
mean(SE__cor$value, na.rm = T)
median(SE__cor$value, na.rm = T)

SW__cor <- Southwest_model %>% 
  filter(as.Date(date) < "2015-01-01") %>%  
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(SW__cor$value, na.rm = T)
mean(SW__cor$value, na.rm = T)
median(SW__cor$value, na.rm = T)

W__cor <- West_model %>% 
  filter(as.Date(date) < "2015-01-01") %>%  
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(W__cor$value, na.rm = T)
mean(W__cor$value, na.rm = T)
median(W__cor$value, na.rm = T)

#### FULL MODEL VARIALBE IMPORTANCE ####
varimp <- function(data){
  varImp(data, scale = FALSE)$importance %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    arrange(desc(Overall)) %>% 
    mutate(Overall = as.character(Overall))
}
## Northeast
varimp(Northeast_els1_full)
varimp(Northeast_els2_full)

## Midwest
varimp(Midwest_els1_full)
varimp(Midwest_els2_full)

## Southeast
varimp(Southeast_els1_full)
varimp(Southeast_els2_full)

## Southwest
varimp(Southwest_els1_full)
varimp(Southwest_els2_full)

## West
varimp(West_els1_full)
varimp(West_els2_full)

#### VSD MODEL VARIABLE IMPORTANCE ####
## Northeast
varimp(Northeast_els1_VSD)
varimp(Northeast_els2_VSD)

## Midwest
varimp(Midwest_els1_VSD)
varimp(Midwest_els2_VSD)

## Southeast
varimp(Southeast_els1_VSD)
varimp(Southeast_els2_VSD)

## Southwest
varimp(Southwest_els1_VSD)
varimp(Southwest_els2_VSD)

## West
varimp(West_els1_VSD)
varimp(West_els2_VSD)


#### SEARCH TERM CORRELATION ENTIRE DATA SET####
# comment out the filter statement in order to get full data set. 
NE__cor <- Northeast_model %>% 
  filter(date<="2014-12-01") %>%
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(NE__cor$value, na.rm = T)
mean(NE__cor$value, na.rm = T)
median(NE__cor$value, na.rm = T)
# write_csv(NE__cor, path = "../Draft 12/Supplemental Tables/S2 Prep Tables/NE__Cor.csv")

MW__cor <- Midwest_model %>% 
  filter(date<="2014-12-01") %>%
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(MW__cor$value, na.rm = T)
mean(MW__cor$value, na.rm = T)
median(MW__cor$value, na.rm = T)
write_csv(MW__cor, path = "../Draft 12/Supplemental Tables/S2 Prep Tables/MW__cor.csv")

SE__cor <- Southeast_model %>% 
  filter(date<="2014-12-01") %>%
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(SE__cor$value, na.rm = T)
mean(SE__cor$value, na.rm = T)
median(SE__cor$value, na.rm = T)
write_csv(SE__cor, path = "../Draft 12/Supplemental Tables/S2 Prep Tables/SE__cor.csv")


SW__cor <- Southwest_model %>% 
  filter(date<="2014-12-01") %>%
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(SW__cor$value, na.rm = T)
mean(SW__cor$value, na.rm = T)
median(SW__cor$value, na.rm = T)
write_csv(SW__cor, path = "../Draft 12/Supplemental Tables/S2 Prep Tables/SW__cor.csv")

W__cor <- West_model %>% 
  filter(date<="2014-12-01") %>%
  summarise_at(vars(-date, -Lyme_Rate), funs(cor(Lyme_Rate,.))) %>% 
  gather() %>% 
  arrange(desc(value))

range(W__cor$value, na.rm = T)
mean(W__cor$value, na.rm = T)
median(W__cor$value, na.rm = T)
write_csv(W__cor, path = "../Draft 12/Supplemental Tables/S2 Prep Tables/W__cor.csv")

#### FIGURES ####

#### REGIONAL RATES GRAPHS ####
#Figure 2 (Regional Rates of Lyme disease)
Region_rates <- read_csv("Region_Rates.csv")
Rate_Bar <-
  Region_rates %>% 
  mutate(year = lubridate::year(date),
         Region = factor(Region,
                         levels = c("Northeast","Midwest", "Southeast", "Southwest", "West")))%>% 
  group_by(Region, year) %>% 
  summarize(total_lyme_count = sum(total_lyme)) %>% 
  mutate(year = as.factor(year)) %>% 
  ggplot(aes(x = year, y = total_lyme_count, fill = Region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("firebrick3", "slateblue4", "slateblue3", "dodgerblue", "dodgerblue3"))+
  ylab("Total Lyme Count") + 
  theme_minimal()+
  theme(#legend.position = "None",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white",
                                     color = "white"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(color = "black", size = .75),
    axis.ticks = element_line(color = "black", size = .75),
    strip.text.x = element_text(size = 13, face = "bold", color = "black"))+
  coord_cartesian(ylim = c(0,21500), expand = FALSE)

Rate_Bar_legend <- get_legend(Rate_Bar)
plot(Rate_Bar_legend )
#162 x 164

#Figure 3A
Region_rates %>%   
  mutate(Region = factor(Region,
                         levels = c("Northeast","Midwest", "Southeast", "Southwest", "West"))) %>% 
  ggplot(aes(x = date, y = Lyme_Rate, color = Region)) + 
  geom_line(size = 1.5) + 
  facet_wrap(~Region, scales = "free_x") + 
  ylab("Lyme Incidence per 100,000") + 
  theme_classic()+ theme(legend.position = "none")+
  scale_x_date(date_labels = "%Y" , date_breaks = '3 years') +
  scale_color_manual(values = c("firebrick3", "slateblue4", "slateblue3", "dodgerblue", "dodgerblue3")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = .75),
        axis.ticks = element_line(color = "black", size = .75),
        strip.text.x = element_text(size = 13, face = "bold", color = "black"))
#911 x 597

#Figure 3B
Region_rates %>% 
  group_by(Region) %>% 
  filter(Region == "Southwest") %>% 
  ggplot(aes(x = date, y = Lyme_Rate)) + 
  geom_line(size = 1.5, color = "dodgerblue") +
  ylab("Lyme Incidence per 100,000") +
  theme_classic() +
  scale_x_date(date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 14, face = "bold", color = "black"))
#619 x 501

#Figure 3C
Region_rates %>% 
  group_by(Region) %>% 
  filter(Region == "West") %>% 
  ggplot(aes(x = date, y = Lyme_Rate)) + 
  geom_line(size = 1.5, color = "dodgerblue3") +
  ylab("Lyme Incidence per 100,000") +
  theme_classic() +
  scale_x_date(date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 14, face = "bold", color = "black"))
#619 x 501

#### FEVER PLOTS ####
# Figure 4 (Fever often negatively correlated with Lyme disease)
feverMW <- Midwest_model %>% 
  select(date, Lyme_Rate, `fever`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  mutate(Region = rep("Midwest", length(search_term)))

feverNE <- Northeast_model %>% 
  select(date, Lyme_Rate, `fever`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  mutate(Region = rep("Northeast", length(search_term)))

feverSE <- Southeast_model %>% 
  select(date, Lyme_Rate, `fever`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  mutate(Region = rep("Southeast", length(search_term)))

feverSW <- Southwest_model %>% 
  select(date, Lyme_Rate, `fever`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  mutate(Region = rep("Southwest", length(search_term)))

feverW <- West_model %>% 
  select(date, Lyme_Rate, `fever`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  mutate(Region = rep("West", length(search_term)))

common_fever <- rbind(feverMW, feverNE, feverSE, feverSW, feverW)

common_fever %>% 
  mutate(Region = factor(Region,
                         levels = c("Northeast", "Midwest", "Southeast", "Southwest", "West"))) %>% 
  ggplot(aes(x=hits, y=Lyme_Rate, color=Region))+
  geom_point(aes(fill=Region),pch=21, color="black", size=2)+ geom_smooth(method = "lm", se = F, size = 1.5, color = "black")+
  facet_wrap(~Region, scales = "free")+
  scale_fill_manual(values = c("firebrick3", "slateblue4", "slateblue3", "dodgerblue", "dodgerblue3"))+
  theme_minimal() +
  ylab("Lyme Incidence Rate per 100,000") + 
  xlab("Propotional Google Hits") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 14, face = "bold", color = "black")) +
  scale_x_continuous(limits = c(0,100), breaks = c(0, 25, 50, 75, 100))
#835 x 623

#### GRAPHS OF MODELS ####
#Figure 5 and 6
#  c("firebrick3", "slateblue4", "slateblue3", "dodgerblue", "dodgerblue3")
#Northeast
Northeast_final %>% 
  mutate(els1 = predict(Northeast_els1_full, newdata = .),
         els2 = predict(Northeast_els2_full, newdata = .),
         VSD_els1 = predict(Northeast_els1_VSD, newdata = .),
         VSD_els2 = predict(Northeast_els2_VSD, newdata = .)) %>%
  ggplot(aes(x = date)) +
  # geom_line(aes(y = els1), color = "red2", size = 2) +
  # geom_line(aes(y = els2), color = "firebrick3", size = 2) +
  # geom_line(aes(y = VSD_els1), color = "red2", size = 2) +
  geom_line(aes(y = VSD_els2), color = "firebrick3", size = 2) +
  geom_line(aes(y = Lyme_Rate), size = 2.5, alpha = .55) +
  geom_vline(aes(xintercept = as.numeric(date[132])), color = "black", linetype = "twodash", size = 1.5)+ #132 is the last day of the training period data set
  # ylab("Lyme Incidence per 100,000")+
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(-.8, 12.05), breaks = c(0, 4, 8, 12))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        # axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"))
#619 x 501

#Midwest
Midwest_final %>% 
  mutate(els1 = predict(Midwest_els1_full, newdata = .),
         els2 = predict(Midwest_els2_full, newdata = .),
         VSD_els1 = predict(Midwest_els1_VSD, newdata = .),
         VSD_els2 = predict(Midwest_els2_VSD, newdata = .)) %>% 
  ggplot(aes(x = date)) +
  # geom_line(aes(y = els1), color = "darkorange1", size = 2) +
  geom_line(aes(y = els2), color = "slateblue4", size = 2) +
  # geom_line(aes(y = VSD_els1), color = "darkorange1", size = 2) +
  # geom_line(aes(y = VSD_els2), color = "slateblue4", size = 2) +
  geom_line(aes(y = Lyme_Rate), size = 2.5, alpha = .55) +
  geom_vline(aes(xintercept = as.numeric(date[132])), color = "black", linetype = "twodash", size = 1.5)+
  # ylab("Lyme Incidence per 100,000")+
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(-.15, 3.05))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        # axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"))
#619 x 501

#Southeast
Southeast_final %>% 
  mutate(els1 = predict(Southeast_els1_full, newdata = .),
         els2 = predict(Southeast_els2_full, newdata = .),
         VSD_els1 = predict(Southeast_els1_VSD, newdata = .),
         VSD_els2 = predict(Southeast_els2_VSD, newdata = .)) %>% 
  ggplot(aes(x = date)) +
  # geom_line(aes(y = els1), color = "gold2", size = 2) +
  geom_line(aes(y = els2), color = "slateblue3", size = 2) +
  # geom_line(aes(y = VSD_els1), color = "gold2", size = 2) +
  # geom_line(aes(y = VSD_els2), color = "slateblue3", size = 2) +
  geom_line(aes(y = Lyme_Rate), size = 2.5, alpha = .55) +
  geom_vline(aes(xintercept = as.numeric(date[132])), color = "black", linetype = "twodash", size = 1.5)+
  # ylab("Lyme Incidence per 100,000")+
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(-.03,1), breaks = c(0, .3, .6, .9))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        # axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"))
#619 x 501

#Southwest
Southwest_final %>% 
  mutate(els1 = predict(Southwest_els1_full, newdata = .),
         els2 = predict(Southwest_els2_full, newdata = .),
         VSD_els1 = predict(Southwest_els1_VSD, newdata = .),
         VSD_els2 = predict(Southwest_els2_VSD, newdata = .)) %>% 
  ggplot(aes(x = date)) +
  # geom_line(aes(y = els1), color = "green3", size = 2) +
  geom_line(aes(y = els2), color = "dodgerblue", size = 2) +
  # geom_line(aes(y = VSD_els1), color = "green3", size = 2) +
  # geom_line(aes(y = VSD_els2), color = "dodgerblue", size = 2) +
  geom_line(aes(y = Lyme_Rate), size = 2.5, alpha = .55) +
  geom_vline(aes(xintercept = as.numeric(date[132])), color = "black", linetype = "twodash", size = 1.5)+
  ylab("Lyme Incidence per 100,000")+
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(-.002,0.061), breaks = c(0, .02, .04, .06))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        # axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"))
#619 x 501

#West
West_final %>% 
  mutate(els1 = predict(West_els1_full, newdata = .),
         els2 = predict(West_els2_full, newdata = .),
         VSD_els1 = predict(West_els1_VSD, newdata = .),
         VSD_els2 = predict(West_els2_VSD, newdata = .)) %>% 
  ggplot(aes(x = date)) +
  # geom_line(aes(y = els1), color = "green", size = 2) +
  geom_line(aes(y = els2), color = "dodgerblue3", size = 2) +
  # geom_line(aes(y = VSD_els1), color = "green", size = 2) +
  # geom_line(aes(y = VSD_els2), color = "dodgerblue3", size = 2) +
  geom_line(aes(y = Lyme_Rate), size = 2.5, alpha = .55) +
  geom_vline(aes(xintercept = as.numeric(date[132])), color = "black", linetype = "twodash", size = 1.5)+
  ylab("Lyme Incidence per 100,000")+
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(-.002,.091), breaks = c(0, .03, .06, .09))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        # axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"))
#619 x 501

#### ONE YEAR PLOTS ####
#Figure 7
Northeast_final %>% 
  filter(date > "2017-01-01" & date < "2018-01-01") %>% 
  mutate(els1 = predict(Northeast_els1_full, newdata = .),
         els2 = predict(Northeast_els2_full, newdata = .),
         VSD_els1 = predict(Northeast_els1_VSD, newdata = .), 
         VSD_els2 = predict(Northeast_els2_VSD, newdata = .)) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y=Lyme_Rate), color = "black", size = 2.5)+ 
  # geom_line(aes(y = els1), color = "red2", size = 2) +
  geom_line(aes(y = els2), color = "firebrick3", size = 2) +
  # geom_line(aes(y = VSD_els1), color = "red2", size = 2) +
  # geom_line(aes(y = VSD_els2), color = "firebrick3", size = 2) +
  ylab("Lyme Incidence per 100,000") + 
  xlab(element_blank())+
  theme_classic()+ 
  # scale_y_continuous(expand = c(0,0), limits = c(0, 10.05), breaks = c(0, 5, 10))+
  scale_x_date(date_labels = "%b-%Y" , date_breaks = '2 month')+
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"))
#816 x 616

Midwest_final %>% 
  filter(date > "2016-01-01" & date < "2017-01-01") %>% 
  mutate(els1 = predict(Midwest_els1_full, newdata = .),
         els2 = predict(Midwest_els2_full, newdata = .),
         VSD_els1 = predict(Midwest_els1_VSD, newdata = .),
         VSD_els2 = predict(Midwest_els2_VSD, newdata = .)) %>% 
  ggplot(aes(x=date, y=Lyme_Rate)) + 
  # geom_line(aes(y = els1), color = "darkorange1", size = 2) +
  geom_line(aes(y = els2), color = "slateblue4", size = 2) +
  # geom_line(aes(y = VSD_els1), color = "darkorange1", size = 2) +
  # geom_line(aes(y = VSD_els2), color = "slateblue4", size = 2) +
  geom_line(size = 3, color = "black") +
  ylab("Lyme Incidence per 100,000") + 
  xlab(element_blank())+
  theme_classic()+ 
  # scale_y_continuous(expand = c(0,0), limits = c(0, 10.05), breaks = c(0, 5, 10))+
  scale_x_date(date_labels = "%b-%Y" , date_breaks = '2 month')+
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 18, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 18, face = "bold", color = "black"))
#816 x 616

#### THREE TOP TERMS ####
# Figure 8

# NorthEast_legend <-
Northeast_model %>% 
  select(date, Lyme_Rate, `july calendar`, `fresh cherry pie`, `bullseye rash`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  ungroup() %>% mutate(hits=hits/9) %>% 
  ggplot(aes(x=date)) +
  geom_line(aes(y = Lyme_Rate), size = 4, color = "black", alpha = 1)+
  geom_line(aes(x=date, y=hits, color = search_term, size = search_term, lty = search_term,))+
  scale_color_manual(values = c("firebrick2", "cornflowerblue", "chartreuse3"))+
  scale_size_manual(values = c(1.5,1.5,1.5))+
  scale_linetype_manual(values = c("solid", "solid", "solid"))+
  ylab("Lyme Incidence per 100,000") + labs(color = "Search Term", 
                                            linetype = "Search Term",
                                            size = "Search Term")+
  xlab(element_blank()) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(-.2, 12), breaks = c(0, 4, 8, 12))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        # legend.position = "none",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        # axis.title.y = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=4)))
#992 x 720
NorthEast_leg <- get_legend(NorthEast_legend)
plot(NorthEast_leg)
#316 x 151

# Midwest_legend <-
Midwest_model %>% 
  select(date, Lyme_Rate, `festivals milwaukee`, `lake beaches`, `kings island discount`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  ungroup() %>% mutate(hits=hits/35) %>% 
  ggplot(aes(x=date)) +
  geom_line(aes(y = Lyme_Rate), size = 4, color = "black", alpha = 1)+
  geom_line(aes(x=date, y=hits, color = search_term, size = search_term, lty = search_term,))+
  scale_color_manual(values = c("firebrick2", "cornflowerblue", "chartreuse3"))+
  scale_size_manual(values = c(1.5,1.5,1.5))+
  scale_linetype_manual(values = c("solid", "solid", "solid"))+
  ylab("Lyme Incidence per 100,000") + labs(color = "Search Term", 
                                            linetype = "Search Term",
                                            size = "Search Term")+
  xlab(element_blank()) +
  theme_classic() +
  # scale_y_continuous(expand = c(0,0), limits = c(-.15, 3.2))+
  # scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  scale_y_continuous(expand = c(0,0), limits = c(-.2, 12), breaks = c(0, 4, 8, 12))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        # axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=4)))
#992 x 720
Midwest_leg <- get_legend(Midwest_legend)
plot(Midwest_leg)
#set the legend title to element_blank() and legend text to element_text(size = 25, face = "bold") to get correct legend stuff
#316 x 151

# SouthEast_legend <-
Southeast_model %>% 
  select(date, Lyme_Rate, `intex pool cover`, `rash`, `swampdogs`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  ungroup() %>% mutate(hits=hits/80) %>% 
  ggplot(aes(x=date)) +
  geom_line(aes(y = Lyme_Rate), size = 4, color = "black", alpha = 1)+
  geom_line(aes(x=date, y=hits, color = search_term, size = search_term, lty = search_term,))+
  scale_color_manual(values = c("firebrick2", "cornflowerblue", "chartreuse3"))+
  scale_size_manual(values = c(1.5,1.5,1.5))+
  scale_linetype_manual(values = c("solid", "solid", "solid"))+
  ylab("Lyme Incidence per 100,000") + labs(color = "Search Term", 
                                            linetype = "Search Term",
                                            size = "Search Term")+
  xlab(element_blank()) +
  theme_classic()+
  # scale_y_continuous(expand = c(0,0), limits = c(-.03,1.15), breaks = c(0, .35, .65, .95))+
  # scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  scale_y_continuous(expand = c(0,0), limits = c(-.2, 12), breaks = c(0, 4, 8, 12))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        # axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=4)))
#992 x 720
SouthEast_leg <- get_legend(SouthEast_legend)
plot(SouthEast_leg)
#316 x 151

# SouthWest_legend <-
Southwest_model %>% 
  select(date, Lyme_Rate, `loans for`, `ca water`, `hotels ca`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  ungroup() %>% mutate(hits=hits/1800) %>% 
  ggplot(aes(x=date)) +
  geom_line(aes(y = Lyme_Rate), size = 4, color = "black", alpha = 1)+
  geom_line(aes(x=date, y=hits, color = search_term, size = search_term, lty = search_term,))+
  scale_color_manual(values = c("firebrick2", "cornflowerblue", "chartreuse3"))+
  scale_size_manual(values = c(1.5,1.5,1.5))+
  scale_linetype_manual(values = c("solid", "solid", "solid"))+
  ylab("Lyme Incidence per 100,000") + labs(color = "Search Term", 
                                            linetype = "Search Term",
                                            size = "Search Term")+
  xlab(element_blank()) +
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(-.002,0.061), breaks = c(0, .02, .04, .06))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  # scale_y_continuous(expand = c(0,0), limits = c(-.2, 12), breaks = c(0, 4, 8, 12))+
  # scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        # legend.position = "none",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 18, face = "bold"),
        # axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=4)))
#992 x 720
SouthWest_leg <- get_legend(SouthWest_legend)
plot(SouthWest_leg)
#316 x 151

# West_legend <-
West_model %>% 
  select(date, Lyme_Rate, `movies in the park`, `concert in the park`, `waterworld denver`) %>% 
  gather(search_term, hits, -date, -Lyme_Rate) %>% 
  ungroup() %>% mutate(hits=hits/900) %>% 
  ggplot(aes(x=date)) +
  geom_line(aes(y = Lyme_Rate), size = 4, color = "black", alpha = 1)+
  geom_line(aes(x=date, y=hits, color = search_term, size = search_term, lty = search_term,))+
  scale_color_manual(values = c("firebrick2", "cornflowerblue", "chartreuse3"))+
  scale_size_manual(values = c(1.5,1.5,1.5))+
  scale_linetype_manual(values = c("solid", "solid", "solid"))+
  ylab("Lyme Incidence per 100,000") + labs(color = "Search Term", 
                                            linetype = "Search Term",
                                            size = "Search Term")+
  xlab(element_blank()) +
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(-.002,.091), breaks = c(0, .03, .06, .09))+
  scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  # scale_y_continuous(expand = c(0,0), limits = c(-.2, 12), breaks = c(0, 4, 8, 12))+
  # scale_x_date(expand = c(0,0), date_labels = "%Y" , date_breaks = '3 years') +
  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        # legend.position = "none",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 18, face = "bold"),
        # axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black"),
        strip.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=4)))
#992 x 720
West_leg <- get_legend(West_legend)
plot(West_leg)
#316 x 151









