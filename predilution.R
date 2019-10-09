
#Ideas: input Sample names and RFU's
#Script will output csv of growth rates, logistic growth and png of plot

# Step 1: Take RFU of cultures and input data into formatted spreadsheet
# Step 2: Run part 1 of the script and use _dilution output to do dilution
# Step 3: Get post dilution RFU's of cultures and put into _postdilution output
# Step 4: Run part 2 of script and calculate growth rate and post dilution avergaes, plot current growth
# Step 4: next time dilution is done, start with _postdilution file and add new rows for the day's data

# TO DO:

# Fix date issues
# 

sample_count <- 6
#Sample names
Sample_treatments <- c("CMA-2-139-1",
                       "CMA-2-139-2",
                       "CMA-2-139-3",
                       "CMA-2-139-4",
                       "CMA-2-139-5",
                       "CMA-2-139-6")

library(dplyr)
library(janitor)
library(ggplot2)
library(matrixStats)

setwd("~/Desktop/FC_Global_Proteome/Methods/FC_GP_culture_dilutions/1_deg")


#Prompt to enter file name
file_name <- 'FC_GP_dilution_1_180719'
RFU_data <- read.csv(paste0(file_name, ".csv"), stringsAsFactors = FALSE)

#sample volume >:(
#Goal: ask and answer
sample_mL <- 3



#Add sample volume to today's entries
RFU_data$vol_sampled[(length(RFU_data$Date)-7):length(RFU_data$Date)] <- sample_mL


#Calculate day number 
RFU_data$Date <- as.Date(RFU_data$Date, origin = "1899-12-30") #Change date to date format 
for (i in 1:length(RFU_data$Date)) { #For every date, calculate day by subtracting date from first date
  RFU_data$Day[i] <- RFU_data$Date[i] - RFU_data$Date[1] 
}


# Get average of RFU triplicate measurements, sd, cv
# Add warning to cv >.05
RFU_data$predil_RFU_avg <- rowMeans( RFU_data[ , c("predil_RFU_1", "predil_RFU_2", "predil_RFU_3")] )
RFU_data$predil_RFU_sd <- apply(RFU_data[, c("predil_RFU_1","predil_RFU_2", "predil_RFU_3")],1,sd)
RFU_data$predil_RFU_cv <- RFU_data$predil_RFU_sd /  RFU_data$predil_RFU_avg

#Check if dilution is needed based on RFU
for (i in 1:length(RFU_data$Date)){
  if (RFU_data$predil_RFU_avg[i] > 60) {
    RFU_data$Dil_needed[i] <- TRUE
  }
  else{ RFU_data$Dil_needed[i] <- FALSE #If not larger than 100, dilution is not needed
  }
  
}

#Ask if dilution will be done. If yes, fill in dilution done = TRUE

#Get average of RFU triplicate measurements after dilution (only if it happened) 
for (i in 1:length(RFU_data$Date)){ #For entire length of table, check if a dilution was done
  if (RFU_data$Dil_done[i] == 'TRUE'){
    RFU_data$postdil_RFU_avg[i] <- mean(c(RFU_data$postdil_RFU_1[i], #If so, get the means of tripliate measurements
                                          RFU_data$postdil_RFU_2[i], 
                                          RFU_data$postdil_RFU_3[i]))
  } else { #If dilution was not done, put NA's in avg column
    RFU_data$postdil_RFU_avg[i] <- RFU_data$predil_RFU_avg[i]
  }
}





#Get average of RFU triplicate measurements after dilution (only if it happened) 
for (i in 1:length(RFU_data$Date)){ #For entire length of table, check if a dilution was done
  if (RFU_data$Dil_done[i] == 'TRUE'){ #if dilution was done, make new presample volume 153
    RFU_data$presample_postdil_vol[i] <- 153
  } 
}
#Calculate standard deviation
RFU_data$postdil_RFU_sd <- apply(RFU_data[, c("postdil_RFU_1","postdil_RFU_2", "postdil_RFU_3")],1,sd)

#Calculate coeffient of variation
RFU_data$postdil_RFU_cv <- RFU_data$postdil_RFU_sd /  RFU_data$postdil_RFU_avg


#Calculate culture volume before the dilution
for (i in (sample_count+1):length(RFU_data$Date)) { 
  if (RFU_data$Dil_done[i] == "FALSE") { #Check if a dolution has been done
    RFU_data$postsample_predil_vol[i] <- RFU_data$postsample_predil_vol[i- sample_count] - RFU_data$vol_sampled[i];
    RFU_data$postsample_postdil_vol[i] <- RFU_data$postsample_predil_vol[i]
  } #If not done, make predil volumne the same as post dil volume
  else {
    RFU_data$presample_postdil_vol[i] <- 153;
    RFU_data$postsample_postdil_vol[i] <- RFU_data$presample_postdil_vol[i] - 3
    RFU_data$postsample_predil_vol[i] <- RFU_data$postsample_postdil_vol[i- sample_count] - RFU_data$vol_sampled[i]
  }
}

#Calculate amount of culture to remove (return RFU to 60)
RFU_data$prop_vol_wanted <- 60/RFU_data$predil_RFU_avg #proportion of volume needed for 60 RFU
RFU_data$vol_keep <- RFU_data$prop_vol_wanted * RFU_data$postsample_predil_vol * (153/RFU_data$postsample_predil_vol)
RFU_data$vol_keep <-  round(RFU_data$vol_keep, digits = 0)

RFU_data$vol_remove <- RFU_data$postsample_predil_vol - RFU_data$vol_keep
RFU_data$vol_remove<- round(RFU_data$vol_remove, digits=0) # Round to no decimal places

#Calculate volume after removal and check that it is equal to desired "keep" volume
RFU_data$postremoval_vol <- RFU_data$postsample_predil_vol - RFU_data$vol_remove
RFU_data$postremoval_vol == RFU_data$vol_keep

#Calculate amount of media to add
RFU_data$media_add <-153- RFU_data$vol_keep # Top back up to 153 mLs to account for postdil sample
RFU_data$media_add <- round(RFU_data$media_add , digits=0)


#------------------------------------------------------

#Table from the current date with removal and addition amounts
#Add checks, cv, u
#Export as new csv
tableout <- select(RFU_data, 
                   Dil_needed, 
                   vol_remove, 
                   media_add)
tableout <- tableout[(length(RFU_data$Date)-5):length(RFU_data$Date), ]
colnames(tableout) <- c('Dilution Needed', 
                        'Volume to Remove', 
                        'Media to Add')
rownames(tableout) <- Sample_treatments

#Export the CSV
write.csv(tableout, 
          file = paste(file_name, ".csv", 
                       sep= "_dilution"))

#-----------------------------------


#Create vector with length of csv
RFU_data$u <- NA

#Calculate growth rate
#If dilution= true, get growth rate from post dil RFU
for (i in (sample_count+1):length(RFU_data$Date)) {
  if (RFU_data$Dil_done[i] == "FALSE"){
    RFU_data$u[i] <- log(RFU_data$predil_RFU_avg[i]/RFU_data$predil_RFU_avg[i-sample_count])/
      (RFU_data$Day[i]-RFU_data$Day[i-sample_count])
  }
  else{
    RFU_data$u[i] <- log(RFU_data$predil_RFU_avg[i]/RFU_data$postdil_RFU_avg[i-sample_count])/
      (RFU_data$Day[i]-RFU_data$Day[i-sample_count])
  }
}

#------------------------------------------------------
#Export the CSV
write.csv(RFU_data, 
          file = paste(file_name, ".csv", 
                       sep= "_postdilution"))

#------------------------------------------------------------
#Calculate growth rate
# Import postdilution data
postdilRFU_data <- read.csv(paste0(file_name, "_postdilution.csv"), stringsAsFactors = FALSE)



#Get average of RFU triplicate measurements after dilution (only if it happened) 
for (i in 1:length(postdilRFU_data$Date)){ #For entire length of table, check if a dilution was done
  if (postdilRFU_data$Dil_done[i] == 'TRUE'){ #if dilution was done, make new presample volume 153
    postdilRFU_data$postdil_RFU_avg <- rowMeans( postdilRFU_data[ , c("postdil_RFU_1", "postdil_RFU_2", "postdil_RFU_3")] )
  } 
}

#Calculate standard deviation
postdilRFU_data$postdil_RFU_sd <- apply(postdilRFU_data[, c("postdil_RFU_1","postdil_RFU_2", "postdil_RFU_3")],1,sd)

#Calculate coeffient of variation
postdilRFU_data$postdil_RFU_cv <- postdilRFU_data$postdil_RFU_sd /  postdilRFU_data$postdil_RFU_avg

#If dilution= true, get growth rate from post dil RFU
for (i in (sample_count+1):length(RFU_data$Date)) {
  if (postdilRFU_data$Dil_done[i] == "FALSE"){
    postdilRFU_data$u[i] <- log(postdilRFU_data$predil_RFU_avg[i]/postdilRFU_data$predil_RFU_avg[i-sample_count])/
      (postdilRFU_data$Day[i]-postdilRFU_data$Day[i-sample_count])
  }
  else{
    postdilRFU_data$u[i] <- log(postdilRFU_data$predil_RFU_avg[i]/postdilRFU_data$postdil_RFU_avg[i-sample_count])/
      (postdilRFU_data$Day[i]-postdilRFU_data$Day[i-sample_count])
  }
}

#Get percent difference in u from dilution to dilution
for (i in 1:(length(RFU_data$Date)-sample_count)) {
  RFU_data$per_dif[i+sample_count] <-  ((RFU_data$u[i+sample_count] - RFU_data$u[i]) / RFU_data$u[i])*100
}
  

#Export the CSV
write.csv(postdilRFU_data, 
          file = paste(file_name, ".csv", 
                       sep= "_postdilution"))



# ----------------------------------------------------------------
#Plot growth rate over time
ggplot(postdilRFU_data, aes(Day, u, colour = Sample)) +
  geom_point() +
  theme_classic() +
  labs(y= 'Growth Rate', x= 'Time (Days)') +
  ylim(0, .5) +
  ggtitle('u 1 deg') +
  theme(text = element_text(size=15))
  

#Plot growth
ggplot(postdilRFU_data, aes(Day, predil_RFU_avg, colour = Sample)) +
  geom_line()
  labs(y= 'log(RFU)', x= 'Time (Days)') +
  ggtitle('RFU 6 deg') +
  theme(text = element_text(size=15))
