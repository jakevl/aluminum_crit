
# File: Aluminum_Criteria_Calculator_R_Code.R
# Updated: 11/15/2019
# Code developed using R version 3.4.1

# Purpose: This script takes data from multiple paired sampling events and calculates  
#          criteria magnitude values consistent with EPA's 2018 final 304(a) aquatic life 
#          criteria recommendations for aluminum.

# Notes for use:
# - You will need to edit capitalized placeholder text in rows # 31, 175, 180, 194-196, 198, 
#   and 205-206 with your specific file paths or column names in order to run the code.
# - Make sure to fill in placeholders as directed, following the examples.
# - To run this script, you will also need "final304a_acute.csv" and "final304a_chronic.csv".
# - These CSV files feed into a function that calculates acute and chronic criteria
#   magnitude values (CMC/CCCs) for sampling events with paired pH, hardness, and DOC data.
# - R is case-sensitive and errors will occur if file names/paths don't match exactly.
# - If you've renamed the two CSV files containing the acute and chronic toxicity datasets, 
#   you will need to edit the file names accordingly in rows # 34 and 46.
# - The following code relies on functions contained in the dplyr and data.table packages.
#   Install these packages if you don't already have them before proceeding.
# - If water chemistry data inputs fall outside accepted model bounds, they will be flagged,
#   and caps will be placed on the inputs and reflected in your output file. 

# Load libraries (install packages if needed) ---------------------------------------------------
library("dplyr")
library("data.table")

# Set working directory where "final304a_acute.csv" and "final304a_chronic.csv" are located. 
# Example: ("C:/Users/MyName/Documents/MyFolder")
setwd("PLACEHOLDER:/WORKING DIRECTORY FILEPATH HERE")

# Read in 304(a) acute and chronic datasets and prepare for processing --------------------------
ac <- read.csv("final304a_acute.csv", stringsAsFactors = FALSE, strip.white = TRUE, skip = 2)
ac$data <- "Acute"
ac <- ac[, c(1, 3:5, 7:8, 11, 16:20)] # remove blank columns
colnm <- c("Species", "Hardness", "pH", "DOC", "EC_LC", "Reference", "Reason_Excluded",  
           "Genus", "Taxa_Group", "DOC_Notes", "Dilution_Water", "data")
names(ac) <- colnm
NArows <- row.names(ac[is.na(ac$pH), ]) # remove blank rows
ac <- ac[!row.names(ac) %in% NArows, ]
rm(NArows)
NOinput <- row.names(ac[ac$Hardness == "-", ]) # remove studies lacking input parameters
ac <- ac[!row.names(ac) %in% NOinput, ]

ch <- read.csv("final304a_chronic.csv", stringsAsFactors = FALSE, strip.white = TRUE, skip = 2)
ch$data <- "Chronic"
ch <- ch[, c(1:4, 7:8, 11, 16:20)] # remove blank columns
names(ch) <- colnm
NArows <- row.names(ch[is.na(ch$pH), ]) # remove blank rows
ch <- ch[!row.names(ch) %in% NArows, ]
rm(NArows)

all <- rbind(ac, ch) # merge acute and chronic datasets
all$EC_LC <- gsub(",", "", all$EC_LC) # remove special characters
fields <- c("Hardness", "pH", "DOC", "EC_LC") # convert fields to numeric format
all[fields] <- sapply(all[fields], as.numeric)

invert <- c("Physa", "Lampsilis", "Lymnaea", "Aeolosoma", "Brachionus", "Ceriodaphnia", 
            "Daphnia", "Stenocypris", "Crangonyx", "Hyalella", "Acroneuria", "Chironomus", 
            "Paratanytarsus", "Nais", "Melanoides")
vert <- c("Hyla", "Rana", "Oncorhynchus", "Salmo", "Salvelinus", "Lepomis", "Hybognathus", 
          "Pimephales", "Micropterus", "Danio", "Poecilia")

all$Grouping <- ifelse(all$Genus %in% invert, "Invertebrate",
  ifelse(all$Genus %in% vert, "Vertebrate",
    ""
  )
)

all <- all[all$Reason_Excluded == "", ] # remove studies excluded from SMAV

# Load function to run 304(a) aluminum MLRs  ----------------------------------------------------
criteria_calc <- function(pH, hardness, DOC) {

  # check inputs against bounds, apply caps if necessary
  Flag <-
    ifelse(pH <= 10.5 &
             pH >= 5 &
             hardness <= 430 &
             DOC <= 12,
           "",
           "Outside MLR model bounds. Caps applied.")
  
  pH <- ifelse(pH > 10.5, 10.5,
    ifelse(pH < 5, 5,
      pH
    )
  )

  hardness <- ifelse(hardness > 430, 430,
                     ifelse(hardness <= 0, 0.01,
                            hardness
                     )
  )
  
  DOC <- ifelse(DOC > 12, 12,
                ifelse(DOC <= 0, 0.08,
                       DOC
                )
  )

  # MLRs
  all$Normalized_Conc <- ifelse(all$Grouping == "Invertebrate",
    # if invertebrate, apply invertebrate MLR
    exp(log(all$EC_LC) - 0.597 * (log(all$DOC) - log(DOC)) -
      8.802 * (all$pH - pH) - 2.089 * (log(all$Hardness) - log(hardness)) +
      0.491 * (all$pH^2 - pH^2) + 0.23 * (all$pH * (log(all$Hardness)) - pH * (log(hardness)))),
    # if not invertebrate, apply vertebrate MLR
    exp(log(all$EC_LC) - 2.209 * (log(all$DOC) - log(DOC)) -
      2.041 * (all$pH - pH) - 1.862 * (log(all$Hardness) - log(hardness)) +
      0.232 * (all$pH * (log(all$Hardness)) - pH * (log(hardness))) +
      0.261 * (all$pH * log(all$DOC) - pH * log(DOC)))
  )

  # calculate SMVs
  all <- all %>%
    group_by(data, Species) %>%
    mutate(Species_Mean_Value_ug_L = exp(mean(log(Normalized_Conc))))

  # calculate GMVs - take geomean of unique SMAVs
  all <- all %>%
    group_by(data, Genus) %>%
    mutate(Genus_Mean_Value_ug_L = exp(mean(log(unique(Species_Mean_Value_ug_L)))))

  # generate ranked GMAV/GMCV table
  summary <- all[, c(16, 8, 12)] # subset to "Genus_Mean_Value_ug_l", "Genus", and "data" 
  summary <- unique(summary)
  summary <- summary %>%
    arrange(data, Genus_Mean_Value_ug_L) %>%
    group_by(data) %>%
    mutate(N = length(Genus_Mean_Value_ug_L)) %>%
    mutate(Rank = c(1:length(Genus_Mean_Value_ug_L))) %>%
    filter(Rank %in% c(1:4)) %>%
    mutate_at(vars(Rank, N), list(~as.numeric)) %>%
    group_by(data) %>%
    mutate(lnGMV = log(Genus_Mean_Value_ug_L)) %>%
    mutate(lnGMV2 = lnGMV^2) %>%
    mutate(P = Rank / (N + 1)) %>%
    mutate(sqrtP = sqrt(P)) %>%
    mutate(sum_lnGMV = sum(lnGMV)) %>%
    mutate(sum_lnGMV2 = sum(lnGMV2)) %>%
    mutate(sum_P = sum(P)) %>%
    mutate(sum_sqrtP = sum(sqrtP)) %>%
    group_by(data) %>%
    mutate(S2 = (sum_lnGMV2 - ((sum_lnGMV^2) / 4)) / (sum_P - ((sum_sqrtP^2) / 4))) %>%
    mutate(L = (sum_lnGMV - (sqrt(S2) * sum_sqrtP)) / 4) %>%
    mutate(A = (sqrt(S2) * sqrt(0.05)) + L) %>%
    mutate(FV = exp(A))

  CCC <- as.numeric(unique(summary[summary$data == "Chronic", "FV"]))
  FAV <- as.numeric(unique(summary[summary$data == "Acute", "FV"]))
  CMC <- FAV / 2
  Final_CMC <- round(CMC, digits = 2 - (1 + trunc(log10((abs(CMC))))))
  Final_CCC <- round(CCC, digits = 2 - (1 + trunc(log10((abs(CCC))))))

  df <- data.frame(pH, hardness, DOC, Flag, FAV, CMC, Final_CMC, CCC, Final_CCC)

  
  ranks <- summary %>% 
    mutate(rowid = row_number()) %>% 
    select(rowid, Genus_Mean_Value_ug_L, Genus, Rank,data)
  ranks <- melt(ranks, id.vars=c("rowid","Rank", "data")) 
  ranks$rowid <- 1
  ranks <- dcast(ranks, rowid ~ Rank+data+variable, value.var="value") 
  print_results <- cbind(df,ranks[,c(4:5,8:9,12:13,16:17,2:3,6:7,10:11,14:15)])
  return(print_results)
}

# Calculate CMC/CCCs ----------------------------------------------------------------------------
# Edit placeholder text to locate water chemistry inputs (pH, hardness, DOC) data file.
# Water chemistry data must be in CSV file format.
# Note: Before converting Excel files to CSV, make sure all significant digits are displayed.
# Example: ("C:/Users/MyName/Documents/MyFolder/SubFolder/FileName.csv")
df <- read.csv("PLACEHOLDER:/FOLDER WITH CHEMISTRY DATA AS CSV/NAME OF DATA FILE.csv")

# Edit placeholder text to specify where you'd like to save your results and the name of 
# your output file. 
# Example: "C:/Users/MyName/Documents/MyFolder/SubFolder/MyOutputFile.csv"
saveFileTo <- "PLACEHOLDER:/FOLDER TO SAVE RESULTS IN/NAME OF NEW OUTPUT FILE.csv"

# In the order specified, edit placeholder text with your chemistry file's pH, hardness, and
# DOC column names. Note that R reads blank spaces in CSV column headers as periods. 
# To double-check your chemistry data's column names, type names(df) into your R console. 
# Example: results <-
#           apply(df[, c(
#             "pH", 
#             "hardness",
#             "DOC"
#           )], 1, function(x)
#             criteria_calc(x["pH"], x["hardness"], x["DOC"]))
results <-
  apply(df[,c(
    "PLACEHOLDER FOR PH COLUMN NAME",
    "PLACEHOLDER FOR HARDNESS COLUMN NAME",
    "PLACEHOLDER FOR DOC COLUMN NAME"
  )], 1, function(x)
    criteria_calc(x["PH COLUMN NAME"], x["HARDNESS COLUMN NAME"], x["DOC COLUMN NAME"]))
results <- rbindlist(results)

# Edit placeholder text to specify column names of fields (e.g., site or sample ID, date) 
# that you would like printed with your results.
combined <-
  cbind(df[, c(
    "PLACEHOLDER FOR SAMPLE ID COLUMN NAME",
    "PLACEHOLDER FOR SAMPLE DATE COLUMN NAME"
  )], results)

# Save your results as a CSV file ---------------------------------------------------------------
write.csv(combined, saveFileTo, row.names = FALSE)
