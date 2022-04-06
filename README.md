# Aluminum criterion calculator
This is an adaptation of the <a href="https://www.epa.gov/wqc/aquatic-life-criteria-aluminum" target="_blank">EPA aluminum criterion calculator R code</a>.  

Changes made to the original calculator primarily included package namespace clarifications to prevent errors in new workspaces.

This repository also includes a draft function wrapper (calcAlCrit.R) for the calculator to:  
1. Make the primary input a data frame format,  
2. Allow flexible naming of the required columns through function arguments,  
3. Store required taxa toxicity data within the function so that external toxicity data does not need to be sourced, and  
4. Include package dependencies and specify namespaces as appropriate if incorporated into a package.  

To run the function version:
```{r}
library(dplyr)
library(data.table)

# Build some test data
pH=sample(seq(6.3, 8.4, 0.1))
Hardness_mgL=sample(seq(25, 550, 25))
DOC_mgL=sample(seq(0.1, 2.2, 0.1))
sampID=paste0(rep('samp', 22), seq(1,22,1))
test_data=data.frame(sampID, pH, Hardness_mgL, DOC_mgL)
head(test_data)

# Source the function (or include in and install from a package)
source('calcAlCrit.R') # source the function from a local copy OR:
devtools::source_url("https://github.com/utah-dwq/wqTools/blob/master/R/calcAlCrit.R?raw=TRUE") # Source function from GitHub

# Apply function to test_data
criteria=calcAlCrit(test_data, 'pH', 'Hardness_mgL', 'DOC_mgL')
```





