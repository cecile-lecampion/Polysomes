# choisir le répertoire de travail
setwd("~/partage/bio-informatique/poly_marwa")

#-------------------------------------------------------------------------------------
# Function to load the data
# Load the csv file and return a one cilumon data frame
# Column is named with parameter colName
# Usage exemple :
#    df <- f_load_polysome_csv("a1.csv", "a1")
#-------------------------------------------------------------------------------------
# The file contains 3 columns:
# 0.0001041666692,0.1534467638,
# 0.000312499993,0.1535461843,
# 0.0005208333605,0.1536204666,
# Only the 2nd column is kept
# Only numerics are kept
# Two first lines and useless lines at the end of the file are ommitted
#-------------------------------------------------------------------------------------
f_load_polysome_csv <- function(csvFile, colName) {
  lines <- readLines(csvFile)
  linesNbToImport <- grep("^\\d+\\.\\d+,\\d+\\.\\d+,", lines, perl = TRUE)
  df <- read.table(text = lines[linesNbToImport], header = FALSE, sep = ",", dec = ".")
  df <- as.data.frame(df$V2)
  colnames(df) <- c(colName)
  return(df)
}
#-------------------------------------------------------------------------------------

######################################################################################
# Partie du script à modifier
######################################################################################

a1 <- f_load_polysome_csv("line1-1.csv", "a1")
a2 <- f_load_polysome_csv("line1-2.csv", "a2")
b1 <- f_load_polysome_csv("line2-1.csv", "b1")
b2 <- f_load_polysome_csv("lien2_2.csv", "b2")
c1 <- f_load_polysome_csv("line3-1.csv", "c1")
c2 <- f_load_polysome_csv("line3-2.csv", "c2")

# Keep the generic names of the data in a variable
NAMES <- c("a1", "a2", "b1", "b2", "c1", "c2")

if (!require(dplyr)) { install.packages('dplyr')}
library(dplyr)

if (!require(purrr)) { install.packages('purrr')}
library(purrr)

# Create a liste and build a data frame with all data
liste <- list(a1, a2, b1, b2, c1, c2)
df <- liste %>% map( ~ .x %>% tibble::rownames_to_column('rn')) %>% reduce(full_join, by = 'rn')

#----------------------------------------------------------------------------------------------------
# Variables used in the script
# Limits for graphs axis
# On X axis
XMIN <- 1000
XMAX <- 3600

# On Y axis
YMIN <- 0
YMAX <- max(df[XMIN:XMAX, 2: ncol(df)], na.rm = TRUE) + (max(df[XMIN:XMAX, 2: ncol(df)], na.rm = TRUE))*0.1

# Area where to search for the monosome peak
ZMONO <- c(2200 : 3000)

# Area for the valley
# Area before the monosome
Z_BMONO <- c(1500:2000)
# Area after the monosome
Z_PMONO <- c(2100:2600)
# Area after the peak in supernatent
Z_PSUR <- c(2800:3200)

# Limits to compute replicats mean
MEAN_START <- 500
MEAN_END <- 5000

# Real names of the data
# Keep here the real names of the data
DATA_REAL_NAME <- c("Line1", "Line2", "Line3")

# Number of replicats
NB_REPLICAT <- 2
#-------------------------------------------------------------------------------------------------------

# Plot the raw data

plot.ts(df[XMIN:XMAX,2:ncol(df)],plot.type=c("single"), type = "l", 
        xlab = "", ylab = "",
        col=rainbow(ncol(df)))


# Identify the position of monosome peak : Find the max between ZMONO and XMAX for each curves

if (!require(splus2R)) { install.packages('splus2R')}
library(splus2R)

mono <-(apply((apply(df[ZMONO, 2:ncol(df)], 2, peaks, span=201, strict=TRUE, endbehavior=0)), 2, which))+2199

# Index of the monosme peak is the minimum of all monosome peak
mono_pic <- min(mono)  

# Align all curves on monosomes : To align all monosome peak one have to remove the n first lines of each column. 
# n is the difference between the monosome peak for the considered curve and the monosome peak as define by varaible mono_pic
# Compute the number of line to remove for each column
N <- mono-mono_pic

#__________________________________________________________________________________________________
# Function that remove the n first lines
remove_n_first_rows <- function(dat, n) if(n <= 0) dat else tail(dat, -n)
#__________________________________________________________________________________________________

# Remove the lines
Liste_mono_align <- mapply(remove_n_first_rows, df[, -1], N)

# Convert the list element to a data frame
Liste_mono_align <- lapply(Liste_mono_align, as.data.frame)
df_mono_align <- Liste_mono_align %>% map( ~ .x %>% tibble::rownames_to_column('rn')) %>% 
  reduce(full_join, by = 'rn')
colnames(df_mono_align) <- c("rn",NAMES)


# Plot the data to check that all monosmes peak are aligned

plot.ts(df_mono_align[XMIN:XMAX,2:ncol(df_mono_align)],plot.type=c("single"), type = "l", 
        xlab = "", ylab = "",
        col=rainbow(ncol(df_mono_align)))

# Minimum of the curves to 0
#Identify the minimum for each curve around the monosome and substract this value to each column
y_min <-apply(df_mono_align[XMIN:mono_pic, -1], 2, min)

# substract
df_zero_align <- as.data.frame(mapply('-', df_mono_align[, -1], y_min))

# Plot the data to check that the minimum of each curve is now 0
plot.ts(df_zero_align[XMIN:XMAX,2:ncol(df_zero_align)],plot.type=c("single"), type = "l", 
        xlab = "", ylab = "",
        col=rainbow(ncol(df_zero_align)))


# Mean of replicats
# The mean is computed in row for the n column corresponding to the number of replicats (NB_REPLICAT)

# Define the number of lines/conditions in the analysis
nb_line <- ncol(df_zero_align)/NB_REPLICAT

# Select the column to porcess
select_col <- c(1:NB_REPLICAT)

# Create an empty data frame (length of the data frame is length of the result that will be storein it)
df_mean <- data.frame(matrix(ncol = 0, nrow = length(c(MEAN_START : MEAN_END))))

#__________________________________________________________________________________________________
# Function that return the row mean of the selected columns   
f_compute_rep_mean <- function (select_col) {
  result <- as.data.frame(rowMeans(df_zero_align[MEAN_START : MEAN_END, select_col]))
  return(cbind(df_mean, result))
}
#__________________________________________________________________________________________________

# Loop over the data frame
for (LineCount in c(1:nb_line)) {
  df_mean <- f_compute_rep_mean(select_col)
  select_col <- select_col+NB_REPLICAT
}
colnames(df_mean) <- DATA_REAL_NAME

# Plot the data
plot.ts(df_mean[XMIN:XMAX, ],plot.type=c("single"), type = "l", 
        ylim = c(0, YMAX-(min(y_min))), 
        xlab = "", ylab = "",
        col=rainbow(ncol(df_mean)))

#smoothing with adjacent averaging method

if (!requireNamespace("forecast", quietly = TRUE)) { install.packages("forecast") }
library(forecast)

df_smoothed <- as.data.frame(apply(df_mean, 2, ma, order = 50, centre = TRUE))

# plot the data
plot.ts(df_smoothed[XMIN:XMAX, ],plot.type=c("single"), type = "l", 
        ylim = c(0, YMAX-(min(y_min))), 
        xlab = "", ylab = "",
        col=rainbow(ncol(df_smoothed)))
legend("topleft", legend = DATA_REAL_NAME,
       col = rainbow(ncol(df_smoothed)), lty = 1, lwd = 2,
       box.lty = 0)

# Check that the intervalls to search for the valley were correctly defined
# If the valley are contained in the intervals then the further analysis will be correct. 
# If not, the intervals must be redefined at the begining and the anlysis must be started again.

plot.ts(df_smoothed[, ],plot.type=c("single"), type = "l", 
        ylim = c(0, YMAX-(min(y_min))), 
        xlab = "", ylab = "",
        col=rainbow(ncol(df_smoothed)))
abline(v=c(Z_BMONO[1], Z_BMONO[length(Z_BMONO)]))
abline(v=c(Z_PMONO[1], Z_PMONO[length(Z_PMONO)]), col= "lightblue")
abline(v=c(Z_PSUR[1], Z_PSUR[length(Z_PSUR)]), col= "grey")

# There is a shift of 1000 in the index value. wich.min() returns indexes that start at the first boundary of the search interval. 
# For real index add (boundary - 1) 
# To identifi index on the smoothed graph, remove 1000.

# Find min before monosme
min_poly <-apply(df_smoothed[Z_BMONO, ], 2, which.min)+1499

# Find min after monosome
min_post_mono <-apply(df_smoothed[Z_PMONO, ], 2, which.min)+2099

# Find min after the peak in supernatent
min_post_sur <-apply(df_smoothed[Z_PSUR, ], 2, which.min)+ 2799

#Compute area under the curve
# add a column for indexes
df_smoothed$X <- c(1:nrow(df_smoothed))

# select the areas
# for polysome area
# create an empty list to store the results
List_auc_poly <- list()

for (i in c(1:nb_line)) {
  index1 <- XMIN
  index2 <- min_poly[[i]]
  df_auc_poly <- dplyr::filter(df_smoothed, X >= index1 & X <= index2) %>% 
    select(all_of(i), 4)
  List_auc_poly[[length(List_auc_poly)+1]]<- df_auc_poly
}

# for monosome area
# create an empty list to store the results
List_auc_mono <- list()

for (i in c(1:nb_line)) {
  index1 <- min_poly[[i]]
  index2 <- min_post_mono[[i]]
  df_auc_mono <- dplyr::filter(df_smoothed, X >= index1 & X <= index2)%>% select(all_of(i), 4)
  List_auc_mono[[length(List_auc_mono)+1]]<- df_auc_mono
}

# for supernatent area
# create an empty list to store the results
List_auc_sur <- list()

for (i in c(1:nb_line)) {
  index1 <- min_post_mono[[i]]
  index2 <- min_post_sur[[i]]
  df_auc_sur <- dplyr::filter(df_smoothed, X >= index1 & X <= index2) %>% select(all_of(i), 4)
  List_auc_sur[[length(List_auc_sur)+1]]<- df_auc_sur
}

# Combine 3 areas in a data frame
df_AUC <- data.frame((matrix(ncol = length(DATA_REAL_NAME), nrow = 3)), 
                     row.names = c("AUC_poly", "AUC_mono", "AUC_sur"))
colnames(df_AUC) <- DATA_REAL_NAME


if (!requireNamespace("DescTools", quietly = TRUE)) { install.packages("DescTools") }
library(DescTools)

# Compute AUC
for (element in c(1:nb_line)) {
  AUC_poly <- AUC(x = List_auc_poly[[element]][,2], y = List_auc_poly[[element]][,1], na.rm = TRUE)
  df_AUC[1, element] <- AUC_poly
}

for (element in c(1:nb_line)) {
  AUC_mono <- AUC(x = List_auc_mono[[element]][,2], y = List_auc_mono[[element]][,1], na.rm = TRUE)
  df_AUC[2, element] <- AUC_mono
}

for (element in c(1:nb_line)) {
  AUC_sur <- AUC(x = List_auc_sur[[element]][,2], y = List_auc_sur[[element]][,1], na.rm = TRUE)
  df_AUC[3, element] <- AUC_sur
}

# Compute percent for further plot
for (line in c(1:nb_line)) {
  nom <- sprintf("%s%d ", "percent Line", line)
  df_AUC[,nom] <- 100*(df_AUC[,line]/sum(df_AUC[,line]))
}

# Plot the results
barplot(as.matrix(df_AUC)[, (nb_line+1) : ncol(df_AUC)], legend = c("Poly", "mono", "surpernatent"), 
        col = c("peachpuff", "thistle1", "lightsteelblue1"), 
        names.arg = DATA_REAL_NAME, width = c(40,40), xlim = c(0,240))

