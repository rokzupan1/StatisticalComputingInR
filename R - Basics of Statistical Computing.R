library(datasets)
head(iris)
install.packages('pacman')
?plot

plot(iris$Species) #Categorical variable
plot(iris$Petal.Length) #Quantitative variable
plot(iris$Species, iris$Petal.Width) # Cat x quant
plot(iris$Petal.Length, iris$Petal.Width) #Quant pair
plot(iris) #Entire data frame

#Plot with options
plot(iris$Petal.Length, iris$Petal.Width,
     col = "#cc0000", # Hex code for datalab.cc red
     pch = 19,        # Use solid circles for points
     main = "Iris: Petal Length vs Petal Width",
     xlab = "Petal Length",
     ylab = "Petal Width")

# PLOT FORMULAS WITH PLOT() ######
plot(cos, 0, 2*pi)
plot(exp, 1, 5)
plot(dnorm,-3,+3)

# Formula plot with options
plot(dnorm, -3, +3,
     col = "#cc0000",
     lwd = 5, # line width
     main = "Standard Normal Distribution", #TITLE
     xlab = "z-scores", 
     ylab = "Density")

library(datasets)
# Load data
?mtcars
head(mtcars)

# BAR CHARTS ###########################
barplot(mtcars$cyl) #doesn't work
cylinders <- table(mtcars$cyl) #Create table
barplot(cylinders) #Bar chart
plot(cylinders) #Default X-Y plot (lines)

# HISTOGRAMS ##########################
#For data that is quantitative, scaled, measured, interval or ratio level.
#Shape, gaps, outliers, symmetry --> what you look for in Histograms
?iris
head(iris)
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Petal.Width)

# HISTOGRAMS BY GROUP ##########################
par(mfrow = c(3,1)) #Put graphs in 3 rows and 1 column

#Histograms for each species using options
hist(iris$Petal.Width [iris$Species == "setosa"],
     xlim = c(0,3),
     breaks = 9,
     main = "Petal Width for Setosa",
     xlab = "",
     col = "red")

hist(iris$Petal.Width [iris$Species == "versicolor"],
     xlim = c(0,3),
     breaks = 9,
     main = "Petal Width for Versicolor",
     xlab = "",
     col = "purple")

hist(iris$Petal.Width [iris$Species == "virginica"],
     xlim = c(0,3),
     breaks = 9,
     main = "Petal Width for Virginica",
     xlab = "",
     col = "blue")

par(mfrow = c(1,1)) #Restore graphic parameter

# SCATTERPLOTS ######################################################
#For visualizing the associtation between two quantitative variables
#Linear, spread, outliers, correlation.

?mtcars
head(cars)

#Good to first check unvariate distributions
hist(mtcars$wt)
hist(mtcars$mpg)

#Basic X-Y plot for two quantitative variables
plot(mtcars$wt, mtcars$mpg,
     pch = 19, #solid circle
     cex = 1.5, #Make 150% size
     col = "#cc0000", #color 
     main = "MPG as a Function of Weight of Cars",
     xlab = "Weight (in 1000 pounds)",
     ylab = "MPG")

#Why overlay plots? Because you get increased
#information density. Use views that complement & support one another.

?lynx
head(lynx)
#histogram
hist(lynx,
     breaks = 14,
     freq = FALSE,
     col = "thistle1",
     main = paste("Histogram of Annual Canadian Lynx",
                  "Trappings, 1821-1934"),
     xlab = "Number of Lynx Trapped")

#add a normal distribution
#curve for normal distribution, density of the normal distribution,
#use the mean and sd of the lynx data
curve(dnorm(x, mean = mean(lynx), sd = sd(lynx)),
      col = "thistle4", #color of curve
      lwd = 2,          #line width of 2 pixels
      add = TRUE)       # superimpose on previous graph

# add two kernel density estimators
lines(density(lynx), col = "blue", lwd = 2)
lines(density(lynx, adjust = 3), col = "purple", lwd = 2)

rug(lynx, lwd = 2, col = "gray")


#after pictures, you want to get some precision: counts for categories,
#quartiles & mean for quantitative variables

summary(iris$Species) #Categorical variable
summary(iris$Sepal.Length) #Quantitative variable
summary(iris)

head(iris)
p_load(psych)
p_help(psych, web = F)

describe(iris)


hist(iris$Petal.Length[iris$Species=="versicolor"],
     main = "Petal Length: Versicolor")

hist(iris$Petal.Length[iris$Species=="virginica"],
     main = "Petal Length: Virginica")

hist(iris$Petal.Length[iris$Species=="setosa"],
     main = "Petal Length: Setosa")

#Short petals only (all Setosa)
hist(iris$Petal.Length[iris$Petal.Length < 2],
     main = "Petal Length <2")

# MULTIPLE SELECTORS ###################################################
hist(iris$Petal.Length[iris$Species == "virginica" &
            iris$Petal.Length < 5.5],
     main = "Petal Length: Short Virginica")

# CREATE SUBSAMPLE #####################################################
i.setosa <- iris[iris$Species == "setosa",] #going in iris data, selecting 
#species and picking just setosa which select rows, and then put comma,
#and leave it blank if you select all the columns.

#Vector 1+ numbers in 1d array
#All same data type, R's basic data type
#Matrix 2 dimensions, same length, same data class, columns not named.

#Data frame - most common, can have vectors of multiple types. All the same
#length, Closest R analogue to spreadsheet, special functions.

#List most flexible, ordered collection of elements, any class length, or 
#structure. list can include lists. 

#Data frame is the most optimal thing.
#Coercion is good in DATA SCIENCE. That is, changing a data object from
#one type to another. (charachter to logical, matrix to data frame, double to
#integer)

#numeric
n1 <- 15 #Double precision by defualt
n1
typeof(n1)

n2 <- 1.5
n2
typeof(n2)

c1 <- "c"
c1
typeof(c1)

c2 <- "a string of text"
c2
typeof(c2)

l1 <- TRUE
l1
typeof(l1)

l2 <- F
l2
typeof(l2)

v1 <- c(1,2,3,4,5)
v1
is.vector(v1)

v2 <- c("a","b","c")
v2
is.vector(v2)

v3 <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
v3
is.vector(v3)

m1 <- matrix(c(T,T,F,F,T,F), nrow = 2)
m1

m2 <- matrix(c("a", "b",
               "c", "d"),
             nrow = 2,
             byrow = T)

#ARRAY, Give data, then dimensions (rows,columns,tables)
a1 <- array(c(1:24), c(4,3,2))
a1

#DATA FRAME - can combine vectors of the same length

vNumeric <- c(1,2,3)
vCharacter <- c("a", "b", "c")
vLogical <-c(T,F,T)

dfa <- cbind(vNumeric,vCharacter,vLogical)
dfa #Matrix of one data type

df <- as.data.frame(cbind(vNumeric,vCharacter,vLogical))
df #Makes a data frame with three different data types

## List #######################################################

o1 <- c(1,2,3)
o2 <- c("a","b","c","d")
o3 <- c(T,F,T,T,F)

list1 <- list (o1,o2,o3)
list1

list2 <- list(o1,o2,o3,list1)
list2

## COERCING TYPES 
## Automatic coercion - goes to "least restrictive" data type

(coerce1 <- c(1,"b",TRUE))
typeof(coerce1)

## Coerce numeric to integer
(coerce2 <- 5)
typeof(coerce2)

(coerce3 <- as.integer(5))
typeof(coerce3)

## Coerce character to numeric
(coerce4 <- c("1","2","3"))
typeof(coerce4)

(coerce5 <- as.numeric(c("1","2","3")))
typeof(coerce5)

#Probably most common you wil do: coerce matrix to data frame
(coerce6 <- matrix(1:9, nrow = 3))
is.matrix(coerce6)

(coerce7 <- as.data.frame(matrix(1:9, nrow = 3)))
is.data.frame(coerce7)

# Clear environment
rm(list = ls())

#FACTORS - Categories & Names, factor is an "attribute" of a vector that 
#specifies the possible values & their order

(x1 <- 1:3)
(y1 <- 1:9)
#combine variables

(df1 <- cbind.data.frame(x1,y1))
typeof(df1$x1)

#it repeats the values for x1 123123123
str(df1)
typeof(x1)

# AS.FACTOR ##########################

(x2 <- as.factor(c(1:3)))
typeof(x2)
(df2 <- cbind.data.frame(x2,y1))
typeof(df2$x2)
str(df2)

rm(list = ls())

#ENTERING DATA

# = is generally poor form in R, <- is used to assign values to a variable
# Assigns number 0 through 10 to x1

x1 <- 0:10
x1

x2 <- 10:0
x2

(x3 <- seq(10))
(x4 <- seq(30,0,by=-3))

x5 <- c(5,4,1,6,7,2,2,3,2,8)
x5 #c=concatenate (or combine or collect)

# SCAN ###################3

x6 <- scan()
x6


x7 <- rep(TRUE,5)
x7

x8 <- rep(c(TRUE, FALSE),5)
x8

x9 <- rep(c(TRUE, FALSE), each = 5)
x9

rm(list = ls())

### Importing data
### CSV, TXT, XLSX, JSON
# R has built in functions for importing data in many formats
# rio - r input output, rio combines all of  R's import functions into one
# simple utility

pacman::p_load(pacman,rio)

#about excel files: better to save in csv files and the import it

#import with RIO
rio_csv <- import("mbb.csv")
head(rio_csv)

rio_txt <- import("mbb.txt")
head(rio_txt)

rio_xlsx <- import("mbb.xlsx")
head(rio_xlsx)

View(rio_csv) #it is sortable and you can access through environment

#R's built-in function for text files (used by rio)
r_txt1 <- read.table("mbb.txt", header = TRUE, sep = "\t")
r_txt1

trends.csv <- read.csv("mbb.csv", header = TRUE)

rm(list = ls())

##HIERARCHIAL CLUSTERING
#LIKE WITH LIKE, but similarity depends on your criteria
#Hierarchical vs set k.
#Measures of distance
#Divisive vs agglomerative

head(mtcars)
cars <- mtcars[,c(1:4,6:7,9:11)]
head(cars)

#COMPUTE AND PLOT CLUSTERS
hc <- cars %>% #Gets cars data
      dist %>% #Compute distance/dissimilarity matrix
      hclust   #Computer hierarcical clusters

plot(hc) #Plot dendrogram

#Add boxes to plot

rect.hclust(hc, k=2, border = "gray")
rect.hclust(hc, k=3, border = "blue")
rect.hclust(hc, k=4, border = "green4")
rect.hclust(hc, k=6, border = "darkred")
rect.hclust(hc, k=7, border = "pink")


#less noise & fewer unheplful variables in data = more meaning
#PRINICPAL COMPONENTS AKA DIMENSIONALITY REDUCTION #################
# PCA - two variables, regression line, perpendicular distance to the 
# regression line, collapse, thats the PC.
# Went from 2D to 1D but maintained the most important information!
# We made analysis & interpreatation easier and more reliable

pacman::p_load(pacman, dplyr, GGally,ggplot2, ggthemes,
               ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny,
               string, tidyr)

head(mtcars)
cars <- mtcars[,c(1:4,6:7,9:11)]
head(cars)

#COMPUTE PCA

pc <- prcomp(cars,
        center = TRUE, #Centers means to 0 (optional)
        scale = TRUE)  #Sets unit variance (heplful)

#other way: to specify variables

pc <- prcomp(~ mpg + cyl + disp + hp + wt + qsec + am +
               gear + carb,
               data = mtcars,
               center = TRUE,
               scale = TRUE)

summary(pc) #get summary stats
plot(pc)  #screenplot for number of components
pc #get standars deviations and rotation
predict(pc) %>% round(2) #see how cases load on PRINCIPAL COMPONENTS

biplot(pc)


## REGRESSION #############################################################
#Out of many variables, one variable -> Out of many scores, one score
#The idea of regression is that u use many variables to predict scores on one
#variable
#There are many versions & adapatitions of regression to make it flexible and
#powerful.

?USJudgeRatings
head(USJudgeRatings)
#We have 6 judges listed by name and we have scores on a number of different
#variables and finishes with are they worthy of RTEN retention.
data <- USJudgeRatings

#Define variable groups
x <- as.matrix(data[-12]) #we created matrix, consist all of the predictor
#variables simultaneously. Read data, read all of the columns except 12.
y <- data[,12] #use all the rows but only read the 12th columns, that is the one
#that has retentionary outcome.

#REGRESSION WITH SIMULTANEOUS ENTRY
# Using variable groups
reg1 <- lm(y ~ x)

# Or specify variables individually
reg1 <- lm(RTEN ~ CONT + INTG + DMNR + DILG + CFMG +
           DECI + PREP + FAMI + ORAL+ WRIT + PHYS,
           data = USJudgeRatings)

# Results
reg1 #Coefficients only
summary(reg1) #Inferential tests

# MORE SUMMARIES ######

anova(reg1)                # Coefficients w/inferential tests
coef(reg1)                 # Coefficients (same as reg1)
confint(reg1)              # CI for coefficients
resid(reg1)                # Residuals case-by-case
hist(residuals(reg1))      # Histogram of residuals

# ADDITIONAL MODELS #######################################################

#Use two additional libraries
p_load(lars, caret)

#Conventional stepwise regression
stepwise <- lars(x,y,type="stepwise")

#Stagewise: Like stepwise but with better generalizability
forward <- lars(x,y, type = "forward.stagewise")

# LAR: Least Angle Regression
lar <- lars(x,y,type = "lar")

# LASSO: Least Absoulute Shrinkage and Selection Operator
lasso <- lars(x,y, type = "lasso")

# Comparsion of R^2 for new models
r2comp <- c(stepwise$R2[6], forward$R2[6],
            lar$R2[6], lasso$R2[6]) %>% 
            round(2)
names(r2comp) <- c("stepwise","forward","lar","lasso")
r2comp

