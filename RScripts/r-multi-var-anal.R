# R STUDIO WINDOW


# USER SESSION

R.Version() # list
R.version # cat information
sessionInfo()
ls()
rm(list = ls()) # clears the whole workspace

# FILE MANAGEMENT

getwd()
dirname('~')
list.dirs() # or dir()
list.files()
dir.create('mydir')
setwd('mydir')
setwd('~') # to come back to home directory
unlink('mydir', recursive = TRUE) # to delete directory
sink() # to sink results

read.csv(obj, path)
write.csv(obj, path)

# IO 

print('hello world!')

name <- 'mk'
sprintf('hello world - %s', name)

i <- 1
sprintf('the number is: %d', i)

i <- 1.5
sprintf('the number is: %f', i)

age = 30
cat('Hi!', name, ', Your age is', age)

# DATA TYPES/STRUCTURES

# c operator
# colon operator
# vectors
# factors
# list
# matrices
# data frames


# OPERTORS
x <- TRUE
y <- FALSE
# AND table
print(x && x) # TRUE
print(x && y) # FALSE
print(y && x) # FALSE
print(y && y) # FALSE

# OR Table
print(x || x) # TRUE
print(x || y) # TRUE
print(y || x) # TRUE
print(y || y) # FALSE

# CONTROL STRUCTURE

# conditional statements

if (TRUE){
  print('true')
} else {
  print('false')
}

x <- c(1, 2, 1, 1, 1, 2, 1, 2)
gender <- ifelse(x == 1, 'male', 'female')

switch(
  'mean', 
  'mean' = mean(1:10),
  'sd' = sd(1:10),
  'max' = max(1:10),
  'min' = min(1:10)
  
)

# loops

for(i in 1:10){
  print(i)
}

for(i in 1:10){
  # print(letters[i])
  print(LETTERS[i])
}

c = 0
while(c <=5){
  print(c)
  c = c+1
}

repeat{
  
  i <- rnorm(1)
  print(i)
  
  if (i <0){
    print('There is negative number!')
    break
  }
}

# FUNCTIONS

# anonymous 

(function(a, b) a+b)(1, 2)

emptyFunction <- function(){
  print('this is empty function')
}
# emptyFunction()

functionWithParams <- function(a, b, c=NULL){
  for (i in c(a, b, c)){
    print(i)
  }
}
# functionWithParams(1, 2)
# functionWithParams(1, 2, 3)

functionWithStructures <- function(a, b, List=NULL){
  for (i in c(a, b, List)){
    if (i == is.list(List)){
      print(dim(List))
    } else {
      print(i)  
    }
    
  }
}
# functionWithStructures(1, 2)
# functionWithStructures(1, 2, list(a='a', b='b'))

functionWithEmptyArgs <- function(...){
  x <- c(...)
  for (i in x){
    print(i)
  }
}
# functionWithEmptyArgs('a', 'b')
# functionWithEmptyArgs(1:5)
# functionWithEmptyArgs(list(a='this is a', b='this is b'))


findSummary <- function(vec, stat){
  out <- switch(
    stat, 
    'mean' = mean(vec),
    'sd' = sd(vec),
    'max' = max(vec),
    'min' = min(vec)
  )
  return(out)
}

# infix functions
# 1+2
# `+`(1, 2)

`%</>%` <- function(a, b){
  if (a > b){
    cat(a, 'is greater than', b)
    # message(a, 'is greater than', b) #prints as warning 
  } else {
    cat(b, 'is greater than', a)
    # message(b, 'is greater than', a)
  }
}

2%</>%1
2%</>%3

# DATA SIMULATIONS

# probabilities 

rnorm(1, 0, 1)

pnorm(1.64)-pnorm(-1.64) # 0.90
pnorm(1.96)-pnorm(-1.96) # probability 0.95
pnorm(2.576)-pnorm(-2.576) # 0.99

-qnorm((1-0.95)/2) # quantiles

# densities - normal distribution
x = seq(-3, 3, by=0.1)
y = dnorm(x, mean(x), sd(x))

# quantile distribution
x <- -1000:1000
qnorm(x, mean(x), sd(x))
plot(x, qnorm(x, mean(x), sd(x)), type='l', col='red')

# probabilities
x <- seq(-3, 3, by=0.1)
y <- pnorm(x)
plot(x, y, type='l', col='red')

# png(file="dnormExample.png")
# scatter.smooth(x, y)
x <- seq(-3:3, by=0.01)
hist(x, freq = F); lines(density(x), col='red')

x <- rnorm(1000)
hist(x, freq = F); lines(density(x), col = 'red') # perfect

# chisq

n <- 1000
x <- rchisq(n = n, df = 3)
hist(x = x, freq = FALSE, xlim = c(0, 15))
lines(x = density(x = x), col = "red")

# binomial

# rbinom
dbinom(3, size = 13, prob = 1 / 6)
# 1 - dbinom(3, size = 13, prob = 1 / 6)
probabilities <- dbinom(x = c(0:10), size = 10, prob = 1 / 6)
data.frame(probabilities)
plot(0:10, probabilities, type = "b", col = 'red')

#pbinom
pbinom(3, size = 13, prob = 1 / 6)
plot(0:10, pbinom(0:10, size = 10, prob = 1 / 6), type = "b")

#qbinom

qbinom(0.8419226, size = 13, prob = 1 / 6)
x <- seq(0, 1, by = 0.1)
y <- qbinom(x, size = 13, prob = 1 / 6)
plot(x, y, type = 'l')

# poisson

dpois(2, 3)
dpois(6, 6)

ppois(2, 3)
ppois(6, 6)

y <- c(.01, .05, .1, .2)
qpois(y, 2)
qpois(y, 6)

rpois(2, 3)
rpois(6, 6)

# t dist
n <- 10000
df <- n - 1
samples <- rt(n, df)
hist(samples, breaks = 'Scott', freq = FALSE)

x <- seq(-4, 4, by = 2)
dt(x, df = 4)

x <- seq(-2, 2, by = 2)
pvalues <- pt(x, df = 5, lower.tail = TRUE)
qt(pvalues, df = 5, lower.tail = FALSE)

# data simulations

sal <- round(runif(30, 10, 100), )
sal
salint <- cut(sal, c(10, 20, 40, 60, 80, 100))
salint
table(salint)

factData <- function(len, cat, type="numerical"){
  if (missing(type)){
    out <- sample(cat, len, replace = TRUE)  
    
  } else if (type=="numerical"){
    out <- sample(1:length(cat), len, replace = TRUE)
    
  }
  return(out)
}

multiVarCatData <- function(len, ...){
  cats <- list(...)
  out <- list()
  for (i in 1:length(cats)){
    out[[i]] <- factData(len, cats[[i]])
  }
  return(out)
}

vecData <- function(len, min=1, max=5){
  out <- runif(len, min, max)
  return(out)
}

makeMonths <- function(n){
  monthslist <- c('january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'septermber', 'october', 'november', 'december')
  if(n <=12){
    out <- monthslist[1:n] 
  } else{
    print('invalid number')
  }
  return(out)
}


multiVarNumData <- function(rows, cols, r=2){
  len <- rows * cols
  out <- matrix(round(runif(len), r), rows, cols)
  return(out)
}

gender <- factData(30, c('male', 'female'), type='numerical')
education <- factData(30, c('primary', 'secondary', 'higher'), type='numerical')
salary <- round(vecData(30), 2)
age <- round(vecData(30), 2)
satisfaction <- round(vecData(30), 0)
hetdataset <- data.frame(satisfaction, gender, education, salary, age)
head(hetdataset)

printCustomer <- function(name){
  print(name)
} 

printRandomNumbers <- function(typeofthedist, n){
  if (typeofthedist == 'normal'){
    out <- rnorm(n)
  } else if(typeofthedist == 'uniform'){
    out <- runif(n, 1, n)
  } else if (typeofthedist == 'binomial'){
    out <- rbinom(n, n, 0.5)
  } else if (typeofthedist == 'poisson'){
    out <- rpois(n, 0.5)
  }
  return(out)
}

#usage

nrd <- printRandomNumbers('normal', 10)
hist(nrd, freq = FALSE); lines(density(nrd), col = 'red') # hist in-built function

sales <- multiNumData(12, 3)
colnames(sales) <- paste('Product.', LETTERS[1:3], sep = '')
rownames(sales) <- makeMonths(12)
sales <- round(sales*100, 2)
head(sales)

# handling missing data
# https://rpubs.com/chibueze99/MissingR
# uni
x <- sample(1:5, 10, replace = T)
x[1] <- NA; x[5] <- NA

x.na.zeros <- ifelse(is.na(x), 0, x)
x.na <- is.na(x)
x.na.means <- ifelse(is.na(x), mean(x, na.rm = T), x)

# t.test(x.na.zeros ~ x.na) # change is significant
# t.test(x.na.means ~ x.na) # change is insignificant


# bi
x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
dimnames(x)[[1]] <- letters[1:8]
apply(x, 2, mean, trim = .2)

x[1, 2] <- NA; x[5, 2] <- NA

for(i in 1:ncol(x)){
  x[,i][is.na(x[, i])] <- mean(x[, i], na.rm = T)
}

# multi

data_set <- matrix(runif(40, 1, 5), 10, 4)
data_frame <- data.frame(data_set)
data_frame <- round(data_frame, 2)
data_frame[sample(1:3, 3, replace = F), sample(1:3, 3, replace = F)] <- NA

data_frame.1 <- data_frame

data_frame.1[, 1][is.na(data_frame.1[, 1])]

missing.x1 <- is.na(data_frame.1$X1)
# missing.x1
tapply(data_frame.1$X1, INDEX = missing.x1, FUN=mean, na.rm=T)

for (i in 1:ncol(data_frame.1)){
  data_frame.1[, i][is.na(data_frame.1)[, i]] <- mean(data_frame.1[,i], na.rm = T)
  }

# list-wise deletion

data_frame.2 <- data_frame

na.omit(data_frame.2)
data_frame.2[complete.cases(data_frame.2), ]

# pair-wise deletion

colnames(data_frame.1)[colSums(is.na(data_frame.1))>0] # to know the cols of missing data
data_frame.pair_wise <- data_frame.1[complete.cases(t(data_frame.1))]

# data_frame.pair_wise

# data coding 

indianpoverty <- read.csv(paste("D:\\Work\\Books\\MultiVarAnal_\\MultiVarAnal\\Data\\indianpoverty_ranked.csv", sep = ""))
names(indianpoverty)
summary(indianpoverty)

help("cut")
indianpoverty$Poverty_interval <- cut(indianpoverty$Poverty, c(0, 10, 20, 30, 40))
indianpoverty$Poverty_ordinal <- as.numeric(indianpoverty$Poverty_interval)
head(indianpoverty)

library(anacor)
data(sleeping)
sleeping_cat <- sleeping
temp_cat <- cut(sleeping$Temperature, c(-20, -1, 7), labels = c("warm", "cold")) 
sleeping_cat$Temperature <- temp_cat
weight_cat <- cut(sleeping$Weight, c(700, 1100, 2200), labels = c("light", "heavy")) 
sleeping_cat$Weight <- weight_cat
price_cat <- cut(sleeping$Price, c(100, 250, 350, 700), 
                 labels = c("cheap", "medium", "expensive"))  
sleeping_cat$Price <- price_cat
sleeping_cat
burtTable(sleeping_cat)


# MULTIPLE CORRELATION

satisfaction <- sample(c(0, 1), 30, replace = TRUE)
age <- round(vecData(30), )
education <- sample(c(0, 1), 30, replace = TRUE)
marital_status <- sample(c(0, 1), 30, replace = TRUE)
mulcordata <- data.frame(satisfaction, age, education, marital_status)

y <- satisfaction
x <- data.frame(age, education, marital_status)

ryx <- matrix(NA, dim(x)[2], 1)
for (i in 1:length(names(x))){
  ryx[i] <- cor(y, x[, i])
}

rxx <- cor(x)

## multiple correlation using matrix multiplication   

r2 <- t(ryx) %*% solve(rxx) %*% ryx
r2

multicor.model <- lm(satisfaction ~ age + education+marital_status, data = mulcordata)
r <- cor(multicor.model$model$satisfaction, multicor.model$fitted.values)

# polycoric correlation
## creating data from cor matrix

library(mvtnorm)
R <- matrix(0, 4, 4)

R[upper.tri(R)] <- runif(6)
diag(R) <- 1
R <- cov2cor(t(R) %*% R)
round(R, 4) # population correlations
data <- rmvnorm(1000, rep(0, 4), R)
round(cor(data), 4)

x1 <- data[,1]
x2 <- data[,2]
y1 <- cut(data[,3], c(-Inf, .75, Inf))
y2 <- cut(data[,4], c(-Inf, -1, .5, 1.5, Inf))
data <- data.frame(x1, x2, y1, y2)
hetcor(data) # Pearson, polychoric, and polyserial correlations, 2-step est.

library(polycor)
# simulation
sat1 <- sample(c(0, 1), 30, replace = TRUE)
sat2 <- sample(c(0, 1), 30, replace = TRUE)
age <- round(vecData(30), )

education <- sample(c(0, 1), 30, replace = TRUE)
marital_status <- sample(c(0, 1), 30, replace = TRUE)
multicordata <- data.frame(sat1, sat2, age, education, marital_status)
head(multicordata)
round(cor(multicordata), 2)

# polychoric correlation

hetcor(multicordata, ML=TRUE)

# polyserial

age_cat <- cut(age, c(20, 40, 60, 80, 100))
multicordata <- data.frame(sat1, sat2, age_cat, education, marital_status)
hetcor(multicordata, ML=TRUE)

# partial 

codeData <- function(dat){
  out <- unclass(as.factor(dat))
  return(out)
}

sanitation <- sample(1:5, 30, replace = TRUE)
social.distance <- sample(1:5, 30, replace = TRUE)
covid19.infection <- sample(c('yes', 'no'), 30, replace = 30)

covid19.dataset <- data.frame(sanitation, social.distance, covid19 = as.numeric(codeData(covid19.infection)))

cor(covid19.dataset)

library(ppcor)
pcor(covid19.dataset)
pcor.test(covid19.dataset$social.distance, covid19.dataset$covid19, covid19.dataset$sanitation)
  
spcor.test(covid19.dataset$social.distance, covid19.dataset$covid19, covid19.dataset$sanitation)
spcor.test(covid19.dataset$sanitation, covid19.dataset$covid19, covid19.dataset$social.distance)

# canonical cor
## data sim
southsales <- multiVarNumData(10, 4)
colnames(southsales) <- paste('city', 1:4)
head(southsales)

westsales <- multiVarNumData(10, 4)
colnames(westsales) <- paste('city', 1:4)
head(westsales)

library(CCA)
out <- cancor(southsales, westsales)
out

# simple regression

y = rnorm(10)
x = rnorm(10)

lm(y~x, data=data.frame(y, x))


icecreamsales <- round(runif(10)*100, )
temperature <- round(runif(10)*100, )
summary(lm(icecreamsales ~ temperature))

# multiple regression

gender <- factData(30, c('male', 'female'), type='numerical')
education <- factData(30, c('primary', 'secondary', 'higher'), type='numerical')
salary <- round(vecData(30), 2)
age <- round(vecData(30), 2)
satisfaction <- round(vecData(30), 0)
hetdataset <- data.frame(satisfaction, gender, education, salary, age)
head(hetdataset)

# summary(lm(satisfaction ~ gender + education + salary + age, data = hetdataset))

mva.out <- lm(satisfaction ~ gender + education + salary + age, data = hetdataset)
summary(mva.out)
confint(mva.out) # conf. int. for coefficients

# manova

# MDS

# non-metric

indianpoverty <- read.csv(paste("D:\\Work\\Books\\MultiVarAnal_\\MultiVarAnal\\Data\\indianpoverty_ranked.csv", sep = ""))
names(indianpoverty)
summary(indianpoverty)

# help("cut")
# indianpoverty$Poverty_interval <- cut(indianpoverty$Poverty, c(0, 10, 20, 30, 40))
# indianpoverty$Poverty_ordinal <- as.numeric(indianpoverty$Poverty_interval)
# head(indianpoverty)

# d <- dist(indianpoverty$Poverty_ordinal)
d <- dist(indianpoverty$Poverty)
fit <- cmdscale(d, eig=TRUE, k=2)

x <- fit$points[, 1]
y <- fit$points[, 2]
plot(x, y, xlab="Dimension 1", ylab="Dimension 2", main="Non-Metric MDS", type="n")
text(x, y, labels = indianpoverty[, 1], cex=.7)













