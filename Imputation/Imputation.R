#' ---
#'title: "Multiple Imputation by Chained Equations"
#'author: Julia van Hees
#' date: April 6th, 2020
#' output: 
#'      html_document:
#'             toc: true
#'             highlight: default
#'             theme: readable
#' ---
#------------- Introduction -------------------------
#' 
#' **Included in this file:**
#' 
#' * Perform imputation using **mice** 
#' * Check variable distribution after imputation
#' <br>
#' <br>
#------------------------ Multiple imputation -----------------------------
#'
#' ## Multiple imputation using **mice** 
#' 
#' First load the data into R:
library(foreign, quietly = TRUE)
library(mice, quietly = TRUE)
part1 <- read.spss("Part1.sav", to.data.frame = TRUE, use.value.labels = TRUE) 
part2 <- read.spss("Part2.sav", to.data.frame = TRUE, use.value.labels = TRUE) 
#' 
#' Start to assess missingness. I start with looking at how many values are missing per variable and express this as percentage. Part 1 first:
cbind(sapply(part1, function(x) sum(is.na(x))), apply(part1, 2, function(col)sum(is.na(col))/length(col)))
#' And part 2 as well:
cbind(sapply(part2, function(x) sum(is.na(x))), apply(part2, 2, function(col)sum(is.na(col))/length(col)))
#' Then, let's visualise the missingness in a plot:
library(rms, verbose = FALSE, quietly = TRUE)
par(mfrow=c(1,1))
na.patterns_part1 <- naclus(part1)
naplot(na.patterns_part1, which=c('na per var'))
#' And plot this as well for part 2:
par(mfrow=c(1,1))
na.patterns_part2 <- naclus(part2)
naplot(na.patterns_part2, which=c('na per var'))
#' I will also want to look at the patterns of missingness. For this, I use the following function:
na.patterns <- naclus(part1)
par(mar=c(3,3,2,1))
plot(na.patterns, ylab="Fraction of NAs in common")
#' Part 2 pattern:
na.patterns <- naclus(part2)
par(mar=c(3,3,2,1))
plot(na.patterns, ylab="Fraction of NAs in common")
#' Some other cool functions for visualizing missingness:
library(VIM, quietly = TRUE)
par(mfrow=c(1,1), mar=c(10,10,10,10))
aggr(part1, sortVars=T, col=c("green", "red"))
#' <br>
#' Now that I have assessed the missigness per variables and its pattern, I am going to impute the missing values using mice:
#' 

#' 
#------------------------ Variable distribution -----------------------------
#'