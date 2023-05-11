# 10/05/2023
## R script with basic commands
###################################################

# get current working directory
getwd()

# set new one, use tabs in R console, be aware that "\" instead of "/"
setwd("C:/Users/giulia.capitoli/Desktop/BANDI/LEZIONI/CORSO SPETTROMETRIA/CODICI-Introduzione") 


#check again
getwd()

# mathematics
3+5
10^3
log(10)
?log

log(10,10)
log(10,2)
log2(10)

# functions and help
?log

# equivalent commands
log10(100) # this is the log with basis 10 of 100
log(100,base=10)

# assignments to objects
# the assignment operator (‘<-’), consists of the two characters ‘<’ (“less than”)
# and ‘-’ (“minus”) occurring strictly side-by-side and it ‘points’ to the object
# receiving the value of the expression. 
# In most contexts the ‘=’ operator can be used as an alternative.
A <- log(10)
B <- (3 + 10) * 666
C<-B + A

# vectors, lists
# To set up a vector named x, consisting of different numbers, use the R command "c"
m_over_z <- c(144,288,287,366,200)
charges <- c(2,4,3,2,1)

# indices
m_over_z[2]
m_over_z[10]

# operations are generally element-wise
# Vectors can be used in arithmetic expressions, in which case the operations are
# performed element by element. 
# Vectors occurring in the same expression need not all be of the same length. 
# ATTENTION: If they are not, the value of the expression is a vector with the same length
# as the longest vector which occurs in the expression. 
# Shorter vectors in the expression are recycled as often as need be (perhaps fractionally)
# until they match the length of the longest vector. In particular a constant is simply repeated.
length(m_over_z)
length(charges)
m_over_z/charges
masses <- m_over_z/charges
charges * 2 

charges_bis<-c(charges,2,5)
length(charges_bis)
m_over_z
charges
m_over_z/charges
m_over_z/charges_bis
charges_bis


# Descriptive values
min(masses) 
max(masses) 
length(masses) 
mean(masses) 
median(masses) 
range(masses) 
sd(masses) 
var(masses) 
sum((masses - mean(masses))^2)/(length(masses)-1) 

masses
sort(masses) 
order(masses) 
sort(m_over_z)
order(m_over_z)

# Generating regular sequences
# R has a number of facilities for generating commonly used sequences of numbers
c(1:30) 
# repeat the number one, 30 times.
rep (1,times=30)
rep(1,30)
#  create five copies of "masses" end-to-end.
rep(masses,time=5)
# repeats each element of "masses" five times before moving on to the next
rep(masses,each=5)
rep(masses,,,5)
# by=value specify a step size
seq (from=1, to=30, by=0.5)
seq (1,30,0.5) 


# Different ways to read a file
my.data <- read.table("ExampleFile.csv")
my.data <- read.table("ExampleFile.csv",fill=T)
my.data <- read.table("ExampleFile.csv",fill=T,sep=",")
my.data <- read.table("ExampleFile.csv",fill=T,sep=",",header=T)

#read.table(file.choose(),header=T) 

# these two are equivalent
my.data <- read.table("ExampleFile.csv",fill=T,sep=",",header=T,comment.char="")
my.data2 <- read.csv("ExampleFile.csv")

# Access to the names of columns and rows
colnames(my.data)
rownames(my.data)
dim(my.data)

# Take a look of the data frame
View(my.data)
# There are some NA in the data set?
summary(my.data)
# whith "head" we can only see the first 6 rows
head(is.na(my.data))

# Categorical variables
my.data[2,3]
my.data[,3]
colnames(my.data)
my.data$X..Proteins
table(my.data$Activation.Type)


# Continuous variables
# Inspection of a column
my.data$Homology.Threshold
mean(my.data$Homology.Threshold)
mean_thre<-mean(my.data$Homology.Threshold,na.rm=TRUE)
median(my.data$Homology.Threshold,na.rm=TRUE)
sd(my.data$Homology.Threshold,na.rm=TRUE)
summary(my.data$Homology.Threshold,na.rm=TRUE) ##237 NA

length(which(is.na(my.data$Homology.Threshold)))

colnames(my.data)[which(is.na(my.data[3,]))]

nomi_colonne<-colnames(my.data)
which(is.na(my.data[3,]))
nomi_colonne[12]
nomi_colonne[33]

y <- my.data$Homology.Threshold[!is.na(my.data$Homology.Threshold)]
length(my.data$Homology.Threshold)
length(y)
length(my.data$Homology.Threshold)-length(y)
hist(my.data$Homology.Threshold)

# Another way to see the distribution of continuous variables
boxplot(my.data$X116.114,my.data$X115.114)
boxplot(my.data$X116.114,my.data$X115.114,ylim=c(0,12))

# Look at correlation
# Take two numeric columns
plot(my.data$X116.114,my.data$X115.114)
plot(my.data$X116.114,my.data$X115.114,xlim=c(0,20),ylim=c(0,20))
?abline
abline(a=1,b=1)
abline(a=1,b=1,col="red")
abline(h=5,col="blue")
abline(v=4,col="green")
# calculate correlation
cor(my.data$X116.114,my.data$X115.114)
cor(my.data$X116.114,my.data$X115.114,na.rm=TRUE)
# Why it doesn't work?
# search for the error type in R
# https://stackoverflow.com/questions/31412514/na-values-not-being-excluded-in-cor
cor(my.data$X116.114,my.data$X115.114,use="complete.obs")

# Mathematical operators to inspect data
# < minor
# <= monor or equal
# > above
# >= above or equal
# == equal
# != different
# & strict intersection (and)
# | not strict intersection (or)
a <- c (1:30)
a <= 15 

# Look at the graph we produced before and try to identify the points 
# above the horizontal line and at the right of the vertical line
# Identify the ID (the rows of the dataset of these samples)
which(my.data$X116.114>=4)
# Who are? Extract their name
my.data$Sequence[which(my.data$X116.114>=4)]
# REMEMBER: it is good to save information in an object in R
id_over4<-which(my.data$X116.114>=4)
# It os more simple then to obtain the information without loosing it
my.data$Sequence[id_over4]

# Horizontal line
id_over5<-which(my.data$X115.114>=5)
my.data$Sequence[id_over5]

# # above the horizontal line and at the right of the vertical line simuntaneously
id_over45<-which(my.data$X116.114>=4 & my.data$X115.114>=5)
my.data$Sequence[id_over45]

my.data$Annotation[which(my.data$X116.114>=4 & my.data$X115.114>=5)]<-"up"
summary(my.data$Annotation)
my.data$Annotation

# Save the information obtain
# The object are string of name, a sequence (vector) of words
# Look at the environment, we store these new object and they are "character" objects
first<-my.data$Sequence[id_over5]
second<-my.data$Sequence[id_over4]
first
second
venn(list(first.vector = first, second.vector = second))
# Error: could not find function "venn"
# We need to install the package
# If you don't have the 'gplots' package, type: install.packages("gplots")
library("gplots")
venn(list(first.vector = first, second.vector = second))
# names shared by the two vectors
intersect(first,second)
# In the visualization plot they are 7 in common, why the "intersect" function
# gives us only 6 name?
# Take a look of the two vector names, we have repeated name!!
# cATITPDEAR is repeated two times!
table(first)
table(second)

# Find the not shared Sequences
# Only two.... Repeated name!!
setdiff(first,second)
setdiff(second,first)

#Attention: some names seems to be very similar,
# R is case sensitive: cSSmQR        cSSMQR

# Some Types of graphs for data visualization
# The type of plots, to be created, depends on the format of your data. 
# The ggplot2 package provides methods for visualizing the following data structures:
# One variable - x: continuous or discrete
# Two variables - x & y: continuous and/or discrete
# Continuous bivariate distribution - x & y (both continuous)
# Continuous function
# Error bar
# Maps
# Three variables

#install.packages("ggplot2")
library(ggplot2)

#website for graph in R
#https://r-graph-gallery.com/index.html 

# Basic scatter plot
library(hrbrthemes)
ggplot(my.data, aes(x = X116.114, y = X115.114, color=Activation.Type)) + 
  geom_point(size=6) +
  theme_ipsum()
# Change the point size, and shape
ggplot(data = my.data, aes(x = X116.114, y = X115.114, color=Activation.Type)) +
  geom_point(size = 2, shape = 23)+xlim(0,12)+ylim(0,12)+
  theme_ipsum()

#geom_boxplot() for boxplot:
ggplot(my.data,aes(x=Activation.Type, y = Delta.M..ppm.))+
  geom_boxplot(col="#F8C070")+ #boxplot
  theme_light() #white background

# Possible layers are:
# For one continuous variable:
# - geom_area() for area plot
# - geom_density() for density plot
# - geom_histogram() for histogram plot
# - geom_smooth() for correlation
# For one discrete variable:
# - geom_bar() for bar plot
names(my.data)

# geom_area(): Create an area plot
ggplot(my.data, aes(x = Homology.Threshold)) + geom_area(stat = "bin")
ggplot(my.data, aes(x = IonScore)) + geom_area(stat = "bin")
# change fill colors by a group variable
table(my.data$Rank)
ggplot(my.data, aes(x = IonScore)) + geom_area(aes(fill = Rank),stat = "bin")
# Fix the errors
# add binwidth=10
# transform "Rank" into a factor
my.data$Rank_fact <- as.factor(my.data$Rank)
ggplot(my.data, aes(x = IonScore)) + geom_area(aes(fill = Rank_fact),stat = "bin",binwidth=10)
# Try another group variable
ggplot(my.data, aes(x = IonScore)) + geom_area(aes(fill = Activation.Type),stat = "bin",binwidth=10)
# smoothed colors and no gray background
# add "alpha" arguments for colors by group
# add the theme for the background
ggplot(my.data, aes(x = IonScore)) + 
  geom_area(aes(fill = Activation.Type),stat = "bin",binwidth=10,alpha=0.6) +
  theme_classic()

# geom_density(): Create a smooth density estimate
ggplot(my.data, aes(x = Homology.Threshold)) + geom_density()
# - geom_density() to create a density plot
# - geom_vline() to add a vertical lines corresponding to group mean values
# - scale_color_manual() to change the color manually by groups
# change line colors by Rank_fact
ggplot(my.data, aes(x = Homology.Threshold)) + geom_density(aes(color = Rank_fact)) 
# Change fill color by Rank_fact
# Use semi-transparent fill: alpha = 0.4
ggplot(my.data, aes(x = Homology.Threshold)) + geom_density(aes(fill = Rank_fact), alpha=0.4)

# Add mean line and Change color manually
# calculate the mean of "Homology.Threshold" by group of "Rank_fact"
table(my.data$Rank_fact)
mu1<-mean(my.data$Homology.Threshold[which(my.data$Rank_fact==1)],na.rm=T)
mu2<-mean(my.data$Homology.Threshold[which(my.data$Rank_fact==2)],na.rm=T)
mu3<-mean(my.data$Homology.Threshold[which(my.data$Rank_fact==3)],na.rm=T)
grp.mean<-c(mu1,mu2,mu3)
group<-as.factor(c(1,2,3))
data_mean<-data.frame(grp.mean,group)
ggplot(my.data, aes(x = Homology.Threshold)) + geom_density(aes(color = Rank_fact)) +
  geom_vline(data=data_mean, aes(xintercept=grp.mean, color=group),
             linetype="dashed") +
  scale_color_manual(values=c("blue", "red","green"))

# geom_histogram(): Histogram
# Basic plot
ggplot(my.data, aes(x = Homology.Threshold)) + geom_histogram()
ggplot(my.data, aes(x = Homology.Threshold)) + geom_histogram(bins=20)
# change line colors by sex
ggplot(my.data, aes(x = Homology.Threshold)) + geom_histogram(aes(color = Rank_fact), fill = "white", position = "dodge",bins=10) 
ggplot(my.data, aes(x = Homology.Threshold)) + geom_histogram(aes(color = Activation.Type), fill = "white", position = "dodge",bins=10) 

# stat_qq(): quantile - quantile plot
ggplot(my.data, aes(x = Homology.Threshold)) + stat_qq()
?stat_qq
# Learn how to solve Error is the 90% of the work in R
ggplot(my.data, aes(sample = Homology.Threshold)) + stat_qq()
ggplot(my.data, aes(sample = X116.114)) + stat_qq()


# One variable: Discrete
# The function geom_bar() can be used to visualize one discrete variable. 
# Basic plot
ggplot(my.data, aes(Rank_fact)) + geom_bar()
# Change color
ggplot(my.data, aes(Rank_fact)) + geom_bar(fill = "steelblue", color ="steelblue") + theme_minimal()

ggplot(my.data, aes(Rank_fact)) + stat_count(aes(fill = Activation.Type))


# geom_smooth(): Add regression line or smoothed conditional mean
# Scatterplot line only
ggplot(data = my.data, aes(x = X116.114, y = X115.114)) +
  geom_point(size = 2, shape = 23)+xlim(0,12)+ylim(0,12)
# Regression line only
ggplot(data = my.data, aes(x = X116.114, y = X115.114)) +
  geom_point(size = 2, shape = 23)+xlim(0,12)+ylim(0,12) + geom_smooth(method = lm)
# Point + regression line
# Remove the confidence interval 
ggplot(data = my.data, aes(x = X116.114, y = X115.114)) +
  geom_point(size = 2, shape = 23)+xlim(0,12)+ylim(0,12) +
  geom_point() + 
  geom_smooth(method = lm, se = FALSE)
# loess method: local regression fitting
ggplot(data = my.data, aes(x = X116.114, y = X115.114)) +
  geom_point(size = 2, shape = 23)+xlim(0,12)+ylim(0,12)+ 
  geom_point() + geom_smooth()
# Change color and shape by groups (Activation.Type)
p<-ggplot(data = my.data, aes(x = X116.114, y = X115.114)) +
  geom_point(size = 2, shape = 23)+xlim(0,12)+ylim(0,12) +
  geom_point(aes(color=Activation.Type, shape=Activation.Type)) + 
  geom_smooth(aes(color=Activation.Type, shape=Activation.Type), 
              method=lm, se=FALSE, fullrange=TRUE)


fig_dir<-getwd()
fig_dir
#to save image for journal
ggplot2::ggsave(filename = paste(fig_dir,"/", "Fig", ".jpeg", sep=""),
                plot = p,
                #device = cairo_ps,
                dpi = 1000,
                width = 17,
                height = 10,
                units = "cm")




# Some programming tips 
# if
if(1==1 & 2>1) {
  print("Yes, we made it")
} else {
  print("R is wrong")
}

# for
x <- seq(0,10,0.1)
y <- rep(NA,length(x))
for (i in x) {
  y[i] <- cos(i)
}
plot(x,y)
plot(x,cos(x))


## Lists, data frames
mylist <- list(mystr=c("A","B"), mynumbers=c(1:10,5), matrix(0,2,3))
mylist$mystr
mylist[["mystr"]]
mylist[[1]]
mylist[[3]]

## Workspace
ls()
rm(list=ls())

