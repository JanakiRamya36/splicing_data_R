## Loading required Libraries
library(tidyverse)
library(pastecs)
library(skimr)
library(corrplot)
library(ggcorrplot)
library(ggplot2)
library(MASS)

### DATA EXPLORATION ###

## Reading the data
splicing_data<-read.csv("C:/Users/DELL/OneDrive - Rutgers University/Desktop/R/Project/splicing_data.csv",header=TRUE)

## SUMMARY STATISTICS

# Using glimpse function to get the basic insights of data
glimpse(splicing_data)  

# Using describe function to get the summary statistics
stat.desc(cbind(splicing_data$SplicingFactor1,splicing_data$SplicingFactor2,splicing_data$SplicingFactor3,splicing_data$SplicingEvent),basic=FALSE,norm=TRUE)  

# Using skim function to get the data summary
skim(splicing_data) 

## CHECKING FOR MISSING VALUES AND OUTLIERS

# Checking missing values
splicing_data %>% map(~sum(is.na(.))) 

# Checking the data type of variables
sapply(splicing_data,typeof)

# Converting the datatype of Subject Id
splicing_data$SubjectID<-as.character(splicing_data$SubjectID)

# Checking Outliers
boxplot(splicing_data[,-1])
# Outliers are noticed for SplicingFactor1 and SplicingFactor2 variables

# Removing outliers from both the variables
splicing_data2<-splicing_data[!splicing_data$SplicingFactor1 %in% boxplot.stats(splicing_data$SplicingFactor1)$out,]
splicing_data2<-splicing_data2[!splicing_data2$SplicingFactor2 %in% boxplot.stats(splicing_data2$SplicingFactor2)$out,]

# Verifying the removal of outliers
boxplot(splicing_data2[,-1])

## CORRELATION MATRIX

# Pair-wise correlations between splicing factors and splicing event
cormatrix<-cor(splicing_data2[,-1],use="complete.obs",method="pearson")
cormatrix

# Creating a corr plot
corrplot(cor(splicing_data2[,-1]))

# Creating a heat map
ggcorrplot(cor(splicing_data2[,-1]),title="Correlation matrix",colors=c("lightyellow","orange","darkblue"))

### DATA VISUALIZATION ###
## Data Distribution
## Using Histogram to check the distribution of data
# Splicing Factor 1
hist_fac1<-ggplot(splicing_data2,aes(SplicingFactor1))+
  geom_histogram(aes(y = ..density..),binwidth=0.4,colour="darkblue",fill="aquamarine")+
  labs(x="Splicing Factor 1",y="Density")
hist_fac1+stat_function(fun=dnorm,args=list(mean=mean(splicing_data2$SplicingFactor1,na.rm=TRUE),
  sd=sd(splicing_data2$SplicingFactor1,na.rm=TRUE)),colour="chocolate",size=1)

# Splicing Factor 2 
hist_fac2<-ggplot(splicing_data2,aes(SplicingFactor2))+
  geom_histogram(aes(y = ..density..),binwidth=0.6,colour="cornsilk2",fill="darkmagenta")+
  labs(x="Splicing Factor 2",y="Density")
hist_fac2+stat_function(fun=dnorm,args=list(mean=mean(splicing_data2$SplicingFactor2,na.rm=TRUE),
  sd=sd(splicing_data2$SplicingFactor2,na.rm=TRUE)),colour="coral",size=1)

# Splicing Factor 3 
hist_fac3<-ggplot(splicing_data2,aes(SplicingFactor3))+
  geom_histogram(aes(y = ..density..),binwidth=0.2,colour="floralwhite",fill="dodgerblue4")+
  labs(x="Splicing Factor 3",y="Density")
hist_fac3+stat_function(fun=dnorm,args=list(mean=mean(splicing_data2$SplicingFactor3,na.rm=TRUE),
  sd=sd(splicing_data2$SplicingFactor3,na.rm=TRUE)),colour="greenyellow",size=1)

# Splicing Event
hist_eve<-ggplot(splicing_data2,aes(SplicingEvent))+
  geom_histogram(aes(y = ..density..),binwidth=0.8,colour="deeppink4",fill="gold1")+
  labs(x="Splicing Event",y="Density")
hist_eve+stat_function(fun=dnorm,args=list(mean=mean(splicing_data2$SplicingEvent,na.rm=TRUE),
  sd=sd(splicing_data2$SplicingEvent,na.rm=TRUE)),colour="darkslategray",size=1)

## Using q-q plots to check the distribution
# Splicing Factor 1
qq_fac1<-ggplot(splicing_data2,aes(sample=SplicingFactor1))+stat_qq()+stat_qq_line()+
  labs(x="Splicing Factor 1",y="Density estimate")
qq_fac1

# Splicing Factor 2
qq_fac2<-ggplot(splicing_data2,aes(sample=SplicingFactor2))+stat_qq()+stat_qq_line()+
  labs(x="Splicing Factor 2",y="Density estimate")
qq_fac2

# Splicing Factor 3
qq_fac3<-ggplot(splicing_data2,aes(sample=SplicingFactor3))+stat_qq()+stat_qq_line()+
  labs(x="Splicing Factor 3",y="Density estimate")
qq_fac3

# Splicing Event
qq_eve<-ggplot(splicing_data2,aes(sample=SplicingEvent))+stat_qq()+stat_qq_line()+
  labs(x="Splicing Event",y="Density estimate")
qq_eve

## Checking the relationships between Splicing Factors and Splicing Event using Scatter plots
# Splicing Factor 1
scatter_fac1<-ggplot(splicing_data2,aes(SplicingFactor1,SplicingEvent))
scatter_fac1+geom_point()+geom_smooth(method="lm",colour="darkblue")+labs(x="Splicing Factor 1",y="Splicing Event")

# Splicing Factor 2
scatter_fac2<-ggplot(splicing_data2,aes(SplicingFactor2,SplicingEvent))
scatter_fac2+geom_point()+geom_smooth(method="lm",colour="darkmagenta")+labs(x="Splicing Factor 2",y="Splicing Event")

# Splicing Factor 3
scatter_fac3<-ggplot(splicing_data2,aes(SplicingFactor3,SplicingEvent))
scatter_fac3+geom_point()+geom_smooth(method="lm",colour="firebrick")+labs(x="Splicing Factor 3",y="Splicing Event")





### PREDICTIVE MODELLING ###

## Splitting the data set
# As the target variable is having a balanced distribution, random sampling can be done 
set.seed(6)
# Using sample function, split the training indices
train_indices<-sample(1:nrow(splicing_data2),size=0.7*nrow(splicing_data2),replace =FALSE,prob=NULL)

# Create training and testing data
train_data<-splicing_data2 %>% slice(train_indices)
test_data<-splicing_data2 %>% slice(-train_indices)

# Checking the number of rows
nrow(train_data)
nrow(test_data)

## Multiple Linear Regression model
eventmodel<-lm(SplicingEvent~SplicingFactor1+SplicingFactor2+SplicingFactor3,data=splicing_data2)
summary(eventmodel)
plot(eventmodel)

# Mean squared error
mse<-mean(eventmodel$residuals^2)
mse

# Root mean squared error
rmse<-sqrt(mse)
rmse

# R2
summary(eventmodel)$r.squared

# Prediction of test data
predictval<-predict(eventmodel,newdata=test_data,interval="confidence")
predictval

# Adding predicted values to the test data as column
test_data<- test_data%>% mutate(predval=predict(eventmodel,newdata=test_data))
head(test_data)
# Calculating residuals
test_data<- test_data%>% mutate(residuals=SplicingEvent-predval)

# Cross-validation correlation and Calculation of R2
test_data %>% summarize(cor=cor(SplicingEvent,predval)) %>%
  mutate(R2=cor^2)

# Calculating difference between the R2
test_data %>% summarize(cor=cor(SplicingEvent,predval)) %>%
  mutate(R2=cor^2, shrinkage=summary(eventmodel)$r.squared-R2)

# Calculating the mean of residuals and mean squared error
test_data %>% summarize(mean=mean(residuals),mse=mean(test_data$residuals^2))

# Residual plot for test data

ggplot(test_data,aes(x=predval,y=residuals))+
  geom_point()+geom_smooth(se=FALSE)+labs(x="Predicted values",y="residuals")

