setwd('')
mycsv='output_7.csv'
library(survey)
mydata = read.csv(mycsv)

#Define agreement column
mydata$correct = 0
mydata$correct[mydata$Strata == mydata$Reference] = 1
mydata$correct[mydata$Strata != mydata$Reference] = 0
mydata$Weights = (1 / mydata$Inclu_Fin)
allclasses = unique(mydata$Strata)


#Define strata indicators
#class 1
mydata$class1 = 0
mydata$class2 = 0
mydata$class3 = 0
mydata$class4 = 0
mydata$class5 = 0
mydata$class6 = 0

mydata$class1[mydata$Strata == 1] = 1
mydata$class2[mydata$Strata == 2] = 1
mydata$class3[mydata$Strata == 3] = 1
mydata$class4[mydata$Strata == 4] = 1
mydata$class5[mydata$Strata == 5] = 1
mydata$class6[mydata$Strata == 6] = 1


mydata$class1ref = 0
mydata$class2ref = 0
mydata$class3ref = 0
mydata$class4ref = 0
mydata$class5ref = 0
mydata$class6ref = 0

mydata$class1ref[mydata$Reference == 1 ] = 1
mydata$class2ref[mydata$Reference == 2 ] = 1
mydata$class3ref[mydata$Reference == 3 ] = 1
mydata$class4ref[mydata$Reference == 4 ] = 1
mydata$class5ref[mydata$Reference == 5 ] = 1
mydata$class6ref[mydata$Reference == 6 ] = 1

mydata$class1cor = 0
mydata$class2cor = 0
mydata$class3cor = 0
mydata$class4cor = 0
mydata$class5cor = 0
mydata$class6cor = 0

mydata$class1cor[mydata$Strata == 1 & mydata$Reference == 1] = 1
mydata$class2cor[mydata$Strata == 2 & mydata$Reference == 2] = 1
mydata$class3cor[mydata$Strata == 3 & mydata$Reference == 3] = 1
mydata$class4cor[mydata$Strata == 4 & mydata$Reference == 4] = 1
mydata$class5cor[mydata$Strata == 5 & mydata$Reference == 5] = 1
mydata$class6cor[mydata$Strata == 6 & mydata$Reference == 6] = 1


#Need sequential units for each PSU
mydata$PSU_ID <- 0
unique_psu <- unique(mydata$Tile)
for (i in 1:length(unique_psu)){
  mydata$PSU_ID[mydata$Tile == unique_psu[i]] <- i
}

#Get inclusion probabilities and weights 
#Note: This is done in the sampling scripts, but can be re-done here if samples are removed. 
unique_classes <- unique(mydata$Strata)
first_1_strata <- length(unique(mydata$Tile[mydata$Strata1 == 1]))
first_2_strata <- length(unique(mydata$Tile[mydata$Strata1 == 2]))

first_1_pop <- mydata$TilePop[mydata$Strata1 == 1 ][1]
first_2_pop <- mydata$TilePop[mydata$Strata1 == 2 ][1]

#First Stage
mydata$Inclu_1[mydata$Strata1 == 1] = first_1_strata / first_1_pop
mydata$Inclu_1[mydata$Strata1 == 2] = first_2_strata / first_2_pop

for (i in unique_classes) {
  
  sample_units <- length(mydata$Strata[mydata$Strata == i])
  sample_pop <- mydata$TotalPix[mydata$Strata == i][1]
  mydata$Inclu_2b[mydata$Strata == i] = sample_units / sample_pop

}

#Weights
mydata$Weights <- 1 / (mydata$Inclu_1 * mydata$Inclu_2)


#Set up the ratio estimator design
mydesign <-
  svydesign(
              id= ~PSU_ID + ID,
              strata = ~Strata1+Strata,
              data = mydata,
              weights = ~Weights ,
              fpc = ~TilePop+TotalPix)


##Accuracy metrics. User's and Producer's accuracy are using equations from (Cochran, 1977 Section 6.11,) and (Wickham, et al. 2013)
#This can be equivelantly produced using the error matrix 

#Users accuracy
options(survey.lonely.psu = "adjust")
u_c1 <- svyratio(~class1cor, ~class1, mydesign)
u_c2 <- svyratio(~class2cor, ~class2, mydesign)
u_c3 <- svyratio(~class3cor, ~class3, mydesign)
u_c4 <- svyratio(~class4cor, ~class4, mydesign)
u_c5 <- svyratio(~class5cor, ~class5, mydesign)
u_c6 <- svyratio(~class6cor, ~class6, mydesign)

#Producer's accuracy
p_c1 <- svyratio(~class1cor, ~class1ref, mydesign)
p_c2 <- svyratio(~class2cor, ~class2ref, mydesign)
p_c3 <- svyratio(~class3cor, ~class3ref, mydesign)
p_c4 <- svyratio(~class4cor, ~class4ref, mydesign)
p_c5 <- svyratio(~class5cor, ~class5ref, mydesign)
p_c6 <- svyratio(~class6cor, ~class6ref, mydesign)

#Overall accuracy
#Create dummy variable for the denominator
mydata$dummy <- 1
ov_call <- svyratio(~correct, ~dummy, mydesign)

cat('Users/Producers Accuracy with Standard Errors in () for Class 1:' ,coef(u_c1),'(',SE(u_c1),')',',',coef(p_c1), '(',SE(p_c1),')')
cat('Users/Producers Accuracy with Standard Errors in () for Class 2:' ,coef(u_c2),'(',SE(u_c2),')',',',coef(p_c2), '(',SE(p_c2),')')
cat('Users/Producers Accuracy with Standard Errors in () for Class 3:' ,coef(u_c3),'(',SE(u_c3),')',',',coef(p_c3), '(',SE(p_c3),')')
cat('Users/Producers Accuracy with Standard Errors in () for Class 4:' ,coef(u_c4),'(',SE(u_c4),')',',',coef(p_c4), '(',SE(p_c4),')')
cat('Users/Producers Accuracy with Standard Errors in () for Class 5:' ,coef(u_c5),'(',SE(u_c5),')',',',coef(p_c5), '(',SE(p_c5),')')
cat('Users/Producers Accuracy with Standard Errors in () for Class 6:' ,coef(u_c6),'(',SE(u_c6),')',',',coef(p_c6), '(',SE(p_c6),')')

#Error matrix
#Todo: Get this automatically
class1_pix <-                
class2_pix <- 
class3_pix <- 
class4_pix <- 
class5_pix <- 
class6_pix <- 1

total_pix = sum(class1_pix,class2_pix,class3_pix,class4_pix,class5_pix,class6_pix)

#Area calculation and error matrix

for (row in seq_along(allclasses)){
  for (col in seq_along(allclasses)){
    weights <- sum(mydata$Weights[mydata$Strata == row & mydata$Reference == col])
    class_prop <- weights / total_pix
    conf_max[row, col] <- class_prop
  }
}

#Show the matrix
conf_max



area_c1 <- svytotal(~class1ref,mydesign)
area_c2 <- svytotal(~class2ref,mydesign)
area_c3 <- svytotal(~class3ref,mydesign)
area_c4 <- svytotal(~class4ref,mydesign)
area_c5 <- svytotal(~class5ref,mydesign)
area_c6 <- svytotal(~class6ref,mydesign)

cat('Estimated Area (Pixels) For Class 1 (SE TODO)' ,area_c1)
cat('Estimated Area (Pixels) For Class 2 (SE TODO)' ,area_c2)
cat('Estimated Area (Pixels) For Class 3 (SE TODO)' ,area_c3)
cat('Estimated Area (Pixels) For Class 4 (SE TODO)' ,area_c4)
cat('Estimated Area (Pixels) For Class 5 (SE TODO)' ,area_c5)
cat('Estimated Area (Pixels) For Class 6 (SE TODO)' ,area_c6)


