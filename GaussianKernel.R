# Copyright (C) 2014 Abhik Seal <abhik1368@gmail.com>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details. 



###################################################################################################################
# This code is created based on the paper (http://bioinformatics.oxfordjournals.org/content/27/21/3036.long)       #
# define a kernel on these profiles, called the Gaussian Interaction Profile (GIP) kernel,                        #   
# and use a simple classifier, (kernel) Regularized Least Squares (RLS),for prediction drug–target interactions.  #
###################################################################################################################

###########################################################################
# Read in the files drug similarity matrix,drug target interactions, and  #
# sequence similarity matrices for calculations.                          #
###########################################################################
RLS <- function(cmp_sim,dti,seqsim){
drugs<-as.matrix(read.csv(cmp_sim,header=T,row.names=1)) ## Compound Similarity
dti<-as.matrix(read.csv(dti,header=T,row.names=1))  # Drug Target adjancency matrix
seqs<-as.matrix(read.csv(seqsim,header=T,row.names=1)) # Normalized similarity matrix 

#########################################################################
# gaussian Internaction profile matrix of drug target interactions      #
# K_GIP<-exp(-sigma*D^2) here D is the distance matrix                  #    
#########################################################################
## Setting the parameter sigma 
t<-as.matrix(rowSums(dti))
d<-as.matrix(rowSums(t(dti)))
s<-0
p<-0
for (i in 1:length(t)){
  s<-sum(s+t[i]/length(t))
}
for (i in 1:length(d)){
  p<-sum(p+d[i]/length(d)
}

Target_Dist<-as.matrix(dist(dti))
K_GIP_Targets<-exp(-(1/s)*Target_Dist^2)

Drug_Dist<-as.matrix(dist(t(dti)))
K_GIP_Drugs<-exp(-(1/p)*Drug_Dist^2)

################################################################
# Final Kernel matrix of drug similarity,target simialrity and #
# the Gaussian interaction profiles                            #
################################################################

Kd<- 0.5*(drugs)+0.5*(K_GIP_Drugs)
Kt<- 0.5*(seqs)+0.5*(K_GIP_Targets)

#################################################################################
# Regularized least squares Y_hat<-K%*%(K+sigmaI)^-1%*%Y                        #  
# K is the kernel matrix for drugs and sequence simialrity                      #
#################################################################################
Y<-dti
d<-solve(Kd+1*diag(nrow(Kd)))
t<-solve(Kt+1*diag(nrow(Kt)))
Yhat<-1/2*(Kd%*%d%*%t(Y))+1/2*t(Kt%*%t%*%Y)
#write.csv(Yhat,"datapredict.csv") # saves the data and also return Yhat (predicted matrix)
return(Yhat)
}




