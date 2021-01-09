# STAT351 FINAL PROJECT YUZHUO KANG
#
#
#### Question 1
#For this question, I use both Kruskal-Wallis test statistics for K-sample comparison and Kendall coefficient of concordance to test the relationship between egg lengths and the host bird. Using Kruskal-Wallis test statistics, $M_i$ is the median of each sample.$H_0$ is $M_1=M_2=...=M_k$ and $H_1$ is at least two of the $M_j$’s differ

library(DescTools);library(readxl)
data1=read_excel("351_problem1_data.xlsx") # load the data
bird1=data1$MDW_PIPIT;bird2=data1$TREE_PIPIT;bird3=data1$HDGE_SPRW
bird4=data1$ROBIN;bird5=data1$PIED_WTAIL;bird6=data1$WREN
bird=c(bird1,bird2,bird3,bird4,bird5,bird6);bird=bird[!is.na(bird)]
g=factor(rep(1:ncol(data1), c(sum(!is.na(data1$MDW_PIPIT)), sum(!is.na(data1$TREE_PIPIT)),sum(!is.na(data1$HDGE_SPRW)),sum(!is.na(data1$ROBIN)),sum(!is.na(data1$PIED_WTAIL)),sum(!is.na(data1$WREN)))))
kruskal.test(bird,g) # Conduct the Kruskal-Wallis rank sum test
```
#The test statistics of Kruskal-Wallis rank sum test is 35.04, which is a chi-squared value because the sample is large. The asymptotic p-value is 1.477e-06, which is very close to 0, so the null hypothesis is rejected. The medians of egg lengths of six birds are not equal.

#Using Kendall coefficient of concordance and Friedman test, $H_0$ is the length of the eggs and the host birds don't have association and $H_1$ is The association exists between the host bird and the length of eggs.
KendallW(data1, correct = TRUE, test = TRUE, na.rm = TRUE)

#The test statistics of Kendall coefficient of concordance and Friedman test Q equals 76.069. W value in the test result represents the association between 6 groups. The asymptotic p-value is 6.016e-11, so it is strictly less than 0.05. The null hypothesis is rejected. We conclude that there exists a relationship between the length of eggs and the host birds. 

#In conclusion, there exists a relationship between the length of eggs and the host birds.

### Question 2
###(a)
data2=read.csv("351_problem2_data.csv",header=TRUE,sep="")
ozone_y=data2$ozone # Response variable
radiation_x1=data2$radiation;temperature_x2=data2$temperature
wind_x3=data2$wind
n=length(ozone_y) # sample size
# Sort the data ascendingly
x1=sort(radiation_x1);x2=sort(temperature_x2);x3=sort(wind_x3)
y1=ozone_y[order(radiation_x1)];y2=ozone_y[order(temperature_x2)]
y3=ozone_y[order(wind_x3)]
# Use the local constant estimate and cross-validation CV to determine the bandwith 
vec_h= seq(12,80, by=0.1);CV=rep(0,length(vec_h))
for ( k in 1:length(vec_h)) {
h=vec_h[k];D_hat_m = rep(0, n)
for ( j in 1:n) {
D_X=x1[-j]; D_Y=y1[-j];x=x1[j]
Kcv=(1-((D_X-x)/h)^2)*(abs((D_X-x)/h)<=1)
wcv=Kcv/sum(Kcv);D_hat_m[j]=sum(wcv*D_Y) # delete j estimate at Xj
}
CV[k]=mean((y1-D_hat_m)^2)
}
(h_CV=vec_h[which(CV==min(CV))]) # the bandwith with smallest error
par(mfrow=c(1,2))
# Make a plot of CV and bandwidth
plot(vec_h,CV,type="l",main="Cross-validation on the bandwidth")
abline(v=h_CV,col="red",xlab="bandwith")
hat_m_1=rep(0,n); h=h_CV # using bandwith from CV
for ( i in 1:n) {# Conduct local constant estimate
x=x1[i]
K1=(1-((x1-x)/h)^2)*(abs((x1-x)/h)<=1)
w=K1/sum(K1);hat_m_1[i]=sum(w*y1)
}
plot(x1,y1,main="Local constant estimate between ozone and radiation",pch=16,cex=1.2,col="purple",xlab="radiation",ylab="ozone")
lines(x1,hat_m_1,col="red")


###(b)

vec_h= seq(5,20, by=0.5);CV=rep(0,length(vec_h))
for ( k in 1:length(vec_h)) {# Use the CV to determine the bandwith
h=vec_h[k];D_hat_m = rep(0, n)
for ( j in 1:n) {
D_X=x2[-j];D_Y=y2[-j];x=x2[j]
Kcv=(1-((D_X-x)/h)^2)*(abs((D_X-x)/h)<=1)
wcv=Kcv/sum(Kcv);D_hat_m[j]=sum(wcv*D_Y) # delete j estimate at Xj
}
CV[k]=mean((y2-D_hat_m)^2)
}
(h_CV=vec_h[which(CV==min(CV))]) # the optimal CV
par(mfrow=c(1,2))
# Make a plot of CV and bandwidth
plot(vec_h,CV,type="l",main="Cross-validation on the bandwidth")
abline(v=h_CV,col="red",xlab="bandwith")
hat_m_2=rep(0,n);h=h_CV
for ( i in 1:n) {# Conduct local constant estimate
x=x2[i]
K2=(1-((x2-x)/h)^2)*(abs((x2-x)/h)<=1)
w=K2/sum(K2);hat_m_2[i]=sum(w*y2)
}
plot(x2,y2,main="Local constant estimate between ozone and temperature",pch=16,cex=1.2,col="orange",xlab="temperature",ylab="ozone")
lines(x2,hat_m_2,col="red")


###(c)

vec_h= seq(2,10, by=0.1); CV=rep(0,length(vec_h))
for ( k in 1:length(vec_h)) {# Use the CV to determine the bandwith
h=vec_h[k];D_hat_m = rep(0, n)
for ( j in 1:n) {
D_X=x3[-j];D_Y=y3[-j];x=x3[j]
Kcv=(1-((D_X-x)/h)^2)*(abs((D_X-x)/h)<=1)
wcv=Kcv/sum(Kcv);D_hat_m[j]=sum(wcv*D_Y) # delete j estimate at Xj
}
CV[k]=mean((y3-D_hat_m)^2)
}
(h_CV=vec_h[which(CV==min(CV))]) # optimal bandwidth
par(mfrow=c(1,2))
# Make a plot of CV and bandwidth
plot(vec_h,CV,type="l",main="Cross-validation on the bandwidth")
abline(v=h_CV,col="red",xlab="bandwith")
hat_m_3=rep(0,n);h=h_CV
for ( i in 1:n) {# Conduct local constant estimate
x=x3[i]
K3=(1-((x3-x)/h)^2)*(abs((x3-x)/h)<=1)
w=K3/sum(K3); hat_m_3[i]=sum(w*y3)
}
plot(x3,y3,main="Local constant estimate between ozone and wind",pch=16,cex=1.2,col="green",xlab="wind",ylab="ozone")
lines(x3,hat_m_3,col="red")


### (d)
#The kernel regression estimation method for these three parts is local constant estimate, and uses Epanechnikov kernel. Kernel regression fit reflects the relationship between ozone and three variables. In contrast with temperature and wind, the relationship between the ozone and radiation may not be linear. From part(a), the optimal bandwidth using the Cross-Validation(CV) is h=51. The scatterplot shows when the radition is between 150 to 200, the ozone level is large. The kernel regression fitted value shows ozone changes significantly as radiation increases in this interval. When the radiation exceeds 200, the association is negative. For part(b), the optimal bandwidth using CV is h=7.The scatterplot shows when the temperature increases, the ozone also increases.The kernel regression fitted value shows there exists a positive linear relationship between the temperature and ozone. For part(c), the optimal bandwith using CV is h=2.9. The scatterplot shows when the wind speed increases, the ozone level declines. The kernel regression fitted value shows there exists a nagative linear relationship between the wind and ozone. 

#The trends shown in kernel regression fit between the ozone and each variable follow the general trend in the scatterplot. It shows that local constant estimate could properly reflect the relationship between the ozone and three variables. 

### Question 3
### Uniform distribution

set.seed(2020); n=100 # set the seed and sample size
X_unif=runif(n,0,1);unif_rank=rank(X_unif,ties.method="average")
up1=cor(X_unif,unif_rank, method="pearson")  # Part a
up2=cor(X_unif[1:(n-1)],unif_rank[2:n], method="pearson") # Part b
up3=cor(X_unif[2:n],unif_rank[1:(n-1)], method="pearson") # Part c
cat("For uniform distribution, the Pearson sample correlation 
coefficient for question a, b, c are", up1,up2,up3)

### Exponential distribution

n=100;X_exp=rexp(n,1);exp_rank=rank(X_exp,ties.method="average")  
ep1=cor(X_exp,exp_rank, method="pearson") # Part a
ep2=cor(X_exp[1:(n-1)],exp_rank[2:n], method="pearson")# Part b
ep3=cor(X_exp[2:n],exp_rank[1:(n-1)], method="pearson")# Part c
cat("For exponential distribution, the Pearson sample correlation 
coefficient for question a, b, c are", ep1,ep2,ep3)

### Normal distribution

n=100;X_norm=rnorm(n,0,1);norm_rank=rank(X_norm,ties.method="average")  
np1=cor(X_norm,norm_rank, method="pearson")  # Part a
np2=cor(X_norm[1:(n-1)],norm_rank[2:n], method="pearson")# Part b
np3=cor(X_norm[2:n],norm_rank[1:(n-1)], method="pearson")# Part c
cat("For normal distribution, the Pearson sample correlation 
coefficient for question a, b, c are", np1,np2,np3)


#For the uniform, exponential, and normal distribution, the variate values and their corresponding ranks have strong positive associations because the pearson "product-moment" sample correlation coefficient $\hat{\rho_1}$ is very close to 1.Their data points are nearly on the line of best fit. Also, the result shows that the variate values and the rank value of the adjacent two variates have very little linear association. The pearson "product-moment" sample correlation coefficient with the subsequent and previous variate's rank $\hat{\rho_2}$ and $\hat{\rho_3}$ are all close to 0. The variate value changes very little when the rank of adjacent variate changes.The variates following uniform, exponential, and normal distribution in the sample are indepedent.The variate following normal distribution have negative association with the the rank of adjacent variables, but that association is positive for the uniform and exponential distribution.

### Question 4

y=c(0,2,-1,4)
xpr=seq(1,4,length=100)
x=1:4
spline1 = smooth.spline(x,y,lambda = 0)
spline2 = smooth.spline(x,y,lambda = 1)
spline3 = smooth.spline(x,y,lambda = 1000)
par(mfrow=c(1,2))
plot(x,y,ylim=c(-5,15),main="lambda=0",xlab="t",pch=16,cex=2)
lines(xpr, predict(spline1, xpr)$y, col=2)
plot(x,y,ylim=c(-5,15),main="lambda=1",xlab="t",pch=16,cex=2)
lines(xpr, predict(spline2, xpr)$y, col=3)
par(mfrow=c(1,1))
plot(x,y,ylim=c(-5,15),main="lambda=1000",xlab="t",pch=16,cex=2)
lines(xpr, predict(spline3, xpr)$y, col=4)




