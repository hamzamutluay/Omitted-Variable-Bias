# **Omitted Variable Bias**

Omitted variable bias occurs in linear regression analysis when one or more relevant independent variables are not included in regression model. 
For omitted variable bias to occur, explanatory variable must be correlated with the omitted variable and, the omitted variable must be a determinant of the dependent variable. As a result of this, the first linear regression model assumptions $E(U|X_{i})=0$ violates. In this situation, ordinary least squares produces biased and inconsistent estimates. 

Suppose that a correctly specified regression model would be

$$y = X_{1}β_{1} + X_{2}β_{2} + ε,$$

where the two parts of X have $K_{1}$ and $K_{2}$ columns, respectively. If we regress y on $X_{1}$ without including $X_{2}$, then the estimator is

$$b_{1} = (X'_{1}X_{1})^{-1}X'_{1}y = β_{1} + (X'_{1}X_{1})^{-1}X'_{1}X_{2}β_{2} + (X'_{1}X_{1})^{-1}X'_{1}ε.$$

Taking the expectation, we see that unless $X_1'X_2 = 0$ or $β_2 = 0$, $b_1$ is biased. The wellknown result is the omitted variable formula:

$$E[b_1|X] = β_1 + P_{1.2}β_{2},$$

where $$P_{1.2} = (X_1'X_1)^{-1}X_1'X_2$$

Each column of the $K_1$ x $K_2$ matrix $P_{1.2}$ is the column of slopes in the least squares regression of the corresponding column of $X_2$ on the columns of $X_1$.

## A Monte Carlo Study-Biased estimator
This do file below produces a Monte Carlo study of omitted variable	bias	in least squares

---
```stata

/*************************************************************************
A Monte Carlo Study of Omitted Variable Bias in Least Squares	   

Finite-Sample Properties of OLS	
							             
14-Jan-2023    		      
								             
Hamza Mutluay			      
  	 			 			   
comments/errors reports: hamzamutluay@gmail.com
****************************************************************************/

cls
clear 		*
se 			mo off, perm
pause 		off, perm


*- Directory --

cd 			"C:\Users\hamza.mutluay\Dropbox\Github\Github Part 4"

se 			se 		100113 // Set the seed to reproduce the same output of simulation

*- Sample size  ------

se obs 		100
loc 		n = _N
*----------------------


*- Define the population parameters -------

loc 		beta_0 = 1
loc 		beta_1 = 1
loc 		beta_2 = 1.5
*-------------------------------------------


loc 		rep	= 5000   // Number of replications


*- Parameters of the variance-covariance matrix of jointly distributed random variables  ----

loc 		sigma12 = 0.3
mat 		vector_mean = (0.5,0.5)

mat 		Covmatrix = (1,`sigma12'\ `sigma12',1)
drawnorm 	x1 x2 , n(`n') cov(Covmatrix) m(vector_mean)

*---------------------------------------------------------


*- Store the estimated beta in a matrix ---

mat 		betahat = J(`rep',2,.)

*-------------------------------------------


forv 		t = 1/`rep' {

g 			u = rnormal(0,1)
	
*- Data Generating Process

g 			y = `beta_0' + `beta_1'*x1 + `beta_2'*x2 + u

*- Estimated model

reg			y x1 

mat 		betahat[`t',1] = r(table)[1,2] 		// betahat_0
mat 		betahat[`t',2] = r(table)[1,1] 		// betahat_1

drop 		y 	u  


}

sca 		sum_betahat_0 = 0
sca 		sum_betahat_1 = 0


forv 		i = 1 / `rep' {

sca 		sum_betahat_0 = sum_betahat_0 + betahat[`i',1]
sca 		sum_betahat_1 = sum_betahat_1 + betahat[`i',2]

}

sca 		bias_betahat_0  = ((sum_betahat_0/`rep') - `beta_0')*100
sca 		bias_betahat_1  = ((sum_betahat_1/`rep') - `beta_1')*100

di 			" Bias of betahat_0 = " bias_betahat_0
di 			" Bias of betahat_1 = " bias_betahat_1

mat			colname betahat = betahat_0 betahat_1  
svmat       betahat, names(col)



*- Graph ---

#delimit 		;

hist 		betahat_1, normal xline(1, lwidth(thick)) scheme(white_piyg)  color(%30) 
					   ti("Sampling Distributions of Beta Estimates", size(*0.8) height(3)) ///
					   subtitle("Biased Estimator", size(*0.7)) 
					   xlab(0.5(0.5)1.5) xti(Estimates, height(3)) 
					   xlabel(0.5 1 "True Beta(`=ustrunescape("\u03b2")'=1)"  1.5 2) 
					   text(1.2 0.7 "E(`=ustrunescape("\u03b2")') `=ustrunescape("\u2260")' `=ustrunescape("\u03b2")'" ) 
					   text(1.24 0.682 "`=ustrunescape("\u0302")'")
;
#delimit cr  				   
   
					   
gr 			export 		biased_estimator.png, replace




```
---

The figure below from the simulation study shows that the OLS estimator is biased in the presence of omitted variable bias.
![combine](https://user-images.githubusercontent.com/101017847/212010487-1be5c8cc-e685-4eea-8c00-5c32a907565b.png)

## A Monte Carlo Study-Inconsistent estimator
This do file below produces a Monte Carlo study of omitted variable bias in least squares

---
```stata

/*************************************************************************
A Monte Carlo Study of Omitted Variable Bias in Least Squares	   

Large-Sample Properties of OLS	
							             
14-Jan-2023    		      
								             
Hamza Mutluay			      
  	 			 			   
comments/errors reports: hamzamutluay@gmail.com
****************************************************************************/

cls
clear 		*
se mo 		off, perm
pause 		off, perm


se se 		100215421  // Set the seed to reproduce the same output of simulation


*- Directory --

cd 			"C:\Users\hamza.mutluay\Dropbox\Github\Github Part 4"


*- Sample Sizes ----------------------------

numlist 	"100(900)10000" // Sample sizes
loc 		numlist `r(numlist)'
loc 		numlist :  subinstr loc numlist " " ",", all
mat			n_all = (`numlist')
loc         dim_col = colsof(n_all)
loc         dim_row = rowsof(n_all)
mat 		li n_all

*-------------------------------------------



*- Define the population parameters -------------------

loc 		beta0 = 2   
loc 		beta1 = 1    
loc 		beta2 = 1.5

*-------------------------------------------


loc 		rep = 2000 // Number of replications


*-------------------------------------------

*- Parameter of the variance-covariance matrix of jointly distributed random variables  

loc 		sigma12 = 0.3 

*-------------------------------------------


loc 		epsilon = 0.05 // Any arbitrarily small positive number


*- Store the estimated beta in a matrix ----------------------

mat 		coef_matrix  = J(`dim_col',3,.)
di  		`dim_row' ,`dim_col'

*-------------------------------------------


forv 		t = 1/`dim_col' {
	

se 			obs `=n_all[`dim_row',`t']'
loc 		n = _N
	
mat 		betahat = J(`rep',2,.)


mat 		vector_mean = (0.5,0.5)
mat 		Covmatrix = (1,`sigma12' \ `sigma12',1)

drawnorm 	x1 x2 , n(`n') cov(Covmatrix) m(vector_mean)

forv 		tt = 1/`rep' {
	

g			u = rnormal(0,1)

** Data Generating Process

g 			y = `beta0' + `beta1'*x1 + `beta2'*x2 + u


** Estimated Model

reg 		y x1 

mat 		betahat[`tt',1] = r(table)[1,2] 	// betahat_0
mat 		betahat[`tt',2] = r(table)[1,1] 	// betahat_1

drop 		y  u


	
}

sca 		sum_betahat_0 = 0
sca 		sum_betahat_1 = 0 


forv 		i = 1/`rep' {

sca 		sum_betahat_0 = sum_betahat_0 + (abs(betahat[`i',1] - `beta0') < `epsilon')
sca 		sum_betahat_1 = sum_betahat_1 + (abs(betahat[`i',2] - `beta1') < `epsilon')
}

mat 		coef_matrix[`t',1] = sum_betahat_0/`rep'
mat 		coef_matrix[`t',2] = sum_betahat_1/`rep'
mat 		coef_matrix[`t',3] = n_all[`dim_row',`t']	

mat drop 	vector_mean Covmatrix
drop 		x1 x2
sca drop	sum_betahat_0 sum_betahat_1
}


mat			colname coef_matrix = betahat_0 betahat_1  samplesize
mat li 		coef_matrix


svmat       coef_matrix, names(col)


*- Graph ---


#delimit 		;

scatter     betahat_1 samplesize,  
			ti("Sampling Distributions of Beta Estimates", size(*0.8) height(3))
			subti("Inconsistent Estimator", height(3) size(*0.8)) 
			ylab(0(0.20)1,nogrid) xla(100(900)10000) 
			yti("Probability", size(*0.9) height(8)) 
			xti("Sample size", size(*0.9) height(8)) 
			graphregion(col(white) lc(black) ) 
			plotregion(col(white)) scheme(white_piyg)


;
#delimit cr  				   

gr 			save 		inconsistent_estimator.gph, replace					   
gr 			export 		inconsistent_estimator.png, replace


```
---

The figure below from the simulation study shows that the OLS estimator is inconsistent in the presence of omitted variable bias.
![combine_large_sample](https://user-images.githubusercontent.com/101017847/212012080-52a53263-51c9-40ee-b100-8c192cff0356.png)


## A Monte Carlo Study with R codes-Inconsistent and biased estimator
```
{r Omitted Variabel  - Biased estimator}

 Define the population parameters
beta0 = 10
beta1 = 1
beta2 = 1.5

# Sample size
n = 100

#Number of replications
replication = 5000


sigma12 = 0.3 

alfahat0 = rep(0,replication)
alfahat1 = rep(0,replication)
x = mvrnorm(n, mu = c(0.5,0.5), Sigma=matrix(c(1,sigma12,sigma12,1),2,2), empirical = TRUE)

x1 = x[,1]
x2 = x[,2]

for (i in (1:replication)) {

epsilon = rnorm(n, mean = 0, sd = 1) 

# DGP
y = beta0 + beta1 * x1 + beta2 * x2 + epsilon
  
# Estimated Model 
model = lm(y ~ x1)
  
alfahat0[i] = summary.lm(model)$coefficients[1,1]
alfahat1[i] = summary.lm(model)$coefficients[2,1]

}

bias_alfa0 = (mean(alfahat0) - beta0)*100
bias_alfa1 = (mean(alfahat1) - beta1)*100

bias_alfa0
bias_alfa1
hist(alfahat1)

m   <- mean(alfahat1)
st  <- sqrt(var(alfahat1))
hist(alfahat1, prob = TRUE, ylim = c(0,4),  col="light coral")
curve(dnorm(x , mean = m, sd = st), col="light coral", lwd=3,
      add = TRUE) 
abline(v = 1, col="red", lwd=3, lty=1)

plotNormalHistogram(alfahat1, prob = TRUE, col="light coral",linecol="blue",border="red",
                    xlim=c(0.8,2), 
                    main = "Sampling Distributions of Beta Estimates-Biased estimator", 
                    xlab="Estimates")
abline(v = 1, col="blue", lwd=3, lty=2)


```

```
{r Omitted Variabel  - Inconsistent estimator}

# Define the population parameters
beta0   = 10
beta1   = 1
beta2   = 1.5


sigma12     = 0.3

replication = 2000

epsilon =  0.05 #Any arbitrarily small positive number


n_all = seq(100, 10000, by = 900) #Sample Sizes

prob1 = rep(0, length(n_all))
prob2 = rep(0, length(n_all))

for (i in (1:length(n_all))) {
  
n = n_all[i]

alfahat0 = rep(0,replication)
alfahat1 = rep(0,replication)

x = mvrnorm(n, mu = c(0.5,0.5), Sigma = matrix(c(1,sigma12,sigma12,1),2,2), empirical = TRUE)

x1 = x[,1]
x2 = x[,2]

for (j in (1:replication)) {

error_term = rnorm(n, mean = 0, sd = 1) 
  
# Data Generating Process
y = beta0 + beta1 * x1 + beta2 * x2 + error_term
    
# Estimated Model 
model = lm(y ~ x1)

alfahat0[j] = summary.lm(model)$coefficients[1,1]
alfahat1[j] = summary.lm(model)$coefficients[2,1]
    
}
  
prob1[i] = mean(abs(alfahat0 - beta0) < epsilon)
prob2[i] = mean(abs(alfahat1 - beta1) < epsilon)
  
}

prob1
prob2

x = seq(100,10000, 900)
plot(y = prob2, x, ylim=c(0,1) ,xlim=c(100,10000), pch=19,
     main="Inconsistent estimator", xlab="Sample size" ,
     ylab="Probability ")


  
```
