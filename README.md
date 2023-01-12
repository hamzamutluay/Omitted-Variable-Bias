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
							             
12-Jan-2023    		      
								             
Hamza Mutluay			      
  	 			 			   
comments/errors reports: hamzamutluay@gmail.com
****************************************************************************/


cls
clear 	 *
se 			 mo off, perm
pause 	 off, perm

se 			 se 	1990123  // Set the seed to reproduce the same output of simulation

*- Directory --

cd 			 "C:\Users\hamza.mutluay\Dropbox\Github\Github Part 4"


*- Sample size  ------

se obs 		100
loc 		  n = _N
*----------------------


*- Define the population parameters -------

loc 		  beta_0 = 1
loc 		  beta_1 = 1
loc 		  beta_2 = 1.5
*-------------------------------------------


loc 		  rep	= 5000   // Number of replications


*- Parameters of the variance-covariance matrix of jointly distributed random variables  ----

loc 		  sigma12 = 0.4  // Cov(x1,x2) = 0.4
loc 		  sigma13 = 0   // Cov(x1,u) = 0
loc 		  sigma23 = 0 // Cov(x2,u) = 0

*---------------------------------------------------------



*- Store the estimated beta in a matrix ---

mat 		  betahat = J(`rep',2,.)

*-------------------------------------------


forv 		t = 1/`rep' {

mat 		  vector_mean = (0.5,0.5,0)
mat 		  Covmatrix = (1,`sigma12',`sigma13' \ `sigma12',1, `sigma23' \ `sigma13',`sigma23',1)

	
drawnorm 	x1 x2 u, n(`n') cov(Covmatrix) m(vector_mean)
	

*- Data Generating Process

g 	      y = `beta_0' + `beta_1'*x1 + `beta_2'*x2 + u

*- Estimated model

reg	      y x1 

mat 		  betahat[`t',1] = r(table)[1,2] 		// betahat_0
mat 		  betahat[`t',2] = r(table)[1,1] 		// betahat_1



drop 		  y 	u  x1 x2 
mat drop  vector_mean Covmatrix


}



sca 		 sum_betahat_0 = 0
sca 		 sum_betahat_1 = 0



forv 		i = 1 / `rep' {

sca 		 sum_betahat_0 = sum_betahat_0 + betahat[`i',1]
sca 		 sum_betahat_1 = sum_betahat_1 + betahat[`i',2]

}

sca 		 bias_betahat_0  = ((sum_betahat_0/`rep') - `beta_0')*100
sca 		 bias_betahat_1  = ((sum_betahat_1/`rep') - `beta_1')*100

di 			 "Bias of betahat_0 =" bias_betahat_0
di 			 "Bias of betahat_1 =" bias_betahat_1

mat			 colname betahat = betahat_0 betahat_1  
svmat    betahat, names(col)



*- Graph ---

#delimit 		;

hist 		 betahat_1, normal xline(1, lwidth(thick)) scheme(white_piyg)  color(%30) 
					   title("Sampling Distributions of Beta Estimates", size(*0.8)) 
					   subtitle("Biased estimator") 
					   xlab(0.5(0.5)1.5) xti(Estimates, height(3)) 
					   xlabel(0 0.5 1 "True Beta(`=ustrunescape("\u03b2")'=1)"  1.5 2) 
					   text(1.2 0.5 "E(`=ustrunescape("\u03b2")') `=ustrunescape("\u2260")' `=ustrunescape("\u03b2")'" ) 
					   text(1.22 0.470 "`=ustrunescape("\u0302")'")
;
#delimit cr  				   
					   
gr 			 save 		biased_estimator.gph, replace					   
gr 			 export 		biased_estimator.png, replace

```
---

The figure below from the simulation study shows that the OLS estimator is biased in the presence of omitted variable bias.
![combine](https://user-images.githubusercontent.com/101017847/212010487-1be5c8cc-e685-4eea-8c00-5c32a907565b.png)

## A Monte Carlo Study-Inconsistent estimator
This do file below produces a Monte Carlo study of omitted variable	bias	in least squares

---
```stata

/*************************************************************************
A Monte Carlo Study of Omitted Variable Bias in Least Squares	   

Large-Sample Properties of OLS	
							             
12-Jan-2023    		      
								             
Hamza Mutluay			      
  	 			 			   
comments/errors reports: hamzamutluay@gmail.com
****************************************************************************/

cls
clear 		*
se mo 		off, perm
pause 		off, perm


se se 		212355421  // Set the seed to reproduce the same output of simulation


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

*- Parameters of the variance-covariance matrix of jointly distributed random variables  

loc 		sigma12 = 0.4  // Cov(x1,x2) = 0.4
loc 		sigma13 = 0    // Cov(x1,u) = 0
loc 		sigma23 = 0   // Cov(x2,u) = 0

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


forv 		tt = 1/`rep' {
	
mat 		vector_mean = (0.5,0.5,0)
mat 		Covmatrix = (1,`sigma12',`sigma13' \ `sigma12',1, `sigma23' \ `sigma13',`sigma23',1)

drawnorm 	x1 x2 u, n(`n') cov(Covmatrix) m(vector_mean)


** Data Generating Process

g 			y = `beta0' + `beta1'*x1 + `beta2'*x2 + u


** Estimated Model

reg 		y x1 

mat 		betahat[`tt',1] = r(table)[1,2] 	// betahat_0
mat 		betahat[`tt',2] = r(table)[1,1] 	// betahat_1

drop 		y x1 x2 u
mat drop 	vector_mean Covmatrix

	
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


sca drop	sum_betahat_0 sum_betahat_1
}


mat			colname coef_matrix = betahat_0 betahat_1  samplesize
mat li 		coef_matrix


svmat       coef_matrix, names(col)


*- Graph ---


#delimit 		;

scatter     betahat_1 samplesize,  
			ti("Sampling Distributions of Beta Estimates", size(*0.8))
			subti("Inconsistent estimator") 
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
