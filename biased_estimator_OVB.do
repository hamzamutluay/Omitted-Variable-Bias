


***********************************
***********************************
***        						***
***   A Monte Carlo Study of    ***
***	  Omitted Variable	Bias	***
***   in Least Squares			***	  https://github.com/hamzamutluay/Omitted-Variable-Bias			 
***   							***
***								***
***	  Finite-Sample 			***
***	  Properties of OLS			***   comments/errors reports: hamzamutluay@gmail.com
***								***
***								***  
***								***
***	  14-Jan-2023    		    ***
***								***
***	  Hamza Mutluay			    *** 
***	  	 			 			***
***********************************
***********************************


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



