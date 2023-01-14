


***********************************
***********************************
***        						***
***   A Monte Carlo Study of    ***
***	  Omitted Variable	Bias	***
***   in Least Squares			***	  https://github.com/hamzamutluay/Omitted-Variable-Bias		 
***   							***
***								***
***	  Large Sample Properties	***
***	          of OLS			***   comments/errors reports: hamzamutluay@gmail.com
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








