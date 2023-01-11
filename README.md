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

