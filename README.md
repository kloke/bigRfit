bigRfit: Rank-based estimation for linear models in a big data setting
----------------------------------------------------------------------

### Update v0.8.1 December 2022 ###
New driver function **bigRreg**
Fits regression model y = 1 alpha + X beta + e 
where 1 is an n x 1 vector of ones and X is an n x p matrix.
The model estimates are returned as alphahat and betahat.
Note: coef is NULL by design as it is reserved for future development.

Uses biglm to obtain projection rather core R qr as was the case in the previous version --- previous version now renamed bigRfit1. 
New driver version of bigRfit is planned.


### Background ###
Rank-based (R) estimation for statistical models is a robust nonparametric alternative to classical estimation procedures such as least squares. R methods have been developed for models ranging from linear models, to linear mixed models, to timeseries, to nonlinear models. Advantages of these R methods over traditional methods such as maximum-likelihood or least squares is that they require fewer assumptions, are robust to gross outliers, and are highly efficient at a wide range of distributions. 
See also [Rfit](https://github.com/kloke/Rfit).

bigRfit 
- uses a partial ranking to save computational time
- coordinate free Newton steps w/ backtracking
- a density estimate of the scale parameter tau
