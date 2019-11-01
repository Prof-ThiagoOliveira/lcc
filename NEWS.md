# lcc 1.0.1

## New argument

*  We included an argument called "interaction " in the lcc function as
   an option to estimate the interaction effect between time and
   method. If interaction is declared FALSE, the lcc function consider
   only the main effects of time and method in the model.

## Bug fixes and minor improvements

* Fixed issue when the wrong name of variable is declared.

* Fixed issue when three or more levels of methods is considered in the
  dataset.

* Fixed issue in lccWrapper function get the correct rho:
  	* previous code:
	```
	if(is.na(rho[[2]])=TRUE)
	```
	* current code:
	```
	if(length(rho)==1){
            return(rho[[1]])
        }else(if(is.na(rho[[2]])){
            return(rho[[1]])
         }else(return(rho[[n.delta]])))
	 ```
# lcc 1.0.2

# Bug fixes and minor improvements

* Fixed issue when the "all.plots" argument is declared FALSE. Now the
list of plot can be split into multiple pages and save them as pdf
using ggsave function. The marrangeGrob function of package gridExtra
was used to solve this problem.

* Fixed issue in lccWrapper function get the correct rho:
  	* previous code:
	```
	if(length(rho)==1){
            return(rho[[1]])
        }else(if(is.na(rho[[2]])){
            return(rho[[1]])
         }else(return(rho[[n.delta]])))
   ```
   	* current code:
    ```
	if(length(rho)==1){
            return(rho[[1]])
        }else(if(sum(is.na(rho[[2]]))!=0){
            return(rho[[1]])
         }else(return(rho[[n.delta]])))
   ```
* We include the parameter 'type' in the lccPlot() function. Now the user can choice among lcc, lpc, or la as plot output.

* Fixed issue when y-axis labels is changed.

* A new parameter called interaction was included in the lcc() function. This parameter allows to estimate or not the interaction effect among predictors variables in the fixed part of the linear predictor.

# lcc 1.0.3

# Bug fixes and minor improvements

* Added a counter for the bootstrap number (lccBootstrap line 90)

* The got (goodness of fittenes) result was changed for CCC between fitted values (mixed effect model) and original observations rather than fitted LCC and sampled CCC. This make more sense because the model goodness of fit have a high impact in the LCC.

* Fixed issue in ciCompute function to get the correct confidence intervals for the LA statistics. Versions before 1.0.3 the CI for the LA is computed using the logit transformation and versions equal or above 1.0.3 the we are using the arcsin(sqrt(x)) rather than logit.

* Updated the summary.lcc output (cleaner and more informative)

# New functions

* print.lcc: print information about the longitudinal concordance correlation.

* fitted.lcc: fitted values from object of class lcc returned by modeling functions.

* plot.lcc: diagnostic plots for conditional error and random effects from the linear mixed-effects fit are obtained.

* print.summary.lcc: information summarizing the fitted longitudinal concordance correlation is printed.

* residuals.lcc: extract the residulas from the model used to estimate the longitudinal concordance correlation function.

* vcov.lcc: Returns the variance-covariance matrix of a fitted lcc model object.

* getVarCov.lcc: returns the variance-covariance matrix of a fitted lcc model object.

* AIC.lcc and BIC.lcc: calculate the Akaike's 'An Information Criterion' or the BIC or SBC (Schwarz's Bayesian criterion) for an object of class lcc.

* ranef.lcc: the estimated random effects 

* coef.lcc: the fixed effects estimated and corresponding random effects estimates are obtained at subject levels.

* logLik: extract fitted log-likelihood
