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
* We include the parameter 'type' in the lccPlot() function. Now the user can choice among lcc, lpc, and la as plot output. This allowed we reduce the number of parameters in the control list. Now we have only the parameters 'ylab' rather than 'LCC_ylab', 'LPC_ylab', and 'LA_ylab'. In the same way,  we replace the parameters 'LCC_scale_y_continuous',
'LPC_scale_y_continuous', and 'LA_scale_y_continuous' by 'scale_y_continuous'.

* Fixed issue when y-axis labels is changed.

* A new parameter called interaction was included in the lcc() function. This parameter allows to estimate or not the interaction effect among predictors variables in the fixed part of the linear predictor.  