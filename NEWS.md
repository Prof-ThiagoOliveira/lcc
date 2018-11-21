# lcc 1.1.0

## New argument

*  We included an argument called "interaction " in the lcc function as an option to estimate the interaction effect between time and method. If interaction is declared FALSE, the lcc function consider only the main effects of time and method in the model.

## Bug fixes and minor improvements

* Fixed issue when the wrong name of variable is declared.  

* Fixed issue when three or more levels of methods is considered in the dataset.

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
