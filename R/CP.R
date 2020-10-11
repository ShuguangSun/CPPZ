
#' @title Calculate conditional power
#' @description This function calculates the conditional power at an interim analysis of a \cr
#' clinical trial
#' @param r Proportion of subjects assigned to active group
#' @param n1 Number of patients / events at first interim analysis (NULL if parameter f provided)
#' @param n2 Number of patients / events at final analysis (NULL if parameter f provided)
#' @param alpha_1s one-sided alpha
#' @param eff_est effect observed/estimated at interim
#' @param eff_planned planned effect (usually alternative hypothesis, but can be modified to obtain \cr
#' conditional power for any assumprion)
#' @param eff_null effect corresponding with null hypothesis (e.g. 1 for hazard rates, 0 for difference) \cr
#' This accomodates for non-inferiority analyses
#' @param SE standard error: has to be provided if type not equal to "HR". If for instance SE =1 \cr
#' then eff_est and eff_planned correspond with z-scores
#' @param type if type="HR", then SE is calculated, if type="general", then SE's have to be provided \cr
#' by user
#' @param pow Power, needed to calculate sample size in second part provided to obtain given power
#' @param f n1/n2 at interim analysis. Only to be provided if n1 and n2 not provided

#' @return a list of three vectors
#' \itemize{
#' \item CP_obs : conditional power, given 1) observed effect at interim, 2) total sample size, \cr
#' 3) Assumed true effect= "eff_est" parameter
#' #' \item CP_planned : conditional power, given 1) observed effect at interim, 2) total sample size, \cr
#' 3) Assumed true effect= "eff_planned" parameter
#' #' \item n2_inc_new : additional patients needed to obtain given power
#' }
#' @references Lan and Wittes. The B-Value: A Tool for Monitoring Data. Biometrics 1988;44:579-585 \cr
#' Mehta CR, Pocock SJ. Adaptive increase in sample size when interim results are promising: A practical \cr
#' guide with examples
#' @export
#' @examples
#' CP(r=0.5,n1=72 ,n2=180,alpha_1s=0.025,eff_est=0.7,eff_planned=0.7 ,eff_null=1,SE=NULL,
#'    type="HR",pow=0.8)
#' CP(r=0.5              ,alpha_1s=0.025,eff_est=1.5,eff_planned=1.96,eff_null=0,SE=1,
#'    type="general",pow=0.8,f=0.5)
#' CP(r=0.5,n1=316,n2=633,alpha_1s=0.025,eff_est=1.5,eff_planned=1.96,eff_null=0,SE=1,
#'    type="general",pow=0.8)


CP<-function(r,n1=NULL,n2=NULL,alpha_1s,eff_est,eff_planned,eff_null,SE=NULL,type,pow=NULL,f=NULL)
{

  if (is.null(n1)==FALSE & is.null(n2)==FALSE){
    n2_inc   <- n2-n1
    f<-n1/n2
  }

  if (is.null(n1)==TRUE  & is.null(n2)==FALSE){
    n1 <- f*n2
    n2_inc   <- n2-n1
  }

  zalpha   <- qnorm(1-alpha_1s)

  # Calculate conditional power for given observed/estimated effect at interim
  #---------------------------------------------------------------------------

  if (type=="HR"){
    z1_obs     <- (-log(eff_est)    -(-log(eff_null)))*(sqrt(n1*(r*(1-r))))
    z1_planned <- (-log(eff_planned)-(-log(eff_null)))*(sqrt(n1*(r*(1-r))))
  }

  if (type=="general"){
    z1_obs     <- (eff_est    -eff_null)/SE
    z1_planned <- (eff_planned-eff_null)/SE
  }

  stat_obs     <- (zalpha-z1_obs*sqrt(f))/(sqrt(1-f))-(z1_obs    *(sqrt(1-f)))/sqrt(f)
  stat_planned <- (zalpha-z1_obs*sqrt(f))/(sqrt(1-f))-(z1_planned*(sqrt(1-f)))/sqrt(f)

  CP_obs     <- 1-pnorm(stat_obs)
  CP_planned <- 1-pnorm(stat_planned)

  if (is.null(pow)==FALSE & is.null(n2)==FALSE){
    zbeta            <- qnorm(pow)
    n2_inc_new_obs   <- (n1/(z1_obs^2)) * (((zalpha*sqrt(n2)-z1_obs*sqrt(n1))/(n2_inc^0.5))+zbeta)^2
  }

  if (is.null(pow)==FALSE & is.null(n2)==FALSE){
    return(list(CP_obs=CP_obs,CP_planned=CP_planned,n2_inc_new=n2_inc_new_obs))
  }

  if (is.null(pow)==TRUE | is.null(n2)==TRUE){
    return(list(CP_obs=CP_obs,CP_planned=CP_planned))
  }
}


