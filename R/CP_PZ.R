
#' @title Calculate promising zone and position of observed effect in promising zone
#' @description This function calculates the promising zone boundaries at an interim analysis of a \cr
#' clinical trial
#' @param r Proportion of subjects assigned to active group
#' @param n1 Number of patients / events at first interim analysis (NULL if parameter f provided)
#' @param n2 Number of patients / events at final analysis
#' @param alpha_1s one-sided alpha
#' @param eff_est effect observed/estimated at interim
#' @param eff_planned planned effect (usually alternative hypothesis, but can be modified to obtain \cr
#' conditional power for any assumprion)
#' @param eff_null effect corresponding with null hypothesis (e.g. 1 for hazard rates, 0 for difference) \cr
#' This accomodates for non-inferiority analyses
#' @param SE standard error: has to be provided if type not equal to "HR". If for instance SE =1 \cr
#' then eff_est and eff_planned correspond with z-scores
#' @param p_c proportion in control group (needs to be provided if type="prop")
#' @param sd_t standard deviation in treatment group (needs to be provided if type="cont")
#' @param sd_c standard deviation in control group (needs to be provided if type="cont")
#' @param type if type="HR", then SE is calculated, if type="general", then SE's have to be provided \cr
#' by user, if type="cont" (continuous) then sd_t and sd_c have to be provided \cr
#' if type="prop" then p_c has to be provided
#' @param pow Power, needed to calculate sample size in second part provided to obtain given power
#' @param f n1/n2 at interim analysis. Only to be provided if n1 and n2 not provided
#' @param max (Maximum sample size)/(original sample size), could be for instance 1.5, 2 or 3
#' @param plot_effect TRUE if plotting boundaries with corresponding effect scale

#' @return a list of three vectors
#' \itemize{
#' \item CP_obs : conditional power, given 1) observed effect at interim, 2) total sample size, \cr
#' 3) Assumed true effect= "eff_est" parameter
#' \item CP_planned : conditional power, given 1) observed effect at interim, 2) total sample size, \cr
#' 3) Assumed true effect= "eff_planned" parameter
#' \item n2_inc_new : additional patients needed to obtain given power
#' \item CP_ll lower boundary of promising zone (z-score, b-value, conditional power, effect scale)
#' \item CP_ul upper boundary of promising zone (z-score, b-value, conditional power, effect scale)
#' }
#' @references Lan and Wittes. The B-Value: A Tool for Monitoring Data. Biometrics 1988;44:579-585 \cr
#' Mehta CR, Pocock SJ. Adaptive increase in sample size when interim results are promising: A practical \cr
#' guide with examples. Statist. Med. 2011;30:3267–3284
#' @export
#' @examples
#' CP_PZ(r=0.5,n1=72,n2=180,alpha_1s=0.025,eff_est=0.7,eff_planned=0.7,eff_null=1,SE=NULL,
#'       type="HR",max=1.5,pow=0.8       ,plot_effect=1)
#' CP_PZ(r=0.5      ,n2=180,alpha_1s=0.025,eff_est=0.7,eff_planned=0.7,eff_null=1,SE=NULL,
#'       type="HR",max=1.5,pow=0.8,f=0.25,plot_effect=1)

CP_PZ<-function(r,n1=NULL,n2=NULL,alpha_1s,eff_est,eff_planned,eff_null,SE=NULL,p_c=NULL,sd_t=NULL,sd_c=NULL,type,pow=NULL,f=NULL,max,plot_effect)
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

  if (type=="prop"){  # pooled estimate (East manual)
    p_t<- p_c+eff_est
    n1_t<-n1*r
    n1_c<-n1*(1-r)
    p  <- (n1_t*p_t + n1_c*p_c)/(n1_t+n1_c)
    SE <- sqrt(p*(1-p)*(1/n1_t + 1/n1_c))
    z1_obs <- -(eff_est-eff_null)/SE # '-' because case lower proportion is better

    p_t<- p_c+eff_planned
    n1_t<-n1*r
    n1_c<-n1*(1-r)
    p  <- (n1_t*p_t + n1_c*p_c)/(n1_t+n1_c)
    SE <- sqrt(p*(1-p)*(1/n1_t + 1/n1_c))
    z1_planned <- -(eff_planned-eff_null)/SE # '-' because case lower proportion is better
  }

  if (type=="cont"){
    n1_t<-n1*r
    n1_c<-n1*(1-r)

    SE<- sqrt(sd_t^2/n1_t+sd_c^2/n1_c)

    z1_obs     <- (eff_est    -eff_null)/SE
    z1_planned <- (eff_planned-eff_null)/SE
  }

  if (type=="general"){
    z1_obs     <- (eff_est    -eff_null)/SE
    z1_planned <- (eff_planned-eff_null)/SE
  }

  stat_obs     <- (zalpha-z1_obs*sqrt(f))/(sqrt(1-f))-(z1_obs    *(sqrt(1-f)))/sqrt(f)
  stat_planned <- (zalpha-z1_obs*sqrt(f))/(sqrt(1-f))-(z1_planned*(sqrt(1-f)))/sqrt(f)

  CP_obs     <- 1-pnorm(stat_obs)
  CP_planned <- 1-pnorm(stat_planned)

  if (is.null(pow)==FALSE & is.null(n1)==FALSE & is.null(n2)==FALSE){
    zbeta            <- qnorm(pow)
    n2_inc_new_obs   <- (n1/(z1_obs^2)) * (((zalpha*sqrt(n2)-z1_obs*sqrt(n1))/(n2_inc^0.5))+zbeta)^2
  }

  # Calculate promising zone (zone where b-value<=qnorm(1-alpha_1s))
  #-----------------------------------------------------------------

  z1           <- seq(0.5,2.05,by=0.001)
  CP           <- 1-pnorm((zalpha -z1*sqrt(f))/ sqrt(1-f)- z1*sqrt(1-f)/sqrt(f ))
  n2_inc_new   <- (n1/z1^2)*((zalpha*sqrt(n2)-z1*sqrt(n1))/sqrt(n2_inc)+zbeta)^2

  n2_new       <- n1+(n2_inc_new)
  n2_real      <- pmin(n2_new,rep(max*n2,length(n2_new)),na.rm=T)
  n2_real_inc  <- n2_real-n1
  b            <- n2_real^(-0.5)*( (sqrt(n2_real_inc/n2_inc))*(zalpha*sqrt(n2)-z1*sqrt(n1))+z1*sqrt(n1))

  Results <- data.frame(cbind(z1,b,CP))
  CP_ll<-Results[min(which(Results$b <= zalpha)),]
  CP_ul<-Results[max(which(Results$b <= zalpha)),]

  plot(CP,b,type="l",xlim=c(0,1),ylim=c(min(b),2.05),lwd=1)
  grid()
  abline(h=zalpha)
  abline(v=CP_obs,col="red")
  abline(v=CP_ll$CP,col="blue")
  abline(v=CP_ul$CP,col="blue")

  text(CP_ll$CP+0.04,min(b)+0.01,"Boundary="      ,font = 2)
  text(CP_ll$CP+0.04,min(b)     ,round(CP_ll$CP,2),font = 2)
  text(CP_ul$CP+0.04,min(b)+0.01,"Boundary="      ,font = 2)
  text(CP_ul$CP+0.04,min(b)     ,round(CP_ul$CP,2),font = 2)

  if (plot_effect==1){

    if (type=="HR"){
      eff <- exp(log(1)-z1/(sqrt(n1*(r*(1-r)))))
    }

    if (type=="prop"){ # to do: grid search, but only for boundaries
    }

    if (type=="general" | type=="cont"){
      eff <- z1*SE+eff_null
    }

    Results <- (cbind(Results,eff))
    CP_ll<-Results[min(which(Results$b <= zalpha)),];rownames(CP_ll)<-NULL
    CP_ul<-Results[max(which(Results$b <= zalpha)),];rownames(CP_ul)<-NULL

    text(CP_ll$CP+0.04,2.01,"Effect="         ,font = 2)
    text(CP_ll$CP+0.04,2   ,round(CP_ll$eff,2),font = 2)
    text(CP_obs  +0.04,2.02,"Observed"        ,font = 2)
    text(CP_obs  +0.04,2.01,"Effect="         ,font = 2)
    text(CP_obs  +0.04,2   ,round(eff_est  ,2),font = 2)
    text(CP_ul$CP+0.04,2.01,"Effect="         ,font = 2)
    text(CP_ul$CP+0.04,2   ,round(CP_ul$eff,2),font = 2)
  }

  return(list(CP_obs=CP_obs,CP_planned=CP_planned,n2_inc_new=n2_inc_new_obs,CP_ll=CP_ll,CP_ul=CP_ul))

}

