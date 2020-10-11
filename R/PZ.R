

#' @title Calculate promising zone
#' @description This function calculates the promising zone boundaries at an interim analysis of a \cr
#' clinical trial
#' @param r Proportion of subjects assigned to active group
#' @param n1 Number of patients / events at first interim analysis (NULL if parameter f provided)
#' @param n2 Number of patients / events at final analysis
#' @param alpha_1s one-sided alpha
#' @param eff_null effect corresponding with null hypothesis (e.g. 1 for hazard rates, 0 for difference) \cr
#' This accomodates for non-inferiority analyses
#' @param SE standard error: has to be provided if type not equal to "HR". If for instance SE =1 \cr
#' then eff_est and eff_planned correspond with z-scores
#' @param type if type="HR", then SE is calculated, if type="general", then SE's have to be provided \cr
#' by user
#' @param pow Power, needed to calculate sample size in second part provided to obtain given power
#' @param f n1/n2 at interim analysis. Only to be provided if n1 and n2 not provided
#' @param max (Maximum sample size)/(original sample size), could be for instance 1.5, 2 or 3
#' @param plot_effect TRUE if plotting boundaries with corresponding effect scale

#' @return a list of three vectors
#' \itemize{
#' \item CP_ll lower boundary of promising zone (z-score, b-value, conditional power, effect scale)
#' \item CP_ul upper boundary of promising zone (z-score, b-value, conditional power, effect scale)
#' }
#' @references Lan and Wittes. The B-Value: A Tool for Monitoring Data. Biometrics 1988;44:579-585 \cr
#' Mehta CR, Pocock SJ. Adaptive increase in sample size when interim results are promising: A practical \cr
#' guide with examples. Statist. Med. 2011;30:3267–3284
#' @export
#' @examples
#' PZ(r=0.5,n2=180,alpha_1s=0.025,type="HR",max=1.5,pow=0.8,f=0.25)
#' PZ(r=0.5,n2=180,alpha_1s=0.025,type="HR",max=1.5,pow=0.9,f=0.25)

PZ<-function(r,n1=NULL,n2=NULL,alpha_1s,eff_null=NULL,SE=NULL,type,pow=NULL,f=NULL,max,plot_effect=TRUE)
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
  zbeta    <- qnorm(pow)

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
  CP_ll<-Results[min(which(Results$b <= zalpha)),];rownames(CP_ll)<-NULL
  CP_ul<-Results[max(which(Results$b <= zalpha)),];rownames(CP_ul)<-NULL

  plot(CP,b,type="l",xlim=c(0,1),ylim=c(min(b),2.05),lwd=1)
  grid()
  abline(h=zalpha)
  abline(v=CP_ll$CP,col="blue")
  abline(v=CP_ul$CP,col="blue")

  text(CP_ll$CP+0.04,min(b)+0.01,"Boundary="      ,font = 2)
  text(CP_ll$CP+0.04,min(b)     ,round(CP_ll$CP,2),font = 2)
  text(CP_ul$CP+0.04,min(b)+0.01,"Boundary="      ,font = 2)
  text(CP_ul$CP+0.04,min(b)     ,round(CP_ul$CP,2),font = 2)

  if (plot_effect==T){

    if (type=="HR"){
      eff <- exp(log(1)-z1/(sqrt(n1*(r*(1-r)))))
    }

    if (type=="general"){
      eff <- z1*SE+eff_null
    }

    Results <- (cbind(Results,eff))
    CP_ll<-Results[min(which(Results$b <= zalpha)),]
    CP_ul<-Results[max(which(Results$b <= zalpha)),]

    text(CP_ll$CP+0.04,2.01,"Effect="         ,font = 2)
    text(CP_ll$CP+0.04,2   ,round(CP_ll$eff,2),font = 2)
    text(CP_ul$CP+0.04,2.01,"Effect="         ,font = 2)
    text(CP_ul$CP+0.04,2   ,round(CP_ul$eff,2),font = 2)
  }

  return(list(CP_ll=CP_ll,CP_ul=CP_ul))
}


