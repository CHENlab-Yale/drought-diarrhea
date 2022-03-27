##################################################################################
####            Associations between long-term drought and diarrhea           ####
####      among children under five in low- and middle-income countries       ####
####                                                                          ####
####   Pin Wang, Ernest Asare, Virginia E. Pitzer, Robert Dubrow, Kai Chen    ####
####                     Yale School of Public Health                         ####
####                             Mar 25, 2022                                 ####
##################################################################################

library(MASS);library(splines)

################# MAIN ANALYSIS ###############################
###Dataset using data.all, including the following variables:
## "country","survey","diarrhea","drought6","drought12","drought18","drought24","age","sex","residence",
## "education","wealth","tmean","rain","month_int","survey_site","year_int", and
## WASH variables: "water_source","time_collect","toilet","place_wash","if_water", "if_soap","any_do"

data.main <- na.omit(data.all[,c("country","survey","diarrhea","drought6","drought12","drought18","drought24","age","sex","residence",
                                 "education","wealth","tmean","rain","month_int","survey_site","year_int")])

# two crude models and main model (all without WASH)
## results used to create Figure 3 and Supplementary Table 3
### replace drought6 with drought12/18/24 for another timescale
m.all1 <- glmmPQL(diarrhea ~ factor(drought6) + 
                  ns(tmean,df=3) + ns(rain,df=3) + ns(month_int,df=3) + ns(year_int,df=3), random= ~1|survey_site, data.all, family=binomial(link="log"))
m.all2 <- glmmPQL(diarrhea ~ factor(drought6) + 
                  ns(tmean,df=3) + ns(rain,df=3) + ns(month_int,df=3) + ns(year_int,df=3), random= ~1|survey_site, data.main, family=binomial(link="log"))
m.all3 <- glmmPQL(diarrhea ~ factor(drought6) + age + factor(sex) + factor(residence) + factor(education) + factor(wealth) + 
                  ns(tmean,df=3) + ns(rain,df=3) + ns(month_int,df=3) + ns(year_int,df=3), random= ~1|survey_site, data.main, family=binomial(link="log"))

################# MEDIATION ANALYSIS ####################
data.med <- na.omit(data.all[,c("country","survey","diarrhea","drought6","drought12","drought18","drought24","age","sex","residence",
                                 "education","wealth","tmean","rain","month_int","water_source","time_collect","toilet","place_wash","if_water",
                                 "if_soap","any_do","survey_site","year_int")]) 

base <- as.formula("diarrhea ~ age + factor(sex) + factor(residence) + factor(education) + factor(wealth)")

# create dataframe to save mediation results
med.all <- data.frame(matrix(ncol=3,nrow=4))
names(med.all) <- c("Effect","6-month-mild","6-month-severe")
med.all$Effect <- c("Total effect","Direct effect","Indirect effect","IE percentage")

# main model and adjusted model
m.med1 <- glmmPQL(update(base, . ~ factor(drought6) + . + 
                           ns(tmean,df=3) + ns(rain,df=3) + ns(month_int,df=3) + ns(year_int,df=3)), random= ~1|survey_site, data.med, family=binomial(link="log"))
m.med2 <- glmmPQL(update(base, . ~ factor(drought6) + . + 
                           factor(water_source) + factor(time_collect) + factor(toilet) + factor(place_wash) + factor(if_water) + factor(if_soap) + factor(any_do) +
                           ns(tmean,df=3) + ns(rain,df=3) + ns(month_int,df=3) + ns(year_int,df=3)), random= ~1|survey_site, data.med, family=binomial(link="log"))

## results used to create Table 2
for (j in 1:2) { # 2 models
  for (k in 1:2) { # 2 drought categories
    TE.ik <- coef(summary(m.med1))[k+1,"Value"]
    med.all[1,k+1] <- TE.ik
    
    DE.ik <- coef(summary(m.med2))[k+1,"Value"]
    med.all[2,k+1] <- DE.ik
    
    IE.ik <- TE.ik-DE.ik
    med.all[3,k+1] <- IE.ik
    
    # Proportion mediated
    med.all[4,k+1] <- IE.ik/TE.ik*100
  }
}

# function to estimate the mediation effects and bootstrap 95% CI's
## only performed for 6-month drought
boot <- list()
for (i in 1:1000) { # 1000 bootstrapping loops
  tryCatch({ # some simulated models don't converge
    data.i <- data.med[sample(1:nrow(data.med), replace=T), ]
    
    m.med.i1 <- glmmPQL(update(base, . ~ factor(drought6) + . + 
                                 ns(tmean,df=3) + ns(rain,df=3) + ns(month_int,df=3) + ns(year_int,df=3)), random= ~1|survey_site, data.i, family=binomial(link="log"))
    # Total effects
    TOT.mild.i <- coef(summary(m.med.i1))[2,"Value"]
    TOT.severe.i <- coef(summary(m.med.i1))[3,"Value"]
    
    m.med.i2 <- glmmPQL(update(base, . ~ factor(drought6) + . + 
                                factor(water_source) + factor(time_collect) + factor(toilet) + factor(place_wash) + factor(if_water) + factor(if_soap) + factor(any_do) +
                                ns(tmean,df=3) + ns(rain,df=3) + ns(month_int,df=3) + ns(year_int,df=3)), random= ~1|survey_site, data.i, family=binomial(link="log"))
    # Direct effects
    DE.mild.i <- coef(summary(m.med.i2))[2,"Value"]
    DE.severe.i <- coef(summary(m.med.i2))[3,"Value"]
    
    # Indirect effects
    IE.mild.i <- TOT.mild.i-DE.mild.i
    IE.severe.i <- TOT.severe.i-DE.severe.i
    
    boot[[i]] <- c(IE.mild.i, DE.mild.i, TOT.mild.i, IE.severe.i, DE.severe.i, TOT.severe.i)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# matrices to store the Indirect effect (IE)
G <- 1000 
IE.mild <- IE.severe <- matrix(rep(0,G), G, 1)
DE.mild <- DE.severe <- matrix(rep(0,G), G, 1)
TOT.mild <- TOT.severe <- matrix(rep(0,G), G, 1)

IE.mild <- unlist(lapply(boot, '[[', 1)) # This returns a vector with the 1st number of each element of the list (IE)
DE.mild <- unlist(lapply(boot, '[[', 2))
TOT.mild <- unlist(lapply(boot, '[[', 3))
IE.severe <- unlist(lapply(boot, '[[', 4)) 
DE.severe <- unlist(lapply(boot, '[[', 5))
TOT.severe <- unlist(lapply(boot, '[[', 6))

# IE estimate and 95%CI
est.mild.IE <- mean(IE.mild)
ci.mild.IE <- quantile(IE.mild, c(0.025, 0.975))
est.severe.IE <- mean(IE.severe)
ci.severe.IE <- quantile(IE.severe, c(0.025, 0.975))