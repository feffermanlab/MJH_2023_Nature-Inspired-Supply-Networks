
library(ggplot2)
library(MuMIn)
library(nlme)
library(emmeans)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Set working directory
setwd()

#Import data
responseStratData <- read.csv("simSummaries_responseStrategies.csv", header = TRUE)
responseStratData$cvFillRate <- responseStratData$sdFillRate / responseStratData$meanFillRate
responseStratData$cvPercentShortfall <- responseStratData$sdPercentShortfall / responseStratData$percentShortfall

#Remove trials 31-50 (mean responses converge following ~30 simulation runs)
responseStratData <- responseStratData[which(responseStratData$simID < 31), ]

#Split data into pre- and post-disruption
#An initial burn-in period (i.e., the first 5 time steps) are removed from the pre-disruption data
responseStratData_preDis <- responseStratData[which(responseStratData$TimeStep < 51 & 
                                                      responseStratData$TimeStep > 5), ]
responseStratData_postDis <- responseStratData[which(responseStratData$TimeStep > 50), ]

#######################
#PreDisruption Analysis
#######################

responseStratData_preDis$ID <- as.factor(responseStratData_preDis$uniqueID)
responseStratData_preDis$StrategyF <- factor(responseStratData_preDis$Strategy, levels = c("None", "inflateOrder", "distPriority", "addSupplier"))
responseStratData_preDis$cvFillRate <- responseStratData_preDis$sdFillRate / responseStratData_preDis$meanFillRate
responseStratData_preDis$cvPercentShortfall <- responseStratData_preDis$sdPercentShortfall / responseStratData_preDis$percentShortfall

####Mean fill rate####
meanFillRateData <- summarySE(data = responseStratData, 
                              measurevar = "meanFillRate", 
                              groupvars = c("Strategy", "TimeStep"))

ggplot(data = meanFillRateData, aes(x = TimeStep, y = meanFillRate, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("Mean order fill rate") +
  geom_line(size = 1.5) + xlim(6,100) + 
  geom_ribbon(aes(ymin = meanFillRate - (1.96*se), ymax = meanFillRate + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                       "Inflate Order", "None")) + 
  scale_color_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                        "Inflate Order", "None")) + 
  theme(legend.position = c(0.8, 0.9), legend.title = element_blank()) +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global models
mFR.M1 <- lme(meanFillRate ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = responseStratData_preDis, na.action = na.fail)

mFR.M0 <- lme(meanFillRate ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_preDis, na.action = na.fail)

#Based on AICc, mFR.M1 is preferred
#Evaluate residuals
R <- resid(mFR.M1, type = "normalized")
Ft <- fitted(mFR.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_preDis$StrategyF)
plot(R ~ responseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dredge(mFR.M1, REML = FALSE)

#Top model strongly favored (weight = 1)
summary(mFR.M1)
intervals(mFR.M1, which = "fixed")

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(mFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 10), adjust = "bonferroni")
emmeans(mFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 25), adjust = "bonferroni")
emmeans(mFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 40), adjust = "bonferroni")

####% shortfall####
percentShortfallData <- summarySE(data = responseStratData, 
                              measurevar = "percentShortfall", 
                              groupvars = c("Strategy", "TimeStep"))

ggplot(data = percentShortfallData, aes(x = TimeStep, y = percentShortfall, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("% shortfall in supply") +
  geom_line(size = 1.5) + xlim(6,100) + 
  geom_ribbon(aes(ymin = percentShortfall - (1.96*se), ymax = percentShortfall + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                       "Inflate Order", "None")) + 
  scale_color_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                        "Inflate Order", "None")) + 
  theme(legend.position = c(0.8, 0.2), legend.title = element_blank()) + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global model
pS.M1 <- lme(percentShortfall ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = responseStratData_preDis, method = "REML", na.action = na.fail)

pS.M0 <- lme(percentShortfall ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_preDis, na.action = na.fail)

#Based on AICc, ps.M1 is preferred
#Evaluate residuals
R <- resid(pS.M1, type = "normalized")
Ft <- fitted(pS.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_preDis$StrategyF)
plot(R ~ responseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dpS.M1 <- dredge(pS.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dpS.M1.Reduced<-subset(dpS.M1, !nested(.))

#Global model strongly favored (weight = 1)
summary(pS.M1)
intervals(pS.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(pS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 10), adjust = "bonferroni")
emmeans(pS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 25), adjust = "bonferroni")
emmeans(pS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 40), adjust = "bonferroni")

#####Price per unit####
priceData <- summarySE(data = responseStratData, 
                                  measurevar = "meanPricePerUnit", 
                                  groupvars = c("Strategy", "TimeStep"))

ggplot(data = priceData, aes(x = TimeStep, y = meanPricePerUnit, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("Mean price per unit") +
  geom_line(size = 1.5) + xlim(6,100) + 
  geom_ribbon(aes(ymin = meanPricePerUnit - (1.96*se), ymax = meanPricePerUnit + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                       "Inflate Order", "None")) + 
  scale_color_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                        "Inflate Order", "None")) + 
  theme(legend.position = c(0.8, 0.2), legend.title = element_blank()) + ylim(15,18)+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global model
mPPU.M1 <- lme(meanPricePerUnit ~ StrategyF * TimeStep, random = ~1|ID, 
             weights = varIdent(form = ~1 | StrategyF), data = responseStratData_preDis, method = "REML", na.action = na.fail)

mPPU.M0 <- lme(meanPricePerUnit ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_preDis, na.action = na.fail)

#Based on AICc, mPPU.M1 is preferred
#Evaluate residuals
R <- resid(mPPU.M1, type = "normalized")
Ft <- fitted(mPPU.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_preDis$StrategyF)
plot(R ~ responseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dmPPU.M1 <- dredge(mPPU.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dmPPU.M1.Reduced<-subset(dmPPU.M1, !nested(.))

#Top model strongly favored (weight = 0.989)
summary(mPPU.M1)
intervals(mPPU.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 10), adjust = "bonferroni")
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 25), adjust = "bonferroni")
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 40), adjust = "bonferroni")

####Change in purchasing strategy####
stratData <- summarySE(data = responseStratData, 
                       measurevar = "percentPricePriority", 
                       groupvars = c("Strategy", "TimeStep"))

ggplot(data = stratData, aes(x = TimeStep, y = percentPricePriority, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("% Price priority") +
  geom_line(size = 1.5) + xlim(6, 100) + 
  geom_ribbon(aes(ymin = percentPricePriority - (1.96*se), ymax = percentPricePriority + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                       "Inflate Order", "None")) + 
  scale_color_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                        "Inflate Order", "None")) + 
  theme(legend.position = c(0.8, 0.2), legend.title = element_blank()) +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global models
pPP.M1 <- lme(percentPricePriority ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = responseStratData_preDis, method = "REML", na.action = na.fail)

pPP.M0 <- lme(percentPricePriority ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_preDis, method = "REML", na.action = na.fail)

#Based on AICc, pPP.M1 is preferred
#Evaluate residuals
R <- resid(pPP.M1, type = "normalized")
Ft <- fitted(pPP.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_preDis$StrategyF)
plot(R ~ responseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dredge(pPP.M1, REML = FALSE)

#Global model strongly favored (weight = 0.992)
summary(pPP.M1)
intervals(pPP.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(pPP.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 10), adjust = "bonferroni")
emmeans(pPP.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 25), adjust = "bonferroni")
emmeans(pPP.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 40), adjust = "bonferroni")

####CV of fill rate####
cvFillRateData <- summarySE(data = responseStratData, 
                              measurevar = "cvFillRate", 
                              groupvars = c("Strategy", "TimeStep"))

ggplot(data = cvFillRateData, aes(x = TimeStep, y = cvFillRate, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("CV of order fill rate") +
  geom_line(size = 1.5) + 
  geom_ribbon(aes(ymin = cvFillRate - (1.96*se), ymax = cvFillRate + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                       "Inflate Order", "None")) + 
  scale_color_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                        "Inflate Order", "None")) + 
  theme(legend.position = c(0.8, 0.25), legend.title = element_blank()) + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global model
cvFR.M1 <- lme(cvFillRate ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = responseStratData_preDis, na.action = na.fail)

cvFR.M0 <- lme(cvFillRate ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_preDis, na.action = na.fail)

#Based on AICc, cvFR.M1 is preferred
#Evaluate residuals
R <- resid(cvFR.M1, type = "normalized")
Ft <- fitted(cvFR.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_preDis$StrategyF)
plot(R ~ responseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dredge(cvFR.M1, REML = FALSE)

#Global model strongly favored (weight = 1)
summary(cvFR.M1)
intervals(cvFR.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 10), adjust = "bonferroni")
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 25), adjust = "bonferroni")
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 40), adjust = "bonferroni")

####CV of % shortfall####
cvPercentShortfallData <- summarySE(data = responseStratData, 
                            measurevar = "cvPercentShortfall", 
                            groupvars = c("Strategy", "TimeStep"))

ggplot(data = cvPercentShortfallData, aes(x = TimeStep, y = cvPercentShortfall, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("CV of % shortfall in supply") +
  geom_line(size = 1.5) + xlim(6, 100) + 
  geom_ribbon(aes(ymin = cvPercentShortfall - (1.96*se), ymax = cvPercentShortfall + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                       "Inflate Order", "None")) + 
  scale_color_manual(values = c("#D55E00", "#CC79A7", "black", "#56B4E9"), labels = c("Add Supplier", "Dist. Priority", 
                                                                                        "Inflate Order", "None")) + 
  theme(legend.position = c(0.8, 0.19), legend.title = element_blank()) + ylim(0.5,1.6)+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global modelS
cvPS.M1 <- lme(cvPercentShortfall ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = responseStratData_preDis, na.action = na.fail)

cvPS.M0 <- lme(cvPercentShortfall ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_preDis, na.action = na.fail)

#Based on AICc, cvPS.M1 is preferred
#Evaluate residuals
R <- resid(cvPS.M1, type = "normalized")
Ft <- fitted(cvPS.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_preDis$StrategyF)
plot(R ~ responseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dredge(cvPS.M1, REML = FALSE)

#Global model strongly favored (weight = 1)
summary(cvPS.M1)
intervals(cvPS.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 10), adjust = "bonferroni")
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 25), adjust = "bonferroni")
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 40), adjust = "bonferroni")

#########################
#Post-disruption analysis
#########################

responseStratData_postDis$ID <- as.factor(responseStratData_postDis$uniqueID)
responseStratData_postDis$StrategyF <- factor(responseStratData_postDis$Strategy, levels = c("None", "inflateOrder", "distPriority", "addSupplier"))
responseStratData_postDis$cvFillRate <- responseStratData_postDis$sdFillRate / responseStratData_postDis$meanFillRate
responseStratData_postDis$cvPercentShortfall <- responseStratData_postDis$sdPercentShortfall / responseStratData_postDis$percentShortfall


####Mean fill rate####

#Fit global models
mFR.M1 <- lme(meanFillRate ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = responseStratData_postDis, na.action = na.fail)

mFR.M0 <- lme(meanFillRate ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_postDis, na.action = na.fail)

#Based on AICc, cvFR.M1 is preferred
#Evaluate residuals
R <- resid(mFR.M1, type = "normalized")
Ft <- fitted(mFR.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_postDis$StrategyF)
plot(R ~ responseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dmFR.M1 <- dredge(mFR.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dmFR.M1.Reduced<-subset(dmFR.M1, !nested(.))

#Global model strongly favored (weight = 1)
summary(mFR.M1)
intervals(mFR.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(mFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 60), adjust = "bonferroni")
emmeans(mFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 75), adjust = "bonferroni")
emmeans(mFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 90), adjust = "bonferroni")

####% shortfall####

#Fit global modelS
pS.M1 <- lme(percentShortfall ~ StrategyF * TimeStep, random = ~1|ID, 
             weights = varIdent(form = ~1 | StrategyF), data = responseStratData_postDis, method = "REML", na.action = na.fail)

pS.M0 <- lme(percentShortfall ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_postDis, na.action = na.fail)

#Based on AICc, pS.M1 is preferred
#Evaluate residuals
R <- resid(pS.M1, type = "normalized")
Ft <- fitted(pS.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_postDis$StrategyF)
plot(R ~ responseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dpS.M1 <- dredge(pS.M1, REML = FALSE)

#Global model strongly favored
summary(pS.M1)
intervals(pS.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(pS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 60), adjust = "bonferroni")
emmeans(pS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 75), adjust = "bonferroni")
emmeans(pS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 90), adjust = "bonferroni")

#####Price per unit####

#Fit global models
mPPU.M1 <- lme(meanPricePerUnit ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = responseStratData_postDis, method = "REML", na.action = na.fail)

mPPU.M0 <- lme(meanPricePerUnit ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_postDis, na.action = na.fail)

#Based on AICc, mPPU.M1 is preferred
#Evaluate residuals
R <- resid(mPPU.M1, type = "normalized")
Ft <- fitted(mPPU.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_postDis$StrategyF)
plot(R ~ responseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dmPPU.M1 <- dredge(mPPU.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dmPPU.M1.Reduced<-subset(dmPPU.M1, !nested(.))

#Global model strongly favored (weight = 1)
summary(mPPU.M1)
intervals(mPPU.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 60), adjust = "bonferroni")
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 75), adjust = "bonferroni")
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 90), adjust = "bonferroni")

####Change in purchasing strategy####

#Fit global models
pPP.M1 <- lme(percentPricePriority ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = responseStratData_postDis, method = "REML", na.action = na.fail)

pPP.M0 <- lme(percentPricePriority ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_postDis, method = "REML", na.action = na.fail)

#Based on AICc, pPP.M1 is preferred
#Evaluate residuals
R <- resid(pPP.M1, type = "normalized")
Ft <- fitted(pPP.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_postDis$StrategyF)
plot(R ~ responseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dpPP.M1 <- dredge(pPP.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dpPP.M1.Reduced<-subset(dpPP.M1, !nested(.))

#Top-ranked model strongly favored (weight = 0.999)
pPP.M2 <- lme(percentPricePriority ~ StrategyF, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = responseStratData_postDis, method = "REML", na.action = na.fail)

summary(pPP.M2)
intervals(pPP.M2)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(pPP.M2, pairwise ~ StrategyF, adjust = "bonferroni")

####CV of fill rate####

#Fit global models
cvFR.M1 <- lme(cvFillRate ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = responseStratData_postDis, na.action = na.fail)

cvFR.M0 <- lme(cvFillRate ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_postDis, na.action = na.fail)

#Based on AICc, cvFR.M1 is preferred
#Evaluate residuals
R <- resid(cvFR.M1, type = "normalized")
Ft <- fitted(cvFR.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_postDis$StrategyF)
plot(R ~ responseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dcvFR.M1 <- dredge(cvFR.M1, REML = FALSE)

#Global model strongly favored (weight = 1)
summary(cvFR.M1)
intervals(cvFR.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 60), adjust = "bonferroni")
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 75), adjust = "bonferroni")
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 90), adjust = "bonferroni")

####CV of % shortfall####

#Fit global modelS
cvPS.M1 <- lme(cvPercentShortfall ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = responseStratData_postDis, na.action = na.fail)

cvPS.M0 <- lme(cvPercentShortfall ~ StrategyF * TimeStep, random = ~1|ID, data = responseStratData_postDis, na.action = na.fail)

#Based on AICc, cvPS.M1 is preferred
#Evaluate residuals
R <- resid(cvPS.M1, type = "normalized")
Ft <- fitted(cvPS.M1)
qqnorm(R)
qqline(R)
plot(R ~ responseStratData_postDis$StrategyF)
plot(R ~ responseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dredge(cvPS.M1, REML = FALSE)

#Global model strongly favored (weight = 0.991)
summary(cvPS.M1)
intervals(cvPS.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 60), adjust = "bonferroni")
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 75), adjust = "bonferroni")
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 90), adjust = "bonferroni")
