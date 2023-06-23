
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
purchaseStratData <- read.csv("simSummaries_priceStrategies.csv", header = TRUE)
purchaseStratData$cvFillRate <- purchaseStratData$sdFillRate / purchaseStratData$meanFillRate
purchaseStratData$cvPercentShortfall <- purchaseStratData$sdPercentShortfall / purchaseStratData$percentShortfall

#Remove trials 31-50 (mean responses converge following ~30 simulation runs)
purchaseStratData <- purchaseStratData[which(purchaseStratData$simID < 31), ]

#Split data into pre- and post-disruption
#An initial burn-in period (i.e., the first 5 time steps) are removed from the pre-disruption data
purchaseStratData_preDis <- purchaseStratData[which(purchaseStratData$TimeStep < 51 & 
                                                      purchaseStratData$TimeStep > 5), ]
purchaseStratData_postDis <- purchaseStratData[which(purchaseStratData$TimeStep > 50), ]

#######################
#PreDisruption Analysis
#######################

purchaseStratData_preDis$ID <- as.factor(purchaseStratData_preDis$uniqueID)
purchaseStratData_preDis$StrategyF <- factor(purchaseStratData_preDis$Strategy, levels = c("Random", "pricePriority", "betHedging", "Mixed"))
purchaseStratData_preDis$cvFillRate <- purchaseStratData_preDis$sdFillRate / purchaseStratData_preDis$meanFillRate
purchaseStratData_preDis$cvPercentShortfall <- purchaseStratData_preDis$sdPercentShortfall / purchaseStratData_preDis$percentShortfall


####Mean fill rate####
meanFillRateData <- summarySE(data = purchaseStratData, 
                              measurevar = "meanFillRate", 
                              groupvars = c("Strategy", "TimeStep"))

ggplot(data = meanFillRateData, aes(x = TimeStep, y = meanFillRate, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("Mean order fill rate") +
  geom_line(size = 1.5) + 
  xlim(6,100) +
  geom_ribbon(aes(ymin = meanFillRate - (1.96*se), ymax = meanFillRate + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                       "Price Priority", "Random")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                        "Price Priority", "Random")) + 
  theme(legend.position = c(0.8, 0.9), legend.title = element_blank()) +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global models
mFR.M0 <- lme(meanFillRate ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_preDis, na.action = na.fail)

mFR.M1 <- lme(meanFillRate ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_preDis, na.action = na.fail)

#Based on AICc, mFR.M1 is preferred
#Evaluate residuals
R <- resid(mFR.M1, type = "normalized")
Ft <- fitted(mFR.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_preDis$StrategyF)
plot(R ~ purchaseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dfR.M1 <- dredge(mFR.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dfR.M1.Reduced<-subset(dfR.M1, !nested(.))

#The global model is strongly supported (weight = 0.983)
summary(mFR.M1)
intervals(mFR.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(mFR.M1, specs=~StrategyF * TimeStep, at = list(TimeStep = 10), data = purchaseStratData_preDis)
contrast(emmeans(mFR.M1, specs=~StrategyF * TimeStep, at = list(TimeStep = 10), data = purchaseStratData_preDis), method = "pairwise", adjust = "bonferroni")

emmeans(mFR.M1, specs=~StrategyF * TimeStep, at = list(TimeStep = 25), data = purchaseStratData_preDis)
contrast(emmeans(mFR.M1, specs=~StrategyF * TimeStep, at = list(TimeStep = 25), data = purchaseStratData_preDis), method = "pairwise", adjust = "bonferroni")

emmeans(mFR.M1, specs=~StrategyF * TimeStep, at = list(TimeStep = 40), data = purchaseStratData_preDis)
contrast(emmeans(mFR.M1, specs=~StrategyF * TimeStep, at = list(TimeStep = 40), data = purchaseStratData_preDis), method = "pairwise", adjust = "bonferroni")

####% shortfall####
percentShortfallData <- summarySE(data = purchaseStratData, 
                              measurevar = "percentShortfall", 
                              groupvars = c("Strategy", "TimeStep"))

ggplot(data = percentShortfallData, aes(x = TimeStep, y = percentShortfall, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("% shortfall in supply") +
  geom_line(size = 1.5) + 
  xlim(6,100) +
  geom_ribbon(aes(ymin = percentShortfall - (1.96*se), ymax = percentShortfall + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                       "Price Priority", "Random")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                        "Price Priority", "Random")) + 
  theme(legend.position = c(0.8, 0.2), legend.title = element_blank()) + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global models
pS.M1 <- lme(percentShortfall ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_preDis, method = "REML", na.action = na.fail)

pS.M0 <- lme(percentShortfall ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_preDis, na.action = na.fail)

#Based on AICc, pS.M1 is preferred
#Evaluate residuals
R <- resid(pS.M1, type = "normalized")
Ft <- fitted(pS.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_preDis$StrategyF)
plot(R ~ purchaseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dpS.M1 <- dredge(pS.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dpS.M1.Reduced<-subset(dpS.M1, !nested(.))

#Top-ranked model strongly favored (weight = 1)
pS.M2 <- lme(percentShortfall ~ StrategyF + TimeStep, random = ~1|ID, 
             weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_preDis, method = "REML", na.action = na.fail)
summary(pS.M2)
intervals(pS.M2)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(pS.M2, pairwise ~ StrategyF, adjust = "bonferroni")

#####Price per unit####
priceData <- summarySE(data = purchaseStratData, 
                                  measurevar = "meanPricePerUnit", 
                                  groupvars = c("Strategy", "TimeStep"))

ggplot(data = priceData, aes(x = TimeStep, y = meanPricePerUnit, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("Mean price per unit") +
  geom_line(size = 1.5) + 
  xlim(6,100) +
  geom_ribbon(aes(ymin = meanPricePerUnit - (1.96*se), ymax = meanPricePerUnit + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                       "Price Priority", "Random")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                        "Price Priority", "Random")) + 
  theme(legend.position = c(0.8, 0.2), legend.title = element_blank()) + ylim(15.3,18.4)+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global models
mPPU.M1 <- lme(meanPricePerUnit ~ StrategyF * TimeStep, random = ~1|ID, 
             weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_preDis, method = "REML", na.action = na.fail)

mPPU.M0 <- lme(meanPricePerUnit ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_preDis, na.action = na.fail)

#Based on AICc, mPPU.M1 is preferred
#Evaluate residuals
R <- resid(mPPU.M1, type = "normalized")
Ft <- fitted(mPPU.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_preDis$StrategyF)
plot(R ~ purchaseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dmPPU.M1 <- dredge(mPPU.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dmPPU.M1.Reduced<-subset(dmPPU.M1, !nested(.))

#Global model strongly favored (weight = 0.989)
summary(mPPU.M1)
intervals(mPPU.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 10), adjust = "bonferroni")
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 25), adjust = "bonferroni")
emmeans(mPPU.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 40), adjust = "bonferroni")

####Change in purchasing strategy####
mixedStrategyData_preDis <- purchaseStratData_preDis[which(purchaseStratData_preDis$Strategy == "Mixed" & 
                                                             purchaseStratData_preDis$TimeStep > 0), ]

mixedStrategyData <- purchaseStratData[which(purchaseStratData$Strategy == "Mixed"), ]

stratData <- summarySE(data = mixedStrategyData, 
                       measurevar = "percentPricePriority", 
                       groupvars = c("TimeStep"))

ggplot(data = stratData, aes(x = TimeStep, y = percentPricePriority, color = "#56B4E9")) + 
  xlab("Time step") + ylab("% Price priority") +
  geom_line(size = 1.5) + 
  xlim(6,100) + 
  geom_ribbon(aes(ymin = percentPricePriority - (1.96*se), ymax = percentPricePriority + (1.96*se), fill = "#56B4E9"), alpha = 0.15) +
  scale_fill_manual(values = c("#56B4E9")) + 
  scale_color_manual(values = c("#56B4E9")) +
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit model
pPP.M1 <- lme(percentPricePriority ~ TimeStep, random = ~1|ID, data = mixedStrategyData_preDis)

#Evaluate residuals
R <- resid(pPP.M1)
Ft <- fitted(pPP.M1)
qqnorm(R)
qqline(R)
plot(R ~ mixedStrategyData_preDis$TimeStep)
plot(R ~ Ft)

summary(pPP.M1)

####CV of fill rate####
cvFillRateData <- summarySE(data = purchaseStratData, 
                              measurevar = "cvFillRate", 
                              groupvars = c("Strategy", "TimeStep"))

ggplot(data = cvFillRateData, aes(x = TimeStep, y = cvFillRate, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("CV of order fill rate") +
  geom_line(size = 1.5) + 
  xlim(6,100) + 
  geom_ribbon(aes(ymin = cvFillRate - (1.96*se), ymax = cvFillRate + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                       "Price Priority", "Random")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                        "Price Priority", "Random")) + 
  theme(legend.position = c(0.8, 0.25), legend.title = element_blank()) + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global models
cvFR.M1 <- lme(cvFillRate ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_preDis, na.action = na.fail)

cvFR.M0 <- lme(cvFillRate ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_preDis, na.action = na.fail)

#Based on AICc, cvFR.M1 is preferred
#Evaluate residuals
R <- resid(cvFR.M1, type = "normalized")
Ft <- fitted(cvFR.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_preDis$StrategyF)
plot(R ~ purchaseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dcvFR.M1 <- dredge(cvFR.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dcvFR.M1.Reduced<-subset(dcvFR.M1, !nested(.))

#Global model strongly favored (weight = 0.997)
summary(cvFR.M1)
intervals(cvFR.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 10), adjust = "bonferroni")
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 25), adjust = "bonferroni")
emmeans(cvFR.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 40), adjust = "bonferroni")

####CV of % shortfall####
cvPercentShortfallData <- summarySE(data = purchaseStratData, 
                            measurevar = "cvPercentShortfall", 
                            groupvars = c("Strategy", "TimeStep"))

ggplot(data = cvPercentShortfallData, aes(x = TimeStep, y = cvPercentShortfall, color = as.factor(Strategy))) + 
  xlab("Time step") + ylab("CV of % shortfall in supply") +
  geom_line(size = 1.5) + xlim(6,100) + 
  geom_ribbon(aes(ymin = cvPercentShortfall - (1.96*se), ymax = cvPercentShortfall + (1.96*se), fill = as.factor(Strategy)), alpha = 0.15) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                       "Price Priority", "Random")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#999999", "#009E73"), labels = c("Bet-hedging", "Mixed", 
                                                                                        "Price Priority", "Random")) + 
  theme(legend.position = c(0.8, 0.19), legend.title = element_blank()) + ylim(0.5,1.6)+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) +
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1))+
  theme(legend.title = element_blank())

#Fit global modelS
cvPS.M1 <- lme(cvPercentShortfall ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_preDis, na.action = na.fail)

cvPS.M0 <- lme(cvPercentShortfall ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_preDis, na.action = na.fail)

#Based on AICc, cvPS.M1 is preferred
#Evaluate residuals
R <- resid(cvPS.M1, type = "normalized")
Ft <- fitted(cvPS.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_preDis$StrategyF)
plot(R ~ purchaseStratData_preDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dcvPS.M1 <- dredge(cvPS.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dcvPS.M1.Reduced<-subset(dcvPS.M1, !nested(.))

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

purchaseStratData_postDis$ID <- as.factor(purchaseStratData_postDis$uniqueID)
purchaseStratData_postDis$StrategyF <- factor(purchaseStratData_postDis$Strategy, levels = c("Random", "pricePriority", "betHedging", "Mixed"))
purchaseStratData_postDis$cvFillRate <- purchaseStratData_postDis$sdFillRate / purchaseStratData_postDis$meanFillRate
purchaseStratData_postDis$cvPercentShortfall <- purchaseStratData_postDis$sdPercentShortfall / purchaseStratData_postDis$percentShortfall

####Mean fill rate####

#Fit global models
mFR.M1 <- lme(meanFillRate ~ StrategyF * TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_postDis, na.action = na.fail)

mFR.M0 <- lme(meanFillRate ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_postDis, na.action = na.fail)

#Based on AICc, mFR.M1 is preferred
#Evaluate residuals
R <- resid(mFR.M1, type = "normalized")
Ft <- fitted(mFR.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_postDis$StrategyF)
plot(R ~ purchaseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dmFR.M1 <- dredge(mFR.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dmFR.M1.Reduced<-subset(dmFR.M1, !nested(.))

#Top ranked model strongly supported (weight = 1)
mFR.M2 <- lme(meanFillRate ~ StrategyF + TimeStep, random = ~1|ID, 
              weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_postDis, na.action = na.fail)
summary(mFR.M2)
intervals(mFR.M2)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(mFR.M2, pairwise ~ StrategyF, adjust = "bonferroni")

####% shortfall####

#Fit global modelS
pS.M1 <- lme(percentShortfall ~ StrategyF * TimeStep, random = ~1|ID, 
             weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_postDis, method = "REML", na.action = na.fail)

pS.M0 <- lme(percentShortfall ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_postDis, na.action = na.fail)

#Based on AICc, pS.M1 is preferred
#Evaluate residuals
R <- resid(pS.M1, type = "normalized")
Ft <- fitted(pS.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_postDis$StrategyF)
plot(R ~ purchaseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dpS.M1 <- dredge(pS.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dpS.M1.Reduced<-subset(dpS.M1, !nested(.))

#Obtain model-averaged estimates from 95% confidence set
pS.modAvg<-model.avg(dpS.M1.Reduced, fit = TRUE, subset = delta <= 5)

#Get MAEs, USEs, and 95% CIs
cbind(coefTable(pS.modAvg, revised.var=TRUE, full=TRUE), confint(pS.modAvg, revised.var=TRUE, full=TRUE))[,-3]

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(pS.modAvg, specs=~StrategyF, data = purchaseStratData_preDis)
contrast(emmeans(pS.modAvg, specs=~StrategyF, data = purchaseStratData_preDis), method = "pairwise", adjust = "bonferroni")

####Price per unit####

#Fit global models
mPPU.M1 <- lme(meanPricePerUnit ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_postDis, method = "REML", na.action = na.fail)

mPPU.M0 <- lme(meanPricePerUnit ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_postDis, na.action = na.fail)

#Based on AICc, mPPU.M1 is preferred
#Evaluate residuals
R <- resid(mPPU.M1, type = "normalized")
Ft <- fitted(mPPU.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_postDis$StrategyF)
plot(R ~ purchaseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dmPPU.M1 <- dredge(mPPU.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dmPPU.M1.Reduced<-subset(dmPPU.M1, !nested(.))

#Obtain model-averaged estimates from 95% confidence set
mPPU.modAvg<-model.avg(dmPPU.M1.Reduced, fit = TRUE, subset = delta <= 5)

#Get MAEs, USEs, and 95% CIs
cbind(coefTable(mPPU.modAvg, revised.var=TRUE, full=TRUE), confint(mPPU.modAvg, revised.var=TRUE, full=TRUE))[,-3]

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(mPPU.modAvg, specs=~StrategyF * TimeStep, at = list(TimeStep = 60), data = purchaseStratData_postDis)
contrast(emmeans(mPPU.modAvg, specs=~StrategyF * TimeStep, at = list(TimeStep = 60), data = purchaseStratData_postDis), method = "pairwise", adjust = "bonferroni")

emmeans(mPPU.modAvg, specs=~StrategyF * TimeStep, at = list(TimeStep = 75), data = purchaseStratData_postDis)
contrast(emmeans(mPPU.modAvg, specs=~StrategyF * TimeStep, at = list(TimeStep = 75), data = purchaseStratData_postDis), method = "pairwise", adjust = "bonferroni")

emmeans(mPPU.modAvg, specs=~StrategyF * TimeStep, at = list(TimeStep = 90), data = purchaseStratData_postDis)
contrast(emmeans(mPPU.modAvg, specs=~StrategyF * TimeStep, at = list(TimeStep = 90), data = purchaseStratData_postDis), method = "pairwise", adjust = "bonferroni")

####Change in purchasing strategy####
mixedStrategyData_postDis <- purchaseStratData_postDis[which(purchaseStratData_postDis$Strategy == "Mixed" & 
                                                             purchaseStratData_postDis$TimeStep > 0), ]
#Fit model
pPP.M1 <- lme(percentPricePriority ~ TimeStep, random = ~1|ID, data = mixedStrategyData_postDis)

#Evaluate residuals
R <- resid(pPP.M1)
Ft <- fitted(pPP.M1)
qqnorm(R)
qqline(R)
plot(R ~ mixedStrategyData_postDis$TimeStep)
plot(R ~ Ft)

summary(pPP.M1)

####CV of fill rate####

#Fit global models
cvFR.M1 <- lme(cvFillRate ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_postDis, na.action = na.fail)

cvFR.M0 <- lme(cvFillRate ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_postDis, na.action = na.fail)

#Based on AICc, cvFR.M1 is preferred
#Evaluate residuals
R <- resid(cvFR.M1, type = "normalized")
Ft <- fitted(cvFR.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_postDis$StrategyF)
plot(R ~ purchaseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dcvFR.M1 <- dredge(cvFR.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dcvFR.M1.Reduced<-subset(dcvFR.M1, !nested(.))

#Top ranked model strongly supported (weight = 1)
cvFR.M2 <- lme(cvFillRate ~ StrategyF + TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_postDis, na.action = na.fail)
summary(cvFR.M2)
intervals(cvFR.M2)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(cvFR.M2, pairwise ~ StrategyF, adjust = "bonferroni")

####CV of % shortfall####

#Fit global modelS
cvPS.M1 <- lme(cvPercentShortfall ~ StrategyF * TimeStep, random = ~1|ID, 
               weights = varIdent(form = ~1 | StrategyF), data = purchaseStratData_postDis, na.action = na.fail)

cvPS.M0 <- lme(cvPercentShortfall ~ StrategyF * TimeStep, random = ~1|ID, data = purchaseStratData_postDis, na.action = na.fail)

#Based on AICc, cvPS.M1 is preferred
#Evaluate residuals
R <- resid(cvPS.M1, type = "normalized")
Ft <- fitted(cvPS.M1)
qqnorm(R)
qqline(R)
plot(R ~ purchaseStratData_postDis$StrategyF)
plot(R ~ purchaseStratData_postDis$TimeStep)
plot(R ~ Ft)

#Run full model set
dcvPS.M1 <- dredge(cvPS.M1, REML = FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dcvPS.M1.Reduced<-subset(dcvPS.M1, !nested(.))

#Global model strongly favored (weight = 1)
summary(cvPS.M1)
intervals(cvPS.M1)

#Get estimated marginal means, 95% CIs, and contrasts
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 60), adjust = "bonferroni")
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 75), adjust = "bonferroni")
emmeans(cvPS.M1, pairwise ~ StrategyF * TimeStep, at = list(TimeStep = 90), adjust = "bonferroni")
