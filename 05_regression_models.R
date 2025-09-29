# Load required packages
library(ggplot2)
library(lme4)
library(mgcv)
library(gamm4)
setwd("/Users/mehrnaz/Documents/Thesis/Data")

#---- 1) different sides

prediction <- expand.grid(
                          location_sl = seq(min(by_bar$location_sl), 
                          max(by_bar$location_sl), length.out = 100),
                          treatment = levels(by_bar$treatment)
                          )

prediction$treatment <- as.factor(prediction$treatment)
str(prediction)


# GAMM model with interaction

by_bar$treatment <- as.factor(by_bar$treatment)
by_bar$image <- as.factor(by_bar$image)

gamm_model_interaction <- gamm4(
        offset_diffside ~ treatment + s(location_sl, by = treatment, k = 8),
        random = ~(1|image),
        data   = by_bar,
        family = gaussian
)

# Extracting predictions and confidence intervals
prediction$gamm_pred_int <- predict(gamm_model_interaction$gam, newdata = prediction, type = "response", se.fit = TRUE)

prediction$gamm_pred_int_fit <- prediction$gamm_pred_int$fit
prediction$gamm_pred_int_se <- prediction$gamm_pred_int$se.fit

prediction$lower <- prediction$gamm_pred_int_fit - 1.96 * prediction$gamm_pred_int_se
prediction$upper <- prediction$gamm_pred_int_fit + 1.96 * prediction$gamm_pred_int_se




summary(gamm_model_interaction$gam)
summary(gamm_model_interaction$mer)





opposit <- ggplot(by_bar, aes(x = location_sl, 
                   y = offset_diffside, 
                   group = treatment,
                   shape = brood, 
                   fill = treatment)) +
        geom_point(aes(color = treatment, group = treatment), alpha = 0.2, size = 1) +
        
        scale_color_manual(values = c("warm" = "tomato", "cool" = "cornflowerblue")) +
        scale_shape_manual(values = c(21:25), name = "Brood")+
        
        geom_line(data = prediction, 
                  aes(x = location_sl, y = gamm_pred_int_fit, color = treatment, group = treatment), 
                  size = 1, inherit.aes = FALSE) +
        
        geom_ribbon(data = prediction, 
                    aes(x = location_sl, ymin = lower, 
                        ymax = upper, fill = treatment, group = treatment), 
                    alpha = 0.2, inherit.aes = FALSE) +
        
        scale_fill_manual(values = c("warm" = "tomato2", "cool" = "#193c98")) +
        
        
        labs(
                #title = "GAM Fit with 95% Confidence Intervals",
                x = "Location on Standard Length",
                y = "Left-Right Offset",
        ) +
        theme_classic(base_size = 16)





#same sides

prediction <- expand.grid(
        location_sl = seq(min(by_bar$location_sl), 
                          max(by_bar$location_sl), length.out = 100),
        treatment = levels(by_bar$treatment))
prediction$treatment <- as.factor(prediction$treatment)
str(prediction)

# Linear Mixed Model (if random effects needed)
lmm_model2 <- lmer(offset_sameside ~ treatment + (1|image), data = by_bar)

summary(lmm_model2)

prediction$lmm_pred <- predict(lmm_model2, newdata = prediction2, re.form = NA)

#GAMM model with interaction

by_bar$treatment <- as.factor(by_bar$treatment)
by_bar$image <- as.factor(by_bar$image)

gamm_model_interaction <- gamm_model_interaction <- gamm4(
        offset_sameside ~ treatment + s(location_sl, by = treatment, k = 8),
        random = ~(1|image),
        data   = by_bar,
        family = gaussian
)

# Extracting predictions and confidence intervals
prediction$gamm_pred_int <- predict(gamm_model_interaction$gam, newdata = prediction, type = "response", se.fit = TRUE)

prediction$gamm_pred_int_fit <- prediction$gamm_pred_int$fit
prediction$gamm_pred_int_se <- prediction$gamm_pred_int$se.fit

prediction$lower <- prediction$gamm_pred_int_fit - 1.96 * prediction$gamm_pred_int_se
prediction$upper <- prediction$gamm_pred_int_fit + 1.96 * prediction$gamm_pred_int_se




summary(gamm_model_interaction$gam)
summary(gamm_model_interaction$mer)



sameside <- ggplot(by_bar, aes(x = location_sl, 
                   y = offset_sameside, 
                   group = treatment,
                   color = treatment,
                   shape = brood,
                   fill = treatment)) +
        geom_point(aes(color = treatment), alpha = 0.2, size = 1) +
        
        scale_color_manual(values = c("warm" = "tomato","cool" = "cornflowerblue")) +
        
        scale_shape_manual(values = c(21:25), name = "Brood")+
        
        geom_line(data = prediction, 
                  aes(x = location_sl, y = gamm_pred_int_fit, color = treatment, group = treatment), 
                  size = 1, inherit.aes = FALSE) +
        
        geom_ribbon(data = prediction, 
                  aes(x = location_sl, ymin = lower, 
                     ymax = upper, group = treatment, fill = treatment), 
                   alpha = 0.2, inherit.aes = FALSE) +
        scale_fill_manual(values = c("warm" = "tomato2", "cool" = "#193c98")) +
        
        labs(
                #title = "Offset vs. anterior-posterior location",
                x = "Location on Standard Length",
                y = "Same Side Offset",
        ) +
        theme_classic(base_size = 16)

#---- 2) GLM

library(ggbeeswarm)

inputDF <- read.csv("by_bar2025-09-17.csv", header = TRUE)
aggregMeta <- aggregate(
        cbind(brood,treatment)         # the columns to aggregate; can be multiple
        ~ image,            # formula: “by this grouping variable”
        data = inputDF,
        FUN  = function(x) x[1]
)
ggplot(inputDF, mapping = aes(x = location_sl, y = offset_diffside, color = treatment))+
        geom_point()+ geom_smooth()
ggplot(inputDF, mapping = aes(x = location_sl, y = offset_sameside, color = treatment))+
        geom_point()+ geom_smooth()
ggplot(inputDF, mapping = aes(x = location_sl, y = offset_diffside, color = image))+
        geom_point()+ geom_smooth(se=F)+
        scale_color_manual(values = ifelse(aggregMeta$treatment=="cool","blue","red"))
ggplot(inputDF, mapping = aes(x = location_sl, y = offset_sameside, color = image))+
        geom_point()+ geom_smooth(se=F)+
        scale_color_manual(values = ifelse(aggregMeta$treatment=="cool","blue","red"))
aggrByRange<-aggregate(
        cbind(offset_diffside, offset_sameside)         # the columns to aggregate; can be multiple
        ~ image,            # formula: “by this grouping variable”
        data = inputDF,
        simplify = TRUE,
        FUN  = function(x) max(x) - min(x)
); colnames(aggrByRange) <- c("image", "offset_diffside_range", "offset_sameside_range")
aggrByDeviation<-aggregate(
        cbind(offset_diffside, offset_sameside)         # the columns to aggregate; can be multiple
        ~ image,            # formula: “by this grouping variable”
        data = inputDF,
        simplify = TRUE,
        FUN  = function(x) mad(x)
); colnames(aggrByDeviation) <- c("image", "offset_diffside_dev", "offset_sameside_dev")
aggrByMedian<-aggregate(
        cbind(offset_diffside, offset_sameside)         # the columns to aggregate; can be multiple
        ~ image,            # formula: “by this grouping variable”
        data = inputDF,
        simplify = TRUE,
        FUN  = function(x) median(x)
); colnames(aggrByMedian) <- c("image", "offset_diffside_med", "offset_sameside_med")
resDF <- data.frame(aggregMeta,aggrByRange[,c(2,3)], aggrByDeviation[,c(2,3)],aggrByMedian[,c(2,3)],check.names = F)
ggplot(resDF,aes(x=treatment, y = offset_diffside_dev))+
        geom_boxplot(width = 0.2, outlier.shape = NA)+
        geom_beeswarm(size = 2,, mapping = aes(color = brood))
ggplot(resDF,aes(x=treatment, y = offset_diffside_med))+
        geom_boxplot(width = 0.2, outlier.shape = NA)+
        geom_beeswarm(size = 2,, mapping = aes(color = brood))
glm1 <- lm(offset_diffside_med ~ treatment + brood, data = resDF)
summary(glm1)
ggplot(resDF,aes(x=treatment, y = offset_sameside_dev))+
        geom_boxplot(width = 0.2, outlier.shape = NA)+
        geom_beeswarm(size = 2,, mapping = aes(color = brood))
ggplot(resDF,aes(x=treatment, y = offset_sameside_med))+
        geom_boxplot(width = 0.2, outlier.shape = NA)+
        geom_beeswarm(size = 2,, mapping = aes(color = brood))
glm2 <- lm(offset_sameside_med ~ treatment + brood, data = resDF)
summary(glm2)
#### binning

library(ggbeeswarm)
library(ggplot2)
inputDF <- read.csv("by_bar2025-09-17.csv", header = TRUE)
inputDF$locationBin<-cut(
        inputDF$location_sl,
        breaks = c(0, 1/3, 2/3, 1), # the breakpoints
        labels = c("head", "mid", "tail"), # the labels
        right = FALSE # whether the intervals should be closed on the right
)
resDFbin<-aggregate(
        cbind(offset_diffside, offset_sameside)         # the columns to aggregate; can be multiple
        ~ image + locationBin + treatment + brood,            # formula: “by this grouping variable”
        data = inputDF,
        simplify = TRUE,
        FUN  = function(x) median(x)
); colnames(resDFbin) <- c("image","locationBin","treatment","brood", "offset_diffside_med", "offset_sameside_med")

ggplot(resDFbin,aes(x=treatment, y = offset_diffside_med))+
        geom_boxplot(width = 0.2, outlier.shape = NA)+
        geom_beeswarm(size = 2.5,cex = 3, mapping = aes(color = brood))+
        facet_wrap(~locationBin) + ylim(0,.15)+
        theme_bw(base_size = 16)


ggplot(resDFbin,aes(x=treatment, y = offset_sameside_med))+
        geom_boxplot(width = 0.2, outlier.shape = NA)+
        geom_beeswarm(size = 1.5,cex = 3, mapping = aes(color = brood))+
        facet_wrap(~locationBin) + ylim(0,.15)


resDFbin$ishead<-resDFbin$locationBin == "head"
resDFbin$ismid<-resDFbin$locationBin == "mid"
resDFbin$istail<-resDFbin$locationBin == "tail"
glm_diffside <- lm(offset_diffside_med ~ treatment+  treatment:ishead + treatment:ismid + treatment:istail  + brood, data = resDFbin)


summary(glm_diffside)
summary_df_diffside  <- as.data.frame(summary(glm_diffside)$coefficients)
colnames(summary_df_diffside) <- c("Estimate", "StdError", "Zvalue", "Pvalue")
summary_df_diffside
drop1(glm_diffside, test = "F")

#tail NA because only one measure BR07 / treatment
glm_sameside <- lm(offset_sameside_med ~ treatment + treatment:ishead + treatment:ismid + treatment:istail + brood, data = resDFbin)
summary(glm_sameside)
summary_df_sameside <- as.data.frame(summary(glm_sameside)$coefficients)
colnames(summary_df_sameside) <- c("Estimate", "StdError", "Zvalue", "Pvalue")
summary_df_sameside
drop1(glm_sameside, test = "F")








GLMoffset_LR <- ggplot(resDFbin, aes(x = treatment,y = offset_diffside_med)) +
        geom_boxplot(aes(fill = factor(treatment)), alpha = 0.2, outliers = F)+                      geom_beeswarm(aes(color = treatment, shape = brood, fill = treatment),
        size = 1.75, cex = 3, alpha = 0.5) +
        facet_wrap(~locationBin) + 
        ylim(0, .15) +
        scale_color_manual(values = c("warm" = "tomato", "cool" = "cornflowerblue")) +
        scale_fill_manual(values = c("warm" = "tomato2", "cool" = "cornflowerblue")) +
        scale_shape_manual(values = 21:25, name = "Brood") +
        theme_classic(base_size = 16,base_line_size = 0.5)



GLMoffset_DV <- ggplot(resDFbin, aes(x = treatment, 
                                    y = offset_sameside_med)) +
        
        geom_boxplot(aes(fill = factor(treatment)), 
                     alpha = 0.2, outliers = F)+        
        geom_beeswarm(aes(color = treatment, 
                          shape = brood, 
                          fill = treatment),
                      size = 1.75, alpha = 0.5) +
        facet_wrap(~locationBin) + 
        ylim(0, .15) +
        
        scale_color_manual(values = c("warm" = "tomato", 
                                      "cool" = "cornflowerblue")) +
        
        scale_fill_manual(values = c("warm" = "tomato2", 
                                     "cool" = "cornflowerblue")) +
        
        scale_shape_manual(values = 21:25, name = "Brood") +
        
        theme_classic(base_size = 16,base_line_size = 0.5)


