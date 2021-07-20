install.packages('lmtest')
require(lmtest)


condSurvival <- tbl_df(read.table("../processed/conditional_survival_analysis.csv", 
                  sep = ",", header = TRUE, 
                  stringsAsFactors = FALSE))

write.table(condSurvival %>% select(-n, -km) %>% filter(treatmentType != "total"), 
            "../tables/table_BC_conditional_survival.csv", sep = ",", col.names = TRUE,
            row.names = FALSE, quote = FALSE)


linTest <- condSurvival %>% mutate(point.estimate = 1 - point.estimate, 
                        X95_conf_low = 1 - X95_conf_low, 
                        X95_conf_high = 1 - X95_conf_high) %>% 
                            filter(treatmentType != "total")


expfits <- foreach(tt = c("3TC", "PI", "NNRTI", "NRTI"), .combine = "rbind")%do%{

    dat <- linTest %>% filter(treatmentType == tt, Year > 1)  %>% select(Year, point.estimate) 
    exp.model <- (lm(log(point.estimate) ~ Year, data= dat))

    timevalues <- seq(1, 10, 0.1)
    Counts.exponential <- exp(predict(exp.model,list(Year=timevalues)))
    fits <- tbl_df(cbind(timevalues, Counts.exponential))
    colnames(fits) <- c("Year", "point.estimate")

    bind_rows(dat %>% mutate(type = "dat"),
              fits %>% mutate(type = "fits")) %>% mutate(treat = tt, r = (summary(exp.model))['r.squared']$r.squared, pval = coef(summary(exp.model))[2, 4]) 
    
}



offset <- 0.05
g.a <- condSurvival %>% mutate(point.estimate = 1 - point.estimate, 
                        X95_conf_low = 1 - X95_conf_low, 
                        X95_conf_high = 1 - X95_conf_high) %>% 
                            filter(treatmentType != "total") %>% 
    mutate(Year = ifelse(treatmentType == "3TC", Year - 1.5 * offset, 
               ifelse(treatmentType == "NNRTI", Year - 0.5 * offset, 
               ifelse(treatmentType == "NRTI", Year + 0.5 * offset, 
               ifelse(treatmentType == "PI", Year + 1.5 * offset, NA))))) %>%
    ggplot() + 
    geom_point(aes(x =  Year, y = point.estimate, col = treatmentType)) + 
    geom_segment(aes(x = Year , xend = Year, y = X95_conf_low, yend = X95_conf_high,
                     col = treatmentType, group = treatmentType)) +       
    scale_color_manual(values = cols) + 
    labs(x = "Year of treatment", y = "Probability of resistance") + 
    theme_classic() + theme(legend.position="bottom") + 
    theme(legend.title = element_blank()) + 
    scale_x_continuous(breaks = seq(0, 10, by = 2), labels = seq(0, 10, by = 2)) +  
    geom_path(aes(x = Year, y = point.estimate, col = treat, group = treat), 
              data = expfits %>% arrange(Year) %>% filter(type == "fits", Year > 1.9)) + 
    geom_path(aes(x = Year, y = point.estimate, col = treat, group = treat), 
              data = expfits %>% arrange(Year) %>% filter(type == "fits", Year < 2.1), lty = "dashed") + 
                  geom_richtext(aes(x = 10, y = 0.1 - ya, 
                            label = paste0("*r<sup>2</sup>* = ", signif(r, 2)), col = treat), show.legend = FALSE, 
data = expfits %>% group_by(treat) %>% slice(n = 1) %>% mutate(ya = 0.0075 * which(treat == unique(expfits$treat))), fill = NA, label.color = NA, hjust = 1)




expfits <- foreach(tt = c("3TC", "PI", "NNRTI", "NRTI"), .combine = "rbind")%do%{

    dat <- linTest %>% filter(treatmentType == tt, Year > 1)  %>% select(Year, point.estimate) 
    exp.model <- (lm(log(point.estimate) ~ Year, data= dat))

    tbl_df(coef(summary(exp.model)) ) %>% mutate(Param = c("Intercept", "Year"), r2 = (summary(exp.model))['r.squared']$r.squared, treat = tt)
    
}


tabToPrint <- expfits %>% mutate(Estimate = ifelse(Param == "Intercept", exp(Estimate), Estimate)) %>% 
    select(treat, r2, Param, Estimate, `Std. Error`,  `t value`, `Pr(>|t|)`) %>% 
    mutate(Param = ifelse(Param == "Intercept", "a", "b"))

write.table(tabToPrint, 
            "../tables/table_BC_exponential_fits.csv", sep = ",", col.names = TRUE,
            row.names = FALSE, quote = FALSE)





