library(readr)
library(dplyr)

standard <- read.csv("outputs/standard_area_glc_glu.csv")
isocordata <- read_tsv("~/R Projects/isocordb/0224/isocordata0224.tsv")
isocorres <- read_tsv("~/R Projects/isocordb/0224/isocordata0224_res.tsv")

model.glc <- lm(conc ~ area.glc,data=standard)
model.glu <- lm(conc ~ area.glu,data=standard)

par(mfrow=c(1,2))
plot(standard$area.glc,standard$conc)
abline(model.glc)
plot(standard$area.glu,standard$conc)
abline(model.glu)

conc.glc <- predict(model.glc,
                    newdata = isocorres %>% 
                      filter(metabolite=="Glc") %>% 
                      filter(isotopologue == 0) %>% 
                      select(corrected_area) %>% 
                      mutate(area.glc = corrected_area) %>% 
                      select(area.glc)) %>% 
  `/`(isocorres %>% 
        filter(metabolite=="Glc") %>% 
        filter(isotopologue == 0) %>% 
        select(isotopologue_fraction) %>% 
        unlist)

conc.glu <- predict(model.glu,
                    newdata = isocorres %>% 
                      filter(metabolite=="Glu") %>% 
                      filter(isotopologue == 0) %>% 
                      select(corrected_area) %>% 
                      mutate(area.glu = corrected_area) %>% 
                      select(area.glu)) %>% 
  `/`(isocorres %>% 
        filter(metabolite=="Glu") %>% 
        filter(isotopologue == 0) %>% 
        select(isotopologue_fraction) %>% 
        unlist)

num.0 <- 6.9e7 # cells
num.72 <- 3.7e8
conc.glu.0 <- 0 # g/mL
conc.glu.72 <- conc.glu*10^-4
conc.glc.0 <- 0.12
# conc.glc.72 <- conc.glc*10^-4
conc.glc.72 <- 0.067 # measured value
t1 <- 0 # hour
t2 <- 72
molarmass.glc <- 180.156 # g/mol
molarmass.glu <- 143.17 

growthrate <- (log(num.72) - log(num.0))/(t2-t1) # per hour
glucoseconsumption <- growthrate*10*(conc.glc.72-conc.glc.0)*(1/molarmass.glc)*1e9/((num.72-num.0)*1e-6) # nmol/1e6 cells/hour
glutamateproduction <- growthrate*10*(conc.glu.72-conc.glu.0)*(1/molarmass.glu)*1e9/((num.72-num.0)*1e-6)
