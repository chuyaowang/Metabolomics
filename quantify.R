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
