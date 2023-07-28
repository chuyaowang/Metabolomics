## Libraries -----

library(dplyr)
library(openxlsx)
library(ggplot2)
library(gslnls)
library(stringr)
library(ggtext)

get_LE <- function(v) {
  LE <- sum(v[2:length(v)])/sum(v)
  return(LE)
}

equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    
    if ((min(x)+d)<0) {
      seqq <- seq(0, max(x)-d, length.out=n)
    }
    else {
      seqq <- seq(0,max(x)-d,length.out=n)
      # was seq(min(x)+d,max(x)-d,length.out=n)
    }
    
    if(seqq[2]-seqq[1] < 10^(-r)) seqq else round(seqq, r)
  }
}

fit_model <- function(LE_median,Time) {
  # Reconstruct the dataframe for group_by to work
  data <- data.frame(LE_median = LE_median,
                     Time = Time)
  
  # set starting values
  a.0 <- -max(data$LE_median) *1.05
  model.0 <- lm(-log(LE_median/a.0 + 1) ~ Time, data=data)
  start <- list(a=a.0, k=coef(model.0)[2])
  
  # fit model
  model <- gsl_nls(fn = LE_median ~ a*exp(-Time*k) - a,
                   data = data,
                   start=start
  )
  
  pred <- data %>% 
    mutate(LE_pred = predict(model,Time))
  r_val <- cor(pred$LE_median,pred$LE_pred)
  
  vals <- c(coef(model)[1],coef(model)[2],r_val)
  
  return(list(model,vals))
}

my_theme <- theme_minimal() +
  theme(text = element_text(size = 10),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.text.x = element_text(margin=margin(t=1)),
        axis.title = element_text(),
        axis.title.y = element_text()) 

## Data directory -----

dir <- "data/20230724/fit_exponential.xlsx"

## Load new data -----

data_og <- read.xlsx(dir) %>%
  group_by(Metabolite, Rep,Time) %>%
  summarize(LE = get_LE(Intensity),.groups = "drop")
data <- data_og %>%
  group_by(Metabolite,Time) %>%
  summarize(LE_median = median(LE),.groups = "drop")

## Zhu paper data for aspartate -----

data_og <- read.xlsx(dir,sheet=2)
data_og <- bind_cols(data_og,str_split_fixed(data_og$Sample.name,"_",n=4)) %>%
  `[`(2:ncol(.))
colnames(data_og) <- c("LE","Metabolite","Part","Group","Time","Rep")
data_og$Time <- str_remove(data_og$Time,"h") %>% as.numeric()
data2 <- data_og %>%
  group_by(Metabolite, Time) %>%
  summarize(LE_median = median(LE),.groups = "drop")

## View data -----

ggplot(data %>% filter(Metabolite=="AKG"),aes(x=Time,y=LE_median)) +
  geom_point() +
  theme_minimal()

## Fit model for all metabolites -----

coefs <- data_test %>% 
  group_by(Metabolite) %>% 
  summarize(a = fit_model(LE_median,Time)[[2]][1],
            k = fit_model(LE_median,Time)[[2]][2],
            R = fit_model(LE_median,Time)[[2]][3])

## Fit model for one metabolite -----

metabolite_name <- "AKG"
. <- data %>% 
  filter(Metabolite == metabolite_name)
model <- fit_model(.$LE_median,.$Time)[[1]]

### Plot -----

pred <- data %>% 
  filter(Metabolite == metabolite_name) %>% 
  mutate(LE_pred = predict(model,Time))
r_val <- cor(pred$LE_median,pred$LE_pred)
lab <- paste("<p>",
             "<i>k</i> = ",
             format(unname(coef(model)[2]), scientific = F, nsmall =0, digits = 2),
             " h<sup class='sup' style='display:inline-block;vertical-align:super;font-size=6pt'>-1</sup>",
             "<br>",
             "<i>R</i> = ",
             format(unname(r_val), scientific = F, nsmall =0, digits = 2),
             "</p>",
             sep="")

lab_x <- 0.75*max(data %>% filter(Metabolite == metabolite_name) %>% select(Time))
lab_y <- quantile(data %>% filter(Metabolite == metabolite_name) %>% select(LE_median) %>% unlist)[2]

p <- ggplot(data %>% filter(Metabolite == metabolite_name),aes(x=Time,y=LE_median)) +
  geom_function(
    fun = function(x) {
      predict(model, newdata = data.frame(Time = x))
    },
    colour = "black"
  ) +
  geom_boxplot(data = data_og %>% filter(Metabolite==metabolite_name), aes(x=Time,y=LE,group=Time),outlier.size=0.5,linewidth=0.2) +
  geom_richtext(label = lab, x = lab_x,y = lab_y,size = 4,hjust = 0,
                fill = NA,label.color = NA) +
  xlab("Time (h)") +
  ylab("LE") +
  my_theme +
  scale_x_continuous(breaks = equal_breaks(n=4,s=0,r=0))
