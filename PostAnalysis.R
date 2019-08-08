library(abind)
library(dplyr)
library(ggplot2)
library(cowplot)

#List of all filenames
filenames <- list.files(pattern = "output", full.names = TRUE)

#Load first file
load(filenames[1])

#Initialize vector for all output
Out <- output$Out

#Time vector
Time <- output$Time

#Harvest parameters from files and remove model runs with Rhat > 1.1
for(i in 2:length(filenames)){
  load(filenames[i])
  Out <- abind(Out, output$Out, along = 1)
  Time <- c(Time, output$Time)
}

#Sample 1000 iterations
# iter <- sort(sample(dim(Out)[1], 1000, replace = FALSE))
# Out <- Out[iter,,]

df <- matrix(NA, nrow = 1, ncol = 5)
Rhat <- matrix(NA, nrow = 7, ncol = 2)

for(u in 1:7){
  rhat <- Out[,17,u]
  tmp <- Out[rhat==1,1:16,u]
  y75 <- apply(X = tmp, MARGIN = 2, FUN = quantile, probs = 0.75, na.rm = TRUE)
  y50 <- apply(X = tmp, MARGIN = 2, FUN = quantile, probs = 0.5, na.rm = TRUE)
  y25 <- apply(X = tmp, MARGIN = 2, FUN = quantile, probs = 0.25, na.rm = TRUE)
  ymax <-  ((y75 - y25) * 1.5) + y75
  ymin <- y25 - ((y75 - y25) * 1.5)
  df <- rbind(df, cbind(cbind(ymin, y25, y50, y75, ymax)))
  Rhat[u,1] <- sum(rhat)
  Rhat[u,2] <- sum(rhat)/length(rhat)
}

df <- df[-1,]
names <- as.factor(rownames(df))
df <- as.data.frame(df)
df$param <- names
rownames(df) <- NULL
df$unit.size <- rep(c(1, 2, 4, 5, 6, 8, 10), each = 16)
df$unit.size <- as.factor(df$unit.size)

df1 <- df
df <- df %>% filter(unit.size != "10000")

Fig1 <- df %>% filter(param == "N.ds.14.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Spatial Unit Size (km) 2014", title = "DS")
#Fig1

Fig2 <- df %>% filter(param == "N_ds.16.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Spatial Unit Size (km) 2016", title = "DS")
#Fig2

Fig3 <- df %>% filter(param == "N.nmix.14.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Spatial Unit Size (km) 2014", title = "Nmix")
#Fig3

Fig4 <- df %>% filter(param == "N.nmix.16.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Spatial Unit Size (km) 2016", title = "Nmix")
#Fig4

Fig5 <- df %>% filter(param == "N.ds.14.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Spatial Unit Size (km) 2014", title = "DS")
#Fig5

Fig6 <- df %>% filter(param == "N_ds.16.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Spatial Unit Size (km) 2016", title = "DS")
#Fig6

Fig7 <- df %>% filter(param == "N.nmix.14.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Spatial Unit Size (km) 2014", title = "Nmix")
#Fig7

Fig8 <- df %>% filter(param == "N.nmix.16.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Spatial Unit Size (km) 2016", title = "Nmix")
#Fig8

Fig9 <- df %>% filter(param == "pi.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Detection Probability", x = "Spatial Unit Size (km)", title = "DS")
#Fig9

Fig10 <- df %>% filter(param == "p.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Detection Probability", x = "Spatial Unit Size (km)", title = "Nmix")
#Fig10

Fig11 <- df %>% filter(param == "pi.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Detection Probability", x = "Spatial Unit Size (km)", title = "DS")
#Fig11

Fig12 <- df %>% filter(param == "p.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Detection Probability", x = "Spatial Unit Size (km)", title = "Nmix")
#Fig12

Fig13 <- df %>% filter(param == "phi.ds.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Availability Probability", x = "Spatial Unit Size (km)", title = "DS")
#Fig13

Fig14 <- df %>% filter(param == "phi.nmix.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Availability Probability", x = "Spatial Unit Size (km)", title = "Nmix")
#Fig14

Fig15 <- df %>% filter(param == "phi.ds.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Availability Probability", x = "Spatial Unit Size (km)", title = "DS")
#Fig15

Fig16 <- df %>% filter(param == "phi.nmix.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Availability Probability", x = "Spatial Unit Size (km)", title = "Nmix")
#Fig16

plot_grid(Fig1, Fig3, Fig2, Fig4, Fig5, Fig7, Fig6, Fig8, 
          Fig9, Fig10, Fig13, Fig14, Fig11, Fig12, Fig15, Fig16,
          ncol = 4)
