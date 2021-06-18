library(abind)
library(dplyr)
library(ggplot2)
library(cowplot)

#List of all filenames
filenames <- list.files(path = "./revision", pattern = "output", full.names = TRUE)

#Load first file
load(filenames[1])

#Initialize vector for all output
Out <- output$Out

#Time vector
Time <- output$Time

#Harvest parameters from files and remove model runs with Rhat > 1.1
for(i in 2:length(filenames)){
  load(filenames[i])
  if(all(output$Out[,17,] == 1)){
    Out <- abind(Out, output$Out, along = 1)
    Time <- c(Time, output$Time)
  }
}

#Sample 1000 iterations
iter <- sort(sample(dim(Out)[1], 1000, replace = FALSE))
Out <- Out[iter,,]

df <- matrix(NA, nrow = 1, ncol = 5)
# Rhat <- matrix(NA, nrow = 7, ncol = 2)

for(u in 1:7){
  # rhat <- Out[,17,u]
  tmp <- Out[,1:16,u]
  y75 <- apply(X = tmp, MARGIN = 2, FUN = quantile, probs = 0.75, na.rm = TRUE)
  y50 <- apply(X = tmp, MARGIN = 2, FUN = quantile, probs = 0.5, na.rm = TRUE)
  y25 <- apply(X = tmp, MARGIN = 2, FUN = quantile, probs = 0.25, na.rm = TRUE)
  ymax <-  ((y75 - y25) * 1.5) + y75
  ymin <- y25 - ((y75 - y25) * 1.5)
  df <- rbind(df, cbind(cbind(ymin, y25, y50, y75, ymax)))
  # Rhat[u,1] <- sum(rhat)
  # Rhat[u,2] <- sum(rhat)/length(rhat)
}

df <- df[-1,]
names <- as.factor(rownames(df))
df <- as.data.frame(df)
df$param <- names
rownames(df) <- NULL
df$unit.size <- rep(c(1, 2, 4, 5, 6, 8, 10), each = 16)
df$unit.size <- as.factor(df$unit.size)

df1 <- df
df <- df %>% filter(unit.size != "10")

Fig1 <- df %>% filter(param == "N.ds.14.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km) 2014", title = "HDS")
#Fig1

Fig2 <- df %>% filter(param == "N.ds.16.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km) 2016", title = "HDS")
#Fig2

Fig3 <- df %>% filter(param == "N.nmix.14.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km) 2014", title = "Nmix")
#Fig3

Fig4 <- df %>% filter(param == "N.nmix.16.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km) 2016", title = "Nmix")
#Fig4

Fig5 <- df %>% filter(param == "N.ds.14.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km) 2014", title = "HDS")
#Fig5

Fig6 <- df %>% filter(param == "N.ds.16.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km) 2016", title = "HDS")
#Fig6

Fig7 <- df %>% filter(param == "N.nmix.14.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km) 2014", title = "Nmix")
#Fig7

Fig8 <- df %>% filter(param == "N.nmix.16.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km) 2016", title = "Nmix")
#Fig8

Fig9 <- df %>% filter(param == "pi.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Detection Probability", x = "Sample Unit Size (km)", title = "HDS")
#Fig9

Fig10 <- df %>% filter(param == "p.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Detection Probability", x = "Sample Unit Size (km)", title = "Nmix")
#Fig10

Fig11 <- df %>% filter(param == "pi.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Detection Probability", x = "Sample Unit Size (km)", title = "HDS")
#Fig11

Fig12 <- df %>% filter(param == "p.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Detection Probability", x = "Sample Unit Size (km)", title = "Nmix")
#Fig12

Fig13 <- df %>% filter(param == "phi.ds.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Availability Probability", x = "Sample Unit Size (km)", title = "HDS")
#Fig13

Fig14 <- df %>% filter(param == "phi.nmix.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Availability Probability", x = "Sample Unit Size (km)", title = "Nmix")
#Fig14

Fig15 <- df %>% filter(param == "phi.ds.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Availability Probability", x = "Sample Unit Size (km)", title = "HDS")
#Fig15

Fig16 <- df %>% filter(param == "phi.nmix.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Availability Probability", x = "Sample Unit Size (km)", title = "Nmix")
#Fig16

plot_grid(Fig1, Fig3, Fig2, Fig4, Fig5, Fig7, Fig6, Fig8, 
          Fig9, Fig10, Fig13, Fig14, Fig11, Fig12, Fig15, Fig16,
          ncol = 4)

#------------------#
#-Individual plots-#
#------------------#

png(filename = "N_ds_14_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.ds.14.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km)")
dev.off()
#Fig1

png(filename = "N_ds_16_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.ds.16.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km)")
dev.off()
#Fig2

png(filename = "N_nmix_14_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.nmix.14.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km)")
dev.off()
#Fig3

png(filename = "N_nmix_16_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.nmix.16.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km)")
dev.off()
#Fig4

png(filename = "N_ds_14_sd.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.ds.14.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km)")
dev.off()
#Fig5

png(filename = "N_ds_16_sd.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.ds.16.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km)")
dev.off()
#Fig6

png(filename = "N_nmix_14_sd.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.nmix.14.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km)")
dev.off()
#Fig7

png(filename = "N_nmix_16_sd.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.nmix.16.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km)")
dev.off()
#Fig8

png(filename = "pi_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "pi.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Detection Probability", x = "Sample Unit Size (km)")
dev.off()
#Fig9

png(filename = "p_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "p.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "Mean Detection Probability", x = "Sample Unit Size (km)")
dev.off()
#Fig10

png(filename = "pi_sd.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "pi.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Detection Probability", x = "Sample Unit Size (km)")
dev.off()
#Fig11

png(filename = "p_sd.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "p.sd") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge") +
  theme_bw() +
  labs(y = "SD Detection Probability", x = "Sample Unit Size (km)")
dev.off()
#Fig12

png(filename = "lam_ds_14_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "lam.ds.14.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "Mean Number of Groups \n (lambda)", x = "Sample Unit Size (km)")
dev.off()
#Fig1

png(filename = "lam_ds_16_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "lam.ds.16.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "Mean Number of Groups \n (lambda)", x = "Sample Unit Size (km)")
dev.off()
#Fig2

png(filename = "lam_nmix_14_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "lam.nmix.14.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "blue") +
  theme_bw() +
  labs(y = "Mean Number of Groups \n (lambda)", x = "Sample Unit Size (km)")
dev.off()
#Fig3

png(filename = "lam_nmix_16_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "lam.nmix.16.mu") %>% ggplot(., aes(x = unit.size)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", size = 0.75, position = "dodge", col = "red") +
  theme_bw() +
  labs(y = "Mean Number of Groups \n (lambda)", x = "Sample Unit Size (km)")
dev.off()
#Fig4

#-Combined types-#

png(filename = "N_14_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.ds.14.mu"|param == "N.nmix.14.mu") %>% 
  ggplot(., aes(x = unit.size, col = param)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", width = 0.75, size = 0.75) +
  scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km)")
dev.off()

png(filename = "N_16_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.ds.16.mu"|param == "N.nmix.16.mu") %>% 
  ggplot(., aes(x = unit.size, col = param)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", width = 0.75, size = 0.75) +
  scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(y = "Mean Abundance (N)", x = "Sample Unit Size (km)")
dev.off()

png(filename = "N_14_sd.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.ds.14.sd"|param == "N.nmix.14.sd") %>% 
  ggplot(., aes(x = unit.size, col = param)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", width = 0.75, size = 0.75) +
  scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km)")
dev.off()

png(filename = "N_16_sd.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "N.ds.16.sd"|param == "N.nmix.16.sd") %>% 
  ggplot(., aes(x = unit.size, col = param)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", width = 0.75, size = 0.75) +
  scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(y = "SD Abundance (N)", x = "Sample Unit Size (km)")
dev.off()

png(filename = "lam_14_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "lam.ds.14.mu"|param == "lam.nmix.14.mu") %>% 
  ggplot(., aes(x = unit.size, col = param)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", width = 0.75, size = 0.75) +
  scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(y = "Mean Number of Groups \n (lambda)", x = "Sample Unit Size (km)")
dev.off()

png(filename = "lam_16_mu.png", res = 600, units="in", width=5, height=5)
df %>% filter(param == "lam.ds.16.mu"|param == "lam.nmix.16.mu") %>% 
  ggplot(., aes(x = unit.size, col = param)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", width = 0.75, size = 0.75) +
  scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(y = "Mean Number of Groups \n (lambda)", x = "Sample Unit Size (km)")
dev.off()

#-----------------#
#-Load 2km output-#
#-----------------#

load(file = "output2km.Rdata")

#--------#
#-Figure-#
#--------#

df <- data.frame(out$summary[1:4, c(3,4,1,6,7)],
                 rep(c("2014", "2016"), 2), 
                 rep(c("HDS", "N-mixture"), each = 2))

colnames(df) <- c("q2.5", "q25", "mean", "q75", "q97.5", "Year", "Model")

Figure3 <- ggplot(df) + 
  geom_errorbar(aes(x = Year, ymin = q2.5, ymax = q97.5, group = Model), 
                width = 0.1, size = 1.25, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = Year, y = mean, group = Model), 
             size = 5, shape = 18, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = Year, y = mean, col = Model), 
             size = 4, shape = 18, position = position_dodge(width = 0.5)) +
  theme_bw() +
  labs(y = "Abundance")

png(filename = "Figure3.png", res = 600, units="in", width=5, height=5)
Figure3
dev.off()

