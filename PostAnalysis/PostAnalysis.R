#----------------------------------------#
#----Script to generate figures 3 - 7----#
#----------------------------------------#


library(abind)
library(dplyr)
library(ggplot2)
library(cowplot)

#List of all filenames
filenames <- list.files(path = "./SimulationOutput", pattern = "output", full.names = TRUE)

#Load first file
load(filenames[1])

#Initialize vector for all output
Out <- output$Out

#Time vector
Time <- output$Time

#Harvest parameters from files and remove model runs with Rhat > 1.1
for(i in 2:length(filenames)){
  load(filenames[i])
  if(max(output$Out[,17,]) < 1.1){
    Out <- abind(Out, output$Out, along = 1)
    Time <- c(Time, output$Time)
  }
}

#Sample 1000 iterations
iter <- sort(sample(dim(Out)[1], 1000, replace = FALSE))
Out <- Out[iter,,]

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
df <- df %>% filter(unit.size != "10")

#----------#
#-Figure 3-#
#----------#

Legend <- ggplotGrob(df %>% filter(param == "N.ds.14.mu"|param == "N.nmix.14.mu") %>% 
                       ggplot(., aes(x = unit.size, col = param)) + 
                       geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                    stat = "identity", width = 0.75, size = 0.75) +
                       scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
                       theme_bw() +
                       theme(legend.title = element_blank(),
                             legend.position = "bottom",
                             text = element_text(size = 12)))$grob[[15]]

Fig3A <- ggplotGrob(df %>% filter(param == "N.ds.14.mu"|param == "N.nmix.14.mu") %>% 
  ggplot(., aes(x = unit.size, col = param)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", width = 0.75, size = 0.75) +
  scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        text = element_text(size = 12),
        plot.margin = unit(c(1,1,0,0),"cm"),
        legend.position = "none") +
  labs(y = "Mean Abundance (N)\n", x = "", title = ""))

Fig3B <- ggplotGrob(df %>% filter(param == "N.ds.14.sd"|param == "N.nmix.14.sd") %>% 
  ggplot(., aes(x = unit.size, col = param)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", width = 0.75, size = 0.75) +
  scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        text = element_text(size = 12),
        plot.margin = unit(c(1,1,0,0),"cm"),
        legend.position = "none") +
  labs(y = "SD Abundance (N)\n", x = ""))

Fig3A$heights <- Fig3B$heights

loc <- Fig3A$layout[grep("ylab-l", Fig3A$layout$name),1:4]

Fig3A <- gtable::gtable_add_grob(Fig3A, grid::textGrob(LETTERS[1], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")

loc <- Fig3B$layout[grep("ylab-l", Fig3B$layout$name),1:4]

Fig3B <- gtable::gtable_add_grob(Fig3B, grid::textGrob(LETTERS[2], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")


png(filename = "Figure3.png", res = 600, units="cm", width=15.6, height=7.8)
gridExtra::grid.arrange(
  gridExtra::gtable_cbind(Fig3A, Fig3B),
  grid::textGrob("Sample Unit Size (km)", 
                 x = 0.5, y = 1.2, 
                 gp = grid::gpar(fontsize = 12, 
                                 fontface = 1)),
  Legend,
  nrow = 3,
  heights = c(1,0.1,0.075))
dev.off()

#----------#
#-Figure 4-#
#----------#

Fig4A <- ggplotGrob(df %>% filter(param == "N.ds.16.mu"|param == "N.nmix.16.mu") %>% 
                      ggplot(., aes(x = unit.size, col = param)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", width = 0.75, size = 0.75) +
                      scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
                      theme_bw() +
                      theme(legend.title = element_blank(),
                            text = element_text(size = 12),
                            plot.margin = unit(c(1,1,0,0),"cm"),
                            legend.position = "none") +
                      labs(y = "Mean Abundance (N)\n", x = "", title = ""))

Fig4B <- ggplotGrob(df %>% filter(param == "N.ds.16.sd"|param == "N.nmix.16.sd") %>% 
                      ggplot(., aes(x = unit.size, col = param)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", width = 0.75, size = 0.75) +
                      scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
                      theme_bw() +
                      theme(legend.title = element_blank(),
                            text = element_text(size = 12),
                            plot.margin = unit(c(1,1,0,0),"cm"),
                            legend.position = "none") +
                      labs(y = "SD Abundance (N)\n", x = ""))

Fig4A$heights <- Fig4B$heights

loc <- Fig4A$layout[grep("ylab-l", Fig4A$layout$name),1:4]

Fig4A <- gtable::gtable_add_grob(Fig4A, grid::textGrob(LETTERS[1], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")

loc <- Fig4B$layout[grep("ylab-l", Fig4B$layout$name),1:4]

Fig4B <- gtable::gtable_add_grob(Fig4B, grid::textGrob(LETTERS[2], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")


png(filename = "Figure4.png", res = 600, units="cm", width=15.6, height=7.8)
gridExtra::grid.arrange(
  gridExtra::gtable_cbind(Fig4A, Fig4B),
  grid::textGrob("Sample Unit Size (km)", 
                 x = 0.5, y = 1.2, 
                 gp = grid::gpar(fontsize = 12, 
                                 fontface = 1)),
  Legend,
  nrow = 3,
  heights = c(1,0.1,0.075))
dev.off()

#----------#
#-Figure 5-#
#----------#

Fig5A <- ggplotGrob(df %>% filter(param == "phi.ds.mu"|param == "phi.nmix.mu") %>% 
                      ggplot(., aes(x = unit.size, col = param)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", width = 0.75, size = 0.75) +
                      scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
                      theme_bw() +
                      theme(legend.title = element_blank(),
                            text = element_text(size = 12),
                            plot.margin = unit(c(1,1,0,0),"cm"),
                            legend.position = "none") +
                      labs(y = "Mean Availability\nProbability", x = "", title = ""))

Fig5B <- ggplotGrob(df %>% filter(param == "phi.ds.sd"|param == "phi.nmix.sd") %>% 
                      ggplot(., aes(x = unit.size, col = param)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", width = 0.75, size = 0.75) +
                      scale_color_manual(values = c("blue", "red"), labels = c("HDS", "N-mixture")) +
                      theme_bw() +
                      theme(legend.title = element_blank(),
                            text = element_text(size = 12),
                            plot.margin = unit(c(1,1,0,0),"cm"),
                            legend.position = "none") +
                      labs(y = "SD Availability\nProbability", x = ""))

Fig5A$heights <- Fig5B$heights

loc <- Fig5A$layout[grep("ylab-l", Fig5A$layout$name),1:4]

Fig5A <- gtable::gtable_add_grob(Fig5A, grid::textGrob(LETTERS[1], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")

loc <- Fig5B$layout[grep("ylab-l", Fig5B$layout$name),1:4]

Fig5B <- gtable::gtable_add_grob(Fig5B, grid::textGrob(LETTERS[2], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")


png(filename = "Figure5.png", res = 600, units="cm", width=15.6, height=7.8)
gridExtra::grid.arrange(
  gridExtra::gtable_cbind(Fig5A, Fig5B),
  grid::textGrob("Sample Unit Size (km)", 
                 x = 0.5, y = 1.2, 
                 gp = grid::gpar(fontsize = 12, 
                                 fontface = 1)),
  Legend,
  nrow = 3,
  heights = c(1,0.175,0.075))
dev.off()

#----------#
#-Figure 6-#
#----------#

Fig6A <- ggplotGrob(df %>% filter(param == "p.mu") %>% 
                      ggplot(., aes(x = unit.size)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", width = 0.75, size = 0.75) +
                      theme_bw() +
                      theme(text = element_text(size = 12),
                            plot.margin = unit(c(1,1,0,0),"cm")) +
                      labs(y = "Mean Detection\nProbability", x = "", title = ""))

Fig6B <- ggplotGrob(df %>% filter(param == "p.sd") %>% 
                      ggplot(., aes(x = unit.size)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", width = 0.75, size = 0.75) +
                      theme_bw() +
                      theme(text = element_text(size = 12),
                            plot.margin = unit(c(1,1,0,0),"cm")) +
                      labs(y = "SD Detection\nProbability", x = ""))

Fig6A$heights <- Fig6B$heights

loc <- Fig6A$layout[grep("ylab-l", Fig6A$layout$name),1:4]

Fig6A <- gtable::gtable_add_grob(Fig6A, grid::textGrob(LETTERS[1], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")

loc <- Fig6B$layout[grep("ylab-l", Fig6B$layout$name),1:4]

Fig6B <- gtable::gtable_add_grob(Fig6B, grid::textGrob(LETTERS[2], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")


png(filename = "Figure6.png", res = 600, units="cm", width=15.6, height=7.8)
gridExtra::grid.arrange(
  gridExtra::gtable_cbind(Fig6A, Fig6B),
  grid::textGrob("Sample Unit Size (km)", 
                 x = 0.5, y = 1.025, 
                 gp = grid::gpar(fontsize = 12, 
                                 fontface = 1)),
  nrow = 2,
  heights = c(1,0.1))
dev.off()


#----------#
#-Figure 7-#
#----------#

Fig7A <- ggplotGrob(df %>% filter(param == "pi.mu") %>% 
                      ggplot(., aes(x = unit.size)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", width = 0.75, size = 0.75) +
                      theme_bw() +
                      theme(text = element_text(size = 12),
                            plot.margin = unit(c(1,1,0,0),"cm")) +
                      labs(y = "Mean Detection\nProbability", x = "", title = ""))

Fig7B <- ggplotGrob(df %>% filter(param == "pi.sd") %>% 
                      ggplot(., aes(x = unit.size)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", width = 0.75, size = 0.75) +
                      theme_bw() +
                      theme(text = element_text(size = 12),
                            plot.margin = unit(c(1,1,0,0),"cm")) +
                      labs(y = "SD Detection\nProbability", x = ""))

Fig7A$heights <- Fig7B$heights

loc <- Fig7A$layout[grep("ylab-l", Fig7A$layout$name),1:4]

Fig7A <- gtable::gtable_add_grob(Fig7A, grid::textGrob(LETTERS[1], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")

loc <- Fig7B$layout[grep("ylab-l", Fig7B$layout$name),1:4]

Fig7B <- gtable::gtable_add_grob(Fig7B, grid::textGrob(LETTERS[2], x = 0.2, y = 1, 
                                                       gp = grid::gpar(fontsize = 12, fontface = 2)),
                                 t = loc$t, l = loc$l, clip = "inherit")


png(filename = "Figure7.png", res = 600, units="cm", width=15.6, height=7.8)
gridExtra::grid.arrange(
  gridExtra::gtable_cbind(Fig7A, Fig7B),
  grid::textGrob("Sample Unit Size (km)", 
                 x = 0.5, y = 1.025, 
                 gp = grid::gpar(fontsize = 12, 
                                 fontface = 1)),
  nrow = 2,
  heights = c(1,0.1))
dev.off()
