#-Libraries-#
library(dplyr)
library(jagsUI)

#-Load data-#

data <- list()
data[[1]] <- read.csv("Aerial2014xy.csv", header = TRUE)
data[[2]] <- read.csv("Aerial2016xy.csv", header = TRUE)
data[[3]] <- read.csv("StartPoint.csv")

#-Data clean-#
data[[2]] <- data[[2]] %>%
  mutate(Distance_C = recode(Distance_C, "A/B" = "A")) %>%
  filter(Distance_C != "") %>%
  droplevels()

#-Extract data-#

#Number of years
nyrs <- 2

#Observations per year
nobs <- NULL
nobs[1] <- dim(data[[1]])[1] #2014
nobs[2] <- dim(data[[2]])[1] #2016

#-Transect breaks-#

start.point <- floor(runif(1, 0, 21000))  

unit.size <- c(1000, 2000, 4000, 5000, 6000, 8000, 10000)

S <- list()

Out <- array(NA, dim = c(1,17,7))
dimnames(Out)[[2]] <- c("N.ds.14.mu", "N.ds.16.mu", "N.nmix.14.mu", "N.nmix.16.mu",
                        "pi.mu", "p.mu", "phi.ds.mu", "phi.nmix.mu",
                        "N.ds.14.sd", "N.ds.16.sd", "N.nmix.14.sd", "N.nmix.16.sd",
                        "pi.sd", "p.sd", "phi.ds.sd", "phi.nmix.sd",
                        "Rhat")

start.time <- Sys.time()
for(u in 1:7){
  
  n.unit <- floor(21000/unit.size[u])
  len.units <- n.unit*unit.size[u]
  extra <- 21000-len.units
  spacing <- floor(extra/n.unit)
  
  breaks <- NULL
  breaks[1] <- start.point
  
  site <- list(rep(NA, nobs[1]), rep(NA, nobs[2]))
  
  for(i in 2:(n.unit+1)){
    if((breaks[i-1] + unit.size[u] + spacing) < 21000){
      breaks[i] <- breaks[i-1] + unit.size[u] + spacing
      for(t in 1:nyrs){
        for(j in 1:nobs[t]){
          start.coord <- as.numeric(data[[3]] %>% 
                                      filter(Transect == data[[t]]$Transect[j] & 
                                               Start_Sout == 1 & 
                                               StudyArea == 1) %>% 
                                      select(y))
          if((breaks[i-1] + start.coord) < data[[t]]$Y[j] &
             (breaks[i] - spacing + start.coord) >= data[[t]]$Y[j]){
            site[[t]][j] <- i - 1 + (data[[t]]$Transect[j]-1)*n.unit
          }#end if
        }#end j
      }#end t
    }else{
      breaks[i] <- unit.size[u] - (21000 - (breaks[i-1] + spacing))
      for(t in 1:nyrs){
        for(j in 1:nobs[t]){
          start.coord <- as.numeric(data[[3]] %>% 
                                      filter(Transect == data[[t]]$Transect[j] & 
                                               Start_Sout == 1 & 
                                               StudyArea == 1) %>% 
                                      select(y))
          if((breaks[i-1] + start.coord) < data[[t]]$Y[j]){
            site[[t]][j] <- i - 1 + (data[[t]]$Transect[j]-1)*n.unit
          }#end if
          if(data[[t]]$Y[j] <= (breaks[i] - spacing + start.coord)){
            site[[t]][j] <- i - 1 + (data[[t]]$Transect[j]-1)*n.unit
          }#end if
        }#end j
      }#end t
    }#end if/else
  }#end i
  obs <- list()
  obs[[1]] <- which(!is.na(site[[1]]))
  obs[[2]] <- which(!is.na(site[[2]]))
  site <- c(site[[1]], site[[2]])
  extract <- which(!is.na(site))
  site <- site[extract]
  year <- c(rep(1, nobs[1]), rep(2, nobs[2]))
  year <- year[extract]
  day <- c(data[[1]]$Survey, data[[2]]$Survey)
  day <- day[extract]
  rep <- c(as.numeric(data[[1]]$Heading_Di), as.numeric(data[[2]]$Heading_Di))
  rep <- rep[extract]
  dclass <- c(as.numeric(data[[1]]$Distance_C), as.numeric(data[[2]]$Distance_C))
  dclass <- dclass[extract]
  mdpt <- c(50, 150, 275, 425, 750)
  nK <- 5
  v <- c(100, 200, 250, 250, 500)
  B <- 1000
  ntransects <- 5
  nsites <- n.unit * ntransects
  ndays <- c(4,3)
  nreps <- max(rep)
  groups <- array(0, dim = c(nyrs, nsites, max(ndays), nreps))
  counts <- array(0, dim = c(nyrs, nsites, max(ndays), nreps))
  nobs[1] <- length(obs[[1]])
  nobs[2] <- length(obs[[2]])
  for(i in obs[[1]]){
    groups[year[i],site[i],day[i],rep[i]] <- groups[year[i],site[i],day[i],rep[i]] + 1
    counts[year[i],site[i],day[i],rep[i]] <- data[[1]]$Group_Size[i] + counts[year[i],site[i],day[i],rep[i]]
  }
  for(i in obs[[2]]){
    groups[year[i + nobs[1]],site[i + nobs[1]],day[i + nobs[1]],rep[i + nobs[1]]] <- groups[year[i + nobs[1]],site[i + nobs[1]],day[i + nobs[1]],rep[i + nobs[1]]] + 1
    counts[year[i + nobs[1]],site[i + nobs[1]],day[i + nobs[1]],rep[i + nobs[1]]] <- data[[2]]$Group_Size[i] + counts[year[i + nobs[1]],site[i + nobs[1]],day[i + nobs[1]],rep[i + nobs[1]]]
  }
  gs <- c(data[[1]]$Group_Size, data[[2]]$Group_Size)
  gs <- gs[extract]
  
  sites <- which(seq(1, n.unit*5) %in% site)
  
  bugs.data <- list(Y = groups, X = groups, gs = gs, nobs = nobs[1] + nobs[2], nyrs = nyrs, nsites = nsites, sites = sites,
                    ndays = ndays, nreps = nreps, nK = nK, mdpt = mdpt, dclass = dclass, v = v, B = B,
                    year = year, site = site, day = day, rep = rep)
  
  Nst <- round(apply(groups, c(1,2,3), max))
  Nst[2,1:nsites,4] <- NA
  Mst <- round(apply(groups, c(1,2), max))
  Nst.x <- round(apply(groups, c(1,2), max)) + 1
  
  inits <- function(){list(Mx = Mst, My = Mst, Nx = Nst, Ny = Nst, gamma = runif(1, 5, 6),
                           p = runif(1, 0.4, 0.5), alphax = runif(nyrs, 0, 2),
                           alphay = runif(nyrs, 0, 2), rx = runif(1, 0, 2),
                           ry = runif(1, 0, 2), rg = runif(1, 0, 2))}
  
  params <- c("mean.Abundx", "mean.Abundy", "alphax", "alphay", "rx", "ry", 
              "p", "gamma", "pi", "phix", "phiy", "rg", "beta")
  
  nb <- 80000
  ni <- 120000
  nt <- 10
  nc <- 3
  na <- 1000
  
  S[[u]] <- jagsUI(bugs.data, inits, params, model = "DS_Nmix_Emm.txt",
                   n.burnin = nb, n.iter = ni, n.thin = nt, 
                   n.chains = nc, n.adapt = na, parallel = TRUE)
  
  print(u)
  
  tryCatch({Out[1,1,u] <- S[[u]]$mean$mean.Abundx[1]}, error = function(e){})
  tryCatch({Out[1,2,u] <- S[[u]]$mean$mean.Abundx[2]}, error = function(e){})
  tryCatch({Out[1,3,u] <- S[[u]]$mean$mean.Abundy[1]}, error = function(e){})
  tryCatch({Out[1,4,u] <- S[[u]]$mean$mean.Abundy[2]}, error = function(e){})
  tryCatch({Out[1,5,u] <- S[[u]]$mean$pi}, error = function(e){})
  tryCatch({Out[1,6,u] <- S[[u]]$mean$p}, error = function(e){})
  tryCatch({Out[1,7,u] <- S[[u]]$mean$phix}, error = function(e){})
  tryCatch({Out[1,8,u] <- S[[u]]$mean$phiy}, error = function(e){})
  tryCatch({Out[1,9,u] <- S[[u]]$sd$mean.Abundx[1]}, error = function(e){})
  tryCatch({Out[1,10,u] <- S[[u]]$sd$mean.Abundx[2]}, error = function(e){})
  tryCatch({Out[1,11,u] <- S[[u]]$sd$mean.Abundy[1]}, error = function(e){})
  tryCatch({Out[1,12,u] <- S[[u]]$sd$mean.Abundy[2]}, error = function(e){})
  tryCatch({Out[1,13,u] <- S[[u]]$sd$pi}, error = function(e){})
  tryCatch({Out[1,14,u] <- S[[u]]$sd$p}, error = function(e){})
  tryCatch({Out[1,15,u] <- S[[u]]$sd$phix}, error = function(e){})
  tryCatch({Out[1,16,u] <- S[[u]]$sd$phiy}, error = function(e){})
  tryCatch({Out[1,17,u] <- all(unlist(S[[u]]$Rhat) < 1.1 | is.na(unlist(S[[u]]$Rhat)))}, error = function(e){})
  
}#end u

end.time <- Sys.time()
Time <- end.time - start.time

#-Save HPCC output-#
ID <- gsub(" ","_",Sys.time())
ID <- gsub(":", "-", ID)
output <- list(Out, Time)
heads <- c("Out", "Time")
output <- setNames(output, nm = heads)
save(output, file = paste("output/output", ID, ".Rdata", sep=""))