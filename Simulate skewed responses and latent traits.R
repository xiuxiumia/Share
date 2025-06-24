tictoc::tic()
set.seed(12345)
{
  # - - - Functions
  gquantile <- function(skewness, threshold){
    quant <- quantile(x = (rgamma(1e6, shape = (2/skewness)^2, rate = 1) - (2/skewness)^2) / sqrt((2/skewness)^2),
                      probs = c(pnorm(-2), pnorm(2)))
    quant <- as.vector(quant)
    seq(quant[1], quant[2], length.out = threshold)
  } #get quantile from gamma distribution
  
  
  # - - - Get x & y threshold
  preout.x <- preout.y <- list()
  predesign.x <- expand.grid(xskew = 1.5, cat = c(3, 5, 7))
  
  for (i in 1:nrow(predesign.x)){
    preout.x[[i]] <- gquantile(predesign.x$xskew[i], predesign.x$cat[i]-1)
  }
  
  preout.x <- do.call(rbind.data.frame, lapply(preout.x, `length<-`, max(lengths(preout.x))))
  names(preout.x) <- paste0("XT", 1:ncol(preout.x), sep="")
  threshold.x <- cbind(predesign.x, preout.x)
  
  # pre matrix for x
  threshold.x.df <- threshold.x
  threshold.x.df[,3:ncol(threshold.x.df)] <- threshold.x.df[,3:ncol(threshold.x.df)] - rowMeans(threshold.x.df[,3:ncol(threshold.x.df)],na.rm=TRUE)
  
  # - - - Design Matrix and Output
  design.pre <- expand.grid(n = 10000,
                            xskew = 1.5,
                            nitems = c(5, 10, 20),
                            cat = c(3, 5, 7))
  
  design <- merge(design.pre, threshold.x.df, by = c("xskew", "cat"))
}
# ------------------------------------------------------------------------------
# - - - Run the Loop
library(dplyr)
library(eRm)
library(moments)
library(energy)
library(IRTest)

# - - - Define Variables
CONDITION <-
  xskew <-
  eskew <-
  cat <-
  n <-
  nitems <- 
  x.true <-
  x.raw <-
  x.pcm <-
  x.dcpcm <- list()

for (j in 1:nrow(design)) {
  CONDITION[[j]] <- j
  xskew[[j]] <- design$xskew[j]
  cat[[j]] <- design$cat[j]
  n[[j]] <- design$n[j]
  nitems[[j]] <- design$nitems[j]
  
  # get x true score
  x <- (rgamma(design$n[j], shape = (2/design$xskew[j])^2, rate = 1) - (2/design$xskew[j])^2) / sqrt((2/design$xskew[j])^2)
  
  x.true[[j]] <- x
  
  # get x sim data
  if(design$cat[j]-1 == 2){
    x.df <- design[j,] %>% select(num_range("XT", 1:2)) %>% as.matrix()
    x.df <- x.df[rep(seq_len(nrow(x.df)), design$nitems[j]), ]
    x.df <- cbind("difficulty" = rep(0,design$nitems[j]), x.df)
  } else if(design$cat[j]-1 == 4) {
    x.df <- design[j,] %>% select(num_range("XT", 1:4)) %>% as.matrix()
    x.df <- x.df[rep(seq_len(nrow(x.df)), design$nitems[j]), ]
    x.df <- cbind("difficulty" = rep(0,design$nitems[j]), x.df)
  } else if(design$cat[j]-1 == 6) {
    x.df <- design[j,] %>% select(num_range("XT", 1:6)) %>% as.matrix()
    x.df <- x.df[rep(seq_len(nrow(x.df)), design$nitems[j]), ]
    x.df <- cbind("difficulty" = rep(0,design$nitems[j]), x.df)
  }
  
  x.df.sim <- matrix(NA, nrow = length(x), ncol = nrow(x.df))
  
  for(k in 1:length(x)) {
    for(g in 1:nrow(x.df)) { # of items
      measure=0
      p <- vector()
      p[1] <- 1
      
      for(h in 2:ncol(x.df)) { # of the categories
        measure <- measure + x[k] - x.df[g,1] - x.df[g,h]
        p[h] <- p[(h-1)] + exp(measure)
      }
      
      U <- runif(1, 0, 1)
      U = U * p[ncol(x.df)]
      
      for(t in 1:ncol(x.df)) { # number of categories
        
        if(U <= p[t]) {
          x.df.sim[k,g] <- (t-1)
          break
        }
      }
    }
  }
  
  # --- get raw scores
  x.raw[[j]] <- rowSums(x.df.sim)
  
  x.df.sim <- as.data.frame(x.df.sim)
  
  # --- get PCM estimates
  x.pcm.mod <- PCM(x.df.sim)
  x.pcm[[j]] <- as.vector(person.parameter(x.pcm.mod)$theta.table[,1])
  
  
  # --- get DC-PCM estimates
  npar <- 10
  dav_x <- vector("list", npar)
  for (p in 1:npar) {
    dav_x[[p]] <- IRTest_Poly(x.df.sim, latent_dist = "DC", h = p, model = "PCM")
  }
  HQ.x <- best_model(dav_x[[1]], dav_x[[2]], dav_x[[3]], dav_x[[4]], dav_x[[5]],
                     dav_x[[6]], dav_x[[7]], dav_x[[8]], dav_x[[9]], dav_x[[10]],
                     criterion = "HQ")$criterion$HQ
  x.dcpcm.mod <- dav_x[[which.min(HQ.x)]]
  x.dcpcm[[j]] <- factor_score(x.dcpcm.mod, ability_method = "EAP")$theta
}

Example_Sim_Result <- list(
  CONDITION = CONDITION,
  xskew = xskew,
  cat = cat,
  n = n,
  nitems = nitems,
  x.true = x.true,
  x.raw = x.raw,
  x.pcm = x.pcm,
  x.dcpcm = x.dcpcm
)

save(Example_Sim_Result, file = paste("Sim_Example1", ".RData", sep = ""))
tictoc::toc() # 2947.937 sec elapsed