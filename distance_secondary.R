library(vroom)
library(foreach)
library(doParallel)
library(doSNOW)
library(ffbase)
library(tidyverse)


stepwise_distance_w_plot <- function(p1, q1, #prices and quantities for curve 1
                                     p2, q2,#prices and quantities for curve 2
                                     l = 1, #l of the distance
                                     w_type = c("none", "max_price", "user_wfun"), #type of w(p)
                                     maxprice, # only for step_max. Maximum price considered.
                                     user_fun, # only for user_fun. Function that integrates up to one.
                                     plot = TRUE){
  
  # Define the two supply functions as R functions
  sfun1  <- stepfun(p1, q1, f = 1, right = T)
  sfun2  <- stepfun(p2, q2, f = 1, right = T)
  
  # Define the weighting function
  # Option 1: W(p) = 1. Last jumps are not considered
  if(w_type == "none"){
    
    # Define a grid of evaluation points
    p <- unique(sort(c(0, p1, p2))) #, Inf)))
    grid <- (p[-1]+p[-length(p)])/2
    
    # Define the weighting function W(p) = 1
    swfun  <- stepfun(-1, c(1, 1), f = 1, right = T)
    
    w <- sapply(1:c(length(p)-1), function(x) integrate(swfun, 
                                                        upper = p[-1][x], 
                                                        lower = p[-length(p)][x])$value)
    
    step_dist <- sum(abs(sfun1(grid)-sfun2(grid))^l*w)
  }
  # Option 2: W(p) = 1 until a p < maxprice
  if(w_type == "max_price"){
    if(missing(maxprice)){maxprice <- 10000}
    
    # Define a grid of evaluation points
    p <- unique(sort(c(0, p1, p2, maxprice)))
    grid <- (p[-1]+p[-length(p)])/2
    
    swfun  <- stepfun(maxprice, c(1, 0), f = 1, right = T)
    
    w <- sapply(1:c(length(p)-1), function(x) integrate(swfun, 
                                                        upper = p[-1][x], 
                                                        lower = p[-length(p)][x])$value)
    
    step_dist <- sum(abs(sfun1(grid)-sfun2(grid))^l*w)
  }
  
  if(w_type == "user_wfun"){
    swfun <- user_fun
    # Define a grid of evaluation points
    p <- unique(sort(c(0, p1, p2, Inf))) # W(t) goes to 0
    grid <- (p[-1]+p[-length(p)])/2
    
    w <- sapply(1:c(length(p)-1), function(x) integrate(swfun, 
                                                        upper = p[-1][x], 
                                                        lower = p[-length(p)][x])$value)
    
    step_dist <- sum(abs(sfun1(grid)-sfun2(grid))^l*w)
  }
  
  if(plot == TRUE & w_type %in% c("none", "max_price", "user_wfun")){
    layout(matrix(c(1, 1, 1, 1, 2, 2), nrow = 3, ncol = 2, byrow = T))
    
    if(missing(maxprice)){maxprice = 0}
    plot(sfun1, xlim = c(0, max(p1, p2, maxprice)), ylim = c(0, max(q1, q2)),
         ylab = "Quantity", xlab = "Price", 
         main = paste0("l", l," Stepwise Distance = ", step_dist))
    plot(sfun2, col = "blue",  xlim = c(0, max(p1, p2, maxprice)), add = TRUE)
    points(c(0, p1), q1, pch = 16)
    points(c(0, p2), q2, pch = 16, col = "blue")
    
    plot(swfun, xlim = c(0, max(p1, p2, maxprice)),
         ylab = "Weight", xlab = "Price", 
         main = paste0("Weigthing function: ", w_type))
  }
  
  return(step_dist)
}

exponential_fun <- function(x) dexp(x, rate = 0.02)


up_offers_all = vroom("secondary_up_bids_offers_all.csv",col_types = list(datetime = col_character()))
up_offers_all= up_offers_all[,2:ncol(up_offers_all)]
up_offers_all_cumsum = t(apply(up_offers_all,1, cumsum))
Qmat = up_offers_all_cumsum[-(1:(2 * 365 * 24)),]

up_price_all = vroom("secondary_up_bids_prices_all.csv",col_types = list(datetime = col_character()))
up_price_all = up_price_all[,2:ncol(up_price_all)]
Pmat = up_price_all[-(1:(2 * 365 * 24)),]
Pmat = as.matrix(Pmat)


n_curves= nrow(Qmat)
inds = combn(seq_len(n_curves), 2)
inds_ff = as.ff(inds)
#ffsave(inds_ff,file="D:/Script/R_script/supply_curves_clustering/Supply curves clustering")
ffload(file="D:/Script/R_script/supply_curves_clustering/Supply curves clustering")


start = 945e6 + 1
end = 950e6
inds_sub= inds_ff[,start:end]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = end - start +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

dist = foreach(i=1:iterations,.combine='c', .options.snow=opts) %dopar% {
                index1 = inds_sub[1,i]
                p1 = Pmat[index1,2:ncol(Pmat)]
                p1 = na.omit(p1)
                p1 = t(p1)
                p1 = c(p1)
                
                
                q1 = Qmat[index1,2:ncol(Qmat)]
                q1 = na.omit(q1)
                q1 = t(q1)
                q1 = c(q1)
                
                min_len = min(length(p1),length(q1))
                p1 = p1[1:min_len]
                q1 = q1[1:min_len]
                p1= p1[-1]
                
                index2 = inds_sub[2,i]
                p2 = Pmat[index2,2:ncol(Pmat)]
                p2 = na.omit(p2) 
                p2 = t(p2)
                p2 = c(p2)
                
                
                q2 = Qmat[index2,2:ncol(Qmat)]
                q2 = na.omit(q2)
                q2 = t(q2)
                q2 = c(q2)
                
                min_len = min(length(p2),length(q2))
                p2 = p2[1:min_len]
                q2 = q2[1:min_len]
                p2= p2[-1]
                
                stepwise_distance_w_plot(p1 = p1, 
                                             q1 = q1,
                                             p2 = p2,
                                             q2 = q2,
                                             l = 2,
                                             w_type = "user_wfun",
                                             user_fun = exponential_fun,
                                             plot = F)
}
close(pb)
stopCluster(cl) 


save_file = paste("dist",floor(end/1e6),".RDS",sep = "")
saveRDS(dist,file = save_file)


path = "distance/1-310"
dist_list310 = list.files(path=path,pattern = ".RDS",full.names = T) %>% 
             map(readRDS) %>%
             unlist(use.names = F)
 

path = "distance"
dist_list_after_310 = list.files(path=path,pattern = ".RDS",full.names = T) %>% 
  map(readRDS) %>%
  unlist(use.names = F) 

dist_list_ff = as.ff(c(dist_list310,dist_list_after_310))
ffsave(dist_list_ff,file="D:/Script/R_script/supply_curves_clustering/Supply curves clustering",add = T)
ffload(file="D:/Script/R_script/supply_curves_clustering/Supply curves clustering")


n_curves = 43848L
mat_dist = matrix(0, ncol = n_curves, nrow = n_curves)
mat_dist = as.ff(mat_dist)

inds_ff_t=as.matrix(inds_ff[,])
inds_ff_t = t(inds_ff_t)




start = 9e8 + 1
end = dim(inds_ff_t)[1]
mat_dist_idx =arrayIndex2vectorIndex(inds_ff_t[start:end,], dim=c(43848,43848))
gc()
set.ff(mat_dist,mat_dist_idx,dist_list_ff[start:end])

ffsave(mat_dist,file="D:/Script/R_script/supply_curves_clustering/Supply curves clustering",add = T)

d = dist_list

attributes(d) <- list(Size = 43848,
                      Labels = as.character(1:43848),
                      Diag = FALSE,
                      Upper = FALSE,
                      method = "user")

class(d) = "dist"

result=hclust(d)
result$merge



# D = vt(mat_dist)
# dim(D)
# sum(D[2,])
# 
# D = as.matrix(D[,])
# d = as.dist(D)
# mat_dist_m = as.matrix(mat_dist[,])
# 
# mat_dist_m[lower.tri(mat_dist_m)] =NA



#mat_dist_t = vt(mat_dist)

#mat_dist_full = mat_dist + mat_dist_t

#mat_dist =as.ffdf(mat_dist)
#write.csv.ffdf(mat_dist, file = "mat_dist.csv")










#dat= matrix(rep(1:16),nrow=4,byrow=F)
#dat = as.ff(dat)

#comb_inds = inds_ff[,1:4]
#inx =arrayIndex2vectorIndex(t(comb_inds), dim=c(4,4))

#set.ff(dat,inx,0)





###########################################################################
# 
# Pmat[23805,]
# Qmat[23805,]
# 
# 
# start=1
# end = 2
# inds_sub= inds_ff[,start:end]
# a = sapply(1:2, function(i){
#                     index1 = inds_sub[1,i]
#                     p1 = Pmat[index1,2:ncol(Pmat)]
#                     p1 = na.omit(p1)
#                     p1 = t(p1)
#                     p1 = c(p1)
#                     p1= p1[-1]
#                     
#                     
#                     q1 = Qmat[index1,2:ncol(Qmat)]
#                     q1 = na.omit(q1)
#                     q1 = t(q1)
#                     q1 = c(q1)
#                     
#                     index2 = inds_sub[2,i]
#                     p2 = Pmat[index2,2:ncol(Pmat)]
#                     p2 = na.omit(p2)
#                     p2 = t(p2)
#                     p2 = c(p2)
#                     p2 = p2[-1]
#                     
#                     q2 = Qmat[index2,2:ncol(Qmat)]
#                     q2 = na.omit(q2)
#                     q2 = t(q2)
#                     q2 = c(q2)
#                     
#                     a = stepwise_distance_w_plot(p1 = p1, 
#                                                  q1 = q1,
#                                                  p2 = p2,
#                                                  q2 = q2,
#                                                  l = 2,
#                                                  w_type = "user_wfun",
#                                                  user_fun = exponential_fun,
#                                                  plot = F)})
# 
# 
# 
# 
# 
# 
# 
# 
# 
# b = sapply(1:2, function(i){
#       index1 = inds_sub[1,i]
#       p1 = Pmat[index1,2:ncol(Pmat)]
#       p1 = na.omit(p1)
#       p1 = t(p1)
#       p1 = c(p1)
#       
#       
#       q1 = Qmat[index1,2:ncol(Qmat)]
#       q1 = na.omit(q1)
#       q1 = t(q1)
#       q1 = c(q1)
#       
#       min_len = min(length(p1),length(q1))
#       p1 = p1[1:min_len]
#       q1 = q1[1:min_len]
#       p1= p1[-1]
#       
#       index2 = inds_sub[2,i]
#       p2 = Pmat[index2,2:ncol(Pmat)]
#       p2 = na.omit(p2)
#       p2 = t(p2)
#       p2 = c(p2)
#       
#       
#       q2 = Qmat[index2,2:ncol(Qmat)]
#       q2 = na.omit(q2)
#       q2 = t(q2)
#       q2 = c(q2)
#       
#       min_len = min(length(p2),length(q2))
#       p2 = p2[1:min_len]
#       q2 = q2[1:min_len]
#       p2= p2[-1]
#       
#       stepwise_distance_w_plot(p1 = p1, 
#                                q1 = q1,
#                                p2 = p2,
#                                q2 = q2,
#                                l = 2,
#                                w_type = "user_wfun",
#                                user_fun = exponential_fun,
#                                plot = F)})
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mat_matrix[lower.tri(mat_matrix)] =a

# 
# for(i in 1:n_curves){
#   print(i)
#   mat_dist[i, i:n_curves] <- sapply(i:n_curves, function(x) 
#     stepwise_distance(p1 = supplies[[curve_1]]$p, 
#                       q1= supplies[[curve_1]]$q, 
#                       p2 = supplies[[x]]$p, 
#                       q2= supplies[[x]]$q,
#                       l = 1))
# }

# 
# a = Pmat[889,2:ncol(Pmat)]
# a = na.omit(a)
# a = t(a)
# a = c(a)
# b = Qmat[889,2:ncol(Qmat)]
# b = na.omit(b)
# b = t(b)
# b = c(b)
# 
# min_len = min(length(a),length(b))
# a = a[1:min_len]
# b = b[1:min_len]
# a= a[-1]
# stepfun(a, b, f = 1, right = T) 



a<-data.frame(site.x=c("A","A","A","B","B","C"),   
                              site.y=c("B","C","D","C","D","D"),Distance=c(67,57,64,60,67,60))
dij2 <- with(a, Distance)
nams <- with(a, unique(c(as.character(site.x), as.character(site.y))))
attributes(dij2) <- with(a, list(Size = length(nams),
                         Labels = nams,
                         Diag = FALSE,
                         Upper = FALSE,
                         method = "user"))
class(dij2) <- "dist"
dif2