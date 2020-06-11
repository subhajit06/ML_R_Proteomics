library(dirichlet)
library(gridExtra)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2); #theme_set(theme_bw())

m = matrix(c(1,1,1,3,3,3,12,12,12,3,12,3,3,3,12,0.8,0.8,0.8),nrow=3)

p = NULL

for(i in 1:6){
  f <- function(v) ddirichlet(v, m[,i])
  
  mesh <- simplex_mesh(.0025) %>% as.data.frame %>% tbl_df
  mesh$f <- mesh %>% apply(1, function(v) f(bary2simp(v)))
  
  str = sprintf("Dirichlet(%0.1f, %0.1f, %0.1f)",m[1,i],m[2,i],m[3,i])
  p[[i]] <- ggplot(mesh, aes(x, y)) +
      geom_raster(aes(fill = f)) +
      coord_equal(xlim = c(0,1), ylim = c(0, 0.85)) + 
      ggtitle(str)
  
}

pdf("~/Desktop/plot_density.pdf",width=10,height=8)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],ncol=3)
dev.off()
