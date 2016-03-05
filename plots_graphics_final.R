image_path <- "/Users/Simon/Studium/MSC/Masterarbeit/write-up/Final/images/"

#METHODS PLOTS

#Loglikelihoods of fitted fire distributions
pdf(paste0(image_path,"loglik_fd.pdf"))
df <- data.frame(unique(years.subs), "lognorm fitted" = ll_lnorm_fit, "lognorm unfitted" = ll_lnorm_unfit, "fitted pearson" = ll_pears, "fitted gamma" = ll_gamma)

plot(df$unique.years.subs., df$lognorm.fitted, type = "n", 
     xlab = "Years with sufficient data", ylab = "log-likelihood",
     ylim = c(min(df[,2:5]- 50), max(df[,2:5] + 50)),
     xaxt = "n",
     cex.axis = 0.8)
axis(side=1, at=round(df$unique.years.subs.), cex.axis = 0.8, las = 2)
points(df$unique.years.subs., df$lognorm.fitted, pch = 20, col = "red")
points(df$unique.years.subs., df$lognorm.unfitted, pch = 21, col = "green")
points(df$unique.years.subs., df$fitted.pearson, pch = 22, col = "blue")
points(df$unique.years.subs., df$fitted.gamma, pch = 23, col = "black")
legend(1978, 500, legend = c("log-normal fitted", "log-normal unfitted", "pearson fitted", "gamma fitted"), col = c("red", "green", "blue", "black"), pch = c(20,21,22,23), text.font = 1, cex = 0.8)
dev.off()

#Weibull distributions fire probabilities
pdf(paste0(image_path,"prob_weib.pdf"))
#plot fire probability curves
plot(Ftl[,6], Ftl[,5], type = "n",
     xlab = "time [y]",
     ylab = "F(t)",
     cex.axis = 0.8,
     xlim = c(0,35),
     ylim = c(0,1)
     )

for (i in 1:5){
  points(Ftl[,6], Ftl[,i], pch = 20, col = i)
  nd <- seq(0,35, by = 0.1)
  lines(nd, predict(mod_probs[[i]], newdata = data.frame(yr = nd)), col = i)  
}
legend(30, 0.2, legend = c("class 1", "class 2", "class 3", "class 4", "class 5"), col = c(1,2,3,4,5), pch = c(20, 20 , 20 , 20 ,20), lty = c(1,1,1,1,1), cex = 0.8)
dev.off()

#Model plots tab, fo, loc, scale
#tab mods
new_data<- data.frame("cumu3" = 0:2500)
preds_mod1 <- predict(tab_mod1, newdata = new_data, se.fit = T)[[1]]
sefits_mod1 <-predict(tab_mod1, newdata = new_data, se.fit = T)[[2]]
preds_mod2 <-predict(tab_mod2, newdata = new_data, se.fit = T)[[1]]
sefits_mod2 <-predict(tab_mod2, newdata = new_data, se.fit = T)[[2]]

df_mod1 <- data.frame("preds" = exp(preds_mod1), "upper" = exp(preds_mod1+ sefits_mod1), "lower" = exp(preds_mod1 - sefits_mod1))
df_mod2 <- data.frame("preds" = (preds_mod2)^2, "upper" = (preds_mod2+ sefits_mod2)^2, "lower" = (preds_mod2 - sefits_mod2)^2)
tab_mods <- list(df_mod1, df_mod2)
pdf(paste0(image_path,"tab_mods.pdf"), 14, 7)
par(mfrow = c(1,2))
plot(total ~ cumu3, pch = 20,
     xlab = "3ycp",
     ylab = "tab [%]",
     cex.axis = 0.8
     )
for (i in c(1,2)){
  coords.y <- c(tab_mods[[i]][,2], rev(tab_mods[[i]][,3]))
  coords.x <- c(0:2500, rev(0:2500))
matlines(tab_mods[[i]], lty = c(1,0,0), col = c(rep(i+2,3)))
polygon(x = coords.x, y = coords.y, col = alpha(col = i + 2, 0.2), border = F)
}

legend(600, 45, 
       legend = c("model 1 (+-SE)", "model 2 (+-SE)", "predictions", "observations"), 
       border = c(0,0,0,0),
       fill = c(alpha(col=3, 0.2), alpha(col=4, 0.2),NA,NA), 
       lty = c(NA,NA,1,NA), 
       pch = c(NA,NA,NA,19), 
       cex = 0.8)

#residual plot tab
plot(density(residuals(tab_mod1)), type = "l", 
     col = 4, ylim = c(0, 0.4), cex.axis = .8, main = "n",
     xlab = "Residuals",
     ylab = "Density", 
     )
lines(density(residuals(tab_mod2)), col = 3)
dev.off()

#loc mods
new_data<- data.frame("cumu3.subs" = 0:2500)
preds_mod1 <- predict(locs_mod1, newdata = new_data, se.fit = T)[[1]]
sefits_mod1 <-predict(locs_mod1, newdata = new_data, se.fit = T)[[2]]
preds_mod2 <-predict(locs_mod2, newdata = new_data, se.fit = T)[[1]]
sefits_mod2 <-predict(locs_mod2, newdata = new_data, se.fit = T)[[2]]

df_mod1 <- data.frame("preds" = preds_mod1, "upper" = preds_mod1+ sefits_mod1, "lower" = preds_mod1 - sefits_mod1)
df_mod2 <- data.frame("preds" = preds_mod2, "upper" = preds_mod2+ sefits_mod2, "lower" = preds_mod2 - sefits_mod2)
locs_mods <- list(df_mod1, df_mod2)

pdf(paste0(image_path,"locs_mods.pdf"), 14, 7)
par(mfrow = c(1,2))
plot(lnmeans ~ cumu3.subs, pch = 20,
     xlab = "3ycp",
     ylab = "loc",
     cex.axis = 0.8
)
for (i in c(1,2)){
  coords.y <- c(locs_mods[[i]][,2], rev(locs_mods[[i]][,3]))
  coords.x <- c(0:2500, rev(0:2500))
  matlines(locs_mods[[i]], lty = c(1,0,0), col = c(rep(i+2,3)))
  polygon(x = coords.x, y = coords.y, col = alpha(col = i+2,0.2), border = F)
}

legend(650, 3, 
       legend = c("model 1 (+-SE)", "model 2 (+-SE)", "predictions", "observations"), 
       border = c(0,0,0,0),
       fill = c(alpha(col=3, 0.2), alpha(col=4, 0.2),NA,NA), 
       lty = c(NA,NA,1,NA), 
       pch = c(NA,NA,NA,19), 
       cex = 0.8)

#residual plot locs
plot(density(residuals(locs_mod1)), type = "l", 
     col = 4, ylim = c(0, 0.4), cex.axis = .8,main = "",
     xlab = "Residuals",
     ylab = "Density"
)
lines(density(residuals(locs_mod2)), col = 3)
dev.off()

#scale mods
new_data<- data.frame("cumu3.subs" = 0:2500)
preds_mod1 <- predict(scales_mod1, newdata = new_data, se.fit = T)[[1]]
sefits_mod1 <-predict(scales_mod1, newdata = new_data, se.fit = T)[[2]]
preds_mod2 <-predict(scales_mod2, newdata = new_data, se.fit = T)[[1]]
sefits_mod2 <-predict(scales_mod2, newdata = new_data, se.fit = T)[[2]]

df_mod1 <- data.frame("preds" = preds_mod1, "upper" = preds_mod1+ sefits_mod1, "lower" = preds_mod1 - sefits_mod1)
df_mod2 <- data.frame("preds" = preds_mod2, "upper" = preds_mod2+ sefits_mod2, "lower" = preds_mod2 - sefits_mod2)
scales_mods <- list(df_mod1, df_mod2)

pdf(paste0(image_path,"scales_mods.pdf"), 14, 7)
par(mfrow = c(1,2))
plot(lnscales ~ cumu3.subs, pch = 20,
     xlab = "3ycp",
     ylab = "scale",
     cex.axis = 0.8
)
for (i in c(1,2)){
  coords.y <- c(scales_mods[[i]][,2], rev(scales_mods[[i]][,3]))
  coords.x <- c(0:2500, rev(0:2500))
  matlines(scales_mods[[i]], lty = c(1,0,0), col = c(rep(i+2,3)))
  polygon(x = coords.x, y = coords.y, col = alpha(col = i + 2, 0.2), border = F)
}

legend(650, 1.99, 
       legend = c("model 1 (+-SE)", "model 2 (+-SE)", "predictions", "observations"), 
       border = c(0,0,0,0),
       fill = c(alpha(col=3, 0.2), alpha(col=4, 0.2),NA,NA), 
       lty = c(NA,NA,1,NA), 
       pch = c(NA,NA,NA,19), 
       cex = 0.8)

#residual plot scales
plot(density(residuals(scales_mod1)), type = "l", 
     col = 4, ylim = c(0, 2), cex.axis = .8, main = "",
     xlab = "Residuals",
     ylab = "Density"
)
lines(density(residuals(scales_mod2)), col = 3)

#fo mod
new_data<- data.frame("pr" = 0:2500)
preds_mod1 <- predict(mod_occ, newdata = new_data, se.fit = T, type = "response")[[1]]
sefits_mod1 <-predict(mod_occ, newdata = new_data, se.fit = T, type = "response")[[2]]

df_mod <- data.frame("preds" = preds_mod1, "upper" = preds_mod1+ sefits_mod1, "lower" = preds_mod1 - sefits_mod1)
dev.off()
pdf(paste0(image_path,"fo_mods.pdf"), 14, 7)
par(mfrow = c(1,2))
plot(oc ~ pr, pch = 20,
     xlab = "3ycp",
     ylab = "p fo",
     cex.axis = 0.8,
     xlim = c(0, 2000)
)

  coords.y <- c(df_mod[,2], rev(df_mod[,3]))
  coords.x <- c(0:2500, rev(0:2500))
  matlines(df_mod, lty = c(1,0,0), col = 3)
  polygon(x = coords.x, y = coords.y, col = alpha(col = 3, 0.2), border = F)

legend(0, 0.95, 
       legend = c("SE range", "predictions", "observations"), 
       border = c(0,0,0),
       lty = c(NA,1,NA),
       col = c(NA,3,1),
       fill = c(alpha(col=3, 0.2), NA , NA), 
       pch = c(NA,NA,19),
       cex = 0.8)

dev.off()
#RESULTS
#Trajectory plots
pdf(paste0(image_path,"pop_traj.pdf"))
matplot(trajectories[[1]][,1], trajectories[[1]][,6], type = "n",
        xlab = "Year",
        ylab = "Metapopulation size",
        cex.axis = 0.8
        )
for (i in 1:3){
coords.y <- c(trajectories[[i]][,3], rev(trajectories[[i]][,5]))
coords.x <- c(trajectories[[i]][,1], rev(trajectories[[i]][,1]))
matlines(trajectories[[i]][,1], trajectories[[i]][,2:6], type = "l", lty = c(3,0,1,0,3), col = c(rep(i + 1, 5)))
polygon(x = coords.x, y = coords.y, col = alpha(col = i + 1, 0.2), border = F)
}
legend(2020, 12000, legend = c("1 SD range low", "1 SD range neutral", "1 SD range high", "mean",  "min/max"), fill = c(alpha(col=2, 0.2), alpha(col=3, 0.2), alpha(col=4, 0.2),0,0), lty = c(0,0,0,1,3), border = c(0,0,0,0,0), cex = 0.8)
dev.off()
