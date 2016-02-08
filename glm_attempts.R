rf <- randomForest(fs ~ (cls19 + cls25 + cls16 + cls33 + cls18 + cls14 + cls11 + tsf), data = raw.data)

varImpPlot(rf)

glm <- glm(fs ~ (tsf + I(tsf^2) + I(tsf^3))^2, data = raw.data, family = binomial("logit"))


require("MASS")
glm.step <- stepAIC(glm)
glm2.step <- stepAIC(glm2)
plot(raw.data$tsf, raw.data$fs)

preds<- predict(glm, newdata = raw.data.new, type = "response")
plot(raw.data.new$tsf, preds)
preds.glm2 <- predict(glm2, newdata = raw.data.new, type = "response")
plot(raw.data$tsf, raw.data$fs)