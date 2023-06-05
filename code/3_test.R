df=iris
set.seed(12334)
df$random1=runif(nrow(df),min=min(df$Sepal.Length),max=max(df$Sepal.Length))
df$mediator=df$Sepal.Length*0.35+df$random1*0.65
df$random2=runif(nrow(df),min=min(df$mediator),max=max(df$mediator))
df$dv=df$mediator*0.35+df$random2*0.65
fit.totaleffect=lm(dv~Sepal.Length,df)
summary(fit.totaleffect)
fit.mediator=lm(mediator~Sepal.Length,df)
summary(fit.mediator)
fit.dv=lm(dv~Sepal.Length+mediator,df)
summary(fit.dv)
library(mediation)
results = mediate(fit.mediator, fit.dv, treat='Sepal.Length', mediator='mediator', boot=T)
summary(results)
