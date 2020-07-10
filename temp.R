

flou <- prepareFlou(df_orig[df_orig$react.ID == 6, -c(1:11)])

df_em <- data.frame(Fluorescence=sort(flou))

set.seed(1234)
emres1 <- emcluster(df_em, init.EM(df_em, nclass = 2), assign.class = TRUE)
# em.bic(df_em, emres1)
em.icl(df_em, emres1)

emres2 <- emcluster(df_em, init.EM(df_em, nclass = 3), assign.class = TRUE)
# em.bic(df_em, emres2)
em.icl(df_em, emres2)

which.max(rank(emres2$Mu))


x <- emclassifier_teigen(flou, 0.85, crit="ICL")


df_em[max(which(x$classification == "neg")),]


x <- teigen(x = flou, Gs=2:3, scale = FALSE, convstyle = "lop") 

plot(x=flou,  y =teigen_dist(flou,df = x$parameters$df[1], mu = x$parameters$mean[1], sigma = c(x$parameters$sigma)[1]))
