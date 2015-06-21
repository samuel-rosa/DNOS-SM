////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//               UNIVERSIDADE FEDERAL RURAL DO RIO DE JANEIRO                 //
//                       INSTITUTO DE AGRONOMIA                               //
//           CURSO DE PÓS-GRADUAÇÃO EM AGRONOMIA-CIÊNCIA DO SOLO              //
//                                                                            //
//            CONTRIBUIÇÃO À CONSTRUÇÃO DE MODELOS DE PREDIÇÃO                //
//                         DE PROPRIEDADES DO SOLO                            //
//                                                                            //
//                           PROJETO DE PESQUISA                              //
//                                                                            //
//                          ALESSANDRO SAMUEL-ROSA                            //
//                                                                            //
//                        Seropédica, março de 2012.                          //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//Descrição:                                                                  //
//Projeto de pesquisa apresentada ao Curso de Pós-Graduação                   //
//em Agronomia-Ciência do Solo, da Universidade Federal                       //
//Rural do Rio de Janeiro (UFRRJ), Rio de Janeiro, como                       //
//requisito parcial para a obtenção do grau de Doutor em                      //
//Agronomia-Ciência do Solo.                                                  //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//Comitê de orientação:                                                       //
//Dra. Lúcia Helena Cunha dos Anjos - Orientador                              //
//Dr. Gustavo de Matos Vasques (Embrapa) - co-orientador                      //
//Dr. Gerard Heuvelink (ISRIC) - co-orientador                                //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//                  Script de trabalho:                                       //
//              Análise dos dados geológicos                                  //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//Descrição:                                                                  //
//                                                                            //
//                                                                            //
//----------------------------------------------------------------------------//
//                                                                            //
//e-mail: alessandrosamuel@yahoo.com.br                                       //
//homepage: soil-scientist.net                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

  
###################################################
### Índice de objetos
###################################################  
# g = 'gstat object'
# ick.cval = 'indicator co-kriging cross-validation'
# ivm = 'indicator variogram model'.
# ivmf = 'indicator variogram model fit'.
# ik = 'indicator kriging'
# ick = 'indicator co-kriging'.
# iv = 'indicator variogram'.
# ikp = 'indicator kriging probability'
# ikpm = 'indicator kriging probability map'
# ickpm = 'indicator co-kriging probability map'
# ickvm = 'indicator co-kriging variance map'
# ickcm = 'indicator co-kriging co-variance map'
# ikcvm = 'indicator kriging cross-validation map'
# ikv = 'indicator kriging variance'
# ik.val = 'indicator kriging validation'
# ik.cval = 'indicator kriging cross-validation'
# icv = 'indicator cross-variogram'
# icvm = 'indicator cross-variogram model'
# pm = 'parent material'
# s = 'sedimentary'
# sed = 'sedimentary'
# sim = 'simulation'
# vul = "vulcânica'
# v = 'vulcânica'
# tmp = 'temporário'
# pts = 'points'


# Configuração inicial da área de trabalho ###########################
options(prompt="> ", continue="+ ", digits=5,
        width=70, show.signif.stars=T)
par(mfrow=c(1,1))
ls()
rm(list=ls())
getwd()

setwd("/home/alessandro/rdata/dnos_sm/geology")
load("/home/alessandro/rdata/dnos_sm/dnos_dados.RData")
load("/home/alessandro/rdata/dnos_sm/geology/pm_models.RData")
ls()

# Carregar os pacotes necessários ####################################
require(sp)
require(gstat)
require(rgdal)
require(maptools)
require(raster)
require(MASS)
require(psych)
require(fitdistrplus)
require(lattice)
require(xtable)

#####################################################
### Configuração inicial das variáveis
#####################################################

#----------------------------------------------------
# Material parental do solo
#----------------------------------------------------
file.show("dnos_rocha.txt")
dnos.rocha <- read.table(file = "dnos_rocha.txt", sep = "\t",
                         dec = ",", header = T, skip = 30)
head(dnos.rocha)
tail(dnos.rocha)
class(dnos.rocha)
str(dnos.rocha)
# Criar variável binária
dnos.rocha$sed <- dnos.rocha$material_parental == "Sedimentar"
dnos.rocha$vul <- dnos.rocha$material_parental == "Vulcânica"
head(dnos.rocha);tail(dnos.rocha)
# Criar objeto espacial
coordinates(dnos.rocha) <- c("xcoord", "ycoord")
str(dnos.rocha)
class(dnos.rocha)
# Georeferenciar com busca do código EPSG
EPSG <- make_EPSG()
(EPSG[grep("SIRGAS 2000 / UTM zone 22S", fixed=T, EPSG$note), ])
proj4string(dnos.rocha)
proj4string(dnos.rocha) <- CRS("+init=epsg:31982")
proj4string(dnos.rocha)
head(coordinates(dnos.rocha))
str(dnos.rocha)
proj4string(dnos.rocha)


# material parental----------------------------------
str(dnos.rocha)
su.1 <- (summary(dnos.rocha$Rocha))
su.2 <- round(summary(
  dnos.rocha$Rocha)/length(dnos.rocha$Rocha)*100, 2)
(desc.rocha <- data.frame(
  rbind(su.1, su.2), row.names=c("n", "%")))
rm(su.1, su.2)

coordinates(dnos.rocha) <- coordinates(dnos.coords)
str(dnos.rocha)

pdf(file="/home/alessandro/rdata/dnos_sm/figures/cal_rocha_pos.pdf",
    width=8, height=8)
print(spplot(dnos.rocha, zcol="material_parental",
             main="Location of soil parent material",
             key.space="bottom", aspect="iso",
             legendEntries=c("Sedimentary", "Volcanic"),
             xlab="E (m)", ylab="N (m)",
             xlim=bbox(dnos.raster)[1,],
             ylim=bbox(dnos.raster)[2,],
             scales = list(draw = TRUE),
             col.regions=c("blue", "red"),
             panel=function(x,y, ...) {
               panel.xyplot(x, y, ...);
               panel.grid(h=-1,v=-1, col="lightgray", lty=2)}))
dev.off()








ls()
str(dnos.rocha)

# Análise variográfica ###############################################
ls()

#-------------------------------------------------
# Variograma omnidirecional - visualização inicial
# para definição do 'cutoff'
#-------------------------------------------------
print(plot(variogram(sed ~ 1, loc=dnos.rocha),
           main="Indicator variogram",
           sub="Soil parent material",
           pch=20, col="blue",
           plot.numbers=T))

#-------------------------------------------------
# Variograma omnidirecional
#-------------------------------------------------
sed.iv <- variogram(sed ~ 1, loc=dnos.rocha,
                    cutoff=1000, width=100)
vul.iv <- variogram(vul ~ 1, loc=dnos.rocha,
                    cutoff=1000, width=100)

#-------------------------------------------------
# Variograma
#-------------------------------------------------
print(plot(sed.iv, main="Indicator variogram",
           sub="Soil parent material", pch=20,
           col="blue", plot.numbers=T,
           xlim=c(0, max(sed.iv$dist)*1.1),
           ylim=c(0, max(sed.iv$gamma)*1.1)))

#-------------------------------------------------
# Superfície variográfica
#-------------------------------------------------
print(plot(variogram(sed~1, dnos.rocha, map=TRUE,
                     cutoff=1000, width=100),
           main="Variogram map",
           sub="Soil parent material",
           col.regions=bpy.colors(64)))

#-------------------------------------------------
# Ajuste visual do modelo - sedimentar
#-------------------------------------------------
print(show.vgms())
(sed.ivm <- vgm(psill=0.12, model="Gau",
                range=800/sqrt(3), nugget=0))
print(plot(sed.iv, main="Indicator variogram",
           sub="Sedimentary rock", pch=20, col="blue",
           xlim=c(0, max(sed.iv$dist)*1.1),
           ylim=c(0, max(sed.iv$gamma)*1.1),
           plot.numbers=T, model=sed.ivm))

#-------------------------------------------------
# Ajuste automático do modelo - sedimentar
#-------------------------------------------------
(sed.ivmf <- fit.variogram(object=sed.iv,
                           model=sed.ivm,
                           warn.if.neg=TRUE))
print(plot(sed.iv, main="Indicator variogram",
           sub="Sedimentary rock", pch=20, col="blue",
           xlim=c(0, max(sed.iv$dist)*1.1),
           ylim=c(0, max(sed.iv$gamma)*1.1),
           plot.numbers=T, model=sed.ivmf))

#-------------------------------------------------
# Ajuste visual do modelo - vulcânica
#-------------------------------------------------
(vul.ivm <- vgm(psill=0.12, model="Gau",
                range=800/sqrt(3), nugget=0))
print(plot(vul.iv, main="Indicator variogram",
           sub="Volcanic rock", pch=20, col="blue",
           xlim=c(0, max(vul.iv$dist)*1.1),
           ylim=c(0, max(vul.iv$gamma)*1.1),
           plot.numbers=T, model=vul.ivm))
str(vul.ivmf)
#-------------------------------------------------
# Ajuste automático do modelo - vulcânica
#-------------------------------------------------
(vul.ivmf <- fit.variogram(object=vul.iv,
                           model=vul.ivm,
                           warn.if.neg=TRUE))
print(plot(vul.iv, main="Indicator variogram",
           sub="Volcanic rock", pch=20, col="blue",
           xlim=c(0, max(vul.iv$dist)*1.1),
           ylim=c(0, max(vul.iv$gamma)*1.1),
           plot.numbers=T, model=vul.ivmf))

###################################################
### Produtos gráficos da análise variográfica
###################################################
ls()
dev.off()

#-------------------------------------------------
# Superfície variográfica
#-------------------------------------------------
#pdf(file="rock_vario_map.pdf", width=8, height=8)
#print(rock.vario.map)
#dev.off()

#-------------------------------------------------
# Modelo de variograma ajustado
#-------------------------------------------------
pdf(file="rock_ivmf.pdf", width=8, height=8)
plot(sed.iv$gamma ~ sed.iv$dist, pch=20, col="blue",
     xlim=c(0, max(sed.iv$dist)*1.1),
     ylim=c(0, max(sed.iv$gamma)*1.1),
     xlab="Separation distance", ylab="Semivariance",
     main="Indicator variogram", cex=1.2,
     sub="Soil parent material")
lines(variogramLine(sed.ivm, maxdist=max(sed.iv$dist)), col="red")
lines(variogramLine(sed.ivmf, maxdist=max(sed.iv$dist)), col="blue")
text(sed.iv$dist, max(sed.iv$gamma), sed.iv$np, pos=3)
text(0, max(sed.iv$gamma)*1.1, "Number of point-pairs per bin", pos=4)
text(1000, 0.04, col="red", "eye fitted", pos=2)
text(1000, 0.03,
     "Gaussian, sill: 0.12, nugget: 0.00, range:461.88", pos=2)
text(1000, 0.02, col="blue", "gstat fitted", pos=2)
text(1000, 0.01,
     "Gaussian, sill: 0.11531, nugget: 0.00, range:368.67", pos=2)
dev.off()
ls()

###################################################
### Predição
###################################################
# Substituir a variável de interesse a cada vez um
# comando for utilizado. Para isso, basta usar o
# comando 'Ctrl + F' e solicitar a substituição. Os
# parâmetros a ser substituídos são:
# 'sed' e 'vul'
# 'sedimentary' e 'volcanic'
# 'lik' - local, e 'gik' - global. Nesse caso é
# preciso alterar o parâmetro 'maxdist'
#--------------------------------------------------
ls()

#-------------------------------------------------
# Krigagem indicatriz
#-------------------------------------------------
time <- proc.time()
# Volcanic - local indicator kriging
vul.lik <- krige(formula = vul ~ 1, locations = dnos.rocha,
                 newdata=dnos.raster.limite, indicators=TRUE,
                 model=vul.ivmf, maxdist=vul.ivmf$range[2]*sqrt(3))
# Volcanic - global indicator kriging
vul.gik <- krige(formula = vul ~ 1, locations = dnos.rocha,
                 newdata=dnos.raster.limite, indicators=TRUE,
                 model=vul.ivmf)
# Sedimentary - local indicator kriging
sed.lik <- krige(formula = sed ~ 1, locations = dnos.rocha,
                 newdata=dnos.raster.limite, indicators=TRUE,
                 model=vul.ivmf, maxdist=sed.ivmf$range[2]*sqrt(3))
# Sedimentary - global indicator kriging
sed.gik <- krige(formula = sed ~ 1, locations = dnos.rocha,
                 newdata=dnos.raster.limite, indicators=TRUE,
                 model=sed.ivmf)
proc.time() - time

#-------------------------------------------------
# Avaliação da krigagem indicatriz
#-------------------------------------------------
ik <- data.frame(vul.lik$var1.pred, vul.gik$var1.pred,
                 sed.lik$var1.pred, sed.gik$var1.pred)
str(ik)
(tmp <- summary(ik))
# Salvar tabela LaTeX
print(xtable(tmp, caption="Descriptive statistics of predicted probabilities. Volcanic and sedimentary parent material. Local and global indicator kriging."), file="pm_ik.tex", caption.placement="top", floating=TRUE, table.placement="ht", booktabs=TRUE)

# Normalização-----------------------------------
sed.lik$var1.pred <- pmin(1, sed.lik$var1.pred)
sed.lik$var1.pred <- pmax(0, sed.lik$var1.pred)
sed.gik$var1.pred <- pmin(1, sed.gik$var1.pred)
sed.gik$var1.pred <- pmax(0, sed.gik$var1.pred)
vul.lik$var1.pred <- pmin(1, vul.lik$var1.pred)
vul.lik$var1.pred <- pmax(0, vul.lik$var1.pred)
vul.gik$var1.pred <- pmin(1, vul.gik$var1.pred)
vul.gik$var1.pred <- pmax(0, vul.gik$var1.pred)

# Estatísticas descritivas
ik <- data.frame(vul.lik$var1.pred, vul.gik$var1.pred,
                 sed.lik$var1.pred, sed.gik$var1.pred)
str(ik)
(tmp <- summary(ik))

# Salvar tabela LaTeX----------------------------
print(xtable(tmp, caption="Renormalized"), file="pm_ik.tex", caption.placement="top", floating=TRUE, table.placement="ht", booktabs=TRUE, append=TRUE)
rm(tmp, ik)

#-------------------------------------------------
# Produtos gráficos
#-------------------------------------------------
# Mapa de probabilidade---------------------------
# Krigagens: vul.lik, vul.gik, sed.lik e sed.gik
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$vul, "green", "blue"))
dev.off()
pdf(file="vul_likp.pdf", width=8, height=8, pointsize=18)
print(spplot(vul.lik, zcol="var1.pred", at=seq(0,1, by=0.1),
             col.regions=rev(heat.colors(10)),
             main="Probability, volcanic parent material",
             xlab="E (m)", ylab="N (m)",
             sub="Local indicator kriging. Green: volcanic. Blue: sedimentary.",
             contour=T, sp.layout=list(pts),
             scales = list(draw = TRUE)))
dev.off()
rm(pts)

# Mapa de pureza----------------------------------
# O mapa de pureza apresenta o maior valor de pro-
# babilidade predito  pela krigagem em uma deter-
# minada célula. Para isso basta usar o valor pre-
# dito de um dos indicadores e aplicar uma função
# lógica. Assim, os valores > 0.50 são mantidos,
# e os valores < 0.5 são subtraídos de 1.00. O
# conceito é de Bierkens & Burrough (The
# indicator approach to categorical soil data.
# Journal of Soil Science, 44:361-368, 1993.
#-------------------------------------------------
ls()
summary(vul.lik$var1.pred)
tmp <- vul.lik
tmp$var1.pred[tmp$var1.pred<0.5] <- 1-tmp$var1.pred[tmp$var1.pred<0.5]
# Estatísticas descritivas da nova variável
summary(tmp$var1.pred)
plotdist(tmp$var1.pred)
# Parâmetros gráficos
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$vul, "red", "blue"))
lim <-list("sp.polygons", dnos.limite, col="black", lwd=3)
dev.off()
# Salvar gráfico em arquivo no formato *.pdf
pdf(file="ik_purity.pdf", width=8, height=8, pointsize=18)
print(spplot(tmp, zcol="var1.pred", at=seq(0.5,1, by=0.001),
             col.regions=gray(seq(0, 1, by=0.002)),
             main="Purity map",
             xlab="E (m)", ylab="N (m)",
             sub="Local indicator kriging\n Red: volcanic. Blue: sedimentary.",
             sp.layout=list(pts, lim), scales = list(draw = TRUE)))
dev.off()
# Remover objetos temporários
rm(pts, lim, tmp)
ls()

# Mapa de variância-------------------------------
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$vul, "red", "blue"))
dev.off()
pdf(file="pm_likv.pdf", width=8, height=8)
print(spplot(vul.lik, zcol="var1.var",
             col.regions=rev(heat.colors(64)),
             main="Variance, parent material",
             xlab="E (m)", ylab="N (m)",
             sub="Local indicator kriging. Red: volcanic. Blue: sedimentary.",
             contour=T, sp.layout=list(pts),
             scales = list(draw = TRUE)))
dev.off()
rm(pts)

###################################################
### Modelo linear de co-regionalização
###################################################
ls()

#-------------------------------------------------
# Construir objeto 'gstat'
#-------------------------------------------------
# Global
pm.gg <- gstat(NULL, id="sed", sed~1, loc = dnos.rocha)
pm.gg <- gstat(pm.gg, "vul", vul~1, loc = dnos.rocha)
summary(pm.gicv <- variogram(pm.gg, cutoff = 1000, width = 100))

# Local
pm.lg <- gstat(NULL, id = "sed", sed~1, loc = dnos.rocha,
               maxdist = 368.97*sqrt(3))
pm.lg <- gstat(pm.lg, id = "vul", vul~1, loc = dnos.rocha,
               maxdist = 368.97*sqrt(3))
summary(pm.licv <- variogram(pm.lg, cutoff = 1000, width = 100))

#-------------------------------------------------
# Ajuste visual do modelo de variograma
#-------------------------------------------------
# Global
(pm.gicvm <- gstat(g = pm.gg, model = vgm(psill = 0.12, model = "Gau",
                                          range = 368.67, nugget = 0),
                   fill.all = TRUE))
plot(pm.gicv, model=pm.gicvm)
# Local
(pm.licvm <- gstat(g = pm.lg, model = vgm(psill = 0.12, model = "Gau",
                                          range = 368.67, nugget = 0), fill.all = TRUE))
plot(pm.licv, model = pm.licvm)

#-------------------------------------------------
# Ajuste automático do modelo de variograma
# lmc - linear model of co-regionalization
#-------------------------------------------------
# Global
(pm.gicvmf <- fit.lmc(v=pm.gicv, g=pm.gicvm, warn.if.neg=TRUE,
                     correct.diagonal=1.01))
dev.off()
pdf(file="pm_gicvmf.pdf", width=8, height=8)
print(plot(pm.gicv, model=pm.gicvmf, xlab="Distance",
           ylab="Semivariance", plot.numbers=T, col="blue",
           pch=20, xlim=c(0, max(pm.gicv$dist)*1.1),
           main="Indicator cross-variogram, Soil parent material",
           sub="Linear model of co-regionalization"))
dev.off()

# Local
(pm.licvmf <- fit.lmc(v=pm.licv, g=pm.licvm, warn.if.neg=TRUE,
                      correct.diagonal=1.01))


###################################################
### Co-kriging
###################################################
ls()

#-------------------------------------------------
# Predição
#-------------------------------------------------
# Local
time <- proc.time()
pm.lick <- predict.gstat(object=pm.licvmf,
                         newdata=dnos.raster.limite,
                         indicators=TRUE, debug.level=-1)
summary(pm.lick)
proc.time() - time

# Global
time <- proc.time()
pm.gick <- predict.gstat(object=pm.gicvmf,
                         newdata=dnos.raster.limite,
                         indicators=TRUE, debug.level=-1)
summary(pm.gick)
(proc.time() - time)

#-------------------------------------------------
# Sumário das predições e erros
#-------------------------------------------------
ick <- data.frame(pm.lick$sed.pred, pm.lick$vul.pred,
                  pm.gick$sed.pred, pm.gick$vul.pred)
(tmp <- summary(ick))
 
# Salvar tabela LaTeX
print(xtable(tmp, caption="Descriptive staistics of predicted probabilities. Sedimentary and volcanic parent material. Local and global indicator co-kriging"), file="pm_ick.tex", caption.placement="top", floating=TRUE, table.placement="ht", booktabs=TRUE)

# Normalização-----------------------------------
pm.lick$sed.pred <- pmin(1, pm.lick$sed.pred)
pm.lick$sed.pred <- pmax(0, pm.lick$sed.pred)
pm.lick$vul.pred <- pmin(1, pm.lick$vul.pred)
pm.lick$vul.pred <- pmax(0, pm.lick$vul.pred)
pm.gick$sed.pred <- pmin(1, pm.gick$sed.pred)
pm.gick$sed.pred <- pmax(0, pm.gick$sed.pred)
pm.gick$vul.pred <- pmin(1, pm.gick$vul.pred)
pm.gick$vul.pred <- pmax(0, pm.gick$vul.pred)

# Salvar tabela LaTeX
ick <- data.frame(pm.lick$sed.pred, pm.lick$vul.pred,
                  pm.gick$sed.pred, pm.gick$vul.pred)
(tmp <- summary(ick))
print(xtable(tmp, caption="Renormalized"), file="pm_ick.tex",
      append=TRUE, caption.placement="top", floating=TRUE,
      table.placement="ht", booktabs=TRUE)
ls()
rm(tmp, ick)

#-------------------------------------------------
# Produtos gráficos
#-------------------------------------------------
# Mapas de probabilidade--------------------------
# Co-krigagens: pm.lick$sed.pred, pm.lick$vul.pred,
#               pm.gick$sed.pred, pm.gick$vul.pred
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$vul, "green", "blue"))
dev.off()
pdf(file="vul_lickp.pdf", width=8, height=8)
print(spplot(pm.lick, zcol="vul.pred", at=seq(0, 1, by=0.1),
         col.regions=rev(heat.colors(10)),
         main="Probability, Volcanic parent material",
         xlab="E (m)", ylab="N (m)",
         sub="Local indicator co-kriging. Green: volcanic. Blue: sedimentary.",
         contour=T, sp.layout=list(pts),
         scales = list(draw = TRUE)))
dev.off()
rm(pts)

# Mapa de pureza----------------------------------
# O mapa de pureza apresenta o maior valor de pro-
# babilidade predito  pela krigagem em uma deter-
# minada célula. Para isso basta usar o valor pre-
# dito de um dos indicadores e aplicar uma função
# lógica. Assim, os valores > 0.50 são mantidos,
# e os valores < 0.5 são subtraídos de 1.00. O
# conceito é de Bierkens & Burrough (The
# indicator approach to categorical soil data.
# Journal of Soil Science, 44:361-368, 1993.
#-------------------------------------------------
ls()
summary(pm.lick$sed.pred)
tmp <- pm.lick
tmp$sed.pred[tmp$sed.pred<0.5] <- 1-tmp$sed.pred[tmp$sed.pred<0.5]
# Estatísticas descritivas da nova variável
summary(tmp$sed.pred)
plotdist(tmp$sed.pred)
# Parâmetros gráficos
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$vul, "red", "blue"))
lim <-list("sp.polygons", dnos.limite, col="black", lwd=3)
dev.off()
# Salvar gráfico em arquivo no formato *.pdf
pdf(file="ick_purity.pdf", width=8, height=8, pointsize=18)
print(spplot(tmp, zcol="sed.pred", at=seq(0.5,1, by=0.001),
             col.regions=gray(seq(0, 1, by=0.002)),
             main="Purity map",
             xlab="E (m)", ylab="N (m)",
             sub="Local indicator co-kriging\n Red: volcanic. Blue: sedimentary.",
             sp.layout=list(pts, lim), scales = list(draw = TRUE)))
dev.off()
# Remover objetos temporários
rm(pts, lim, tmp)
ls()



# Mapa de variância-------------------------------
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$sed, "blue", "red"))
dev.off()
pdf(file="pm_lickvm.pdf", width=8, height=8)
print(spplot(pm.lick, zcol="sed.var", xlab="E (m)", ylab="N (m)",
         col.regions=rev(heat.colors(64)),
         main="Variance, Soil parent material",
         sub="Local indicator co-kriging. Red: volcanic. Blue: sedimentary",
         contour=T, sp.layout=list(pts),
         scales = list(draw = TRUE)))
dev.off()
rm(pts)

#-------------------------------------------------
# Mapa de covariância
#-------------------------------------------------
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$sed, "blue", "green"))
dev.off()
pdf(file="pm_gickcm.pdf", width=8, height=8)
print(spplot(pm.gick, zcol="cov.sed.vul",
         col.regions=rev(heat.colors(64)),
         main="Covariance, Soil parent material",
         xlab="E (m)", ylab="N (m)",
         sub="Global indicator co-kriging. Green: volcanic. Blue: sedimentary",
         contour=T, sp.layout=list(pts),
         scales = list(draw = TRUE)))
dev.off()
rm(pts)

###################################################
### Validação
################################################### 
ls()
str(dnos.validation)

#-------------------------------------------------
# Krigagem indicatriz ordinária. A validação das
# duas variáveis ('sed' e 'vul') é feita aqui.
# Basta mudar a variável nas linhas de comando.
#-------------------------------------------------
# Validação cruzada
#vik.cval <- krige.cv(vul ~ 1, nfold=5, nmax=50,
#                  loc=dnos.rocha, model=vul.ivmf)
#summary(vik.cval@data)

# Normalizar as probabilidades
#vik.cval$var1.pred <- pmin(1, vik.cval$var1.pred)
#vik.cval$var1.pred <- pmax(0, vik.cval$var1.pred)
#summary(vik.cval@data)

# Gráfico de posição espacial
#dev.off()
#pdf(file="vul_ikcvm.pdf", width=8, height=8)
#plot(coordinates(vik.cval), asp=1, cex=0.2+2*vik.cval$var1.pred,
#     col=ifelse(vik.cval$observed, "green", "red"),
#     xlab="E (km)", ylab="N (km)",
#     main="Probability of TRUE indicator, Volcanic parent material",
#     sub="Actual indicator: red/green = FALSE/TRUE")
#grid(col="gray")
#dev.off()

# Comparar indicadores observados com as probabilidades preditas
#dev.off()
#pdf(file="vul_ikcval.pdf", width=8, height=3)
#plot(vik.cval$observed ~ vik.cval$var1.pred, pch="|", cex=1.5,
#     xlab="Predicted probability of TRUE indicator",
#     xlim=c(0,1), yaxt="n", ylab="",
#     col=ifelse(vik.cval$observed, "green", "red"),
#     main="Cross-validation: IK, p(Volcanic parent material)",
#     sub="Green/Red: indicator TRUE/FALSE")
#abline(v=sum(vik.cval$observed)/length(vik.cval$observed),
#       lty=2, col="darkblue")
#dev.off()

# Validação externa

slik.val <- krige(formula=sed ~ 1, loc=dnos.rocha,
                  newdata=dnos.validation, indicators=TRUE,
                  model=sed.ivmf, maxdist=sed.ivmf$range[2]*sqrt(3))
summary(slik.val)

# Histograma de probabilidades
#print(histogram(vik.val$var1.pred, col="lightgray",
#                main="Histogram of probabilities",
#                xlab="Volcanic parent material"))

# Normalização das probabilidades
slik.val$var1.pred <- pmin(1, slik.val$var1.pred)
slik.val$var1.pred <- pmax(0, slik.val$var1.pred)
summary(slik.val@data)

#print(histogram(vik.val$var1.pred, col="lightgray",
#                main="Histogram of probabilities",
#                xlab="Volcanic parent material"))

# Comparar as observações (indicadores) e as predições (probabilidades)
# usando um boxplot e um gráfico de linhas
summary(tmp <- (dnos.validation$sed == TRUE))
#sum(tmp)/length(tmp)
#sum(dnos.rocha$sed)/length(dnos.rocha$sed)
tmp <- data.frame(id = as.numeric(row.names(slik.val@data)),
                  ised = tmp, pred = slik.val$var1.pred, 
                  tran = dnos.validation$transecto)
coordinates(tmp) <- coordinates(slik.val)
#str(compare)
#order(compare$pred)
tmp <- tmp[order(tmp$pred),]
#str(tmp)
#tmp[1, ]
#tmp[dim(tmp)[1], ]
#head(row.names(tmp@data))
#row.names(tmp@data) <- 1:dim(tmp)[1]
#head(row.names(tmp@data))
#tmp@data
#boxplot(tmp@data$pred~tmp@data$ised,
#      main="Validation: Indicator kriging, Sedimentary parent material",
#        col="lightgray", xlab="Indicator",
#        ylab="Predicted probability")

dev.off()
pdf(file="sed_likval.pdf", width=8, height=3)
plot(tmp$ised ~ tmp$pred, pch="|", cex=1.5,
     xlab="Predicted probability of TRUE indicator",
     xlim=c(0,1), yaxt="n", ylab="",
     col=ifelse(tmp$ised, "green", "red"),
     main="Validation: Local IK, p(Sedimentary parent material)",
     sub="Green/Red: indicator TRUE/FALSE")
abline(v=sum(tmp$ised)/length(tmp$ised), lty=2, col="darkblue")
text(tmp$pred, tmp$ised, tmp$id, pos=1)
text(tmp$pred, tmp$ised, tmp$id, pos=3)
dev.off()
ls()
rm(tmp)

#-------------------------------------------------
# Co-krigagem indicatriz
#-------------------------------------------------
# Apenas duas variáveis estão sob análise. Como a
# correlação entre elas é perfeita, a validação
# cruzada é inútil. Isso porque a krigagem é o
# BLUP, ou seja, o valor predito é igual ao valor
# observado em cada observação. A validação
# cruzada será interessante quando haver mais de
# duas variáveis sendo preditas.
#-------------------------------------------------

# Validação cruzada
#pm.ick.cval <- gstat.cv(object=pm.icvmf, nmax=50,
#                        verbose=F, nfold=5)
#str(pm.ick.cval)
#summary(pm.ick.cval@data)
#sqrt(mean(pm.ick.cval$residual^2))
#mean(pm.ick.cval$residual)
#mean(pm.ick.cval$residual^2/pm.ick.cval$ltpb.var)

# Gráfico de posição espacial
#dev.off()
#pdf(file="pm_ickcvm.pdf", width=8, height=8)
#plot(coordinates(pm.ick.cval), asp=1,
#     col=ifelse(pm.ick.cval$observed, "green", "red"),
#     cex=0.2+2*pm.ick.cval$sed.pred,
#     xlab="E (km)", ylab="N (km)",
#     main="Probability of TRUE indicator, Sedimentary parent material",
#     sub="Actual indicator: red/green = FALSE/TRUE")
#grid(col="gray")
#dev.off()

# Comparar indicadores observados com as probabilidades preditas
#plot(pm.ick.cval$observed ~ pm.ick.cval$sed.pred, pch="|", cex=1.5,
#     xlab="Predicted probability of TRUE indicator",
#     xlim=c(0,1), yaxt="n", ylab="",
#     col=ifelse(pm.ick.cval$observed, "green", "red"),
#     main="Cross-validation: IK, p(Sedimentary parent material)",
#     sub="Green/Red: indicator TRUE/FALSE")
#abline(v=sum(pm.ick.cval$observed)/length(pm.ick.cval$observed),
#       lty=2, col="darkblue")

# Validação externa-------------------------------
pm.gick.val <- predict.gstat(object=pm.gicvmf,
                             newdata=dnos.validation,
                             indicators=TRUE)
summary(pm.gick.val)
# Normalização das probabilidades
pm.gick.val$sed.pred <- pmin(1, pm.gick.val$sed.pred)
pm.gick.val$sed.pred <- pmax(0, pm.gick.val$sed.pred)
pm.gick.val$vul.pred <- pmin(1, pm.gick.val$vul.pred)
pm.gick.val$vul.pred <- pmax(0, pm.gick.val$vul.pred)
summary(pm.gick.val)

# Comparar as observações (indicadores) e as predições (probabilidades)
# usando um boxplot e um gráfico de linhas. A validação das
# duas variáveis ('sed' e 'vul') é feita aqui.
# Basta mudar a variável nas linhas de comando.
summary(tmp <- (dnos.validation$sed == TRUE))
#sum(tmp)/length(tmp)
#sum(dnos.rocha$sed)/length(dnos.rocha$sed)
tmp <- data.frame(id = as.numeric(row.names(pm.gick.val@data)),
                  ised = tmp, pred = pm.gick.val$sed.pred, 
                  tran = dnos.validation$transecto)
coordinates(tmp) <- coordinates(pm.gick.val)
#str(compare)
#order(compare$pred)
tmp <- tmp[order(tmp$pred),]
#str(tmp)
#tmp[1, ]
#tmp[dim(tmp)[1], ]
#head(row.names(tmp@data))
#row.names(tmp@data) <- 1:dim(tmp)[1]
#head(row.names(tmp@data))
#tmp@data
#boxplot(tmp@data$pred~tmp@data$ised,
#        main="Validation: Indicator co-kriging, Sedimentary parent #material",
#        col="lightgray", xlab="Indicator",
#        ylab="Predicted probability")

dev.off()
pdf(file="sed_gickval.pdf", width=8, height=3)
plot(tmp$ised ~ tmp$pred, pch="|", cex=1.5,
     xlab="Predicted probability of TRUE indicator",
     xlim=c(0,1), yaxt="n", ylab="",
     col=ifelse(tmp$ised, "green", "red"),
     main="Validation: Global ICo-K, p(Sedimentary parent material)",
     sub="Green/Red: indicator TRUE/FALSE")
abline(v=sum(tmp$ised)/length(tmp$ised), lty=2, col="darkblue")
text(tmp$pred, tmp$ised, tmp$id, pos=1)
text(tmp$pred, tmp$ised, tmp$id, pos=3)
dev.off()
ls()
rm(tmp)

###################################################
### Simulação geoestatística condicional
################################################### 
ls()

# Simulação com krigagem ordinária-----------------
sik.sim <- krige(formula=sed~1, loc=dnos.rocha,
                  dnos.raster.limite, indicators=TRUE,
                  model=vgm(psill=0.11531, model="Gau",
                            range=sed.ivmf$range[2]*sqrt(3),
                            nugget=0.00000001),
                  nsim=4, nmax=40, debug.level=-1)



# Estatísticas descritivas-------------------------
#summary(slik.sim)
#str(slik.sim)

# Normalização das probabilidades
#slik.sim$sim1 <- pmin(1, slik.sim$sim1)
#slik.sim$sim1 <- pmax(0, slik.sim$sim1)
#slik.sim$sim2 <- pmin(1, slik.sim$sim2)
#slik.sim$sim2 <- pmax(0, slik.sim$sim2)
#slik.sim$sim3 <- pmin(1, slik.sim$sim3)
#slik.sim$sim3 <- pmax(0, slik.sim$sim3)
#slik.sim$sim4 <- pmin(1, slik.sim$sim4)
#slik.sim$sim4 <- pmax(0, slik.sim$sim4)

# Estatísticas descritivas-------------------------
#summary(slik.sim)

# Mapas simulados----------------------------------
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$sed, "blue", "red"))
dev.off()
pdf(file="sik_sim.pdf", width=16, height=16)
print(spplot(sik.sim, zcol=1:4,
             main="Conditional simulations, local OK",
             sub="Sedimentary parent material: Red: volcanic. Blue: sedimentary",
             sp.layout=list(pts)))
dev.off()
rm(pts)


# Simulação com co-krigagem------------------------
ls()
lick.sim <- predict.gstat(object=pm.gicvmf, newdata=dnos.raster.limite,
                          nsim=4, indicators=TRUE,
                          debug.level=-1)

summary(lick.sim@data)

# Mapas simulados----------------------------------
pts <- list("sp.points", dnos.rocha, pch=20,
            col=ifelse(dnos.rocha$sed, "blue", "red"))
dev.off()
pdf(file="lick_sim.pdf", width=16, height=16)
print(spplot(lick.sim, zcol=1:4,
             main="Conditional simulations, local Co-K",
             sub="Sedimentary parent material: Red: volcanic. Blue: sedimentary",
             sp.layout=list(pts)))
dev.off()
rm(pts)

###################################################
### Produtos gráficos diversos
################################################### 
pdf(file="utopia.pdf", width=8, height=8)
par(mgp=c(2,0,0))
plot(seq(1, 10, by=1), seq(1, 10, by=1), type="n",
     asp=1, yaxt="n", xaxt="n", xlab="", ylab="")
par(cex=1.5)
title(ylab="Visually appealing", xlab="Does the job",
      main="Kriging vs. Co-Kriging\n Local vs. Global")
par(cex=1.2)
text(1.5, 10, "Global\n Indicator\n Kriging", pos=1)
text(5, 10, "Global\n Indicator\n Co-Kriging", pos=1)
text(5, 1, "Local\n Indicator\n Kriging", pos=3)
text(9.5, 1, "Local\n Indicator\n Co-Kriging", pos=3)
text(9.5, 9.5, "Utopia?", pos=1)
lines(c(1.5, 5, 5, 9.5), c(8.5, 8.5, 2.5, 2.5), lwd=2)
arrows(1.5, 8.5, 9.5, 2.5, lwd=2, col="red")
dev.off()





# Salvar os objetos--------------------------------------
ls()
rm()
save(pm.gg, pm.gicv, pm.gicvm, pm.gicvmf,
     pm.lg, pm.licv, pm.licvm, pm.licvmf,
     pm.gick, 
     pm.lick, 
     sed.iv, sed.ivm, sed.ivmf,
     vul.iv, vul.ivm, vul.ivmf,
     sed.gik,
     sed.lik,
     vul.gik,
     vul.lik,
     sik.sim, lick.sim,
     file="pm_models.RData")

rm(list=ls())





