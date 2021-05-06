library(vegan)
library(janitor)

datos <- read.csv("hieracium1819.csv", header=T, sep = ";")

row.names(datos) <- datos[,1]
datos <- datos[,-1]

datos <- clean_names(datos)
colnames(datos) <- abbreviate(colnames(datos),minlength=6)

# escalaBB <- c("0", "r", "+", "1", "2", "3", "4", "5")
# recodifico según 

datos[datos == "r" ] <- 0.01
datos[datos == "+" ] <- 0.1
datos[datos == 1 ] <- 2.5
datos[datos == 1 ] <- 2.5
datos[datos == 2 ] <- 17.5
datos[datos == 3 ] <- 37.5
datos[datos == 4 ] <- 62.5
datos[datos == 5 ] <- 87.5
datos[] <- lapply(datos, as.numeric)

datos.t <- t(datos)

# 1 - riqueza -------------------------------------------------------------

riqueza <- apply(datos>0,MARGIN = 1,sum) # "1" es por fila
riqueza

barplot(riqueza)
barplot(riqueza,cex.names = .7, col = 'darkgreen',las = 2)

datos.z <- read.csv("zonas_tiempos.csv", header=T, sep = ";")
datos.z$zona <- as.factor(datos.z$zona)

specpool(x = datos,pool = datos.z$zona)



# 2 - diversidad ----------------------------------------------------------

diversity(datos,index="shannon", MARGIN = 1) #spp x sitio
diversity(datos.t,index="shannon", MARGIN = 2) # sito x spp

datos.agg <- aggregate(datos, by = list(datos.z$zona), FUN='sum')
row.names(datos.agg) <- datos.agg[,1]
datos.agg <- datos.agg[,-1]

diversity(datos.agg,MARGIN = 1, index = "shannon")
diversity(datos.agg,MARGIN = 1, index = "simpson")

# 3 - curvas RA -----------------------------------------------------------

curva.ra <- radfit(datos.t)
plot(curva.ra)

curva.amb <- radfit(datos.agg)
plot(curva.amb)

summary(curva.amb)

# 4 - PCA -----------------------------------------------------------------

pca <- rda(datos.t)

pca
eigenvals(pca)

decostand(datos.agg,method="max",digits=2)

pca <- rda(datos.t, scale = T)

screeplot(pca)

plot(pca,display = "sites",type="text",scaling=1)
text(pca,display="species",labels=colnames(datos.t),col="red",cex=0.6)
scores(pca,choices=1:2,"sites", scaling=1)
scores(pca,choices=1:2,"species",scaling=1)

plot(pca,display = "sites",type="text",scaling=2)
text(pca,display="species",labels=colnames(datos.t),col="red",cex=0.6)
scores(pca,choices=1:2,"sites", scaling=2)
scores(pca,choices=1:2,"species",scaling=2)
