####        INSTITUTO TECNOLOGICO VALE - DESENVOLVIMENTO SUSTENTAVEL            ####
#### P0S-GRADUACAO EM USO SUSTENTAVEL DE RECURSOS NATURAIS EM REGIOES TROPICAIS ####

####          AULA PRATICA DE MODELOS DE DISTRIBUIÇÃO DE ESPECIES (SDM)         ####
####                   PADRONIZACAO DOS DADOS AMBIENTAIS                        ####

#### Scripts by Jeronymo Dalapicolla ####

####### OBJETIVOS: PREPARAR AS VARIÁVEIS AMBIENTAIS PARA INPUT DOS MODELOS
#######               A. Cortar as rasters usando uma mascara
#######               B. Reprojetar as variaveis
#######               C. Converter os rasters para o formato '.asc'



############ 1. INSTALACAO DOS PACOTES E CARREGAR AS FUNCOES AUXILIARES #############
# Carregar os pacotes necessarios:
library(rgdal)
library(raster)

#delimitar uma projecao espacial para lat/long e para UTM:
longlat_WGS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM_proj = CRS("+proj=utm +zone=22 +south +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

####################### 2. CARREGAR OS ARQUIVOS NECESSARIOS #########################
####### FAZER PARA O FUTURO E O PRESENTE
#carregar um raster para representar a resolução que você deseja
exemplo_res = raster("~/wc2.1_2.5m_MIROC/wc2.1_2.5m_bioc_MIROC6_ssp585_2081-2100.tif")
crs(exemplo_res) = longlat_WGS

#carregar todas os arquivos rasters de uma pasta, selecionando pelo formato. Recursive = TRUE incluirá todas as subpastas, As camadas devem estar na mesma resolução. Mude o p"pattern" para o formato dos rasters
camadas = list.files(path="~/wc2.1_2.5m_MIROC/", pattern =".tif", full.names=TRUE)
camadas = stack(camadas)
#definir a projecao:
crs(camadas) = longlat_WGS
#Verificar
camadas

#carregar a mascara para representar a area de estudo:
mascara = shapefile("~/Crowned/mascaras/South_America.shp")
#define a projection, the same one for all!!!!
crs(mascara) = longlat_WGS
plot(mascara)





################# 3. ALTERAR A RESOLUCAO DAS VARIAVEIS AMBIENTAIS ###################
#reduzir para a area de estudo a camada usada como exemplo para a resolucao. Isso reduzira o tempo das analises:
exemplo_res_red = crop(exemplo_res, extent(mascara), snap="out") #cria uma area retangular
plot(exemplo_res_red)

#Mudar a resolucao das variaveis ambientais. Essa etapa pode levar muito tempo dependendo da area amostrada.
camadas_res = resample(camadas, exemplo_res_red, method="bilinear", bylayer=TRUE, progress='text', snap="out")
camadas_res
plot(camadas_res[[1]])





############ 4. CORTAR AS VARIAVEIS AMBIENTAIS PARA A FORMA DA MASCARA ##############
#Clip the resampled layers with the study area:
camadas_res_mas = mask(camadas_res, mascara, bylayer=TRUE) #exactamente da forma da mascara
plot(camadas_res_mas[[1]])

plot(camadas_res)



######## 5. CONVERTER E SALVAR AS VARIAVEIS AMBIENTAIS PARA O FORMATO ASC ###########
#Salvar como .ascii em uma nova pasta
writeRaster(camadas_res, paste0("~/Crowned/EnvironmentalData/Future/ssp585/2080_2100/", paste0("bio",".asc")), driver='ascii', bylayer=TRUE)
