library(biomod2);
library(raster);
library(RColorBrewer);
library(dismo);

setwd("~/Crowned/Localities/")
#Get data points
crowned <- read.csv(file = "griseoaur_corrigido.csv", header = T);
crowned <- cbind(crowned, rep.int(1, length(nrow(crowned)))); #Adds another column indicating these are presence points
colnames(crowned) <- c("X", "Y", "Response");

#Get environmental variables
setwd("~/Crowned/Ms/Crowned/");
envtList <- list.files(pattern = ".asc");
envt.st <- stack(envtList);
crs(envt.st) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#Get projection variables
setwd("~/Crowned/EnvironmentalData/Future/ssp585/2080_2100/")
projectionList <- list.files(pattern = ".asc");
proj.st <- stack(projectionList);
crs(proj.st) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

names(proj.st)=names(envt.st)

#Setting up data file for Biomod2
bmData <- BIOMOD_FormatingData(resp.var = crowned[,3],
                               resp.xy = crowned[,1:2], 
                               resp.name = "GriseoAur",
                               expl.var = envt.st,
                               PA.nb.rep=1
);

#Download Maxent na pasta correta
utils::download.file(url = "https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar", 
                     destfile = paste0(system.file("java", package = "dismo"), 
                                       "/maxent.jar"), mode = "wb")

#Setting up Maxent run
myBiomodOption <- Print_Default_ModelingOptions();
myBiomodOption@MAXENT.Phillips$path_to_maxent.jar = paste(system.file(package="dismo"), "/java", sep='');
myBiomodOption@MAXENT.Phillips$memory_allocated = 2048; #Allocates 2048 MB/2 GB of memory to modeling
myBiomodOption@MAXENT.Phillips$maximumiterations = 10000;
myBiomodOption@MAXENT.Phillips$threshold = F;
myBiomodOption@MAXENT.Phillips$hinge = F;
myBiomodOption@MAXENT.Phillips$visible = F;
myBiomodOption@MAXENT.Phillips$beta_lqp = .95;

#Running Maxent
setwd("~/Crowned/TestModelRun/R/")
myMaxentModel <- BIOMOD_Modeling(data=bmData,
                                    models=c('MAXENT.Phillips'),
                                    models.options=myBiomodOption,
                                    NbRunEval=10,
                                    do.full.models = F,
                                    DataSplit=50,
                                    models.eval.meth = c('KAPPA','TSS','ROC'),
                                    SaveObj = T
);

#Ensemble of all models--combines model runs using a user-selected evaluation metric
myMaxentEnsemble <- BIOMOD_EnsembleModeling( modeling.output = myMaxentModel,
                                   chosen.models = 'all',
                                   em.by = 'all',
                                   eval.metric = c('TSS'),
                                   eval.metric.quality.threshold = NULL,
                                   models.eval.meth = c('TSS','ROC','KAPPA'),
                                   prob.median = TRUE )

#Projecting your model to the present
myBiomodProjPres <- BIOMOD_Projection(modeling.output = myMaxentModel,
                                    new.env = envt.st,
                                    proj.name = 'Present',
                                    selected.models = 'all',
                                    compress = 'gzip',
                                    clamping.mask = T,
                                    output.format = '.grd',
                                    do.stack=T
);

mod_projPres <- get_predictions(myBiomodProjPres);
presentResult <- calc(mod_projPres,fun = median); #Choose whatever descriptive statistic you'd like
plot(presentResult, main = "Presente - Griseotyrannus Aurantioatrocristatus");
writeRaster(presentResult, filename = "griseoaurPresent", format = "GTiff", overwrite = T);
#Projecting the ensemble model in the present
myBiomodProjPresEnsemble <- BIOMOD_EnsembleForecasting(myMaxentEnsemble,
                            projection.output = myBiomodProjPres,
                            selected.models = 'all',
                            compress = 'gzip'
);
mod_projPresEnsemble <- get_predictions(myBiomodProjPresEnsemble);
presentEnsembleResult <- mod_projPresEnsemble[[2]] #This is the median model ensemble
plot(presentEnsembleResult, main = "Presente - Griseotyrannus Aurantioatrocristatus");
writeRaster(presentEnsembleResult, filename = "griseoaurPresent", format = "GTiff", overwrite = T);

#Projecting your model to the future
myBiomodProj2100 <- BIOMOD_Projection(modeling.output = myMaxentModel,
                                    new.env = proj.st,
                                    proj.name = 'In2100',
                                    selected.models = 'all',
                                    compress = 'gzip',
                                    clamping.mask = T,
                                    output.format = '.grd',
                                    do.stack=T
);

mod_proj2100 <- get_predictions(myBiomodProj2100);
result2100 <- calc(mod_proj2100,fun = median); #Choose whatever descriptive statistic you'd like
plot(result2100, main = "Em 2100 - Griseotyrannus Aurantioatrocristatus");

writeRaster(result2100, filename = "griseoaur2100", format = "GTiff", overwrite = T);

#Projecting the ensemble model in 2100
myBiomodProj2100Ensemble <- BIOMOD_EnsembleForecasting(myMaxentEnsemble,
                                                       projection.output = myBiomodProj2100,
                                                       selected.models = 'all',
                                                       compress = 'gzip'
);
mod_proj2100Ensemble <- get_predictions(myBiomodProj2100Ensemble);
ensembleResult2100 <- mod_proj2100Ensemble[[2]] #This is the median model ensemble
plot(ensembleResult2100, main = "Em 2100 - Griseotyrannus Aurantioatrocristatus");
writeRaster(ensembleResult2100, filename = "griseoaur2100Ensemble", format = "GTiff", overwrite = T);


#Evaluating models
## Variable response curves
response.plot2(models = BIOMOD_LoadModels(myMaxentModel, models='MAXENT.Phillips'),
               Data = get_formal_data(myMaxentModel,'expl.var'),
               show.variables= get_formal_data(myMaxentModel,'expl.var.names'),
               do.bivariate = FALSE,
               fixed.var.metric = 'median',
               col = brewer.pal(10, "Spectral"),
               legend = TRUE,
               data_species = get_formal_data(myMaxentModel,'resp.var')
);

##Varible contributions; only for Maxent, not possible for other models
forSetup <- read.csv(file = paste("~/Crowned/TestModelRun/R/Crowned/models/1604549077/GriseoAur_PA1_RUN1_MAXENT.Phillips_outputs/maxentResults.csv", sep = ""), header = T)#Choose the appropriate model folder with the seed of the analysis you want
variableContributions <- matrix(data = NA, nrow = length(forSetup[, grep('.contribution', names(forSetup))]), ncol = 10);
rownames(variableContributions) <- names(forSetup[, grep('.contribution', names(forSetup))])
colnames(variableContributions) <- c("Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8", "Run9", "Run10")
variablePermutationImportance <- matrix(data = NA, nrow = length(forSetup[, grep('.permutation.importance', names(forSetup))]), ncol = 10);
colnames(variablePermutationImportance) <- c("Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8", "Run9", "Run10")
count <- 1;
while (count <= 10){
  temporary <- read.csv(file = paste("~/Crowned/TestModelRun/R/Crowned/models/1604549077/GriseoAur_PA1_RUN", count, "_MAXENT.Phillips_outputs/maxentResults.csv", sep = ""), header = T);
  variableContributions[,count] <- unlist(unname(temporary[, grep('.contribution', names(temporary))]))
  variablePermutationImportance[,count] <- unlist(unname(temporary[, grep('.permutation.importance', names(temporary))]))
  count <- count + 1;
}
write.csv(variableContributions, "VariableContributions.csv", quote = F);
write.csv(variablePermutationImportance, "VariablePermutationImportance.csv", quote = F);

##Calculate MESS for 2100
mess2100 <- mess(proj.st, rasterToPoints(envt.st)[,-1:-2],-9999);
writeRaster(mess2100, filename = "GriseoAur2100MESS", format = "GTiff", overwrite = T);

##Create dataset for evaluation
###"Cutoff" gives threshold to optimize evaluation metric, "Sensitivity" and "Specificity" are based on this threshold
myMaxentModelEval <- get_evaluations(myMaxentModel, as.data.frame = F); 
write.csv(myMaxentModelEval["TSS",,,,],"TSSEvaluation.csv", quote = F);#Results for True Skill Score
write.csv(myMaxentModelEval["ROC",,,,],"TSSEvaluation.csv", quote = F);#Results for Area Under the Curve
write.csv(myMaxentModelEval["Kappa",,,,],"TSSEvaluation.csv", quote = F) #Results for Cohen's Kappa
