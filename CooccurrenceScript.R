##A script for analyzing cooccurrence in conversational data
##install required packages
packages <- c("ggplot2", "dplyr", "tidytext", "stringr", "devtools", "magrittr", "MASS", "cluster", "Matrix", "plotly")
install.packages(setdiff(packages, rownames(installed.packages())))
library(devtools)
devtools::install_github("TimoSchuer/ExmaraldaR", upgrade = "never")
devtools::install_github("TimoSchuer/Cooccurrence", ref = "master", upgrade = "never")
library(dplyr)
library(ggplot2)
library(tidytext)
library(stringr)
library(Cooccurrence)
library(ExmaraldaR)
library(magrittr)
library(tidyr)
library(MASS)
library(cluster)
library(Matrix)
library(umap)
library(plotly)

# Read in transcription ---------------------------------------------------
##single file
path <- rstudioapi::selectFile()
##uses see vignette/ help file of the ExmaraldaR package for description of the parameters
#path has to be an partitur editor transcription .exb
Corpus <- ExmaraldaR::read_exb_file(path, readAnn = TRUE, annotation = "linear", addMetaData = FALSE)

##multiple files
##files have to be stored in a single folder
path <- rstudioapi::selectDirectory()
Corpus <- ExmaraldaR::read_exb_dir(path, addMetaData = FALSE, readAnn = TRUE, annotation = "linear")


# get insights in the data ------------------------------------------------

##for this to work the Variable name has to be a column and the variant and all dependent factors each have a column

Variable <- "COLUMN NAME" # insert the column name in the quotation marks
Variant <- c("COLUMN NAME VARIANT", "COLUMN NAME FACTOR EG STRESS") # insert the the name of the column where the variant is defined. You can also add factors like stress if you want to check if stressed and unstress variants differ

counts <- Cooccurrence::countVars(Corpus, Variable,Variante)

##plot frequencies for overview
counts %>% ggplot(aes(x=Variable, y=n, fill= Variante)) +geom_col()+ facet_grid(cols =  vars(Variable))
# as the latter option tends to get quite cpnfusing with many variables it can be useful to get one plot per Variable
for (k in unique(counts$Variable )) {
  g <- counts %>% dplyr::filter(Variable==k) %>% ggplot(aes(x=Variable, y=n, fill= Variante))+ geom_col(position = "dodge")
  plot(g)
}

# Recoding ----------------------------------------------------------------
#after inspection it can be useful to recode the data
#while annotation errors should corrected directly in the transcript,
#it can be useful to coerce variants only virtually and leave the original data untouched as this often leads to information loss

##some useful snippets for recoding for further useful functions see: https://statisticsglobe.com/replace-values-in-data-frame-conditionally-in-r
# which() evaluates conditions to select which values should be replaced
# the first condition Corpus$V=="CH" selects a variable so replace "V" with the column that contains the Variablename
# the second condition specifies which values in the column where the variant is defined should be replaced; the values are separated with | which means LOGICAL OR.
# lastly you select the column where the replacement should occur. In this case "Variant"
#if you are not familiar with these please look for "subsetting data.frames r". There are a lot of tutorials online

Corpus[which(Corpus$V=="CH"& str_detect(Corpus$Variante,"çx|ç")),"Variante"] <- "Fri"

#second it can be useful to exculde some very rare and unusual variants
#you either just filter them out
Corpus <- Corpus %>% filter(!(V=="CH" & Variante=="n"))
# or juse the snippe above to replace the value with NA
Corpus[which(Corpus$V=="CH"& str_detect(Corpus$Variante,"çx|ç")),c("Variable","Variante")] <- NA



# calculate the cooccurrence matrix ---------------------------------------
CoocMat <- calcCooc(Corpus,format= "matrix")


# visualize the similiarities in cooccurrence patterns --------------------

# first step ist to calculate a dissimilarity matrix
distMat <- daisy(as.matrix(CoocMat),metric = "manhattan")
#second plot it via non metric multidimensional scaling
#it is essential to curate the data carefully as this method is sensibel to objects that don't differ at all. This is mostly the case when there are very few samples of a variant
nMDS <- MASS::isoMDS(distMat, k=2)
nMDS[["points"]] %>% as.data.frame() %>% ggplot(aes(x=V1, y=V2, label= rownames(.))) +geom_jitter() + geom_text() + labs(subtitle = str_glue("Stress=", nMDS[["stress"]]/100))

##a very similiar method is umap LINK EINFÜGEN. As a nmds already solves the case it can be jused to see how good the coocurrences are represented or if there are too many variants that have few occurrences.
# in the latter case you should rather curate your annotations than just simlply use umap
#anyway especially if you have lots of variants this method is more robust.
umap_conf <- umap.defaults
umap_conf$n_neighbors <- 5# you can play arround with this parameter. low values focus on local structeres, high value on global structures
umap_conf$input <- "dist"
umap.cooc <- distMat %>% as.matrix() %>% umap(config= umap_conf, input= "dist")

data <- umap.cooc$layout %>% as.data.frame() %>% mutate(Var= rownames(.)) %>%
  separate(Var, into = c("Variable", "Variante"), sep="_")
g <- ggplot(data, aes(x=V1, y= V2, label= rownames(data)), shape= factor(Variable))+
  geom_point(aes(colour= Variable)) + geom_text()
g
