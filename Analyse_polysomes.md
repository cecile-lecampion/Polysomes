# Analyse des données de polysomes

Ce document décrit l’utilisation du script `script_poly_final.R`

Le script `script_poly_final.R` permet d’analyser les données issue de l’analyse de polysome prfiling telle que décrite dans :

> An Easy Method for Plant Polysome Profiling. Lecampion, C., Floris, M., Fantino, J.R., Robaglia, C., Laloi, C. J. Vis. Exp. (114), e54231,doi:10.3791/54231 (2016).



:point_right: *Le script `polysome_profiling.Rmd`permet de générer les résultats sous la forme d’un rapport au format `.html`. Les paramètres à modifier sont les mêmes que ceux décrit ci dessous.* L’exécution du script est globale (pas de pas à pas) et se lance avec le bouton `Knit` présent dans le cadre en haut à gauche.



Les données issues de la spectroscopie sont contenues dans des fichiers `.csv` qui ont la dorme suivante :

0.0001041666692,0.1534467638,

0.000312499993,0.1535461843,

0.0005208333605,0.1536204666,

Le script permet le chargement de ces données et leur analyse quelque soit le nombre d’échantillons et de réplicats

## Le script pas à pas

### Choix du répertoire de travail

Il s’agit du répertoire où se trouve les données et où les résultats seront enregistrés.

```R
# choisir le répertoire de travail
setwd("PATH_TO_DIRECTORY")
```



### Le chargement des données

- Fonction pour le chargement

On utilise une fonction pour ne conserver que les données utiles du fichier `.csv` : On ne conserve que la seconde colonne et les éléments numériques. Les deux premières lignes et les lignes à la fin sont éliminées. 

```R
#-------------------------------------------------------------------------------------
# Function to load the data
# Load the csv file and return a one cilumon data frame
# Column is named with parameter colName
# Usage exemple :
#    df <- f_load_polysome_csv("a1.csv", "a1")
#-------------------------------------------------------------------------------------
# The file contains 3 columns:
# 0.0001041666692,0.1534467638,
# 0.000312499993,0.1535461843,
# 0.0005208333605,0.1536204666,
# Only the 2nd column is kept
# Only numerics are kept
# Two first lines and useless lines at the end of the file are ommitted
#-------------------------------------------------------------------------------------
f_load_polysome_csv <- function(csvFile, colName) {
  lines <- readLines(csvFile)
  linesNbToImport <- grep("^\\d+\\.\\d+,\\d+\\.\\d+,", lines, perl = TRUE)
  df <- read.table(text = lines[linesNbToImport], header = FALSE, sep = ",", dec = ".")
  df <- as.data.frame(df$V2)
  colnames(df) <- c(colName)
  return(df)
}
#-------------------------------------------------------------------------------------

```

Cette fonction prend en paramètres le nom du fichier en `.csv ` et le nom que l’on souhaite donner à la colonne qui constitura l’objet créé.

- Chargement des données

On donne aux objets contenant les données un nom générique de format alphanumérique du type a1, b1... avec la lettre représentant une lignée ou une condition et le chiffre représentant le numéro de réplicat.

```R
a1 <- f_load_polysome_csv("line1-1.csv", "a1")
a2 <- f_load_polysome_csv("line1-2.csv", "a2")
b1 <- f_load_polysome_csv("line2-1.csv", "b1")
b2 <- f_load_polysome_csv("lien2_2.csv", "b2")
c1 <- f_load_polysome_csv("line3-1.csv", "c1")
c2 <- f_load_polysome_csv("line3-2.csv", "c2")

# Keep the generic names of the data in a variable
NAMES <- c("a1", "a2", "b1", "b2", "c1", "c2")
```

S’il y a plus de lignées/réplicats is suffira de rajouter des lignes construites de la même façon et de compléter la variable `NAMES`

Par exemple :

```R
a3 <- f_load_polysome_csv("line&-3.csv", "a3")

NAMES <- c("a1", "a2", "a3", "b1", "b2", "c1", "c2")
```

Les données sont ensuite associées en un seul data frame en vue du traitement.

:warning: Ne pas oublié de mettre les nom de fichier et de colonnes entre guillemet

### Personnalisation des variables

Ces variables permettent de définir des zones pour le tracé des courbes et pour identifier les minimum et les maximum sur les courbes. La valeur de ces varaibles est arbitraire et dépend uniquement des données. Il est recommandé de faire un premier plot avec les valeur par défaut pui d’ajuster au fur et à mesure.

- Les limites sur les axes

Le début de la courbe n’a pas vraiment d’intérêt et peut être omit, à la fin la montée due à la chlorophylle induit des DO élevées qui écrasent le plot si elles sont conservées.

Sur l’axe des ordonnées le YMAX est calculé en fonction des données.

```R
# Limits for graphs axis
# On X axis
XMIN <- 1000		# Une valeur numérique
XMAX <- 3600		# Une valeur numérique

# On Y axis
YMIN <- 0			# Une valeur numérique

# Cette valeur est calculée
YMAX <- max(df[XMIN:XMAX, 2: ncol(df)], na.rm = TRUE) + (max(df[XMIN:XMAX, 2: ncol(df)], na.rm = TRUE))*0.1
```

- Zone de recherche du monosome et des minimum

Pour une série de données, le monosome est toujours situé dans la même zone qui peut être définie par un vecteur.

Les minimums avant, après le monosome et après le pic dans le surnageant permettrons de séparer les fractions.

```R
# Area where to search for the monosome peak
ZMONO <- c(2200 : 3000)		# Un vecteur numérique définissant un intervalle

# Area for the valley
# Area before the monosome
Z_BMONO <- c(1500:2000)			# Un vecteur numérique définissant un intervalle
# Area after the monosome
Z_PMONO <- c(2100:2600)			# Un vecteur numérique définissant un intervalle
# Area after the peak in supernatent
Z_PSUR <- c(2800:3200)			# Un vecteur numérique définissant un intervalle
```

- Bornes pour le calcul de la moyenne

La moyenne sera représentée sur le graph final, les valeurs non significative au début du jeu de données ne sont pas utilisées.

```R
# Limits to compute replicats mean
MEAN_START <- 500			# Une valeur numérique
MEAN_END <- 5000			# Une valeur numérique
```

- Noms réels des donénes

Le nom réel de chaque échantillon est conservé dans une variable de façon à être utilisé dans la représentation finale des données.

```R
# Real names of the data
# Keep here the real names of the data
DATA_REAL_NAME <- c("Line1", "Line2", "Line3")			# Un vecteur de chaine de charactère
```

:warning: Attention à ne pas oublié les guillemets autour de chaque nom de lignée

- Nombre de réplicat

```R
# Number of replicats
NB_REPLICAT <- 2			# Une valeur numérique
```

## Exécution du script

Le reste du script ne doit pas être modifié (à partir de la ligne 87).

Le script peut être excuté ligne par ligne en utilisant le bouton `Run` ou d’un bloc en utilisant le bouton `Source`. Ces deux boutons sont dans le cadre supérieur gauche.

Les graphs générés apparaissent dans le cadre en bas à droite :

- Les données brutes
- Les données alignées sur le monosome
- Les données avec le minimum à 0 (sur l’axe Y)
- La moyenne
- La moyenne smoothée par la méthode de l’adjacent averaging
- La moyenne smoothée par la méthode de l’adjacent averaging avec les bornes définie au début : `Z_BMONO`, `Z_PMONO`, `Z_PSUR` de façon à vérifier leur pertinence
- La proportion relative des fractions polysome, monosome et surnageant sous la forme d’un Barplot représentant, en pourcentage, les aires sous la courbe de chaque fraction.

 Chacun de ces graphiques peut être enregistré en utilisant le bouton `Export` du volet en bas à droite puis la boite de dialogue qui permet de choisir le répertoire de sauvegarde et le format.