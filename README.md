# ArcticPlantDynamics
This repository contains data and code for the article: García Criado et al. (2025). Plant diversity dynamics over space and time in a warming Arctic. Nature.

## Authors
Mariana García Criado, Isla H. Myers-Smith, Anne D. Bjorkman, Sarah C. Elmendorf, Signe Normand, Peter Aastrup, Rien Aerts, Juha M. Alatalo, Lander Baeten, Robert G. Björk, Mats P. Björkman, Noémie Boulanger-Lapointe, Ethan E. Butler, Elisabeth J. Cooper, J. Hans C. Cornelissen, Gergana N. Daskalova, Belen Fadrique, Bruce C. Forbes, Greg H. R. Henry, Robert D. Hollister, Toke Thomas Høye, Ida Bomholt Dyrholm Jacobsen, Annika K. Jägerbrand, Ingibjörg S. Jónsdóttir, Elina Kaarlejärvi, Olga Khitun, Kari Klanderud, Tiina H. M. Kolari, Simone I. Lang, Nicolas Lecomte, Jonathan Lenoir, Petr Macek, Julie Messier, Anders Michelsen, Ulf Molau, Robert Muscarella, Marie-Louise Nielsen, Matteo Petit Bon, Eric Post, Katrine Raundrup, Riikka Rinnan, Christian Rixen, Ingvild Ryde, Josep M. Serra-Diaz, Gabriela Schaepman-Strub, Niels M. Schmidt, Franziska Schrodt, Sofie Sjögersten, Manuel J. Steinbauer, Lærke Stewart, Beate Strandberg, Anne Tolvanen, Craig E. Tweedie and Mark Vellend. 

Contact: Mariana García Criado, mariana.garcia.criado@gmail.com

## Data use guidelines
Data output from this manuscript are publicly available using a Creative Commons Attribution 4.0 International copyright (CC BY 4.0; see license.txt). Data are fully public but should be appropriately referenced by citing the paper. Although not mandatory, we additionally suggest that data users contact and collaborate with data contributors if this dataset will contribute a substantial proportion of observations used in a particular paper or analysis. DOI for this dataset is XXX and can also be found at XXXZENODOLINKXXX.

## Data availability & access
The data and code for this manuscript will be mantained at this GitHub repository (https://github.com/marianagarciacriado/ArcticPlantDynamics). 

## Citation
García Criado, M., Myers-Smith, I.H., Bjorkman, A.D., Elmendorf, S.C., Normand, S., Aastrup, P., Aerts, R., Alatalo, J.M., Baeten, L., Björk, R.G., Björkman, M.P., Boulanger-Lapointe, N., Butler, E.E., Cooper, E.J., Cornelissen, J.H.C., Daskalova, G.N., Fadrique, B., Forbes, B.C., Henry, G.H.R., Hollister, R.D., Høye, T.T., Jacobsen, I.B.D., Jägerbrand, A.K., Jónsdóttir, I.S., Kaarlejärvi, E., Khitun, O., Klanderud, K., Kolari, T.H.M., Lang, S.I., Lecomte, N., Lenoir, J., Macek, P., Messier, J., Michelsen, A., Molau, U., Muscarella, R., Nielsen, M.-L., Petit Bon, M., Post, E., Raundrup, K., Rinnan, R., Rixen, C., Ryde, I., Serra-Diaz, J.M., Schaepman-Strub, G., Schmidt, N.M., Schrodt, F., Sjögersten, S., Steinbauer, M.J., Stewart, L., Strandberg, B., Tolvanen, A., Tweedie, C.E. and Vellend, M. 

## Data
Climate data from CHELSA can be accessed at https://chelsa-climate.org/ and snow data is available at https://springernature.figshare.com/collections/ARCLIM_bioclimatic_indices_for_the_terrestrial_Arctic/6216368 

Plant composition data is included here in the 'data' folder. Script 00_data_prep shows the process of generating this input file, but the original raw files including the whole of the ITEX+ database are not provided in this repository. The full ITEX+ database will be available as part of an upcoming data paper, which will be available here: https://github.com/annebj/ITEX30_VegComp

When a raw data file is not included in this repository, this is specified in the script as "this file is not available in the repo as it contains raw data". All summarised data files are provided in this repository as input files and enable the reproducibility of the figures and analyses of this manuscript.

## Scripts
All the analyses undertaken for this manuscript are split between multiple R scripts within the `scripts`folder.
They can be followed in a sequential order (i.e., 1 to 8), with both 0_ scripts showing the process of generating summarised input files.

## Figures
The figures generated in R are stored in the `figures` folder.

## Model outputs
Full model outputs for statistical analyses are stored in the `models` folder.

## Software requirements
R version 4.2.0. or greater.

R packages: `betapart, brms, broom, corrplot, cowplot, ggdist, ggOceanMaps, ggeffects, ggnewscale, ggpubr, ggrepel, ggtern, ggtext, modelr, paletteer, randomcoloR, Taxonstand, tidyverse, vegan, viridis`
