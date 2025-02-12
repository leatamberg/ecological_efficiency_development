# The ecological efficiency of human development in 20th century socialist and capitalist economies

Dylan Sullivan (1,2\*†), Lea A. Tamberg (3\*†), Jason Hickel (1,4,5)

(1) Institute of Environmental Science and Technology (ICTA-UAB), Autonomous University of Barcelona, Barcelona, Spain

(2) School of Social Sciences, Macquarie University, Sydney, Australia

(3) Faculty of Geosciences and Environment, University of Lausanne, Chavannes-près-Renens, Switzerland

(4) ICREA, Barcelona, Spain

(5) International Inequalities Institute, London School of Economics and Political Science, London, United Kingdom.

(\*) Corresponding authors. E-mails: [lea.tamberg\@unil.ch](mailto:lea.tamberg@unil.ch), [dylan.sullivan\@mq.edu.au](mailto:dylan.sullivan@mq.edu.au)

(†) These authors contributed equally.

## Abstract

A key challenge for the 21st century is to achieve good social outcomes for all with sustainable levels of resource use. Existing research indicates this can be realised by organising production around human wellbeing, including through socialist policies such as public finance, universal services and job guarantees. We investigate this question empirically by comparing the performance of socialist and capitalist economies in the late 20th century across 19 human development indicators, as a function of aggregate production and resource use. We find that socialist economies performed better than capitalist economies at any given level of GDP, material footprint and, less significantly, CO2 emissions per capita, for 18 of the 19 indicators covering poverty, life expectancy, child mortality, healthcare, education, essential services, income equality and gender equality, except for physical integrity rights. These results indicate socialist policy may help countries to achieve good lives for all within ecological limits, particularly when paired with robust democratic participation.

## Contents of this repository

*   [code](code): data preparation, analysis notebooks, utility functions
*   [data](data): raw data from public sources, our main dataset with time series for each country, the aggregated datasets for each analysis
*   [figures](figures): empty, expected by scripts to store generated plots and tables
*   [renv](renv): contains the configuration of the exact environment (R and package versions) needed to reproduce the results

## Reproducing the results

When opening the project for the first time, `renv` will bootstrap itself and then inform you that one or more packages recorded in the lockfile are not installed. In order to install the required packages in the local environment, execute `renv::restore()` in the console (see the [renv documentation](https://rstudio.github.io/renv/articles/renv.html) for more details).

To reproduce the analysis, you can either re-build the main dataset ([data_all.csv](data/data_all.csv)) from the raw data by executing the [data preparation notebook](code/prepare_datasets.Rmd) or use the uploaded version. The analysis for [GDP](code/GDP_efficiency.Rmd), [material footprint](code/MF_efficiency.Rmd), and [CO<sub>2</sub> emissions](code/CO2_efficiency.Rmd) all have their own notebook. At the end of each notebook, the corresponding plots and tables of the paper and supplementary information are generated.
