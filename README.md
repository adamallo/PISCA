# PISCA
Plugin for BEAST 1.8.x in order to use somatic copy number alterations (SCA) to perform phylogenetic estimation.

##Installation
In order to install the compiled version of PISCA, execute the install.sh included in the package, indicating the root directory of your BEAST 1.8.X installation as the only argument (i.e., ./install.sh ~/bin/beast1.8.3/).

##Usage
**WARNING:** So far PISCA is only compatible with the java version of BEAST, and therefore it is necessary to include the argument *--beagle_off* (or deselect BEAGLE in the GUI).

PISCA does not include a modified version of Beauti. Therefore, you will have to put together the input xml manually (not recommended) use an script or manually modify one made using Beauti. The script [generate_genotypes.pl](https://github.com/adamallo/scripts_singlecrypt/blob/master/generate_genotypes.pl) could be a good starting point.

##Included models
* Cenancestor likelihood:
* CNV substitution matrix:
* Ascertainment correction bias:
* Strict molecular clock:
* Random local clock:

