# PISCA
Plugin for BEAST 1.8.x in order to use somatic copy number alterations (SCA) to perform phylogenetic estimation.

##Installation
In order to install the compiled version of PISCA, execute the install.sh included in the package, indicating the root directory of your BEAST 1.8.X installation as the only argument. Example:

```
./install.sh ~/bin/beast1.8.3/
```

##Usage
**WARNING:** So far PISCA is only compatible with the java version of BEAST, and therefore it is necessary to include the argument *--beagle_off* (or deselect BEAGLE in the GUI).

PISCA does not include a modified version of Beauti. Therefore, you will have to put together the input xml manually (not recommended) use an script or manually modify one made using Beauti. The script [generate_genotypes.pl](https://github.com/adamallo/scripts_singlecrypt/blob/master/generate_genotypes.pl) could be a good starting point.

##Included models
The models are presented in the order they should appear in the xml file. The different clock models are different alternatives (only one can be used at a time) and are ordered by complexity.

####CNV substitution matrix:
It requires alignments with CNV states recoded as alphanumerical characters using a general data type like the following:

```
<generalDataType id="cnv">
        <state code="@"/> <!-- Genotype: 0,0 ; Beast State: 0 -->
        <state code="A"/> <!-- Genotype: 0,1 ; Beast State: 1 -->
        <state code="B"/> <!-- Genotype: 0,2 ; Beast State: 2 -->
        <state code="C"/> <!-- Genotype: 0,3 ; Beast State: 3 -->
        <state code="D"/> <!-- Genotype: 0,4 ; Beast State: 4 -->
        <state code="E"/> <!-- Genotype: 0,5 ; Beast State: 5 -->
        <state code="F"/> <!-- Genotype: 0,6 ; Beast State: 6 -->
        <state code="G"/> <!-- Genotype: 1,0 ; Beast State: 7 -->
        <state code="H"/> <!-- Genotype: 1,1 ; Beast State: 8 -->
        <state code="I"/> <!-- Genotype: 1,2 ; Beast State: 9 -->
        <state code="J"/> <!-- Genotype: 1,3 ; Beast State: 10 -->
        <state code="K"/> <!-- Genotype: 1,4 ; Beast State: 11 -->
        <state code="L"/> <!-- Genotype: 1,5 ; Beast State: 12 -->
        <state code="M"/> <!-- Genotype: 2,0 ; Beast State: 13 -->
        <state code="N"/> <!-- Genotype: 2,1 ; Beast State: 14 -->
        <state code="O"/> <!-- Genotype: 2,2 ; Beast State: 15 -->
        <state code="P"/> <!-- Genotype: 2,3 ; Beast State: 16 -->
        <state code="Q"/> <!-- Genotype: 2,4 ; Beast State: 17 -->
        <state code="R"/> <!-- Genotype: 3,0 ; Beast State: 18 -->
        <state code="S"/> <!-- Genotype: 3,1 ; Beast State: 19 -->
        <state code="T"/> <!-- Genotype: 3,2 ; Beast State: 20 -->
        <state code="U"/> <!-- Genotype: 3,3 ; Beast State: 21 -->
        <state code="V"/> <!-- Genotype: 4,0 ; Beast State: 22 -->
        <state code="W"/> <!-- Genotype: 4,1 ; Beast State: 23 -->
        <state code="X"/> <!-- Genotype: 4,2 ; Beast State: 24 -->
        <state code="Y"/> <!-- Genotype: 5,0 ; Beast State: 25 -->
        <state code="Z"/> <!-- Genotype: 5,1 ; Beast State: 26 -->
        <state code="["/> <!-- Genotype: 6,0 ; Beast State: 27 -->
        <ambiguity code="-" states="@ABCDEFGHIJKLMNOPQRSTUVWXYZ["/>
        <ambiguity code="?" states="@ABCDEFGHIJKLMNOPQRSTUVWXYZ["/>
</generalDataType>
```
**The data type must be specified before the alignments**

####Ascertainment correction bias:
Corrects for the absence of a specific invariable state (in principle the wild-type) in the sequences. **Important: There cannot be invariable positions with the sate that is being corrected for in the alignment.** If some are being obtained, but most are lost and therefore they need to be corrected, those instances have to be trimmed from the alignment. The block structure is almost equivalent to the usual CharacterPatterns, with the addition of the parameter _state_ that indicates the missing invariable state that has to be corrected for.
 
```
<ascertainedCharacterPatterns id="patterns">
        <alignment idref="alignment"/>
        <state code='H'/>
</ascertainedCharacterPatterns>
```
####Strict molecular clock:
####Random local clock:
####Cenancestor likelihood:
