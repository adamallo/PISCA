# PISCA
Phylogenetic Inference using Somatic Chromosomal Alterations. Plugin for BEAST 1.8.x.

## Citation
If you use this software, please cite:

* [Martinez P, Mallo D, Paulson TG, Li X, Sanchez CA, Reid BJ, Graham TA, Kuhner MK and Maley CC (2018) Evolution of Barrett's Esophagus through space and time at single-crypt and whole-biopsy levels. Nat. Commun 9: 794.](https://www.nature.com/articles/s41467-017-02621-x)
* [Kostadinov RL, Kuhner MK, Li X, et al. 2013. NSAIDs Modulate Clonal Evolution in Barrett’s Esophagus. PLoS Genet. 9.](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003553)
* [Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7 Molecular Biology And Evolution 29: 1969-1973](https://www.ncbi.nlm.nih.gov/pubmed/22367748)


If you use the random local model, please cite:

* [Drummond AJ & Suchard MA (2010) Bayesian random local clocks, or one rate to rule them all. BMC Biology 8, 114](https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-8-114)

If you use the ascertainment bias correction, please cite:

* [Alekseyenko AV, Lee C, Suchard MA (2008) Wagner and Dollo: a stochastic duet by composing two parsiminious solos. Systematic Biology. 57, 772-784](https://www.ncbi.nlm.nih.gov/pubmed/18853363)

## Known problems and limitations
**WARNING:** So far PISCA is only compatible with the java version of BEAST, and therefore it is necessary to include the argument *-beagle_off* (or deselect BEAGLE in the GUI). This makes the current version of PISCA incompatible with BEAST v1.10.X.

**WARNING:** PISCA is compatible with the Metropolis-coupled MCMC algorithm implemented in BEAST. However, there is a bug in BEAST v1.8.X that prevents plugins from using the MC3 algorithm. In order to circumvent this problem, you can download and compile [my modified version of BEAST v1.8.3 that incorporates PISCA](https://github.com/adamallo/beast-mcmc) or patch your BEAST v1.8.4 beast.jar (located in BEAST_ROOT/lib) with the class BeastMain.class distributed with PISCA (patch/BeastMain.class). In order to do so, you need to run something similar to this:

```
cd $BEAST_ROOT/lib
mkdir -p dr/app/beast
cp ~/Downloads/PISCAv1.0.X/patch/BeastMain.class dr/app/beast/
jar uf beast.jar dr/app/beast/BeastMain.class
rm -rf dr
```

## Installation
### Compiled version
In order to install the compiled version of PISCA, execute the install.sh included in the package, indicating the root directory of your BEAST 1.8.X installation as the only argument. 
Example:

```
./install.sh ~/bin/beast1.8.3/
```
### Sources
In order to install PISCA from its sources, execute ant test-install. This requires to have functional copies of Java JDK and Ant and having specified the root directory of your BEAST 1.8.X installation in the file beast_sdk.properties.
Example:

```
echo "beast.root=/Users/your/path/to/beast-mcmc/here" > beast_sdk.properties
ant test-install
```

## Usage
**WARNING:** So far PISCA is only compatible with the java version of BEAST, and therefore it is necessary to include the argument *-beagle_off* (or deselect BEAGLE in the GUI). This makes the current version of PISCA incompatible with BEAST v1.10.X.

PISCA does not include a modified version of Beauti. Therefore, you will have to put together the input xml manually (not recommended) use an script or manually modify one made using Beauti. The script [generate_genotypes.pl](https://github.com/adamallo/scripts_singlecrypt/blob/master/generate_genotypes.pl) could be a good starting point.

Below you can find the models implemented by PISCA and XML examples of their implementation.

## Included models
The models are presented in the order they should appear in the xml file. The different clock models are different alternatives (only one can be used at a time) and are ordered by complexity.

#### General data type for the CNV substitution model
The CNV substitution model requires the usage of a specific generalDatatype. This has to be **specified before the alignments** and look like:

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

#### Ascertainment correction bias:
Corrects for the absence of a specific invariable state (in principle the wild-type) in the sequences. **Important: There cannot be invariable positions with the sate that is being corrected for in the alignment.** If some are present in the data (but most are lost and therefore they need to be corrected for) they must be trimmed from the alignment.

The block structure is almost equivalent to the usual CharacterPatterns, with the addition of the parameter _state_ that indicates the missing invariable state that has to be corrected for (H is wildtype in the generalDataType used by our CNV substitution matrix).
 
```
<ascertainedCharacterPatterns id="patterns">
        <alignment idref="alignment"/>
        <state code='H'/>
</ascertainedCharacterPatterns>
```
#### Strict molecular clock:
The strict molecular clock has the same format and parameters as the original strict clock, with a name change (strictClockCenancestorBranchRates). If it is going to be estimated, an appropriate prior and operators must be included. Below you can see an example of a fixed rate of 1 (not estimated). 

```
<strictClockCenancestorBranchRates id="branchRates">
        <rate>
                <parameter id="clock.rate" value="1"/>
        </rate>
</strictClockCenancestorBranchRates>
```

#### Random local clock:
The random local clock and associated statistics have been adapted to the cenancestor-tree likelihood. However, their only change is their name (from randomLocalClockModel to randomLocalClockModelCenancestor), since their parameters, priors and operators are the same.

Example:

```
<randomLocalClockModelCenancestor id="branchRates" ratesAreMultipliers="false">
        <treeModel idref="treeModel"/>
        <rates>
                <parameter id="localClock.relativeRates"/>
        </rates>
        <rateIndicator>
                <parameter id="localClock.changes"/>
        </rateIndicator>
        <clockRate>
                <parameter id="clock.rate" value="1.0" lower="0.0"/>
        </clockRate>
</randomLocalClockModelCenancestor>

<sumStatistic id="rateChanges" name="rateChangeCount" elementwise="true">
        <parameter idref="localClock.changes"/>
</sumStatistic>

<rateStatisticCenancestor id="meanRate" name="meanRate" mode="mean" internal="true" external="true">
        <treeModel idref="treeModel"/>
        <randomLocalClockModelCenancestor idref="branchRates"/>
</rateStatisticCenancestor>
<rateStatisticCenancestor id="coefficientOfVariation" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">
        <treeModel idref="treeModel"/>
        <randomLocalClockModelCenancestor idref="branchRates"/>
</rateStatisticCenancestor>
<rateCovarianceStatistic id="covariance" name="covariance">
        <treeModel idref="treeModel"/>
        <randomLocalClockModelCenancestor idref="branchRates"/>
</rateCovarianceStatistic>
```
Also, you can get the relative rates, the indicator variables and the rate logged in the trees as traits.

Example:

```
<logTree id="treeFileLog" logEvery="2000" nexusFormat="true" fileName="rlc_tests.trees" sortTranslationTable="true">
                <treeModel idref="treeModel"/>
                <trait name="rate" tag="rate">
                        <randomLocalClockModelCenancestor idref="branchRates"/>
                </trait>
                <trait name="rates" tag="relRates">
                        <randomLocalClockModelCenancestor idref="branchRates"/>
                </trait>
                <trait name="rateIndicator" tag="indicator">
                        <randomLocalClockModelCenancestor idref="branchRates"/>
                </trait>
                <posterior idref="posterior"/>
<!--            <rateStatisticCenancestor idref="meanRate"/>-->
        </logTree>

```

#### Binary LOH substitution matrix:
TBI

#### CNV substitution matrix:
It requires alignments with CNV states recoded as alphanumerical characters using the general data type indicated at the beginning of this README. It also requires a frequency model that indicates the state (or partial likelihood) of the estate of the cenancestor of the tree. This has been designed to include in the model the last common ancestor with a normal (wildtype) genome, for which the frequency model should be:

```
<frequencyModel id="frequencies">
        <dataType idref="cnv"/>
        <frequencies>
                <parameter id="cnv.frequencies" value="0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"/>
        </frequencies>
</frequencyModel>
```

The CNV model itself must be specified after the frequencyModel.

Example of the model:

```
<CNVModel id="cnv_subsmodel">
        <frequencies>
                <frequencyModel idref="frequencies"/>
        </frequencies>
        <gain_rate>
                <parameter id="cnv.gain" value="1" lower="0"/>
        </gain_rate>
        <loss_rate>
                <parameter id="cnv.loss" value="1" lower="0"/>
        </loss_rate>
        <conversion_rate>
                <parameter id="cnv.conversion" value="1" lower="0"/>
        </conversion_rate>
</CNVModel>
```

The model depends on three parameters, namely: _Gain_, _Loss_ and _Conversion_ rate. However, only two parameters are free and the substitution matrix is normalized by the sum of the three rates. **We advice to fix one of the parameters and estimate the other two, which would be relative to the sum of the three**. The absolute rates will then be obtained dividing the estimated rate between the sum of the all of them (e.g., given estimated parameters _pLoss_, _pConversion_ and fixed _pGain_=1; _RelativeGainRate_= _pGain_/(_pLoss_+_pConversion_+_pGain_); _AbsoluteGainRate_=_RelativeGainRate_ * ClockRate).

You need to add operators and priors for the rates you want to estimate.

Example operators:

```
<scaleOperator scaleFactor="0.25" weight="0.25">
	<parameter idref="cnv.loss"/>
</scaleOperator>
<scaleOperator scaleFactor="0.25" weight="0.25">
	<parameter idref="cnv.conversion"/>
</scaleOperator>
```

Example priors:

```
<exponentialPrior mean="1.0" offset="0.0">
	<parameter idref="cnv.loss"/>
</exponentialPrior>
<exponentialPrior mean="1.0" offset="0.0">
	<parameter idref="cnv.conversion"/>
</exponentialPrior>
```


### Cenancestor Tree Likelihood

The cenancestor tree likelihood implements a likelihood calculation that takes into account the extra branch from the most recent common ancestor to the last common ancestor with normal genome. Moreover, it includes two parameters to determine the cenancestor height and the branch length of the extra branch. The model must operate on the second, and the prior can be set in either of those (although the height is usually easier). **Important**: This tree likelihood model can only use cenancestor-aware clock rate models (i.e., clock rates models implemented in PISCA).

Example model:

```
<cenancestorTreeLikelihood id="treeLikelihood" useAmbiguities="false">
        <patterns idref="patterns"/>
        <treeModel idref="treeModel"/>
        <siteModel idref="siteModel"/>
        <cenancestorHeight>
                <parameter id="luca_height" value="5.0" lower="5.0" upper="60.0"/> <!-- Without the value it does not add the bounds -->
        </cenancestorHeight>
        <cenancestorBranch>
                <parameter id="luca_branch" value="1" upper="55.0" lower="0.0"/>
                <!-- Value 1 as a safe starting value -->
        </cenancestorBranch>
        <strictClockCenancestorBranchRates idref="branchRates"/>
</cenancestorTreeLikelihood>
```

Example operators:

```
<scaleOperator scaleFactor="0.2" weight="1.0"> <!-- We operate the branch since it is relative to the root. Operating luca_height is error prone, since it depends on the root -->
	<parameter idref="luca_branch"/>
</scaleOperator>
```

Example priors:

```
<uniformPrior lower="5.0" upper="60.0">
	<parameter idref="luca_height"/>
</uniformPrior>
```

## Examples

In the examples folder you can find two examle xml files, one with an strict molecular clock and another one with a random local clock.
