
- [HigginsChenLab/methylCIPHER](#higginschenlabmethylcipher)
  - [Installation](#installation)
  - [Calculating Epigenetic Clocks and
    Predictors](#calculating-epigenetic-clocks-and-predictors)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# HigginsChenLab/methylCIPHER

<!-- badges: start -->
<!-- badges: end -->

The goal of methylCIPHER is to allow users to easily calculate their
choice of CpG clocks using simple commands, from a single source. CpG
epigenetic clocks are currently found in many places, and some require
users to send data to external portals, which is not advisable when
working with protected or restricted data. The current package allows
you to calculate reported epigenetic clocks, or where not precisely
disclosed, our best estimates–all performed locally on your own machine
in R Studio. We would like to acknowledge the authors of the original
clocks and their valuable contributions to aging research, and the
requisite citations for their clocks can be found at the bottom of the
current page. Please do not forget to cite them in your work!

## Installation

You can install the released version of methylCIPHER and its imported
packages from [Github](https://Cgithub.com/HigginsChenLab/methylCIPHER)
with:

``` r
devtools::install_github("danbelsky/DunedinPoAm38")
devtools::install_github("danbelsky/DunedinPACE")
devtools::install_github("HigginsChenLab/methylCIPHER")
```

## Calculating Epigenetic Clocks and Predictors

### Running single “clock” calculations

The current package contains a large number of currently available CpG
clocks or CpG based predictors. While we strove to be inclusive of such
published CpG-based epigenetic clocks to our knowledge, if you find we
are missing a clock, please contact us and we will do our best to
promptly include it, if possible. You can do so by raising an issue on
this repo or emailing us directly at
\<a\[dot\]higginschen\[at\]yale\[dot\]edu\>.

In order to calculate a CpG clock, you simply need to use the
appropriate function, typically named “calc\[ClockNameHere\]”. For
example:

``` r
library(methylCIPHER)
calcPhenoAge(exampleBetas, examplePheno, imputation = F)
```

<div class="kable-table">

| name              | geo_accession | gender |  age | group | sample | PhenoAge |
|:------------------|:--------------|:-------|-----:|------:|-------:|---------:|
| 7786915023_R02C02 | GSM1343050    | M      | 57.9 |     1 |      1 | 52.29315 |
| 7786915135_R04C02 | GSM1343051    | M      | 42.0 |     1 |      2 | 41.05867 |
| 7471147149_R06C01 | GSM1343052    | M      | 47.4 |     1 |      3 | 43.54460 |
| 7786915035_R05C01 | GSM1343053    | M      | 49.3 |     1 |      4 | 43.96697 |
| 7786923035_R01C01 | GSM1343054    | M      | 52.5 |     1 |      5 | 40.35242 |

</div>

Alternatively, if you would just like to receive a vector with the clock
values to use, rather than appending it to an existing phenotype/
demographic dataframe, simply use:

``` r
calcPhenoAge(exampleBetas, imputation = F)
#> 7786915023_R02C02 7786915135_R04C02 7471147149_R06C01 7786915035_R05C01 
#>          52.29315          41.05867          43.54460          43.96697 
#> 7786923035_R01C01 
#>          40.35242
```

### Categories of Epigenetic Clocks

Due to the abundance of epigenetic clocks are overlapping reasons for
use, it is important to keep track of which clocks are most related to
each other. This will allow you to steer clear of multiple testing and
collinearity problems if you use all available clocks for your analysis.

Another important note is that with a few exceptions, the list of
following clocks were trained and validated almost exclusively in blood.
This can lead to a number of observed effects, which include
unreasonable shifts in age (intercept)As bespoke clocks are developed
for use in additional human tissues, we will include these in their own
section below.

``` r
suppressMessages(getClockInfo())
```

<div class="kable-table">

| Clock Name            | 1st Author   | Year |     PMID | Trained Phenotype   | \# of CpGs | Cohort Trained                                               | Tissues Derived                                                                                                  | Age Range Trained                                                                 | Array Type Trained |
|:----------------------|:-------------|-----:|---------:|:--------------------|-----------:|:-------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|:-------------------|
| AdaptAge              | Ying         | 2024 | 38243142 | Chronological Age   |       1000 | London Life Sciences Prospective Population (LOLIPOP)        | Blood                                                                                                            | NA                                                                                | 450K               |
| Alcohol               | McCartney    | 2018 | 30257690 | Clinical Phenotype  |        450 | Generation Scotland: The Scottish Family Health Study \[GS\] | Blood                                                                                                            | 18–98                                                                             | 450K               |
| BMI                   | McCartney    | 2018 | 30257690 | Clinical Phenotype  |       1109 | Generation Scotland: The Scottish Family Health Study \[GS\] | Blood                                                                                                            | 18-98                                                                             | 450K               |
| Bocklandt             | Bocklandt    | 2011 | 21731603 | Chronological Age   |         88 | See Misc                                                     | Saliva                                                                                                           | 21-55                                                                             | 27K                |
| Bohlin                | Bohlin       | 2016 | 27717397 | Gestational Age     |        251 | MoBa1                                                        | NA                                                                                                               | NA                                                                                | 450K               |
| CausAge               | Ying         | 2024 | 38243142 | Chronological Age   |        586 | London Life Sciences Prospective Population (LOLIPOP)        | Blood                                                                                                            | NA                                                                                | 450K               |
| DNAmClockCortical     | Shireby      | 2020 | 33300551 | Chronological Age   |        347 | Multiple cohorts                                             | Brain Tissue                                                                                                     | 1-108                                                                             | 450K               |
| DNAmTL                | Lu           | 2019 |  3142238 | Telomere Length     |        140 | WHI+ JHS- Women’s Health Initiative & Jackson Heart Study    | Blood                                                                                                            | 50.2, 66.5, 80.2 (WHI min, median, max) - 22.2, 56.6, 93.1 (JHS min, median, max) | 450K and EPICv1    |
| DamAge                | Ying         | 2024 | 38243142 | Chronological Age   |       1090 | London Life Sciences Prospective Population (LOLIPOP)        | Blood                                                                                                            | NA                                                                                | 450K               |
| DunedinPoAm38         | Belsky       | 2020 | 32367804 | Pace of Aging       |         47 | Dunedin Study                                                | Blood                                                                                                            | 26-38                                                                             | 450K and EPICv1    |
| DunedinPACE           | Belsky       | 2022 | 35029144 | Pace of Aging       |        173 | Dunedin Study                                                | Blood                                                                                                            | 26-45                                                                             | 450K and EPICv1    |
| EpiTOC1               | Yang         | 2016 | 27716309 | Mitotic Divisions   |        385 | UCSD and WCH (GSE40279)                                      | Blood                                                                                                            | See Misc                                                                          | 450K               |
| EpiTOC2               | Teschendorff | 2020 | 32580750 | Mitotic Divisions   |        163 | UCSD and WCH (GSE40279)                                      | Blood                                                                                                            | 19-101                                                                            | 450K               |
| Garagnani             | Garagnani    | 2012 | 23061750 | Chronological Age   |          1 | See Misc                                                     | Blood                                                                                                            | 42–83 and 9–52                                                                    | 450K               |
| GrimAge1              | Lu           | 2019 | 30669119 | Mortality           |       1030 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| GrimAge2              | Lu           | 2022 | 36516495 | Mortality           |       1030 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | 40 (min), 59 (25th), 66.1(mean), 73(75th), 92 (max)                               | 450K               |
| HRSInCHPhenoAge       | Higgins-Chen | 2022 | 36277076 | Mortality via proxy |        959 | HRS and InCHIANTI                                            | Blood                                                                                                            | 21-100                                                                            | 450K               |
| Hannum                | Hannum       | 2013 | 23177740 | Chronological Age   |         71 | UCSD and WCH (GSE40279)                                      | Blood                                                                                                            | 19-101                                                                            | 450K               |
| Horvath1              | Horvath      | 2013 | 24138928 | Chronological Age   |        353 | Multiple cohorts                                             | 51 healthy tissues and cell types                                                                                | ~3 (min) to 100 (max)                                                             | 450K and 27K       |
| Horvath2              | Horvath      | 2018 | 30048243 | Chronological Age   |        391 | Multiple cohorts                                             | Human fibroblasts, keratinocytes, buccal cells, endothelial cells, lymphoblastoid cells, skin, blood, and saliva | -0.3(min) to 94 (max)                                                             | 450K and EPICv1    |
| Knight                | Knight       | 2016 | 27717399 | Gestational Age     |        148 | Multiple cohorts                                             | umbilical cord blood or blood spots                                                                              | neonates                                                                          | 450K and 27K       |
| LeeControl            | Lee          | 2019 | 31235674 | Gestational Age     |        546 | Multiple cohorts                                             | Placenta                                                                                                         | 5 to 42 weeks gestation                                                           | 450K and EPICv1    |
| LeeRefinedRobust      | Lee          | 2019 | 31235674 | Gestational Age     |        395 | Multiple cohorts                                             | Placenta                                                                                                         | 5 to 42 weeks gestation                                                           | 450K and EPICv1    |
| LeeRobust             | Lee          | 2019 | 31235674 | Gestational Age     |        558 | Multiple cohorts                                             | Placenta                                                                                                         | 5 to 42 weeks gestation                                                           | 450K and EPICv1    |
| Lin                   | Weidner      | 2014 | 24490752 | Chronological Age   |         99 | HNR study                                                    | Blood                                                                                                            | 0-78                                                                              | 27K                |
| Mayne                 | Mayne        | 2017 | 27894195 | Gestational Age     |         62 | Multiple cohorts                                             | Placenta                                                                                                         | 8-42 weeks gestation                                                              | 450K and 27K       |
| MiAge                 | Youn         | 2018 | 29160179 | Mitotic Divisions   |        268 | NA                                                           | 8 cancer types and adjacent tissues                                                                              | NA                                                                                | 450K and EPICv1    |
| PedBE                 | McEwan       | 2019 | 31611402 | Chronological Age   |         94 | Multiple cohorts                                             | Buccal                                                                                                           | 0-20                                                                              | 450K and EPICv1    |
| PhenoAge              | Levine       | 2018 | 29676998 | Mortality via proxy |        513 | InCHIANTI                                                    | Blood                                                                                                            | 21-100                                                                            | 450K               |
| Smoking               | McCartney    | 2018 | 30257690 | Clinical Phenotype  |        233 | Generation Scotland: The Scottish Family Health Study \[GS\] | Blood                                                                                                            | 18–98                                                                             | 450K               |
| VidalBralo            | Vidal-Bralo  | 2016 | 27471517 | Chronological Age   |          8 | Multiple cohorts                                             | Blood                                                                                                            | 20-78                                                                             | 27K                |
| Weidner               | Weidner      | 2014 | 24490752 | Chronological Age   |          3 | Multiple cohorts                                             | Blood                                                                                                            | 0-78                                                                              | 450K               |
| Zhang2019             | Zhang        | 2019 | 31443728 | Chronological Age   |        514 | Multiple cohorts                                             | Blood and saliva                                                                                                 | 2-104                                                                             | 450K               |
| Zhang                 | Zhang        | 2017 | 28303888 | Mortality           |         10 | ESTHER                                                       | Blood                                                                                                            | 50-75                                                                             | 450K               |
| StocZ                 | Tong         | 2024 | 38724732 | Simulated Age       |        514 | Simulated dataset                                            | Blood                                                                                                            | 45-83                                                                             | 450K and EPICv1    |
| StocH                 | Tong         | 2024 | 38724732 | Simulated Age       |        353 | Simulated dataset                                            | Blood                                                                                                            | 45-83                                                                             | 450K and EPICv1    |
| StocP                 | Tong         | 2024 | 38724732 | Simulated Age       |        513 | Simulated dataset                                            | Blood                                                                                                            | 45-83                                                                             | 450K and EPICv1    |
| HypoClock             | Teschendorff | 2020 | 32580750 | Mitotic Divisions   |        678 | See Misc                                                     | See Misc                                                                                                         | See Misc                                                                          | 450K               |
| Retroelement-Age 450K | Ndhlovu      | 2024 | 38106164 | Chronological Age   |         NA | TruDiagnostic BioBank                                        | Blood                                                                                                            | 12-100                                                                            | EPICv1             |
| PCHorvath1            | Higgins-Chen | 2022 | 36277076 | Chronological Age   |      78464 | Multiple cohorts                                             | Multiple                                                                                                         | -0.5 - 105                                                                        | 450K               |
| PCHorvath2            | Higgins-Chen | 2022 | 36277076 | Chronological Age   |      78464 | Multiple cohorts                                             | Skin and Blood                                                                                                   | -0.3 - 101                                                                        | 450K               |
| PCHannum              | Higgins-Chen | 2022 | 36277076 | Chronological Age   |      78464 | UCSD and WCH (GSE40279)                                      | Blood                                                                                                            | 19 - 101                                                                          | 450K               |
| PCPhenoAge            | Higgins-Chen | 2022 | 36277076 | Mortality via proxy |      78464 | HRS and InCHIANTI                                            | Blood                                                                                                            | 21 - 101                                                                          | 450K               |
| PCDNAmTL              | Higgins-Chen | 2022 | 36277076 | Telomere Length     |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCGrimAge             | Higgins-Chen | 2022 | 36277076 | Mortality           |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCPACKYRS             | Higgins-Chen | 2022 | 36277076 | Clinical Phenotype  |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCADM                 | Higgins-Chen | 2022 | 36277076 | Protein             |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCB2M                 | Higgins-Chen | 2022 | 36277076 | Protein             |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCCystatinC           | Higgins-Chen | 2022 | 36277076 | Protein             |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCGDF15               | Higgins-Chen | 2022 | 36277076 | Protein             |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCLeptin              | Higgins-Chen | 2022 | 36277076 | Protein             |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCPAI1                | Higgins-Chen | 2022 | 36277076 | Protein             |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| PCTIMP1               | Higgins-Chen | 2022 | 36277076 | Protein             |      78464 | FHS- Framingham heart study                                  | Blood                                                                                                            | 24 - 92                                                                           | 450K               |
| DNAmPACKYRS           | Lu           | 2019 |  6366976 | Clinical Phenotype  |        172 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| DNAmADM               | Lu           | 2019 |  6366976 | Protein             |        186 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| DNAmB2M               | Lu           | 2019 |  6366976 | Protein             |         91 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| DNAmCystatinC         | Lu           | 2019 |  6366976 | Protein             |         87 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| DNAmGDF15             | Lu           | 2019 |  6366976 | Protein             |        137 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| DNAmPAI1              | Lu           | 2019 |  6366976 | Protein             |        211 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| DNAmlogA1C            | Lu           | 2022 | 36516495 | Clinical Biomarker  |         86 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | 40 (min), 59 (25th), 66.1(mean), 73(75th), 92 (max)                               | 450K               |
| DNAmlogCRP            | Lu           | 2022 | 36516495 | Clinical Biomarker  |        132 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | 40 (min), 59 (25th), 66.1(mean), 73(75th), 92 (max)                               | 450K               |
| DNAmLeptin            | Lu           | 2019 |  6366976 | Protein             |        187 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| DNAmTIMP1             | Lu           | 2019 |  6366976 | Protein             |         42 | FHS- Framingham heart study Offspring Cohort                 | Blood                                                                                                            | ~40(min) to ~90 (max)                                                             | 450K               |
| SystemsAge            | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Blood                 | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Brain                 | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Inflammation          | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Heart                 | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Hormone               | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Immune                | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Kidney                | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Liver                 | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Metabolic             | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Lung                  | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| MusculoSkeletal       | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |
| Age_prediction        | Sehgal       | 2023 | 37503069 | Mortality           |     125175 | HRS and FHS                                                  | Blood                                                                                                            | 24 - 100                                                                          | 450K               |

</div>

#### Running A User-Defined List of Epigenetic Clocks

The user is welcome to specify a vector of clocks that they would like
to calculate, rather than running each individual clock calculation. In
this case, you will need to choose from the following options:

``` r
clockOptions()
#>  [1] "calcAlcoholMcCartney"            "calcBMIMcCartney"               
#>  [3] "calcBocklandt"                   "calcBohlin"                     
#>  [5] "calcClockCategory"               "calcDNAmClockCortical"          
#>  [7] "calcDNAmTL"                      "calcDunedinPoAm38"              
#>  [9] "calcEpiTOC"                      "calcEpiTOC2"                    
#> [11] "calcGaragnani"                   "calcGrimAgeV1"                  
#> [13] "calcGrimAgeV2"                   "calcHannum"                     
#> [15] "calcHorvath1"                    "calcHorvath2"                   
#> [17] "calcHRSInChPhenoAge"             "calcHypoClock"                  
#> [19] "calcKnight"                      "calcLeeControl"                 
#> [21] "calcLeeRefinedRobust"            "calcLeeRobust"                  
#> [23] "calcLin"                         "calcMayne"                      
#> [25] "calcMiAge"                       "calcPCClocks"                   
#> [27] "calcPEDBE"                       "calcPhenoAge"                   
#> [29] "calcSmokingMcCartney"            "calcSystemsAge"                 
#> [31] "calcVidalBralo"                  "calcWeidner"                    
#> [33] "calcZhang"                       "calcZhang2019"                  
#> [35] "prcPhenoAge::calcPRCPhenoAge"    "prcPhenoAge::calcnonPRCPhenoAge"
#> [37] "DunedinPoAm38::PoAmProjector"
```

To do so, here is an example:

``` r
userClocks <- c("calcSmokingMcCartney","calcPhenoAge","calcEpiTOC2")
calcUserClocks(userClocks, exampleBetas, examplePheno, imputation = F)
```

<div class="kable-table">

| name              | geo_accession | gender |  age | group | sample | Smoking_McCartney | PhenoAge |  epiTOC2 |
|:------------------|:--------------|:-------|-----:|------:|-------:|------------------:|---------:|---------:|
| 7786915023_R02C02 | GSM1343050    | M      | 57.9 |     1 |      1 |          3.993508 | 52.29315 | 5012.412 |
| 7786915135_R04C02 | GSM1343051    | M      | 42.0 |     1 |      2 |          4.501657 | 41.05867 | 4622.625 |
| 7471147149_R06C01 | GSM1343052    | M      | 47.4 |     1 |      3 |          3.173744 | 43.54460 | 2956.300 |
| 7786915035_R05C01 | GSM1343053    | M      | 49.3 |     1 |      4 |          3.216788 | 43.96697 | 3446.410 |
| 7786923035_R01C01 | GSM1343054    | M      | 52.5 |     1 |      5 |          4.414541 | 40.35242 | 3245.157 |

</div>

### Missing beta values

Of course, all of the CpG clocks work best when you have all of the
necessary probes’ beta values for each sample. However, sometimes after
preprocessing, beta values will be removed for a variety of reasons. For
each CpG clock, you have the option to impute missing values for CpGs
that were removed across all samples. In this case, you will need to
impute using a vector of your choice (e.g. mean methylation values
across CpGs from an independent tissue-matched dataset). However, by
default, imputation will not be performed and the portion of the clock
that is reliant upon those CpGs will not be considered. To check quickly
whether this is the case for your data and clock(s) of interest, we have
created the following helper function:

``` r
getClockProbes(exampleBetas)
```

<div class="kable-table">

| Clock             | Total.Probes | Present.Probes | Percent.Present |
|:------------------|-------------:|---------------:|:----------------|
| Alcohol           |          450 |            450 | 100%            |
| BMI               |         1109 |           1109 | 100%            |
| Bocklandt         |            1 |              1 | 100%            |
| Bohlin            |          251 |              8 | 3%              |
| DNAmClockCortical |          347 |             33 | 10%             |
| DNAmTL            |          140 |             11 | 8%              |
| EpiToc2           |          163 |            163 | 100%            |
| EpiToc            |          385 |            385 | 100%            |
| Garagnani         |            1 |              1 | 100%            |
| GrimAge1          |         1139 |             83 | 7%              |
| GrimAge2          |         1362 |            111 | 8%              |
| HRSInCHPhenoAge   |          959 |            959 | 100%            |
| Hannum            |           71 |             71 | 100%            |
| Horvath1          |          353 |            353 | 100%            |
| Horvath2          |          391 |            390 | 100%            |
| Knight            |          148 |             16 | 11%             |
| LeeControl        |          546 |             13 | 2%              |
| LeeRefinedRobust  |          395 |              9 | 2%              |
| LeeRobust         |          558 |              9 | 2%              |
| Lin               |           99 |             39 | 39%             |
| Mayne             |           62 |              5 | 8%              |
| MiAge             |          268 |              4 | 1%              |
| PCClocks          |        78464 |           2976 | 4%              |
| PEDBE             |           94 |             14 | 15%             |
| PhenoAge          |          513 |            513 | 100%            |
| Smoking           |          233 |            233 | 100%            |
| SystemsAge        |       125175 |           4024 | 3%              |
| VidalBralo        |            8 |              5 | 62%             |
| Weidner           |            3 |              3 | 100%            |
| Zhang2019         |          514 |            131 | 25%             |
| Zhang             |           10 |             10 | 100%            |
| hypoClock         |          678 |            678 | 100%            |

</div>

Please note that this will not count columns of NAs for named CpGs as
missing! If you want to check for this you can run the following line of
code to find the column numbers that are all NAs. If you get “named
integer(0)” then you don’t have any. We recommend that you remove any
identified columns from your beta matrix entirely to avoid errors, and
then rerun the code producing the table above.

``` r
which(apply(exampleBetas, 2, function(x)all(is.na(x))))
```

In the case that you have CpGs missing from only some samples, we
encourage you to be aware of this early on. Run the following line and
check that it is 0.

``` r
sum(is.na(betaMatrix))
```

If this does not end up being 0, you might consider running mean
imputation within your data so that NA values for single/ few samples at
least have mean values rather than being ignored.
