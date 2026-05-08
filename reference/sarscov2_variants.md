# SARS-CoV-2 VOC/VOI Curated Variant Metadata

A pre-built named list containing curated biological and epidemiological
metadata for 12 SARS-CoV-2 Variants of Concern (VOC) and Variants of
Interest (VOI): Alpha, Beta, Delta, Epsilon, Eta, Gamma, Iota, Kappa,
Lambda, Omicron, Theta, Zeta. Includes mutation profiles, nomenclature,
temporal detection records, defining SNPs, and a fully citable reference
table with 21 verified references.

## Usage

``` r
data(sarscov2_variants)
```

## Format

A named list of 13 elements (12 metadata fields plus a structured
References object). See Details for per-element descriptions.

## Source

Compiled from 21 peer-reviewed and surveillance sources; see
`sarscov2_variants$References$data` for the full reference table with
DOIs and URLs. Built from `SARS_CoV_2_VOC_VOI.xlsx` via
[`get_variants`](https://vadimtyuryaev.github.io/ViralEntropR/reference/get_variants.md).

## Details

The object is available automatically after
[`library(ViralEntropR)`](https://github.com/vadimtyuryaev/ViralEntropR)
via lazy loading. It is produced by
[`get_variants`](https://vadimtyuryaev.github.io/ViralEntropR/reference/get_variants.md)
from the curated Excel workbook `SARS_CoV_2_VOC_VOI.xlsx` (not bundled
with the package; see `data-raw/sarscov2_variants.R` for the
reproducible build script).

**List elements:**

- `WHO_Label`:

  List of 12. WHO variant label strings (e.g. `"Alpha"`, `"Omicron"`).

- `Pango_Lineage`:

  List of 12. Pango lineage designations (e.g. `"B.1.1.7"`).

- `GISAID_Clade`:

  Character(12). GISAID clade strings (e.g. `"GRY"`, `"GR/484A"`).

- `Nextstrain_Clade`:

  Character(12). Nextstrain clade strings (e.g. `"20I/501Y.V1"`).

- `Country_First_Detected`:

  List of 12. Country of first documented detection per variant.

- `Date_Earliest_Sample`:

  Character(12). Month-Year of earliest documented sample (e.g.
  `"Sep-2020"`).

- `Date_First_Detected`:

  Character(12). Month-Year of world-level first detection.

- `Date_First_Detected_US`:

  Character(12). Month-Year of first US detection.

- `Spike_Mutations`:

  List of 12. Character vectors of Spike protein mutations per variant,
  read from the Excel workbook.

- `Mutation_Sites`:

  List of 12. Integer vectors of amino acid positions extracted from
  `Spike_Mutations`.

- `Defining_SNPs`:

  List of 12. Canonical defining SNP strings (`NA` for variants not
  characterised in the reference set).

- `Defining_SNP_Sites`:

  List of 12. Integer positions extracted from `Defining_SNPs`.

- `References`:

  Named list with three elements: `$data` – data frame of 21 verified
  references (columns: ID, Authors, Year, Title, Journal_Source,
  Volume_Issue_Pages, DOI, URL, Type, Variants_Covered, Data_Field,
  Citation); `$display(variant = NULL)` – renders an interactive
  [`datatable`](https://rdrr.io/pkg/DT/man/datatable.html), optionally
  filtered by WHO label (requires the DT package); `$cite(variant)` –
  returns formatted citation strings for a given WHO label, suitable for
  manuscript use.

## References

The 21 curated references cover peer-reviewed articles, CDC MMWR
government reports, and outbreak.info surveillance database records.
Full provenance is available via `sarscov2_variants$References$data` or
the interactive display.

## Examples

``` r
# Object is available immediately after library(ViralEntropR)

# --- Basic access ---------------------------------------------------------

# All 12 WHO labels
unlist(sarscov2_variants$WHO_Label)
#>  [1] "Alpha"   "Beta"    "Delta"   "Epsilon" "Eta"     "Gamma"   "Iota"   
#>  [8] "Kappa"   "Lambda"  "Omicron" "Theta"   "Zeta"   

# Pango lineage for Alpha
sarscov2_variants$Pango_Lineage[[
  which(unlist(sarscov2_variants$WHO_Label) == "Alpha")
]]
#> [1] "B.1.1.7"

# World detection dates for all variants
data.frame(
  Variant  = unlist(sarscov2_variants$WHO_Label),
  Detected = sarscov2_variants$Date_First_Detected,
  Country  = unlist(sarscov2_variants$Country_First_Detected)
)
#>    Variant Detected         Country
#> 1    Alpha Dec-2020  United Kingdom
#> 2     Beta Oct-2020    South Africa
#> 3    Delta Dec-2020           India
#> 4  Epsilon Jul-2020             USA
#> 5      Eta Nov-2020         Nigeria
#> 6    Gamma Dec-2020          Brazil
#> 7     Iota Nov-2020             USA
#> 8    Kappa Oct-2020           India
#> 9   Lambda Dec-2020            Peru
#> 10 Omicron Nov-2021    South Africa
#> 11   Theta Feb-2021 The Philippines
#> 12    Zeta Apr-2020          Brazil

# --- Mutation sites -------------------------------------------------------

# Defining SNPs and sites for Delta
idx <- which(unlist(sarscov2_variants$WHO_Label) == "Delta")
sarscov2_variants$Defining_SNPs[[idx]]
#> [1] "T19R"  "L452R" "T478K" "P681R" "D950N"
sarscov2_variants$Defining_SNP_Sites[[idx]]
#> [1]  19 452 478 681 950

# All Spike mutation sites for Omicron
sarscov2_variants$Mutation_Sites[[
  which(unlist(sarscov2_variants$WHO_Label) == "Omicron")
]]
#>  [1]  67  95 142 145 211 212 339 346 371 373 375 417 440 446 477 478 484 493 496
#> [20] 498 501 505 547 614 655 679 681 764 796 856 954 969 981

# --- Nomenclature ---------------------------------------------------------

data.frame(
  WHO        = unlist(sarscov2_variants$WHO_Label),
  PANGO      = unlist(sarscov2_variants$Pango_Lineage),
  GISAID     = sarscov2_variants$GISAID_Clade,
  Nextstrain = sarscov2_variants$Nextstrain_Clade
)
#>        WHO           PANGO      GISAID  Nextstrain
#> 1    Alpha         B.1.1.7         GRY 20I/501Y.V1
#> 2     Beta         B.1.351  GH/501Y.V2 20H/501Y.V2
#> 3    Delta       B.1.617.2   G/478K.V1  21A/S:478K
#> 4  Epsilon B.1.427/B.1.429  GH/452R.V1         21C
#> 5      Eta         B.1.525   G/484K.V3  21D/S:484K
#> 6    Gamma             P.1  GR/501Y.V3 20J/501Y.V3
#> 7     Iota         B.1.526  GH/253G.V1         21F
#> 8    Kappa       B.1.617.1   G/452R.V3  21B/S:154K
#> 9   Lambda            C.37  GR/452Q.V1         21G
#> 10 Omicron       B.1.1.529     GR/484A         21K
#> 11   Theta             P.3 GR/1092K.V1  20B/S:265C
#> 12    Zeta             P.2  GR/484K.V2  20B/S:484K

# --- References -----------------------------------------------------------

# Full reference data frame (21 rows)
sarscov2_variants$References$data[,
  c("ID", "Authors", "Year", "Journal_Source")
]
#>    ID                                    Authors Year
#> 1   1                   Ghosh N, Nandi S, Saha I 2022
#> 2   2                        Abulsoud AI, et al. 2023
#> 3   3           Aleem A, Akbar Samad AB, Vaqar S 2023
#> 4   4 Zella D, Giovanetti M, Benedetti F, et al. 2021
#> 5   5  Wink PL, Volpato FCZ, Monteiro FL, et al. 2022
#> 6   6       Chatterjee S, Bhattacharya M, et al. 2023
#> 7   7       Washington NL, Gangavarapu K, et al. 2021
#> 8   8  Long SW, Olsen RJ, Christensen PA, et al. 2021
#> 9   9         Zhang W, Davis BD, Chen SS, et al. 2021
#> 10 10       Gangavarapu K, Abdel Latif A, et al. 2023
#> 11 11                     Chen C, Shi Q, Dong XP 2021
#> 12 12  Firestone MJ, Lorentz AJ, Meyer S, et al. 2021
#> 13 13                 CDC COVID-19 Response Team 2021
#> 14 14       Gangavarapu K, Abdel Latif A, et al. 2023
#> 15 15       Gangavarapu K, Abdel Latif A, et al. 2023
#> 16 16       Gangavarapu K, Abdel Latif A, et al. 2023
#> 17 17     Kannan SR, Spratt AN, Sharma K, et al. 2022
#> 18 18        Dhawan M, Saied AA, Mitra S, et al. 2022
#> 19 19           Tao K, Tzou PL, Nouhin J, et al. 2021
#> 20 20                                 Brookes AJ 1999
#> 21 21            Salleh MZ, Derrick JP, Deris ZZ 2021
#>                                  Journal_Source
#> 1              International Immunopharmacology
#> 2                 Biomedicine & Pharmacotherapy
#> 3  StatPearls [Internet], StatPearls Publishing
#> 4                   Journal of Medical Virology
#> 5     Infection Control & Hospital Epidemiology
#> 6                                       Viruses
#> 7                                          Cell
#> 8             The American Journal of Pathology
#> 9                                          JAMA
#> 10                                outbreak.info
#> 11                                     Zoonoses
#> 12                    MMWR Morb Mortal Wkly Rep
#> 13                    MMWR Morb Mortal Wkly Rep
#> 14                                outbreak.info
#> 15                                outbreak.info
#> 16                                outbreak.info
#> 17                      Journal of Autoimmunity
#> 18                Biomedicine & Pharmacotherapy
#> 19                      Nature Reviews Genetics
#> 20                                         Gene
#> 21  International Journal of Molecular Sciences

# Formatted citation strings for Gamma
sarscov2_variants$References$cite("Gamma")
#> [1] "Ghosh N, Nandi S, Saha I (2022). A review on evolution of emerging SARS-CoV-2 variants based on spike glycoprotein. International Immunopharmacology. 105:108565. doi:10.1016/j.intimp.2022.108565"                                                                                                                                                                      
#> [2] "Abulsoud AI, et al. (2023). Mutations in SARS-CoV-2: Insights on structure, variants, vaccines, and biomedical interventions. Biomedicine & Pharmacotherapy. 157:113977. doi:10.1016/j.biopha.2022.113977"                                                                                                                                                               
#> [3] "Aleem A, Akbar Samad AB, Vaqar S (2023). Emerging Variants of SARS-CoV-2 and Novel Therapeutics Against Coronavirus (COVID-19). StatPearls [Internet], StatPearls Publishing."                                                                                                                                                                                           
#> [4] "Zella D, Giovanetti M, Benedetti F, et al. (2021). The variants question: What is the problem?. Journal of Medical Virology. 93(12):6479-6485. doi:10.1002/jmv.27196"                                                                                                                                                                                                    
#> [5] "Long SW, Olsen RJ, Christensen PA, et al. (2021). Sequence Analysis of 20,453 Severe Acute Respiratory Syndrome Coronavirus 2 Genomes from the Houston Metropolitan Area Identifies the Emergence and Widespread Distribution of Multiple Isolates of All Major Variants of Concern. The American Journal of Pathology. 191(6):983-992. doi:10.1016/j.ajpath.2021.03.004"
#> [6] "Firestone MJ, Lorentz AJ, Meyer S, et al. (2021). First Identified Cases of SARS-CoV-2 Variant P.1 in the United States - Minnesota, January 2021. MMWR Morb Mortal Wkly Rep. 70(10):346-347. doi:10.15585/mmwr.mm7010e1"                                                                                                                                                
#> [7] "Tao K, Tzou PL, Nouhin J, et al. (2021). The biological and clinical significance of emerging SARS-CoV-2 variants. Nature Reviews Genetics. 22(12):757-773. doi:10.1038/s41576-021-00408-x"                                                                                                                                                                              
#> [8] "Salleh MZ, Derrick JP, Deris ZZ (2021). Structural Evaluation of the Spike Glycoprotein Variants on SARS-CoV-2 Transmission and Immune Evasion. International Journal of Molecular Sciences. 22(14):7425. doi:10.3390/ijms22147425"                                                                                                                                      

if (FALSE) { # \dontrun{
# Interactive DT table (requires the DT package)
sarscov2_variants$References$display()

# Filtered to Omicron references only
sarscov2_variants$References$display(variant = "Omicron")
} # }
```
