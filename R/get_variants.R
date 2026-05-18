#' @title Parse SARS-CoV-2 VOC/VOI Variant Metadata from Excel
#'
#' @description Reads a structured Excel workbook of SARS-CoV-2 Variants of
#'   Concern (VOC) and Variants of Interest (VOI) and returns a named list
#'   containing mutation profiles, nomenclature, temporal detection metadata,
#'   defining SNPs, and a fully citable reference table with interactive
#'   display support.
#'
#' @details
#' Intended to be called once during data preparation via
#' \code{data-raw/sarscov2_variants.R}:
#'
#' \preformatted{
#'   variants_dat  <- readxl::read_excel("SARS_CoV_2_VOC_VOI.xlsx")
#'   variants_list <- get_variants(variants_dat)
#'   saveRDS(variants_list,
#'           file = "inst/extdata/sarscov2_variants.rds")
#' }
#'
#' The saved object can then be loaded anywhere in the package or vignettes:
#'
#' \preformatted{
#'   voc_data <- readRDS(system.file("extdata", "sarscov2_variants.rds",
#'                                   package = "ViralEntropR"))
#' }
#'
#' \strong{Column order of variants in the Excel workbook:}
#' Alpha, Beta, Epsilon, Eta, Iota, Kappa, Delta, Lambda, Gamma, Zeta,
#' Theta, Omicron.
#'
#' \strong{Returned list elements:}
#' \describe{
#'   \item{\code{WHO_Label}}{List. WHO variant label strings (e.g. \code{"Alpha"}).}
#'   \item{\code{Pango_Lineage}}{List. Pango lineage designations.}
#'   \item{\code{GISAID_Clade}}{Character vector. GISAID clade strings.}
#'   \item{\code{Nextstrain_Clade}}{Character vector. Nextstrain clade strings.}
#'   \item{\code{Country_First_Detected}}{List. Country of first documented
#'     detection per variant.}
#'   \item{\code{Date_Earliest_Sample}}{Character vector. Month-Year of the
#'     earliest documented sample per variant.}
#'   \item{\code{Date_First_Detected}}{Character vector. Month-Year of
#'     world-level first detection.}
#'   \item{\code{Date_First_Detected_US}}{Character vector. Month-Year of
#'     first US detection.}
#'   \item{\code{Spike_Mutations}}{List. Spike protein mutation strings read
#'     from the Excel workbook.}
#'   \item{\code{Mutation_Sites}}{List. Integer site positions extracted from
#'     \code{Spike_Mutations}.}
#'   \item{\code{Defining_SNPs}}{List. Canonical defining SNP strings per
#'     variant (\code{NA} where not characterised).}
#'   \item{\code{Defining_SNP_Sites}}{List. Integer positions extracted from
#'     \code{Defining_SNPs}.}
#'   \item{\code{References}}{Named list with three elements:
#'     \code{$data} — data frame of 21 verified references;
#'     \code{$display(variant = NULL)} — interactive \code{\link[DT]{datatable}},
#'     optionally filtered by WHO label;
#'     \code{$cite(variant)} — character vector of formatted citation strings
#'     suitable for manuscript use.}
#' }
#'
#' @param tibble A \code{tibble} or \code{data.frame} produced by
#'   \code{\link[readxl]{read_excel}} from \code{SARS_CoV_2_VOC_VOI.xlsx}.
#'   Column 1 must be the \code{Variant} field column; columns 2 onward are
#'   one column per variant in the order listed above.
#' @param check Logical. If \code{TRUE}, validates that all internal vectors
#'   have equal length and stops with an informative message on mismatch.
#'   Default \code{FALSE}.
#'
#' @return A named list with 13 elements; see Details for full descriptions.
#'
#' @keywords internal
get_variants = function(tibble, check = FALSE) {
  
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Package 'readxl' is required. Install with: install.packages('readxl')",
         call. = FALSE)
  if (!requireNamespace("DT", quietly = TRUE))
    stop("Package 'DT' is required. Install with: install.packages('DT')",
         call. = FALSE)
  
  n  <- dim(tibble)[1]
  m  <- dim(tibble)[2]
  mm <- m - 1L
  
  # ── Temporal metadata ────────────────────────────────────────────────────────
  # Order follows variant column order in the Excel sheet (see header above).
  
  # Alphabetical order: Alpha, Beta, Delta, Epsilon, Eta, Gamma, Iota, Kappa,
  #                     Lambda, Omicron, Theta, Zeta
  
  earliest_documented_samples <- c(
    "Sep-2020",  # Alpha   -- ref 2
    "May-2020",  # Beta    -- ref 2
    "Oct-2020",  # Delta   -- ref 2
    "Mar-2020",  # Epsilon -- ref 2
    "Nov-2020",  # Eta     -- ref 3
    "Nov-2020",  # Gamma   -- ref 2
    "Nov-2020",  # Iota    -- ref 2
    "Oct-2020",  # Kappa   -- ref 2
    "Aug-2020",  # Lambda  -- ref 5
    "Nov-2021",  # Omicron -- ref 6
    "Jan-2021",  # Theta   -- ref 2
    "Apr-2020"   # Zeta    -- ref 2
  )
  
  detection_date <- c(
    "Dec-2020",  # Alpha   -- ref 3
    "Oct-2020",  # Beta    -- ref 3
    "Dec-2020",  # Delta   -- ref 3
    "Jul-2020",  # Epsilon -- ref 9
    "Nov-2020",  # Eta     -- ref 3
    "Dec-2020",  # Gamma   -- ref 3
    "Nov-2020",  # Iota    -- ref 3
    "Oct-2020",  # Kappa   -- ref 4
    "Dec-2020",  # Lambda  -- ref 11
    "Nov-2021",  # Omicron -- ref 3
    "Feb-2021",  # Theta   -- ref 3
    "Apr-2020"   # Zeta    -- ref 3
  )
  
  US_detection_date <- c(
    "Dec-2020",  # Alpha   -- ref 7
    "Jan-2021",  # Beta    -- ref 8
    "Mar-2021",  # Delta   -- ref 3
    "Jul-2020",  # Epsilon -- ref 9  ---> in 1 of 1247 samples, then October 2020
    "Nov-2020",  # Eta     -- ref 3
    "Jan-2021",  # Gamma   -- ref 12
    "Nov-2020",  # Iota    -- ref 3
    "Feb-2021",  # Kappa   -- ref 10
    "Feb-2021",  # Lambda  -- ref 11
    "Dec-2021",  # Omicron -- ref 13
    "Jan-2021",  # Theta   -- ref 15
    "Oct-2020"   # Zeta    -- ref 14
  )
  
  # -- Defining SNPs ------------------------------------------------------------
  
  defining_snps <- list(
    c("N501Y", "A570D", "P681H", "T716I", "S982A", "D1118H"),                   # Alpha   -- ref 4
    c("D80A", "D215G", "K417N", "N501Y", "E484K", "A701V"),                     # Beta    -- ref 4
    c("T19R", "L452R", "T478K", "P681R", "D950N"),                              # Delta   -- ref 4
    NA,                                                                          # Epsilon -- not in ref 4
    c("Q52R", "E484K", "Q677H", "F888L"),                                        # Eta     -- ref 4
    c("L18F", "T20N", "P26S", "D138Y", "R190S", "K417T",
      "E484K", "N501Y", "H655Y", "T1027I"),                                      # Gamma   -- ref 4
    NA,                                                                          # Iota    -- not in ref 4
    c("E154K", "L452R", "E484Q", "P681R"),                                       # Kappa   -- ref 4
    NA,                                                                          # Lambda  -- not in ref 4
    c("A76V", "T95I", "L212I", "G339D", "S371L","S373P", "S375F",
      "K417N", "N440K", "G446S", "S477N","T478K", "E484A", "Q493R",
      "G496S", "Q498R", "N501Y","Y505H", "T547K", "D614G", "H655Y",
      "N679K", "P681H","N764K", "D796Y", "N856K", "Q954H", "N969K", "L981F"),   # Omicron -- refs 6, 17, 18
    NA,                                                                          # Theta   -- not in ref 4
    NA                                                                           # Zeta    -- not in ref 4
  )
  
  # -- Nomenclature -------------------------------------------------------------
  # Note: GISAID clades are not unique across variants.
  
  GISAID_clade <- c(
    "GRY",         # Alpha
    "GH/501Y.V2",  # Beta
    "G/478K.V1",   # Delta
    "GH/452R.V1",  # Epsilon
    "G/484K.V3",   # Eta
    "GR/501Y.V3",  # Gamma
    "GH/253G.V1",  # Iota
    "G/452R.V3",   # Kappa
    "GR/452Q.V1",  # Lambda
    "GR/484A",     # Omicron
    "GR/1092K.V1", # Theta
    "GR/484K.V2"   # Zeta
  )
  
  Nextstrain_clade <- c(
    "20I/501Y.V1",  # Alpha
    "20H/501Y.V2",  # Beta
    "21A/S:478K",   # Delta
    "21C",          # Epsilon
    "21D/S:484K",   # Eta
    "20J/501Y.V3",  # Gamma
    "21F",          # Iota
    "21B/S:154K",   # Kappa
    "21G",          # Lambda
    "21K",          # Omicron
    "20B/S:265C",   # Theta
    "20B/S:484K"    # Zeta
  )
  
  # -- WHO label, Pango lineage, country -- read from Excel ---------------------
  
  who_label     <- vector(length = mm, mode = "list")
  pango_lineage <- vector(length = mm, mode = "list")
  country       <- vector(length = mm, mode = "list")
  
  for (i in seq_len(mm)) {
    who_label[[i]]     <- colnames(tibble)[i + 1L]
    pango_lineage[[i]] <- as.character(tibble[1L, i + 1L])
    country[[i]]       <- as.character(tibble[2L, i + 1L])
  }
  
  # -- Mutations -- read from Excel (rows 3 onward) -----------------------------
  
  tibble_temp <- tibble[3:n, ]
  mutations   <- vector(length = mm, mode = "list")
  
  for (i in seq_len(mm)) {
    mutations[[i]] <- c(tibble_temp[!is.na(tibble_temp[, i + 1L]), ]$Variant)
  }
  
  # -- Site extraction ----------------------------------------------------------
  
  get_sites <- function(mutation) {
    as.numeric(gsub("[^0-9]", "", mutation))
  }
  
  sites <- lapply(mutations, get_sites)
  
  # -- Reorder Excel-read vectors to alphabetical order ----------------------
  # Excel column order: Alpha, Beta, Epsilon, Eta, Iota, Kappa,
  #                     Delta, Lambda, Gamma, Zeta, Theta, Omicron
  # Alphabetical order: Alpha, Beta, Delta, Epsilon, Eta, Gamma,
  #                     Iota, Kappa, Lambda, Omicron, Theta, Zeta
  # Reorder index (position in Excel for each alphabetical slot):
  alpha_order <- c("Alpha", "Beta", "Delta", "Epsilon", "Eta", "Gamma",
                   "Iota", "Kappa", "Lambda", "Omicron", "Theta", "Zeta")
  excel_labels <- unlist(who_label)
  reorder_idx  <- match(alpha_order, excel_labels)
  
  who_label     <- who_label[reorder_idx]
  pango_lineage <- pango_lineage[reorder_idx]
  country       <- country[reorder_idx]
  mutations     <- mutations[reorder_idx]
  sites         <- sites[reorder_idx]
  
  defining_snps_sites <- lapply(
    defining_snps,
    function(x) if (length(x) == 1L && is.na(x)) NA else get_sites(x)
  )
  
  # -- Integrity check (optional) -----------------------------------------------
  
  if (isTRUE(check)) {
    vec_list <- list(
      who_label, pango_lineage, GISAID_clade, Nextstrain_clade,
      country, earliest_documented_samples, detection_date,
      US_detection_date, mutations, sites, defining_snps, defining_snps_sites
    )
    if (length(unique(sapply(vec_list, length))) != 1L)
      stop("Vectors have unequal lengths -- check Excel sheet dimensions.", call. = FALSE)
  }
  
  # -- Reference table ----------------------------------------------------------
  #
  # Type vocabulary:
  #   "Peer-reviewed article" -- journal article with DOI
  #   "Government report"     -- CDC MMWR public health release
  #   "Surveillance database" -- outbreak.info live database records (backed by
  #                             Gangavarapu et al. 2023, Nature Methods)
  
  refs_data <- data.frame(
    
    ID = 1:21,
    
    Authors = c(
      "Ghosh N, Nandi S, Saha I",                     #  1
      "Abulsoud AI, et al.",                          #  2
      "Aleem A, Akbar Samad AB, Vaqar S",             #  3
      "Zella D, Giovanetti M, Benedetti F, et al.",   #  4
      "Wink PL, Volpato FCZ, Monteiro FL, et al.",    #  5
      "Chatterjee S, Bhattacharya M, et al.",         #  6
      "Washington NL, Gangavarapu K, et al.",         #  7
      "Long SW, Olsen RJ, Christensen PA, et al.",    #  8
      "Zhang W, Davis BD, Chen SS, et al.",           #  9
      "Gangavarapu K, Abdel Latif A, et al.",         # 10
      "Chen C, Shi Q, Dong XP",                       # 11
      "Firestone MJ, Lorentz AJ, Meyer S, et al.",    # 12
      "CDC COVID-19 Response Team",                   # 13
      "Gangavarapu K, Abdel Latif A, et al.",         # 14
      "Gangavarapu K, Abdel Latif A, et al.",         # 15
      "Gangavarapu K, Abdel Latif A, et al.",         # 16
      "Kannan SR, Spratt AN, Sharma K, et al.",       # 17
      "Dhawan M, Saied AA, Mitra S, et al.",          # 18
      "Tao K, Tzou PL, Nouhin J, et al.",             # 19
      "Brookes AJ",                                   # 20
      "Salleh MZ, Derrick JP, Deris ZZ"              # 21
    ),
    
    Year = c(
      2022L,  #  1 -- Ghosh
      2023L,  #  2 -- Abulsoud
      2023L,  #  3 -- Aleem (last updated May 2023; original PMID: 34033342)
      2021L,  #  4 -- Zella
      2022L,  #  5 -- Wink (epub Sep 2021; published in volume Dec 2022)
      2023L,  #  6 -- Chatterjee
      2021L,  #  7 -- Washington
      2021L,  #  8 -- Long
      2021L,  #  9 -- Zhang
      2023L,  # 10 -- Gangavarapu et al. (outbreak.info, Nature Methods 2023)
      2021L,  # 11 -- Chen/Cao
      2021L,  # 12 -- Firestone
      2021L,  # 13 -- CDC MMWR
      2023L,  # 14 -- Gangavarapu (outbreak.info Zeta)
      2023L,  # 15 -- Gangavarapu (outbreak.info Theta)
      2023L,  # 16 -- Gangavarapu (outbreak.info Lambda)
      2022L,  # 17 -- Kannan (epub Dec 2021; published 2022)
      2022L,  # 18 -- Dhawan
      2021L,  # 19 -- Tao
      1999L,  # 20 -- Brookes
      2021L   # 21 -- Salleh
    ),
    
    Title = c(
      "A review on evolution of emerging SARS-CoV-2 variants based on spike glycoprotein",
      "Mutations in SARS-CoV-2: Insights on structure, variants, vaccines, and biomedical interventions",
      "Emerging Variants of SARS-CoV-2 and Novel Therapeutics Against Coronavirus (COVID-19)",
      "The variants question: What is the problem?",
      "First identification of SARS-CoV-2 lambda (C.37) variant in Southern Brazil",
      "A Detailed Overview of SARS-CoV-2 Omicron: Its Sub-Variants, Mutations and Pathophysiology, Clinical Characteristics, Immunological Landscape, Immune Escape, and Therapies",
      "Emergence and rapid transmission of SARS-CoV-2 B.1.1.7 in the United States",
      "Sequence Analysis of 20,453 Severe Acute Respiratory Syndrome Coronavirus 2 Genomes from the Houston Metropolitan Area Identifies the Emergence and Widespread Distribution of Multiple Isolates of All Major Variants of Concern",
      "Emergence of a Novel SARS-CoV-2 Variant in Southern California",
      "Kappa (B.1.617.1) Variant Report",
      "SARS-CoV-2 Lambda Variant: Spatiotemporal Distribution and Potential Public Health Impact",
      "First Identified Cases of SARS-CoV-2 Variant P.1 in the United States - Minnesota, January 2021",
      "SARS-CoV-2 B.1.1.529 (Omicron) Variant - United States, December 1-8, 2021",
      "Zeta (P.2) Variant Report",
      "Theta (P.3) Variant Report",
      "Lambda (C.37) Variant Report",
      "Omicron SARS-CoV-2 variant: Unique features and their impact on pre-existing antibodies",
      "Omicron variant (B.1.1.529) and its sublineages: What do we know so far amid the emergence of recombinant variants of SARS-CoV-2?",
      "The biological and clinical significance of emerging SARS-CoV-2 variants",
      "The essence of SNPs",
      "Structural Evaluation of the Spike Glycoprotein Variants on SARS-CoV-2 Transmission and Immune Evasion"
    ),
    
    Journal_Source = c(
      "International Immunopharmacology",              #  1
      "Biomedicine & Pharmacotherapy",                 #  2
      "StatPearls [Internet], StatPearls Publishing",  #  3
      "Journal of Medical Virology",                   #  4
      "Infection Control & Hospital Epidemiology",     #  5
      "Viruses",                                       #  6
      "Cell",                                          #  7
      "The American Journal of Pathology",             #  8
      "JAMA",                                          #  9
      "outbreak.info",                                 # 10
      "Zoonoses",                                      # 11
      "MMWR Morb Mortal Wkly Rep",                     # 12
      "MMWR Morb Mortal Wkly Rep",                     # 13
      "outbreak.info",                                 # 14
      "outbreak.info",                                 # 15
      "outbreak.info",                                 # 16
      "Journal of Autoimmunity",                       # 17
      "Biomedicine & Pharmacotherapy",                 # 18
      "Nature Reviews Genetics",                       # 19
      "Gene",                                          # 20
      "International Journal of Molecular Sciences"    # 21
    ),
    
    Volume_Issue_Pages = c(
      "105:108565",            #  1 -- Ghosh
      "157:113977",            #  2 -- Abulsoud
      NA_character_,           #  3 -- StatPearls (no volume/pages)
      "93(12):6479-6485",      #  4 -- Zella
      "43(12):1996-1997",      #  5 -- Wink
      "15(1):167",             #  6 -- Chatterjee
      "184(10):2587-2594.e7",  #  7 -- Washington
      "191(6):983-992",        #  8 -- Long
      "325(13):1324-1326",     #  9 -- Zhang
      NA_character_,           # 10 -- outbreak.info (database record)
      "1(1)",                  # 11 -- Chen/Cao, Zoonoses vol 1 issue 1
      "70(10):346-347",        # 12 -- Firestone MMWR
      "70(50):1731-1734",      # 13 -- CDC MMWR
      NA_character_,           # 14 -- outbreak.info
      NA_character_,           # 15 -- outbreak.info
      NA_character_,           # 16 -- outbreak.info
      "126:102779",            # 17 -- Kannan
      "154:113522",            # 18 -- Dhawan
      "22(12):757-773",        # 19 -- Tao
      "234(2):177-186",        # 20 -- Brookes
      "22(14):7425"            # 21 -- Salleh
    ),
    
    DOI = c(
      "10.1016/j.intimp.2022.108565",   #  1 -- Ghosh
      "10.1016/j.biopha.2022.113977",   #  2 -- Abulsoud
      NA_character_,                    #  3 -- StatPearls (PMID: 34033342; no DOI)
      "10.1002/jmv.27196",              #  4 -- Zella
      "10.1017/ice.2021.390",           #  5 -- Wink
      "10.3390/v15010167",              #  6 -- Chatterjee
      "10.1016/j.cell.2021.03.052",     #  7 -- Washington
      "10.1016/j.ajpath.2021.03.004",   #  8 -- Long
      "10.1001/jama.2021.1612",         #  9 -- Zhang
      "10.1038/s41592-023-01769-3",     # 10 -- Gangavarapu (outbreak.info)
      "10.15212/ZOONOSES-2021-0009",    # 11 -- Chen/Cao
      "10.15585/mmwr.mm7010e1",         # 12 -- Firestone
      "10.15585/mmwr.mm7050e1",         # 13 -- CDC MMWR
      "10.1038/s41592-023-01769-3",     # 14 -- Gangavarapu (outbreak.info Zeta)
      "10.1038/s41592-023-01769-3",     # 15 -- Gangavarapu (outbreak.info Theta)
      "10.1038/s41592-023-01769-3",     # 16 -- Gangavarapu (outbreak.info Lambda)
      "10.1016/j.jaut.2021.102779",     # 17 -- Kannan
      "10.1016/j.biopha.2022.113522",   # 18 -- Dhawan
      "10.1038/s41576-021-00408-x",     # 19 -- Tao
      "10.1016/S0378-1119(99)00219-X",  # 20 -- Brookes
      "10.3390/ijms22147425"            # 21 -- Salleh
    ),
    
    URL = c(
      "https://www.sciencedirect.com/science/article/pii/S1567576922000492",              #  1
      "https://www.sciencedirect.com/science/article/pii/S075333222201366X",              #  2
      "https://pubmed.ncbi.nlm.nih.gov/34033342/",                                        #  3
      "https://onlinelibrary.wiley.com/doi/full/10.1002/jmv.27196",                       #  4
      "https://pmc.ncbi.nlm.nih.gov/articles/PMC8564022/",                                #  5
      "https://www.mdpi.com/1999-4915/15/1/167",                                          #  6
      "https://www.sciencedirect.com/science/article/pii/S0092867421003834",              #  7
      "https://www.sciencedirect.com/science/article/pii/S0002944021001085",              #  8
      "https://jamanetwork.com/journals/jama/fullarticle/2776543",                        #  9
      "https://outbreak.info/situation-reports/kappa",                                    # 10
      "https://www.scienceopen.com/hosted-document?doi=10.15212/ZOONOSES-2021-0009",      # 11
      "https://www.cdc.gov/mmwr/volumes/70/wr/mm7010e1.htm",                              # 12
      "https://www.cdc.gov/mmwr/volumes/70/wr/mm7050e1.htm",                              # 13
      "https://outbreak.info/situation-reports/Zeta",                                     # 14
      "https://outbreak.info/situation-reports/Theta",                                    # 15
      "https://outbreak.info/situation-reports/lambda",                                   # 16
      "https://www.sciencedirect.com/science/article/pii/S0896841121001876",              # 17
      "https://www.sciencedirect.com/science/article/pii/S0753332222009118",              # 18
      "https://www.nature.com/articles/s41576-021-00408-x",                               # 19
      "https://www.sciencedirect.com/science/article/pii/S037811199900219X",               # 20
      "https://www.mdpi.com/1422-0067/22/14/7425"                                          # 21
    ),
    
    Type = c(
      "Peer-reviewed article",  #  1
      "Peer-reviewed article",  #  2
      "Peer-reviewed article",  #  3 
      "Peer-reviewed article",  #  4
      "Peer-reviewed article",  #  5
      "Peer-reviewed article",  #  6
      "Peer-reviewed article",  #  7
      "Peer-reviewed article",  #  8
      "Peer-reviewed article",  #  9
      "Surveillance database",  # 10
      "Peer-reviewed article",  # 11
      "Government report",      # 12
      "Government report",      # 13
      "Surveillance database",  # 14
      "Surveillance database",  # 15
      "Surveillance database",  # 16
      "Peer-reviewed article",  # 17
      "Peer-reviewed article",  # 18
      "Peer-reviewed article",  # 19
      "Peer-reviewed article",  # 20
      "Peer-reviewed article"   # 21
    ),
    
    Variants_Covered = c(
      "Alpha, Beta, Epsilon, Eta, Iota, Kappa, Delta, Lambda, Gamma, Zeta, Theta, Omicron",    #  1
      "Alpha, Beta, Epsilon, Eta, Iota, Kappa, Delta, Lambda, Gamma, Zeta, Theta, Omicron",    #  2
      "Alpha, Beta, Epsilon, Eta, Iota, Kappa, Delta, Lambda, Gamma, Zeta, Theta, Omicron",    #  3
      "Alpha, Beta, Epsilon, Eta, Iota, Kappa, Delta, Lambda, Gamma, Zeta, Theta",             #  4
      "Lambda",                                                                                #  5
      "Omicron",                                                                               #  6
      "Alpha",                                                                                 #  7
      "Alpha, Beta, Epsilon, Gamma, Zeta",                                                     #  8
      "Epsilon",                                                                               #  9
      "Kappa",                                                                                 # 10
      "Lambda",                                                                                # 11
      "Gamma",                                                                                 # 12
      "Omicron",                                                                               # 13
      "Zeta",                                                                                  # 14
      "Theta",                                                                                 # 15
      "Lambda",                                                                                # 16
      "Omicron",                                                                               # 17
      "Omicron",                                                                               # 18
      "Alpha, Beta, Delta, Gamma",                                                             # 19
       NA_character_,                                                                          # 20
      "Alpha, Beta, Delta, Epsilon, Eta, Gamma, Iota, Kappa, Lambda, Theta, Zeta"              # 21
    ),
    
    Data_Field = c(
      "Mutations, Sites, GISAID_clade, Nextstrain_clade, M-Y_Earliest_Documented_Sample",         #  1 -- Ghosh: phylogenetic analysis of 77,681 genomes; GISAID clade assignments; variant emergence timeline
      "General review of SARS-CoV-2 mutations and variant biology",                               #  2 -- Abulsoud: broad review; does not provide the specific dated records used in this list
      "M-Y_Detection, M-Y_Earliest_Documented_Sample (Eta), M-Y_Detection_US (Eta, Iota, Delta)", #  3 -- Aleem StatPearls: archived review; primary source for world detection dates and US detection dates lacking dedicated primary studies
      "Defining_SNPs, Defining_SNPs_Sites, M-Y_Detection (Kappa)",                                #  4 -- Zella: VOC/VOI mutational characterization; also source for Kappa world detection date (Oct-2020)
      "Lambda emergence context",                                                                 #  5 -- Wink: first Lambda identification in Southern Brazil mid-2021; supports Lambda spread narrative
      "Defining_SNPs, Defining_SNPs_Sites, M-Y_Earliest_Documented_Sample (Omicron)",             #  6 -- Chatterjee: detailed Omicron overview; source for Omicron earliest documented sample (Nov-2021); Omicron world detection date now sourced from ref 3
      "M-Y_Detection_US (Alpha)",                                                                 #  7 -- Washington: B.1.1.7 US emergence and logistic growth analysis
      "M-Y_Detection_US (Beta)",                                                                  #  8 -- Long: B.1.351 identified in Houston (replaces SC DHEC press release)
      "M-Y_Detection_US, M-Y_Detection (Epsilon)",                                                #  9 -- Zhang: Epsilon (B.1.427/B.1.429) first identified in Southern California; US origin means US detection = world detection (Jul-2020)
      "M-Y_Detection_US (Kappa)",                                                                 # 10 -- outbreak.info Kappa live report
      "M-Y_Detection, M-Y_Detection_US (Lambda)",                                                 # 11 -- Chen/Cao: spatiotemporal distribution; source for Lambda world detection (Dec-2020) and first documented US presence (Feb-2021)
      "M-Y_Detection_US (Gamma)",                                                                 # 12 -- Firestone: first P.1 cases in Minnesota, January 2021
      "M-Y_Detection_US (Omicron)",                                                               # 13 -- CDC MMWR: Omicron US December 1-8, 2021
      "M-Y_Detection_US (Zeta)",                                                                  # 14 -- outbreak.info Zeta live report
      "M-Y_Detection_US (Theta)",                                                                 # 15 -- outbreak.info Theta live report
      "M-Y_Earliest_Documented_Sample (Lambda)",                                                  # 16 -- outbreak.info Lambda: earliest GISAID record Aug-2020 (Peru); Feb-2021 first appearance in US databases
      "Defining_SNPs, Defining_SNPs_Sites (Omicron)",                                             # 17 -- Kannan: 46 high-prevalence Omicron mutations; structural antibody impact
      "Defining_SNPs, Defining_SNPs_Sites (Omicron)",                                             # 18 -- Dhawan: Omicron sublineages and recombinant variants
      "General evolutionary context",                                                             # 19 -- Tao: biological and clinical significance of SARS-CoV-2 VOCs
      "SNP methodology",                                                                          # 20 -- Brookes: foundational SNP definition and population genetics methodology
      "GISAID_Clade, Nextstrain_Clade"                                                            # 21 -- Salleh: structural evaluation of spike variants; source for GISAID/Nextstrain clade assignments
    ),
    
    stringsAsFactors = FALSE
  )
  
  # Assemble formatted citation string per row for manuscript use.
  # Format: Authors (Year). Title. Journal. Vol(Issue):Pages. doi:DOI
  # Components absent for StatPearls, outbreak.info, etc. are omitted gracefully.
  refs_data$Citation <- mapply(function(auth, yr, ttl, jrnl, vip, doi) {
    parts <- c(
      sprintf("%s (%d).", auth, yr),
      sprintf("%s.", ttl),
      if (!is.na(jrnl)) sprintf("%s.", jrnl),
      if (!is.na(vip))  sprintf("%s.", vip),
      if (!is.na(doi))  sprintf("doi:%s", doi)
    )
    paste(parts, collapse = " ")
  },
  refs_data$Authors, refs_data$Year, refs_data$Title,
  refs_data$Journal_Source, refs_data$Volume_Issue_Pages, refs_data$DOI,
  SIMPLIFY = TRUE)
  
  # -- Reference accessor functions ---------------------------------------------
  
  # display(variant = NULL)
  #   Renders an interactive, filterable DT table.
  #   Supply a WHO label (e.g. "Alpha") to show only references for that variant.
  display_refs <- function(variant = NULL) {
    
    if (!requireNamespace("DT", quietly = TRUE))
      stop("Package 'DT' is required to display the reference table.\n",
           "Install it with: install.packages('DT')",
           call. = FALSE)
    
    df <- refs_data
    
    if (!is.null(variant)) {
      keep <- grepl(variant, df$Variants_Covered, fixed = TRUE)
      if (!any(keep))
        stop(sprintf("No references found for variant '%s'.", variant), 
             call. = FALSE)
      df <- df[keep, ]
    }
    
    df$Title_Link <- sprintf(
      '<a href="%s" target="_blank">%s</a>', df$URL, df$Title
    )
    
    display_cols <- c("ID", "Authors", "Year", "Title_Link",
                      "Journal_Source", "Volume_Issue_Pages",
                      "Type", "Variants_Covered", "Data_Field")
    col_names    <- c("ID", "Authors", "Year", "Title",
                      "Journal / Source", "Vol(Issue):Pages",
                      "Type", "Variants Covered", "Data Field")
    
    caption_text <- if (!is.null(variant))
      sprintf("References for variant: %s", variant)
    else
      "ViralEntropR - SARS-CoV-2 VOC/VOI Reference Table (21 verified references)"
    
    DT::datatable(
      df[, display_cols],
      colnames  = col_names,
      escape    = FALSE,
      rownames  = FALSE,
      filter    = "top",
      caption   = caption_text,
      options   = list(
        paging  = FALSE,       # show all rows -- no pagination (fixes Rmd truncation)
        scrollX = TRUE,
        dom     = "fti",       # filter + table + info; no pagination controls
        columnDefs = list(
          list(width = "240px", targets = 3L),  # Title
          list(width = "180px", targets = 7L),  # Variants Covered
          list(width = "200px", targets = 8L)   # Data Field
        )
      )
    )
  }
  
  # cite(variant)
  #   Returns a character vector of formatted citation strings for a given
  #   WHO variant label -- suitable for copy-pasting into a manuscript.
  #   Example: variants_list$References$cite("Omicron")
  cite_refs <- function(variant) {
    keep <- grepl(variant, refs_data$Variants_Covered, fixed = TRUE)
    if (!any(keep))
      stop(sprintf("No references found for variant '%s'.", variant), call. = FALSE)
    refs_data$Citation[keep]
  }
  
  # -- Assemble and return ------------------------------------------------------
  #
  # All existing list elements are preserved exactly.
  # References is a named list with three elements:
  #   $data    -- the full reference data frame (21 rows x 11 columns)
  #   $display -- function: renders interactive DT, optionally filtered by variant
  #   $cite    -- function: returns formatted citation strings for a variant
  
  list(
    WHO_Label              = who_label,
    Pango_Lineage          = pango_lineage,
    GISAID_Clade           = GISAID_clade,
    Nextstrain_Clade       = Nextstrain_clade,
    Country_First_Detected = country,
    Date_Earliest_Sample   = earliest_documented_samples,
    Date_First_Detected    = detection_date,
    Date_First_Detected_US = US_detection_date,
    Spike_Mutations        = mutations,
    Mutation_Sites         = sites,
    Defining_SNPs          = defining_snps,
    Defining_SNP_Sites     = defining_snps_sites,
    References             = list(
      data    = refs_data,
      display = display_refs,
      cite    = cite_refs
    )
  )
}