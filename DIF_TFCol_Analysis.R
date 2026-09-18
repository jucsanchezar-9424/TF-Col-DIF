# =============================================================================
# Differential Item Functioning of the TF-Col Financial Toxicity Measure
# Across Clinical and Socioeconomic Subgroups in Colombian Oncology
#
# Author:  Julio Cesar Sanchez-Arevalo · jucsanchezar@unal.edu.co
#          Universidad Nacional de Colombia / Fundacion Universitaria de
#          Ciencias de la Salud, Bogota, Colombia
#
# Companion to the instrument development and validation study:
#   Sanchez-Arevalo JC, Sierra-Matamoros FA, Perez Sosa NM, Cuellar Rivera DI,
#   Sanchez Pedraza R. Development and psychometric validation of the TF-Col.
#   Value in Health Regional Issues 2026;101695.
#   https://doi.org/10.1016/j.vhri.2026.101695
#
# METHOD
#   Item-level DIF by ordinal logistic regression with iterative purification,
#   as implemented in lordif (Zumbo 1999; Choi, Gibbons & Crane 2011).
#   Items classified by Nagelkerke dR2 (Jodoin & Gierl 2001).
#
# TWO THRESHOLDS, ONE CRITERION -- read this before the results
#   0.020  internal purification threshold used by lordif to decide which items
#          leave the anchor set. Items reaching it are POTENTIAL SIGNALS.
#   0.035  prespecified criterion for MEANINGFUL DIF. This is the one that
#          classifies an item. No item in this study reached it.
#   Four items reached 0.020: IF07 and IPS08 (respondent type/site), IF10 and
#   IPS10 (insurance regime). They are reported transparently as signals, never
#   as DIF.
#
# MATCHING VARIABLE
#   lordif does NOT accept an external matching variable. It calibrates a graded
#   response model on the items supplied and conditions on the EAP trait
#   estimate, re-derived from the current anchor set at each purification
#   iteration. The total sum score is computed below only as a descriptive
#   summary and for the correlation reported in Methods.
#
# COMPARISON TIERS
#   Primary     (n >= 100 in both groups): respondent type/site, sex, insurance.
#   Exploratory (n <  100 in one group):   education, stratum, employment.
#   The split reflects attainable precision and was NOT preregistered.
#
# DATA
#   The participant-level dataset is not distributed: it contains sensitive
#   patient information governed by the confidentiality conditions of the
#   participating institutions. De-identified data may be requested from the
#   corresponding author with institutional approval.
#
#   ONE SHEET. The analysis reads a single sheet, `ANALISIS`, of the workbook
#   BDTFCol_ANALISIS.xlsx. Point the script at it with an environment variable:
#       Sys.setenv(TFCOL_DATA = "/full/path/BDTFCol_ANALISIS.xlsx")
#       Sys.setenv(TFCOL_OUT  = "./output")
#
#   IF YOU IMPORT THE SHEET BY HAND (RStudio's Import Dataset, or your own
#   read_excel call), assign it to an object named `datos` before sourcing this
#   file. Section 2 detects it and skips the file read entirely, so no path has
#   to resolve on your machine:
#       datos <- readxl::read_excel("BDTFCol_ANALISIS.xlsx", sheet = "ANALISIS")
#       source("DIF_TFCol_Analysis.R")
#
# OUTPUTS
#   Table1_Sample.csv          Sample characteristics (Table 1)
#   Table2_DIF.csv             Item-level DIF, all comparisons (Table 2)
#   TableS1_DTF.csv            Differential test functioning (Online Resource 6)
#   TableS2_Impact.csv         Trait distribution by group (Online Resource 7)
#   DomainDIF_sensitivity.csv  DIF repeated within each domain
#   ESM_2_ICC_RespondentType.pdf   ICC, items flagged at 0.020 (Online Resource 2)
#   ESM_3_ICC_InsuranceRegime.pdf  ICC, items flagged at 0.020 (Online Resource 3)
#   ESM_4_TCC.pdf              Test characteristic curves (Online Resource 4)
#   ESM_5_DTF.pdf              DTF vs measurement precision (Online Resource 5)
#   SessionInfo.txt            R session information
#
# REFERENCES
#   Zumbo BD (1999). A Handbook on the Theory and Methods of DIF. Ottawa:
#     Directorate of Human Resources Research and Evaluation.
#   Choi SW, Gibbons LE, Crane PK (2011). lordif. J Stat Softw 39(8):1-30.
#   Jodoin MG, Gierl MJ (2001). Appl Meas Educ 14(4):329-349.
#   Olsson I et al. (2026). Qual Life Res 35(1):24.
# =============================================================================


# =============================================================================
# SECTION 0: CONFIGURATION
# =============================================================================

# The analysis reads ONE sheet. BDTFCol_ANALISIS.xlsx contains a single sheet,
# `ANALISIS`: sheet V4 of the source workbook (348 records, 20 items, grouping
# variables, summed scores) with two columns appended from sheet theta_TRI --
# theta_EAP and SE_theta -- merged on the participant key (V4 `Record ID` <->
# theta_TRI `ID`). Nothing was recoded, renamed, reordered or imputed.
ruta_datos <- Sys.getenv("TFCOL_DATA",  unset = "BDTFCol_ANALISIS.xlsx")
hoja_datos <- Sys.getenv("TFCOL_SHEET", unset = "ANALISIS")
dir_output <- Sys.getenv("TFCOL_OUT",   unset = ".")

UMBRAL_PURIFICACION <- 0.020   # lordif internal threshold -> potential signal
UMBRAL_MODERADO     <- 0.035   # Jodoin & Gierl: meaningful DIF (the criterion)
UMBRAL_GRANDE       <- 0.070   # Jodoin & Gierl: large DIF
ALPHA_LRT           <- 0.05
MIN_CASILLA         <- 5       # minimum cell frequency per response category
MIN_N_GRUPO         <- 100     # n per group separating primary from exploratory
SEED                <- 2025

GRID  <- seq(-3, 3, length.out = 121)      # reporting grid for DTF
FINE  <- seq(-4.5, 4.5, length.out = 3001) # fine grid for curve inversion
N_CAT <- 4                                 # thresholds per item (5 options)

set.seed(SEED)


# =============================================================================
# SECTION 1: PACKAGES
# =============================================================================

paquetes <- c("lordif", "mirt", "readxl", "dplyr", "MASS")
faltan <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]
if (length(faltan))
  stop("Missing packages: ", paste(faltan, collapse = ", "),
       "\nInstall them with: install.packages(c(\"",
       paste(faltan, collapse = "\", \""), "\"))")
suppressPackageStartupMessages({
  library(lordif); library(readxl); library(dplyr)
})
# MASS::polr is called with its namespace prefix so that MASS::select does not
# mask dplyr::select.

cat(R.version$version.string, "\n")
cat("lordif", as.character(packageVersion("lordif")), "\n\n")

dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# SECTION 2: DATA LOADING AND PREPARATION
# =============================================================================

# If `datos` was already imported by hand, use it; otherwise read the workbook.
if (exists("datos") && is.data.frame(datos)) {
  cat("Using the data frame `datos` already in the environment.\n")
} else {
  if (!file.exists(ruta_datos))
    stop("Data file not found: '", ruta_datos, "'.\n",
         "Either set TFCOL_DATA to its full path, or import the sheet by hand ",
         "into an object named `datos` and source this script again.")
  cat("Loading:", ruta_datos, "| sheet:", hoja_datos, "\n")
  datos <- read_excel(ruta_datos, sheet = hoja_datos)
}
datos <- as.data.frame(datos)
cat("Records loaded:", nrow(datos), "| columns:", ncol(datos), "\n")

# The sheet must be the analytic one, whatever route it arrived by. V3 of the
# source workbook has the same 348 rows and the same 20 items and is easy to
# grab by mistake; it does not carry theta_EAP and SE_theta, so this stops.
requeridas <- c("theta_EAP", "SE_theta", "PUNTAJE_TF")
faltan_col <- requeridas[!requeridas %in% names(datos)]
if (length(faltan_col))
  stop("Columns missing from the data: ", paste(faltan_col, collapse = ", "),
       "\nThis is not the ANALISIS sheet. Sheet V4 alone lacks theta_EAP and ",
       "SE_theta; see the DATA note in the header.")
if (nrow(datos) != 348)
  warning("Expected 348 records, found ", nrow(datos),
          ". Check that the imported sheet is ANALISIS.")
cat("\n")

# --- 2.1 Locale-safe accessors -----------------------------------------------
# Column names and category labels in the workbook contain accented characters.
# Matching them against literal UTF-8 strings fails SILENTLY when the session
# locale is not UTF-8: the subgroup comes back empty and the script runs to
# completion with meaningless output. Everything below matches ASCII substrings.

col <- function(patron) {
  i <- grep(patron, names(datos), fixed = TRUE)[1]
  if (is.na(i)) stop("Column not found: ", patron)
  datos[[i]]
}
has <- function(x, patron) grepl(patron, as.character(x), fixed = TRUE)

# --- 2.2 Item response matrix (N x 20) ---------------------------------------
items_IF   <- sprintf("imf%02d_IF",   1:10)
items_IPS  <- sprintf("imps%02d_Ips", 1:10)
items_cols <- c(items_IF, items_IPS)
etiquetas  <- c(sprintf("IF%02d", 1:10), sprintf("IPS%02d", 1:10))

resp <- as.data.frame(lapply(datos[, items_cols], as.integer))

cat("Item range check (expected 0-4):",
    min(sapply(resp, min, na.rm = TRUE)), "to",
    max(sapply(resp, max, na.rm = TRUE)),
    "| missing item responses:", sum(is.na(resp)), "\n")

# --- 2.3 Descriptive total score ---------------------------------------------
# Descriptive only. NOT the conditioning variable -- see header.
score_total <- rowSums(resp, na.rm = TRUE)

# --- 2.4 Grouping variables ---------------------------------------------------
tipo_pob <- col("Tipo de poblaci")
sexo     <- col("5_sexo_paciente")
reg_raw  <- col("26_regimen_now")
edu_raw  <- col("28_niveleducativo")
est_raw  <- col("25_estrato_now")
trabajo  <- col("29_trabajo_now")
ciudad   <- col("1_ciudad")
dx_cat   <- col("11_dx_cat_paciente")

# Insurance regime: pragmatic grouping into subsidized vs non-subsidized. The
# non-subsidized category pools contributory with the exception and special
# regimes because those two together comprise only 11 respondents and cannot be
# analysed separately. No claim is made that their coverage or copayment
# structures are functionally equivalent; the heterogeneity is a limitation.
regimen_bin <- ifelse(has(reg_raw, "Subsidiado"), "Subsidiado",
  ifelse(has(reg_raw, "Contributivo") | has(reg_raw, "Excepci") |
         has(reg_raw, "Especial"), "No-subsidiado", NA))
n_exc_esp <- sum(has(reg_raw, "Excepci") | has(reg_raw, "Especial"), na.rm = TRUE)

educ_bin <- ifelse(has(edu_raw, "Sin educaci") | has(edu_raw, "Primaria") |
                   has(edu_raw, "Secundaria")  | has(edu_raw, "cnico o tecn"),
                   "Basica-Media",
  ifelse(has(edu_raw, "Pregrado") | has(edu_raw, "Posgrado"), "Superior", NA))

estrato_bin <- ifelse(has(est_raw, "Estrato 1") | has(est_raw, "Estrato 2") |
                      has(est_raw, "Sin estrato"), "Bajo",
  ifelse(has(est_raw, "Estrato 3") | has(est_raw, "Estrato 4") |
         has(est_raw, "Estrato 5") | has(est_raw, "Estrato 6"), "Medio-Alto", NA))

# Guard against the silent-locale failure described in 2.1
stopifnot(sum(!is.na(regimen_bin)) == nrow(datos),
          sum(!is.na(estrato_bin)) == nrow(datos),
          sum(tipo_pob %in% c("Adulto", "Infantil")) == nrow(datos))

# --- 2.5 The six comparisons --------------------------------------------------
COMPARACIONES <- list(
  list(clave = "respondent", etiqueta = "Respondent type / site", nivel = "Primary",
       g = tipo_pob,    ref = "Adulto",       foc = "Infantil",
       ref_en = "Adult patient",      foc_en = "Caregiver"),
  list(clave = "sex",        etiqueta = "Sex", nivel = "Primary",
       g = sexo,        ref = "mujer",        foc = "hombre",
       ref_en = "Female",             foc_en = "Male"),
  list(clave = "insurance",  etiqueta = "Health insurance regime", nivel = "Primary",
       g = regimen_bin, ref = "Subsidiado",   foc = "No-subsidiado",
       ref_en = "Subsidized",         foc_en = "Non-subsidized"),
  list(clave = "education",  etiqueta = "Educational attainment", nivel = "Exploratory",
       g = educ_bin,    ref = "Basica-Media", foc = "Superior",
       ref_en = "Basic-intermediate", foc_en = "Higher"),
  list(clave = "stratum",    etiqueta = "Socioeconomic stratum", nivel = "Exploratory",
       g = estrato_bin, ref = "Bajo",         foc = "Medio-Alto",
       ref_en = "Low",                foc_en = "Middle-high"),
  list(clave = "employment", etiqueta = "Employment status", nivel = "Exploratory",
       g = trabajo,     ref = "Si",           foc = "No",
       ref_en = "Employed",           foc_en = "Not employed")
)

cat("\nComparison sizes:\n")
for (cmp in COMPARACIONES) {
  n1 <- sum(cmp$g == cmp$ref, na.rm = TRUE)
  n2 <- sum(cmp$g == cmp$foc, na.rm = TRUE)
  nivel <- if (min(n1, n2) >= MIN_N_GRUPO) "Primary" else "Exploratory"
  if (nivel != cmp$nivel)
    warning("Tier mismatch for ", cmp$etiqueta, ": computed ", nivel)
  cat(sprintf("  %-24s %-12s %s = %3d | %s = %3d\n",
              cmp$etiqueta, cmp$nivel, cmp$ref_en, n1, cmp$foc_en, n2))
}
cat("\nException/special regimes pooled into non-subsidized: n =", n_exc_esp, "\n")


# =============================================================================
# SECTION 3: ASSUMPTION CHECKS
# =============================================================================

cat("\n=== ASSUMPTION CHECKS ===\n")

# Unidimensionality. The graded response model that lordif calibrates assumes
# essential unidimensionality. The validation study reported ECV = 0.694 and
# omega_h = 0.700, supporting a dominant general factor. Both sit close to
# rather than comfortably above conventional benchmarks, which is why the
# domain-level sensitivity analysis in Section 6 is run.
cat("Unidimensionality: ECV = 0.694, omega_h = 0.700 (validation study).\n")
cat("Domain-level sensitivity analysis run in Section 6 as a check.\n\n")

revisar_casillas <- function(resp_m, g, ref, foc, etiqueta) {
  for (nivel in c(ref, foc)) {
    sub <- resp_m[!is.na(g) & g == nivel, , drop = FALSE]
    malos <- names(sub)[sapply(sub, function(y) {
      f <- table(factor(y, levels = 0:4)); any(f > 0 & f < MIN_CASILLA)
    })]
    if (length(malos))
      cat(sprintf("  %-24s %-14s cells < %d in: %s\n",
                  etiqueta, nivel, MIN_CASILLA, paste(malos, collapse = ", ")))
  }
}
cat("Cell frequency check (flagging categories with fewer than",
    MIN_CASILLA, "observations):\n")
for (cmp in COMPARACIONES)
  revisar_casillas(resp, cmp$g, cmp$ref, cmp$foc, cmp$etiqueta)


# =============================================================================
# SECTION 4: ITEM-LEVEL DIF (Table 2)
# =============================================================================

correr_dif <- function(resp_m, g, ref, foc, etiqueta) {
  keep <- g %in% c(ref, foc)
  set.seed(SEED)
  # lordif reports its purification iterations by printing them, and stores no
  # iteration counter in the returned object. Capturing the console trace is the
  # only route to the count that does not rely on undocumented internals; it
  # also keeps the log tidy.
  traza <- utils::capture.output(
    fit <- lordif(resp.data = resp_m[keep, , drop = FALSE],
                  group     = factor(g[keep], levels = c(ref, foc)),
                  criterion = "R2", pseudo.R2 = "Nagelkerke",
                  alpha     = ALPHA_LRT, maxIter = 10, minCell = MIN_CASILLA)
  )
  fit$keep   <- keep
  fit$n_iter <- sum(grepl("Iteration:", traza, fixed = TRUE))
  fit
}

extraer <- function(fit, cmp) {
  st  <- fit$stats
  d13 <- st[, "pseudo13.Nagelkerke"]
  data.frame(
    Comparison   = cmp$etiqueta,
    Tier         = cmp$nivel,
    Groups       = paste(cmp$ref_en, "vs", cmp$foc_en),
    n_ref        = sum(cmp$g[fit$keep] == cmp$ref),
    n_foc        = sum(cmp$g[fit$keep] == cmp$foc),
    Item         = etiquetas,
    dR2_uniform  = round(st[, "pseudo12.Nagelkerke"], 4),
    dR2_nonunif  = round(st[, "pseudo23.Nagelkerke"], 4),
    dR2_total    = round(d13, 4),
    p_uniform    = round(st[, "chi12"], 4),
    p_nonuniform = round(st[, "chi23"], 4),
    p_total      = round(st[, "chi13"], 4),
    Signal_020   = d13 >= UMBRAL_PURIFICACION,   # potential signal, NOT DIF
    DIF_035      = d13 >= UMBRAL_MODERADO,       # the criterion
    Category     = ifelse(d13 < UMBRAL_MODERADO, "A",
                   ifelse(d13 < UMBRAL_GRANDE,   "B", "C")),
    stringsAsFactors = FALSE)
}

cat("\n\n=== ITEM-LEVEL DIF ===\n")
ajustes <- list(); tabla2 <- NULL
for (cmp in COMPARACIONES) {
  cat("\n", strrep("-", 62), "\n", cmp$etiqueta, " (", cmp$nivel, ")\n", sep = "")
  fit <- correr_dif(resp, cmp$g, cmp$ref, cmp$foc, cmp$etiqueta)
  ajustes[[cmp$clave]] <- fit
  tb <- extraer(fit, cmp)
  tabla2 <- rbind(tabla2, tb)
  cat("  purification iterations:", fit$n_iter,
      "| max dR2:", sprintf("%.4f", max(tb$dR2_total)),
      "(", tb$Item[which.max(tb$dR2_total)], ")\n")
  cat("  potential signals (>= 0.020):",
      if (any(tb$Signal_020)) paste(tb$Item[tb$Signal_020], collapse = ", ") else "none", "\n")
  cat("  meaningful DIF   (>= 0.035):",
      if (any(tb$DIF_035)) paste(tb$Item[tb$DIF_035], collapse = ", ") else "none", "\n")
}

# --- Consistency checks -------------------------------------------------------
# Model 2 nests Model 1, so Nagelkerke R2 cannot decrease when the interaction
# term is added: dR2_total must always be >= dR2_uniform. A violation means the
# values were extracted or transcribed wrongly, not that the model misbehaved.
violaciones <- subset(tabla2, dR2_total < dR2_uniform - 1e-9)
if (nrow(violaciones)) {
  warning("dR2_total < dR2_uniform in ", nrow(violaciones), " cell(s).")
  print(violaciones[, c("Comparison", "Item", "dR2_uniform", "dR2_total")])
} else {
  cat("\nOK: dR2_total >= dR2_uniform in all", nrow(tabla2), "cells.\n")
}

cat("\nEffect-size summary across all", nrow(tabla2), "item-by-comparison tests:\n")
cat(sprintf("  median %.4f | mean %.4f | max %.4f | < 0.005: %d | >= 0.020: %d | >= 0.035: %d\n",
            median(tabla2$dR2_total), mean(tabla2$dR2_total), max(tabla2$dR2_total),
            sum(tabla2$dR2_total < 0.005), sum(tabla2$Signal_020), sum(tabla2$DIF_035)))

senales <- tabla2[tabla2$Signal_020, c("Comparison", "Item", "dR2_total", "Category")]
cat("\nPotential signals reported in the paper:\n"); print(senales, row.names = FALSE)


# =============================================================================
# SECTION 5: DOMAIN-LEVEL SENSITIVITY ANALYSIS
# =============================================================================
# The TF-Col has a correlated two-factor structure. Because lordif calibrates
# the graded response model on the item set it is given, passing one domain at a
# time yields a domain-specific trait estimate. This checks that no signal was
# masked by conditioning on a single trait estimate over all 20 items.

cat("\n\n=== DOMAIN-LEVEL SENSITIVITY ===\n")
DOMINIOS <- list(
  list(nombre = "Financial Impact",    cols = items_IF,  labs = etiquetas[1:10]),
  list(nombre = "Psychosocial Impact", cols = items_IPS, labs = etiquetas[11:20])
)

tabla_dom <- NULL
for (dom in DOMINIOS) {
  resp_dom <- as.data.frame(lapply(datos[, dom$cols], as.integer))
  for (cmp in COMPARACIONES) {
    keep <- cmp$g %in% c(cmp$ref, cmp$foc)
    set.seed(SEED)
    fit <- lordif(resp.data = resp_dom[keep, , drop = FALSE],
                  group     = factor(cmp$g[keep], levels = c(cmp$ref, cmp$foc)),
                  criterion = "R2", pseudo.R2 = "Nagelkerke",
                  alpha = ALPHA_LRT, maxIter = 10, minCell = MIN_CASILLA)
    d13 <- fit$stats[, "pseudo13.Nagelkerke"]
    tabla_dom <- rbind(tabla_dom, data.frame(
      Domain = dom$nombre, Comparison = cmp$etiqueta, Tier = cmp$nivel,
      max_dR2 = round(max(d13), 4), max_item = dom$labs[which.max(d13)],
      n_signals_020 = sum(d13 >= UMBRAL_PURIFICACION),
      items_020 = paste(dom$labs[d13 >= UMBRAL_PURIFICACION], collapse = ", "),
      n_dif_035 = sum(d13 >= UMBRAL_MODERADO),
      stringsAsFactors = FALSE))
    cat(sprintf("  %-20s %-24s max dR2 %.4f (%s) | >=0.020: %d | >=0.035: %d\n",
                dom$nombre, cmp$etiqueta, max(d13), dom$labs[which.max(d13)],
                sum(d13 >= UMBRAL_PURIFICACION), sum(d13 >= UMBRAL_MODERADO)))
  }
}
cat("\nDomain-level items reaching the 0.035 criterion:", sum(tabla_dom$n_dif_035), "\n")


# =============================================================================
# SECTION 6: DIFFERENTIAL TEST FUNCTIONING
# =============================================================================
# Two models, answering different questions.
#
#  (1) Signal-restricted. Only the items flagged at 0.020 in a given comparison
#      take group-specific parameters; the rest are constrained equal. This is
#      the analogue here of the naive-vs-purified contrast of Olsson et al.
#      (2026). It is defined only for comparisons that have flagged items; where
#      there are none it is zero by definition, which is a property of the
#      definition and not a finding.
#
#  (2) All-items-free. Every item may differ between groups. A CONSERVATIVE
#      STRESS TEST of possible cumulative distortion -- deliberately not called
#      an upper bound, since item differences can cancel and the procedure also
#      absorbs sampling variability.
#
# One common metric throughout: a single graded response calibration of the full
# analytic sample, no per-comparison rescaling. SE(theta) from the validation
# calibration (column SE_theta) is mapped linearly onto that metric before any
# comparison against it.

cat("\n\n=== DIFFERENTIAL TEST FUNCTIONING ===\n")

# --- 6.1 Common metric --------------------------------------------------------
set.seed(SEED)
grm_full   <- mirt::mirt(resp, 1, itemtype = "graded", verbose = FALSE)
theta_full <- as.numeric(mirt::fscores(grm_full, method = "EAP")[, 1])

pendiente  <- sd(theta_full) / sd(as.numeric(datos[["theta_EAP"]]))
se_mapeado <- as.numeric(datos[["SE_theta"]]) * pendiente
REFERENCIA <- unname(quantile(se_mapeado, 0.25))
cat(sprintf("Metric link slope: %.5f | reference value (25th pct of SE): %.3f theta\n",
            pendiente, REFERENCIA))
cat(sprintf("Agreement with the validation calibration: r = %.5f\n",
            cor(theta_full, as.numeric(datos[["theta_EAP"]]))))
cat("The reference value gauges practical smallness. It is not an inferential test.\n")

# --- 6.2 Estimation helpers ---------------------------------------------------
# Item parameters conditioning on a FIXED theta: a proportional-odds fit of the
# item on theta gives the discrimination directly and the thresholds as zeta/a.
parametros <- function(y, theta) {
  yf <- factor(y)
  if (nlevels(yf) < 2) return(list(a = 0, b = rep(NA_real_, N_CAT)))
  fit <- suppressWarnings(MASS::polr(yf ~ theta, method = "logistic", Hess = FALSE))
  a <- unname(coef(fit)[1])
  b <- if (a != 0) unname(fit$zeta) / a else unname(fit$zeta)
  length(b) <- N_CAT
  list(a = a, b = b)
}

puntaje_esperado <- function(theta, pars) {
  es <- numeric(length(theta))
  for (p in pars)
    for (k in seq_len(N_CAT))
      if (!is.na(p$b[k])) es <- es + plogis(p$a * (theta - p$b[k]))
  es
}

# libres = indices allowed to take group-specific parameters
modelo_dtf <- function(resp_m, grp, theta, libres) {
  comun <- lapply(seq_len(ncol(resp_m)), function(j) parametros(resp_m[, j], theta))
  pars  <- list(ref = comun, foc = comun)
  for (j in libres) {
    pars$ref[[j]] <- parametros(resp_m[grp == 0, j], theta[grp == 0])
    pars$foc[[j]] <- parametros(resp_m[grp == 1, j], theta[grp == 1])
  }
  es_ref <- puntaje_esperado(GRID, pars$ref)
  es_foc <- puntaje_esperado(GRID, pars$foc)

  fino <- puntaje_esperado(FINE, pars$ref); o <- order(fino)
  d_theta <- approx(fino[o], FINE[o], xout = es_foc, rule = 2)$y - GRID

  dens <- rowSums(exp(-0.5 * ((outer(GRID, theta, "-")) / 0.3)^2))
  dens <- dens / sum(dens)
  lo <- unname(quantile(theta, .025)); hi <- unname(quantile(theta, .975))
  w  <- GRID >= lo & GRID <= hi

  list(es_ref = es_ref, es_foc = es_foc, d_score = es_foc - es_ref,
       d_theta = d_theta, lo = lo, hi = hi,
       max_score = max(abs((es_foc - es_ref)[w])),
       max_theta = max(abs(d_theta[w])),
       media_theta = sum(d_theta * dens),
       pct_sobre = mean(abs(d_theta[w]) > REFERENCIA))
}

# --- 6.3 Run both models per comparison ---------------------------------------
dtf <- list(); tablaS1 <- NULL
for (cmp in COMPARACIONES) {
  keep  <- cmp$g %in% c(cmp$ref, cmp$foc)
  resp_m <- as.matrix(resp[keep, , drop = FALSE])
  grp   <- as.integer(cmp$g[keep] == cmp$foc)
  theta <- theta_full[keep]                       # common metric, not rescaled

  senal <- tabla2$Item[tabla2$Comparison == cmp$etiqueta & tabla2$Signal_020]
  libres <- match(senal, etiquetas)

  restringido <- if (length(libres)) modelo_dtf(resp_m, grp, theta, libres) else NULL
  libre       <- modelo_dtf(resp_m, grp, theta, seq_len(ncol(resp_m)))
  dtf[[cmp$clave]] <- list(cmp = cmp, restringido = restringido, libre = libre,
                           senal = senal,
                           n_ref = sum(grp == 0), n_foc = sum(grp == 1))

  tablaS1 <- rbind(tablaS1, data.frame(
    Comparison = cmp$etiqueta, Tier = cmp$nivel,
    Groups = paste(cmp$ref_en, "vs", cmp$foc_en),
    Items_flagged_020 = if (length(senal)) paste(senal, collapse = ", ") else "none",
    Restricted_max_theta  = if (is.null(restringido)) "not defined" else sprintf("%.3f", restringido$max_theta),
    Restricted_mean_theta = if (is.null(restringido)) "not defined" else sprintf("%+.3f", restringido$media_theta),
    Stress_max_points = sprintf("%.2f", libre$max_score),
    Stress_max_theta  = sprintf("%.3f", libre$max_theta),
    Stress_mean_theta = sprintf("%+.3f", libre$media_theta),
    Stress_pct_above  = sprintf("%.0f%%", 100 * libre$pct_sobre),
    stringsAsFactors = FALSE))

  cat(sprintf("  %-24s flagged: %-12s restricted %s | stress %.3f theta\n",
              cmp$etiqueta, if (length(senal)) paste(senal, collapse = ",") else "-",
              if (is.null(restringido)) "     n/a  " else sprintf("%.3f theta", restringido$max_theta),
              libre$max_theta))
}


# =============================================================================
# SECTION 7: FIGURES
# =============================================================================
# cairo_pdf embeds fonts, so the Greek theta survives on any viewer and in any
# locale; pdf() with plotmath is the fallback where cairo is unavailable.

# cairo_pdf embeds the Greek theta reliably, so it is preferred. But
# capabilities("cairo") can report TRUE on a machine where the cairo library
# fails to load at run time (observed on macOS), in which case cairo_pdf warns
# and may not open a device at all. Checking whether a device actually appeared
# is the only dependable test, so the fallback is driven by that rather than by
# capabilities() alone.
abrir_pdf <- function(archivo, w, h) {
  pdf_base <- function()
    grDevices::pdf(archivo, width = w, height = h, family = "serif",
                   useDingbats = FALSE, encoding = "ISOLatin1")
  if (!capabilities("cairo")) { pdf_base(); return(invisible("pdf")) }

  n_dispositivos <- length(grDevices::dev.list())
  aviso <- NULL
  abierto <- withCallingHandlers(
    tryCatch({
      grDevices::cairo_pdf(archivo, width = w, height = h, family = "serif")
      length(grDevices::dev.list()) > n_dispositivos
    }, error = function(e) { aviso <<- conditionMessage(e); FALSE }),
    warning = function(w) { aviso <<- conditionMessage(w); invokeRestart("muffleWarning") })

  if (!abierto) {
    message("cairo_pdf unavailable (", aviso, "); using the base pdf device for ",
            basename(archivo), " -- confirm that theta renders in the output.")
    pdf_base()
    return(invisible("pdf"))
  }
  if (!is.null(aviso))
    message("cairo_pdf warned while opening ", basename(archivo), " (", aviso,
            ") -- confirm that theta renders in the output.")
  invisible("cairo")
}
lab_theta <- expression(paste("Financial toxicity (", theta, ")"))
lab_sesgo <- expression(paste("Bias in ", theta, " units (focal - reference)"))

AZUL <- "#2A78D6"; NARANJA <- "#EB6834"
TINTA <- "#0B0B0B"; TINTA2 <- "#52514E"; REJILLA <- "#E3E3DF"
FUERA <- "#F4F4F1"; BANDA <- "#DFE7F2"

# --- 7.1 ICC for the items flagged at 0.020 (Online Resources 2 and 3) --------
# Exactly the four items the paper names: IF07 and IPS08 for respondent
# type/site, IF10 and IPS10 for insurance regime.
graficar_icc <- function(clave, archivo) {
  fit   <- ajustes[[clave]]
  d13   <- fit$stats[, "pseudo13.Nagelkerke"]
  cuales <- which(d13 >= UMBRAL_PURIFICACION)
  if (!length(cuales)) { message("No flagged items for ", clave); return(invisible()) }
  abrir_pdf(file.path(dir_output, archivo), 11, 5)
  for (i in cuales) {
    plot(fit, labels = etiquetas, item = i)
    title(sub = sprintf("%s | dR2 = %.3f (purification threshold %.3f; criterion %.3f)",
                        etiquetas[i], d13[i], UMBRAL_PURIFICACION, UMBRAL_MODERADO),
          cex.sub = 0.85, col.sub = "gray40")
  }
  dev.off()
  cat("  saved:", archivo, "-", paste(etiquetas[cuales], collapse = ", "), "\n")
}
cat("\n\n=== FIGURES ===\n")
graficar_icc("respondent", "ESM_2_ICC_RespondentType.pdf")
graficar_icc("insurance",  "ESM_3_ICC_InsuranceRegime.pdf")

# --- 7.2 Panel helper ---------------------------------------------------------
panel <- function(r, ylim) {
  plot(NA, xlim = c(-3, 3), ylim = ylim, xlab = "", ylab = "", axes = FALSE)
  rect(-3, ylim[1], r$libre$lo, ylim[2], col = FUERA, border = NA)
  rect(r$libre$hi, ylim[1], 3, ylim[2], col = FUERA, border = NA)
  abline(h = axTicks(2), col = REJILLA, lwd = .5)
  abline(v = -3:3, col = REJILLA, lwd = .5)
  axis(1, col = TINTA2, col.axis = TINTA2, cex.axis = .8, lwd = .6)
  axis(2, col = TINTA2, col.axis = TINTA2, cex.axis = .8, lwd = .6, las = 1)
  box(bty = "l", col = TINTA2, lwd = .6)
  mtext(r$cmp$etiqueta, side = 3, line = 1.35, adj = 0, font = 2, cex = .72, col = TINTA)
  mtext(sprintf("reference %s (n = %d)  |  focal %s (n = %d)",
                r$cmp$ref_en, r$n_ref, r$cmp$foc_en, r$n_foc),
        side = 3, line = .35, adj = 0, cex = .56, col = TINTA2)
}

# --- 7.3 Test characteristic curves (Online Resource 4) -----------------------
abrir_pdf(file.path(dir_output, "ESM_4_TCC.pdf"), 9.2, 6.1)
par(mfrow = c(2, 3), mar = c(3.2, 3.4, 3, .8), oma = c(3, 1.2, 2.2, .6),
    mgp = c(2, .6, 0), tcl = -.25)
for (r in dtf) {
  panel(r, c(0, 80))
  lines(GRID, r$libre$es_ref, col = AZUL,    lwd = 2)
  lines(GRID, r$libre$es_foc, col = NARANJA, lwd = 2, lty = 2)
  mtext(lab_theta, side = 1, line = 2, cex = .62, col = TINTA2)
  mtext("Expected TF-Col total score", side = 2, line = 2.2, cex = .62, col = TINTA2)
  text(2.95, 3, sprintf("max |DTF| %.1f pts", r$libre$max_score),
       adj = 1, cex = .68, col = TINTA2)
}
mtext("Test characteristic curves by subgroup (all-items-free sensitivity model)",
      outer = TRUE, side = 3, line = .5, adj = 0, font = 2, cex = .82, col = TINTA)
par(fig = c(0, 1, 0, 1), oma = rep(0, 4), mar = rep(0, 4), new = TRUE)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("bottom", horiz = TRUE, bty = "n", cex = .8, seg.len = 2.6,
       legend = c("Reference group", "Focal group"),
       col = c(AZUL, NARANJA), lty = 1:2, lwd = 2, text.col = TINTA)
dev.off()
cat("  saved: ESM_4_TCC.pdf\n")

# --- 7.4 DTF against measurement precision (Online Resource 5) ----------------
# y-limits contain every curve, so nothing is clipped at the tails.
yr <- max(sapply(dtf, function(r) max(abs(r$libre$d_theta)))) * 1.08
abrir_pdf(file.path(dir_output, "ESM_5_DTF.pdf"), 9.2, 6.1)
par(mfrow = c(2, 3), mar = c(3.2, 3.6, 3, .8), oma = c(4.2, 1.2, 2.2, .6),
    mgp = c(2, .6, 0), tcl = -.25)
for (r in dtf) {
  panel(r, c(-yr, yr))
  rect(-3, -REFERENCIA, 3, REFERENCIA, col = BANDA, border = NA)
  abline(h = 0, col = TINTA2, lwd = .6)
  lines(GRID, r$libre$d_theta, col = AZUL, lwd = 2)
  if (!is.null(r$restringido))
    lines(GRID, r$restringido$d_theta, col = NARANJA, lwd = 2, lty = 2)
  mtext(lab_theta, side = 1, line = 2,   cex = .62, col = TINTA2)
  mtext(lab_sesgo, side = 2, line = 2.4, cex = .62, col = TINTA2)
  text(0, -yr * .93, sprintf("stress-test mean %+.3f", r$libre$media_theta),
       adj = .5, cex = .68, col = TINTA2)
}
mtext("Differential test functioning relative to measurement precision",
      outer = TRUE, side = 3, line = .5, adj = 0, font = 2, cex = .82, col = TINTA)
par(fig = c(0, 1, 0, 1), oma = rep(0, 4), mar = rep(0, 4), new = TRUE)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("bottom", ncol = 2, bty = "n", cex = .72, seg.len = 2.4,
       legend = as.expression(c(
         "All-items-free sensitivity model",
         "Signal-restricted DTF (flagged items only)",
         bquote("Within measurement precision (" * phantom() %+-%
                .(sprintf("%.2f", REFERENCIA)) ~ theta * ")"),
         "Outside central 95% of the trait distribution")),
       col = c(AZUL, NARANJA, BANDA, FUERA), lty = c(1, 2, NA, NA),
       lwd = c(2, 2, NA, NA), pch = c(NA, NA, 15, 15),
       pt.cex = c(NA, NA, 1.8, 1.8), text.col = TINTA)
dev.off()
cat("  saved: ESM_5_DTF.pdf\n")


# =============================================================================
# SECTION 8: TABLE 1 — SAMPLE CHARACTERISTICS
# =============================================================================

cat("\n\n=== TABLE 1 ===\n")
N  <- nrow(datos)
ad <- tipo_pob == "Adulto"
cu <- tipo_pob == "Infantil"

npct <- function(m) {
  m <- m & !is.na(m)
  c(sprintf("%d (%.1f)", sum(m), 100 * sum(m) / N),
    sprintf("%d (%.1f)", sum(m & ad), 100 * sum(m & ad) / sum(ad)),
    sprintf("%d (%.1f)", sum(m & cu), 100 * sum(m & cu) / sum(cu)))
}
msd <- function(x) { x <- as.numeric(x); sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)) }
fila <- function(etq, v) data.frame(Characteristic = etq, Total = v[1],
                                    Adult_patients = v[2], Caregivers = v[3],
                                    stringsAsFactors = FALSE)

edad_pac <- as.numeric(col("6_edad_paciente"))
edad_cui <- as.numeric(col("8_edad_cuidador"))

tabla1 <- rbind(
  fila("N", c(N, sum(ad), sum(cu))),
  fila("Age, adult patients, mean (SD)",     c(msd(edad_pac[ad]), msd(edad_pac[ad]), "-")),
  fila("Age, caregivers, mean (SD)",         c(msd(edad_cui[cu]), "-", msd(edad_cui[cu]))),
  fila("Age, paediatric patients, mean (SD)",c(msd(edad_pac[cu]), "-", msd(edad_pac[cu]))),
  fila("Sex, female",       npct(sexo == "mujer")),
  fila("Sex, male",         npct(sexo == "hombre")),
  fila("Sex, not recorded", npct(is.na(sexo))),
  fila("INC, Bogota",       npct(has(ciudad, "INC"))),
  fila("HOMI, Bogota",      npct(has(ciudad, "HOMI"))),
  fila("Manizales",         npct(has(ciudad, "MANIZALES"))),
  fila("Pasto",             npct(has(ciudad, "PASTO"))),
  fila("Leukemia",          npct(has(dx_cat, "Leucemia"))),
  fila("Thyroid cancer",    npct(has(dx_cat, "tiroides"))),
  fila("Breast cancer",     npct(has(dx_cat, "de mama"))),
  fila("Lymphoma / sarcoma",npct(has(dx_cat, "Linfoma") | has(dx_cat, "Sarcoma"))),
  fila("Other tumours",     npct(!(has(dx_cat, "Leucemia") | has(dx_cat, "tiroides") |
                                   has(dx_cat, "de mama")  | has(dx_cat, "Linfoma")  |
                                   has(dx_cat, "Sarcoma")))),
  fila("Subsidized regime",     npct(regimen_bin == "Subsidiado")),
  fila("Non-subsidized regime", npct(regimen_bin == "No-subsidiado")),
  fila("Education, basic-intermediate", npct(educ_bin == "Basica-Media")),
  fila("Education, higher",             npct(educ_bin == "Superior")),
  fila("Education, not recorded",       npct(is.na(educ_bin))),
  fila("Stratum, low",         npct(estrato_bin == "Bajo")),
  fila("Stratum, middle-high", npct(estrato_bin == "Medio-Alto")),
  fila("Currently employed",     npct(trabajo == "Si")),
  fila("Not currently employed", npct(trabajo == "No")),
  fila("TF-Col total (0-80), mean (SD)",
       c(msd(col("PUNTAJE_TF")), msd(col("PUNTAJE_TF")[ad]), msd(col("PUNTAJE_TF")[cu]))),
  fila("Financial Impact (0-40), mean (SD)",
       c(msd(col("PUNTAJE_IF")), msd(col("PUNTAJE_IF")[ad]), msd(col("PUNTAJE_IF")[cu]))),
  fila("Psychosocial Impact (0-40), mean (SD)",
       c(msd(col("PUNTAJE_Ips")), msd(col("PUNTAJE_Ips")[ad]), msd(col("PUNTAJE_Ips")[cu])))
)
print(tabla1, row.names = FALSE)

# Itemisation of "Other tumours" for the Table 1 footnote. Two spellings of the
# same diagnosis coexist in the source data; merging them is what reconciles the
# footnote with the cell total.
otros <- dx_cat[!(has(dx_cat, "Leucemia") | has(dx_cat, "tiroides") |
                  has(dx_cat, "de mama")  | has(dx_cat, "Linfoma")  |
                  has(dx_cat, "Sarcoma"))]
otros <- ifelse(has(otros, "piel"), "Cancer de piel", as.character(otros))
cat("\nOther tumours itemisation (sums to", length(otros), "):\n")
print(sort(table(otros), decreasing = TRUE))


# =============================================================================
# SECTION 8.1: TRAIT DISTRIBUTION BY GROUP — IMPACT, NOT BIAS
# =============================================================================
# Groups can differ in the amount of financial toxicity they actually experience.
# That is IMPACT: a genuine difference in the construct, not a defect of the
# instrument. DIF is the separate question of whether two respondents at the
# SAME trait level answer an item differently, which is why every model in
# Sections 4-6 conditions on theta. The two are routinely confused, so the
# distribution is reported explicitly rather than left implicit.
#
# Theta is the single full-sample calibration of Section 6.1, so this table is
# on the same metric as the DTF results. The summed score is shown alongside it
# because it is the interpretable one for readers.

cat("\n\n=== TRAIT DISTRIBUTION BY GROUP (IMPACT) ===\n")

resumen_grupo <- function(x, m) {
  x <- as.numeric(x[m])
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

puntaje_tf <- as.numeric(col("PUNTAJE_TF"))

impacto <- do.call(rbind, lapply(COMPARACIONES, function(cmp) {
  filas <- lapply(c(cmp$ref, cmp$foc), function(nivel) {
    m <- !is.na(cmp$g) & cmp$g == nivel
    data.frame(
      Comparison = cmp$etiqueta,
      Tier       = cmp$nivel,
      Group      = if (nivel == cmp$ref) cmp$ref_en else cmp$foc_en,
      n          = sum(m),
      Theta      = resumen_grupo(theta_full, m),
      TF_score   = resumen_grupo(puntaje_tf, m),
      stringsAsFactors = FALSE)
  })
  salida <- do.call(rbind, filas)
  # Standardised mean difference on the trait metric, focal minus reference.
  a <- theta_full[!is.na(cmp$g) & cmp$g == cmp$ref]
  b <- theta_full[!is.na(cmp$g) & cmp$g == cmp$foc]
  s <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) /
            (length(a) + length(b) - 2))
  salida$Theta_SMD <- c("", sprintf("%+.2f", (mean(b) - mean(a)) / s))
  salida
}))

print(impacto, row.names = FALSE)
cat("\nTheta is the full-sample calibration of Section 6.1; SMD is the focal-minus-\n",
    "reference standardised difference. These are group differences in the trait,\n",
    "that is impact, and carry no implication of item bias.\n", sep = "")


# =============================================================================
# SECTION 9: EXPORT
# =============================================================================

cat("\n\n=== EXPORT ===\n")
guardar <- function(x, archivo) {
  write.csv(x, file.path(dir_output, archivo), row.names = FALSE)
  cat("  ", archivo, "\n", sep = "")
}
guardar(tabla1,    "Table1_Sample.csv")
guardar(tabla2,    "Table2_DIF.csv")
guardar(tablaS1,   "TableS1_DTF.csv")
guardar(tabla_dom, "DomainDIF_sensitivity.csv")
guardar(impacto,   "TableS2_Impact.csv")

writeLines(capture.output(sessionInfo()), file.path(dir_output, "SessionInfo.txt"))
cat("   SessionInfo.txt\n")

cat("\n", strrep("=", 62), "\n", sep = "")
cat("ANALYSIS COMPLETE — outputs in ", normalizePath(dir_output), "\n", sep = "")
cat("Items reaching the 0.035 criterion for meaningful DIF: ",
    sum(tabla2$DIF_035), "\n", sep = "")
cat("Items reaching the 0.020 purification threshold: ",
    sum(tabla2$Signal_020), " (",
    paste(unique(tabla2$Item[tabla2$Signal_020]), collapse = ", "), ")\n", sep = "")
cat(strrep("=", 62), "\n", sep = "")

