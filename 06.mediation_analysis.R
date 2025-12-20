# ============================================================
# UKB-anchored pairwise alignment (CLEAN + ROBUST READER)
# Outputs:
# 1) UKB_pairwise_reviewer_table.csv
# 2) UKB_aligned_scores_ALLtargets_except_CLHLS.csv
# 3) each pair folder: alpha0_psi0.csv, BrainVital8_aligned_long.csv, result.rds, README.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(sirt)
  library(tibble)
})

# -----------------------------
# 0) Paths + settings
# -----------------------------
paths <- c(
  UKB   = "./input/UKB_BrainVital8_domains_for_alignment.csv",
  SHARE = "./input/SHARE_BrainVital8_domains_for_alignment.csv",
  MHAS  = "./input/MHAS_BrainVital8_domains_for_alignment.csv",
  KLOSA = "./input/KLoSA_BrainVital8_domains_for_alignment.csv",
  HRS   = "./input/HRS_BrainVital8_domains_for_alignment.csv",
  ELSA  = "./input/ELSA_BrainVital8_PlanA_domains_for_alignment.csv",
  CHARLS= "./input/CHARLS_BrainVital8_scores.csv",
  NHANES= "./input/NHANES_BrainVital8_domains_for_alignment_2ndOrder.csv",
  LASI  = "./input/LASI_BrainVital8_domains_for_alignment_SEM.csv",
  CLHLS = "./input/CLHLS_BrainVital8_domains_for_alignment.csv"
)

out_dir <- "./output/UKB_pairwise_alignment_clean2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

domains <- c(
  "respiratory_comp","glucose_comp","sleep_comp","grip_comp",
  "activity_comp","SDOH_comp","emotional_comp","smoke_comp"
)

psi0_cut <- 10
scale_within_group <- TRUE
min_complete_domains <- 6

# -----------------------------
# 1) ID + BrainVital8 mappings (optional)
# -----------------------------
id_map <- c(
  UKB    = "eid",
  SHARE  = "mergeid",
  MHAS   = "rahhidnp",
  KLOSA  = "id",
  HRS    = "hhidpn",         # 
  ELSA   = "idauniq",        # 
  CHARLS = "id",
  NHANES = "seqn",
  LASI   = "ID1",
  CLHLS  = "id"
)

brain_map <- c(
  UKB    = "BrainVital8",
  SHARE  = "BrainVital8",
  MHAS   = "BrainVital8",
  KLOSA  = "BrainVital8",
  HRS    = "BrainVital8",
  ELSA   = "BrainVital8",
  CHARLS = "BrainVital8",
  NHANES = "BrainVital8",
  LASI   = "BrainVital8",
  CLHLS  = "BrainVital8"
)

# -----------------------------
# 2) Helpers: choose columns
# -----------------------------
pick_brain <- function(df) {
  nms <- names(df)
  hit <- nms[str_detect(tolower(nms), "brainvital8")]
  if (length(hit) == 0) {    # 
    hit2 <- nms[str_detect(tolower(nms), "bv8|brain_vital")]
    if (length(hit2) == 0) {
      stop("❌ Cannot find BrainVital8 column. First 80 cols: ",
           paste(head(nms, 80), collapse = ", "))
    }
    return(hit2[1])
  }
  hit[1]
}

pick_id <- function(df) {
  cands <- c("eid","mergeid","id","idauniq","pid","id_tmp","rahhidnp","seqn","id1","hhidpn")
  nms_low <- tolower(names(df))
  hit <- intersect(nms_low, tolower(cands))
  if (length(hit) == 0) return(NULL)
  names(df)[match(hit[1], nms_low)]
}

# -----------------------------
# 3) Read one cohort (robust)
# -----------------------------
read_one <- function(cohort, path, domains) {
  df <- read_csv(path, show_col_types = FALSE)
  
  # find id
  idc <- if (cohort %in% names(id_map) && id_map[[cohort]] %in% names(df)) id_map[[cohort]] else pick_id(df)
  if (is.null(idc) || !(idc %in% names(df))) {
    # No ID allowed
    idc <- NULL
  }
  
  # find BrainVital8
  bc <- if (cohort %in% names(brain_map) && brain_map[[cohort]] %in% names(df)) brain_map[[cohort]] else pick_brain(df)
  if (!(bc %in% names(df))) stop("❌ ", cohort, " cannot find BrainVital8 column.")
  
  # domains check
  miss <- setdiff(domains, names(df))
  if (length(miss) > 0) {
    stop("❌ ", cohort, " missing domain cols: ", paste(miss, collapse = ", "),
         "\nFirst 80 cols: ", paste(head(names(df), 80), collapse = ", "))
  }
  
  out <- df %>%
    transmute(
      cohort = cohort,
      id = if (!is.null(idc)) as.character(.data[[idc]]) else as.character(row_number()),
      across(all_of(domains), ~ as.numeric(.x)),
      BrainVital8_raw = as.numeric(.data[[bc]])
    )
  
  # drop rows with too many missing domains
  out <- out %>%
    mutate(n_non_na = rowSums(!is.na(across(all_of(domains))))) %>%
    filter(n_non_na >= min_complete_domains) %>%
    select(-n_non_na)
  
  out
}

# -----------------------------
# 4) Within-cohort z-score standardization
# -----------------------------
scale_domains_within <- function(dat, domain_vec) {
  dat %>%
    group_by(cohort) %>%
    mutate(across(all_of(domain_vec), ~{
      x <- as.numeric(.x)
      s <- sd(x, na.rm = TRUE)
      if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
      (x - mean(x, na.rm = TRUE)) / s
    })) %>%
    ungroup()
}

# -----------------------------
# 5) Extract effect sizes (R2 / sqrtU2) from summary text (stable)
# -----------------------------
extract_effects_from_summary <- function(aln) {
  txt <- capture.output(suppressWarnings(try(print(summary(aln)), silent = TRUE)))
  i0 <- which(str_detect(txt, "^Effect Sizes of Approximate Invariance"))
  if (length(i0) == 0) {
    return(tibble(
      R2_loadings = NA_real_, R2_intercepts = NA_real_,
      sqrtU2_loadings = NA_real_, sqrtU2_intercepts = NA_real_
    ))
  }
  i0 <- i0[1]
  sub <- txt[(i0+1):min(length(txt), i0+12)]
  
  parse_row <- function(label) {
    line <- sub[str_detect(sub, paste0("^", label, "\\b"))]
    if (length(line) == 0) return(c(NA_real_, NA_real_))
    line <- line[1]
    nums <- str_extract_all(line, "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?")[[1]]
    # loadings & intercepts
    if (length(nums) < 2) return(c(NA_real_, NA_real_))
    as.numeric(nums[1:2])
  }
  
  r2 <- parse_row("R2")
  su <- parse_row("sqrtU2")
  
  tibble(
    R2_loadings = r2[1],
    R2_intercepts = r2[2],
    sqrtU2_loadings = su[1],
    sqrtU2_intercepts = su[2]
  )
}

# -----------------------------
# 6) Run one pair: UKB vs target
# -----------------------------
run_pair <- function(dat_all, target, out_dir, domains,
                     psi0_cut = 10, scale_within_group = TRUE, verbose = FALSE) {
  
  pair_dir <- file.path(out_dir, paste0("UKB_vs_", target))
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  
  dat_pair <- dat_all %>% filter(cohort %in% c("UKB", target))
  if (nrow(dat_pair) == 0) {
    return(list(
      summary = tibble(
        target = target, ok = FALSE, degenerate = NA,
        psi0_min = NA_real_, psi0_max = NA_real_,
        R2_loadings = NA_real_, R2_intercepts = NA_real_,
        sqrtU2_loadings = NA_real_, sqrtU2_intercepts = NA_real_,
        msg = "NO DATA"
      ),
      aligned = NULL
    ))
  }
  
  if (scale_within_group) dat_pair <- scale_domains_within(dat_pair, domains)
  
  X <- dat_pair %>% select(all_of(domains))
  grp <- dat_pair$cohort
  
  ok <- TRUE
  msg <- "OK"
  cfg <- NULL
  aln <- NULL
  
  tmp <- try({
    cfg <- sirt::invariance_alignment_cfa_config(dat = X, group = grp, model = "2PM", verbose = verbose)
    aln <- sirt::invariance.alignment(lambda = cfg$lambda, nu = cfg$nu, optimizer = "nlminb")
  }, silent = TRUE)
  
  if (inherits(tmp, "try-error") || is.null(aln)) {
    ok <- FALSE
    msg <- paste0("FAILED: ", as.character(tmp))
  }
  
  if (!ok) {
    return(list(
      summary = tibble(
        target = target, ok = FALSE, degenerate = NA,
        psi0_min = NA_real_, psi0_max = NA_real_,
        R2_loadings = NA_real_, R2_intercepts = NA_real_,
        sqrtU2_loadings = NA_real_, sqrtU2_intercepts = NA_real_,
        msg = msg
      ),
      aligned = NULL
    ))
  }
  
  # psi0 range + degeneracy
  psi <- aln$pars[, "psi0"]
  psi0_min <- suppressWarnings(min(psi, na.rm = TRUE))
  psi0_max <- suppressWarnings(max(psi, na.rm = TRUE))
  degenerate <- isTRUE(is.finite(psi0_max) && psi0_max > psi0_cut)
  
  # effect sizes from summary
  eff <- extract_effects_from_summary(aln)
  
  # alpha0/psi0 table
  link_tbl <- as.data.frame(aln$pars) %>%
    rownames_to_column("cohort") %>%
    select(cohort, alpha0, psi0)
  
  # aligned score
  dat_aligned <- dat_pair %>%
    left_join(link_tbl, by = "cohort") %>%
    mutate(BrainVital8_aligned = alpha0 + psi0 * BrainVital8_raw)
  
  # save outputs
  write_csv(link_tbl, file.path(pair_dir, "alpha0_psi0.csv"))
  write_csv(dat_aligned, file.path(pair_dir, "BrainVital8_aligned_long.csv"))
  saveRDS(list(cfg = cfg, aln = aln, link_tbl = link_tbl, dat_aligned = dat_aligned),
          file.path(pair_dir, "result.rds"))
  
  cat(
    "Pairwise alignment: UKB vs ", target, "\n",
    "scale_within_group: ", scale_within_group, "\n",
    "psi0_cut(degenerate): ", psi0_cut, "\n",
    "psi0_min: ", psi0_min, "\n",
    "psi0_max: ", psi0_max, "\n",
    "degenerate: ", degenerate, "\n",
    sep = "",
    file = file.path(pair_dir, "README.txt")
  )
  
  sum_row <- tibble(
    target = target,
    ok = TRUE,
    degenerate = degenerate,
    psi0_min = psi0_min,
    psi0_max = psi0_max,
    R2_loadings = eff$R2_loadings,
    R2_intercepts = eff$R2_intercepts,
    sqrtU2_loadings = eff$sqrtU2_loadings,
    sqrtU2_intercepts = eff$sqrtU2_intercepts,
    msg = ifelse(degenerate, paste0("DEGENERATE: psi0_max>", psi0_cut), "OK")
  )
  
  list(summary = sum_row, aligned = dat_aligned)
}

# -----------------------------
# 7) RUN ALL
# -----------------------------
message("Reading all cohorts...")
dat_all <- imap_dfr(paths, ~ read_one(.y, .x, domains))

stopifnot("UKB" %in% unique(dat_all$cohort))
targets <- setdiff(sort(unique(dat_all$cohort)), "UKB")

message("Targets: ", paste(targets, collapse = ", "))

res_list <- setNames(vector("list", length(targets)), targets)
for (tg in targets) {
  message("\n--- UKB vs ", tg, " ---")
  res_list[[tg]] <- run_pair(dat_all, tg, out_dir, domains,
                             psi0_cut = psi0_cut,
                             scale_within_group = scale_within_group,
                             verbose = FALSE)
}

# reviewer table
review_tbl <- map_dfr(res_list, "summary") %>%
  arrange(degenerate, desc(psi0_max))

write_csv(review_tbl, file.path(out_dir, "UKB_pairwise_reviewer_table.csv"))
print(review_tbl)

# aligned score export (exclude CLHLS per your requirement)
aligned_all <- map_dfr(res_list, function(x) x$aligned)
aligned_no_clhls <- aligned_all %>% filter(cohort != "CLHLS")

write_csv(aligned_no_clhls, file.path(out_dir, "UKB_aligned_scores_ALLtargets_except_CLHLS.csv"))

message("\n✅ DONE. Outputs in: ", out_dir)
