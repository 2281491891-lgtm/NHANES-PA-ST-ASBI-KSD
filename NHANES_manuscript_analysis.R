# NHANES manuscript and supplementary analyses
# Input: the already cleaned, fixed NHANES analysis data frame.
# No raw-data reading, sample cleaning, automatic analysis, or figure saving.

FONT_FAMILY <- "Arial"
RCS_ANNOTATION_SIZE <- getOption("nhanes.rcs.annotation_size", 18)
FIGURE_LABEL_SIZE <- getOption("nhanes.figure.label_size", 18)
ST_480_GROUP <- "S4"
ASBI_CONSTANT <- .1
RCS_MODEL <- "M3"
RUN_THRESHOLD <- TRUE
EXPECTED_N <- 19977L
EXPECTED_CASES <- 1713L
DEMO_COVARS <- c("sex", "age_group", "race", "education", "marital", "pir")

FULL_COVARS <- c(DEMO_COVARS, "diabetes", "hyperlipid", "hypertension", "bmi",
                 "waist", "drink", "smoke", "water_g", "uric_g")

KNOT_PROBS <- c(.05,.35,.65,.95)

DISPLAY_PROBS <- c(.01,.99)

SPLINE_NAMES <- c("rcs_linear","rcs_nonlinear1","rcs_nonlinear2")

# 1. Input validation and common functions

validate_analysis_input <- function(d) {
  continuous <- c("SEQN","age","bmi","waist","water","uric_acid","PAQ","sedentary",
    "WTMEC2YR","WTMEC12YR","SDMVPSU","SDMVSTRA","PAQ_work_vig","PAQ_work_mod",
    "PAQ_transport","PAQ_leisure_vig","PAQ_leisure_mod")
  categorical <- setdiff(FULL_COVARS,c("bmi","waist"))
  required <- c(continuous,categorical,"kidney_stone")
  missing <- setdiff(required,names(d))
  if(length(missing)) stop("Missing cleaned variables: ",paste(missing,collapse=", "))
  stopifnot(is.data.frame(d),nrow(d)==EXPECTED_N,!anyDuplicated(d$SEQN),
    !anyNA(d[required]),all(vapply(d[continuous],is.numeric,logical(1))),
    all(vapply(d[continuous],function(x)all(is.finite(x)),logical(1))),
    all(vapply(d[categorical],is.factor,logical(1))),
    all(vapply(d[categorical],function(x)!is.ordered(x),logical(1))))
  y <- as.integer(as.character(d$kidney_stone))
  stopifnot(all(y %in% 0:1),sum(y)==EXPECTED_CASES,
    all(d$WTMEC12YR>0),all(d$PAQ>=0),all(d$sedentary>=0),
    max(abs(d$WTMEC12YR-d$WTMEC2YR/6))<1e-12,
    max(abs(d$PAQ-rowSums(d[c("PAQ_work_vig","PAQ_work_mod","PAQ_transport",
      "PAQ_leisure_vig","PAQ_leisure_mod")])/7))<1e-8)
  invisible(TRUE)
}

check_packages <- function() {
  needed <- c("dplyr", "survey", "Hmisc", "ggplot2")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) stop("Please install: ", paste(missing, collapse=", "))
  invisible(TRUE)
}

rcs_weighted_median <- function(x, d, keep = rep(TRUE, nrow(d))) {
  stopifnot(length(x) == nrow(d), length(keep) == nrow(d), !anyNA(keep),
    all(is.finite(x)), all(is.finite(d$WTMEC12YR)), all(d$WTMEC12YR > 0), any(keep))
  qdata <- data.frame(SDMVPSU=d$SDMVPSU, SDMVSTRA=d$SDMVSTRA,
    WTMEC12YR=d$WTMEC12YR, rcs_reference_value=x)
  qdesign <- survey::svydesign(ids=~SDMVPSU, strata=~SDMVSTRA,
    weights=~WTMEC12YR, nest=TRUE, data=qdata)
  qdesign <- qdesign[keep, ]
  as.numeric(coef(survey::svyquantile(~rcs_reference_value, qdesign,
    quantiles=.5, ci=FALSE, qrule="math")))
}

bind_rows_base <- function(x) {
  if (!length(x)) return(data.frame())
  out <- do.call(rbind, x); rownames(out) <- NULL; out
}

p_text <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))

p_annotate <- function(p) ifelse(p < .001, "< 0.001", paste0("= ", sprintf("%.3f", p)))

formula_text <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = " ")

asbi_value <- function(d, constant) log((d$PAQ / 60 + constant) / (d$sedentary / 60 + constant))

model_covars <- function(model, exposure) {
  covars <- switch(model, M1 = character(), M2 = DEMO_COVARS, M3 = FULL_COVARS,
                   stop("未知模型: ", model))
  if (model == "M3" && exposure == "PA") covars <- c(covars, "ST_group")
  if (model == "M3" && exposure == "ST") covars <- c(covars, "PA_group")
  covars
}

make_design <- function(d) survey::svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~WTMEC12YR, nest = TRUE, data = d)

fit_survey <- function(f, design) {
  fit <- survey::svyglm(f, design = design, family = quasibinomial(), na.action = na.fail,
    control = glm.control(epsilon = 1e-10, maxit = 100), model = TRUE)
  stopifnot(isTRUE(fit$converged), nobs(fit) == nrow(design$variables),
    identical(rownames(fit$model), rownames(design$variables)),
    fit$rank == length(coef(fit)), !anyNA(coef(fit)), all(is.finite(vcov(fit))), fit$df.residual > 0)
  fit
}

wald_test <- function(fit, term, model, exposure, test_label = "Overall") {
  tt <- survey::regTermTest(fit, reformulate(term), method = "Wald")
  data.frame(model = model, exposure = exposure, test = test_label,
    F = as.numeric(tt$Ftest), numerator_df = as.numeric(tt$df),
    denominator_df = as.numeric(tt$ddf), P = as.numeric(tt$p))
}

category_descriptives <- function(d, design, variable, categories) {
  bind_rows_base(lapply(categories, function(g) {
    ii <- as.character(d[[variable]]) == g
    ds <- design[ii, ]
    cp <- survey::svyciprop(~kidney_stone, ds, method = "logit", df = survey::degf(ds))
    ci <- as.numeric(confint(cp))
    data.frame(category = g, n = sum(ii), cases = sum(d$kidney_stone[ii]),
      weighted_population_pct = 100 * sum(d$WTMEC12YR[ii]) / sum(d$WTMEC12YR),
      crude_prevalence_pct = 100 * as.numeric(coef(cp)),
      crude_low_pct = 100 * ci[1], crude_high_pct = 100 * ci[2])
  }))
}

category_estimates <- function(fit, variable, categories, reference, model, exposure, role) {
  ref <- categories == reference
  idx <- match(paste0(variable, categories), names(coef(fit)))
  stopifnot(identical(is.na(idx), ref))
  b <- unname(coef(fit)[idx]); se <- unname(sqrt(diag(vcov(fit)))[idx])
  b[ref] <- 0
  critical <- qt(.975, fit$df.residual)
  ans <- data.frame(model = model, exposure = exposure, role = role, category = categories,
    reference = reference, log_OR = b, SE = se, OR = exp(b),
    OR_low = exp(b - critical * se), OR_high = exp(b + critical * se),
    P = 2 * pt(-abs(b / se), df = fit$df.residual), denominator_df = fit$df.residual)
  ans$OR_95CI <- ifelse(ref, "1.000 (Reference)",
                       sprintf("%.3f (%.3f-%.3f)", ans$OR, ans$OR_low, ans$OR_high))
  ans$P_display <- p_text(ans$P)
  ans
}

baseline_table <- function(d, design) {
  continuous <- c("age", "bmi", "waist", "water", "uric_acid", "PAQ", "sedentary", "ASBI")
  categorical <- unique(c(setdiff(FULL_COVARS, c("bmi", "waist")), "PA_group", "ST_group", "ASBI_group"))
  ans <- list()
  for (group in c("All", "No_history", "History")) {
    ii <- if (group == "All") rep(TRUE, nrow(d)) else d$kidney_stone == as.integer(group == "History")
    dd <- d[ii, , drop = FALSE]; ds <- design[ii, ]
    for (v in continuous) {
      q <- as.numeric(coef(survey::svyquantile(reformulate(v), ds, quantiles = c(.25, .5, .75), ci = FALSE)))
      ans[[paste(group, v)]] <- data.frame(group = group, variable = v, level = "",
        n = nrow(dd), weighted_pct = NA_real_, median = q[2], q25 = q[1], q75 = q[3],
        display = sprintf("%.2f (%.2f-%.2f)", q[2], q[1], q[3]))
    }
    for (v in categorical) for (g in levels(dd[[v]])) {
      n <- sum(dd[[v]] == g); pct <- 100 * sum(dd$WTMEC12YR[dd[[v]] == g]) / sum(dd$WTMEC12YR)
      ans[[paste(group, v, g)]] <- data.frame(group = group, variable = v, level = g,
        n = n, weighted_pct = pct, median = NA_real_, q25 = NA_real_, q75 = NA_real_,
        display = sprintf("%d (%.2f%%)", n, pct))
    }
  }
  bind_rows_base(ans)
}

prepare_main_data <- function(data) {
  stopifnot(ST_480_GROUP %in% c("S3", "S4"), RCS_MODEL %in% c("M1", "M2", "M3"),
            identical(ASBI_CONSTANT, .1))
  d <- as.data.frame(data)
  d$kidney_stone <- as.integer(as.character(d$kidney_stone))
  d$PA_group <- cut(d$PAQ, c(-Inf, 20, 200, 600, Inf), right = FALSE, labels = paste0("P", 1:4))
  st <- ifelse(d$sedentary < 240, 1L, ifelse(d$sedentary < 360, 2L,
    ifelse(if (ST_480_GROUP == "S3") d$sedentary <= 480 else d$sedentary < 480, 3L, 4L)))
  d$ST_group <- factor(paste0("S", st), levels = paste0("S", 1:4))
  # 同步原变量名，避免内存中出现两套矛盾的480分钟分组。
  d$sedentary_g <- factor(st, levels = 1:4)
  # PA/60为MET-h/day，ST/60为h/day；ASBI不是单纯活动时长/久坐时长之比。
  d$ASBI <- asbi_value(d, ASBI_CONSTANT)
  raw_cuts <- as.numeric(coef(survey::svyquantile(~ASBI, make_design(d),
    quantiles=c(.25,.5,.75),ci=FALSE,qrule="math")))
  cuts <- round(raw_cuts,2)
  stopifnot(all(diff(cuts)>0),identical(cuts,c(-2.14,-.47,.91)))
  d$ASBI_group <- relevel(cut(d$ASBI,c(-Inf,cuts,Inf),right=FALSE,
    labels=paste0("Q",1:4)),ref="Q3")
  if("ASBI_g" %in% names(d)) d$ASBI_g <- factor(sub("Q","",as.character(d$ASBI_group)),levels=1:4)
  d$exclude_P1S1 <- d$PA_group=="P1" & d$ST_group=="S1"
  joint_levels <- unlist(lapply(paste0("P", 1:4), function(p) paste(p, paste0("S", 1:4), sep = "_")))
  d$joint_group <- relevel(factor(paste(d$PA_group, d$ST_group, sep = "_"),
                                 levels = joint_levels), ref = "P4_S1")
  required <- c("SEQN", "kidney_stone", "PAQ", "sedentary", "ASBI", "PA_group", "ST_group",
                "ASBI_group", "joint_group", "WTMEC12YR", "SDMVPSU", "SDMVSTRA", FULL_COVARS)
  stopifnot(!anyNA(d[required]), !anyDuplicated(d$SEQN),
    nrow(d) == EXPECTED_N, sum(d$kidney_stone) == EXPECTED_CASES,
    all(d$kidney_stone %in% 0:1), all(is.finite(d$ASBI)),
    all(is.finite(d$PAQ) & d$PAQ >= 0), all(is.finite(d$sedentary) & d$sedentary >= 0),
    all(is.finite(d$WTMEC12YR) & d$WTMEC12YR > 0),
    max(abs(d$WTMEC12YR - d$WTMEC2YR / 6)) < 1e-12,
    all(as.character(d$ST_group[d$sedentary == 480]) == ST_480_GROUP),
    all(as.character(d$sedentary_g) == sub("S", "", as.character(d$ST_group))),
    all(table(d$joint_group) > 0))
  # 保留同一批参与者；不因某个模型较少协变量而扩大样本。
  rownames(d) <- as.character(d$SEQN)
  d
}

joint_category_label <- function(x) {
  gsub("^(P[1-4])\\s*[_/&]\\s*(S[1-4])$", "\\1&\\2", as.character(x), perl = TRUE)
}

joint_display_table <- function(tab) {
  for (nm in intersect(c("category", "reference"), names(tab))) {
    if (any(grepl("^P[1-4]\\s*[_/&]\\s*S[1-4]$", as.character(tab[[nm]]), perl = TRUE))) {
      tab[[nm]] <- joint_category_label(tab[[nm]])
    }
  }
  tab
}

analysis_notes <- function() data.frame(item=c("Sample","Survey","Models","ASBI","ASBI categories",
 "RCS","Joint","Threshold","Mediation","Output"),details=c(
 "NHANES 2007-2018. Fixed complete-coded sample, 19977 participants, 1713 KSD histories.",
 "SDMVSTRA, SDMVPSU, nest=TRUE, WTMEC2YR/6. Domain subsets retain the full survey design.",
 paste("M1 unadjusted. M2:",paste(DEMO_COVARS,collapse=", "),". M3:",paste(FULL_COVARS,collapse=", "),". PA/ST mutually adjust the other category. ASBI does not adjust PA/ST."),
 "ln[(PA/60+c)/(ST/60+c)], c=0.1, PA numerical MET-h/day, ST numerical h/day. ASBI is a score without a physical unit. Offsets are defined on the specified numerical scales.",
 "Survey-weighted 25th/50th/75th quantiles, qrule=math, rounded to 2 decimals BEFORE grouping: -2.14, -0.47, 0.91. Left-closed/right-open. Q3 reference.",
 "Survey weighted. Original unweighted 5/35/65/95 knots and 1-99% plotting range. Weighted-median OR reference. Original annotations retained, except boundary minima are not annotated. Pointwise contrast CIs use the full covariance.",
  "16 joint categories of PA and ST formed by crossing P1-P4 with S1-S4. Labels P1&S1 through P4&S4. P4&S1 reference. The heatmap shows nominal P values; the forest shows BH q values across all 15 comparisons. ORs and 95% CIs are covariate-adjusted model estimates; CIs are pointwise, not multiplicity-adjusted. The reference has no P or q value.",
 "Piecewise slopes compared with a design-based Wald test. Conditional on a minimum selected in this sample. Boundary minima are flagged/skipped.",
 "Continuous total PA, six separate mediator models. Fixed +100 MET-min/day contrast from weighted median. Survey JKn variance. Signed indirect/total percentages may be negative and are not causal prevented fractions.",
 "Every figure is a standalone ggplot. No automatic image or data export by default. export_manuscript_tables() saves numeric tables only. No Word or combined figure generation."))

# 2. Main analyses

fit_threshold_effect <- function(dat, x, covariates, threshold, exposure, method,
                                 increment, units, source, constant = NA_real_) {
  # A fixed split conditional on the selected RCS minimum; this is not an
  # independent change-point test and does not propagate threshold-selection uncertainty.
  stopifnot(length(x) == nrow(dat), all(is.finite(x)), is.finite(threshold), increment > 0,
            threshold > min(x), threshold < max(x), !anyNA(dat[c("kidney_stone", covariates)]))
  d <- dat
  d$threshold_x <- (x - threshold) / increment
  d$threshold_below <- pmin(d$threshold_x, 0)
  d$threshold_above <- pmax(d$threshold_x, 0)
  f_linear <- reformulate(c("threshold_x", covariates), "kidney_stone")
  f_piece <- reformulate(c("threshold_below", "threshold_above", covariates), "kidney_stone")
  f_hinge <- reformulate(c("threshold_x", "threshold_above", covariates), "kidney_stone")
  weighted <- identical(method, "Survey-weighted")
  if (weighted) {
    design <- survey::svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                                 weights = ~WTMEC12YR, nest = TRUE, data = d)
    fit_model <- function(f) survey::svyglm(f, design = design, family = quasibinomial(),
                                            control = glm.control(epsilon = 1e-10, maxit = 100))
  } else {
    fit_model <- function(f) glm(f, data = d, family = binomial(),
                                  control = glm.control(epsilon = 1e-10, maxit = 100))
  }
  linear <- fit_model(f_linear); piece <- fit_model(f_piece); hinge <- fit_model(f_hinge)
  for (fit in list(linear, piece, hinge)) {
    stopifnot(isTRUE(fit$converged), nobs(fit) == nrow(d),
              !anyNA(coef(fit)), fit$rank == length(coef(fit)))
  }
  extract <- function(fit, term, label) {
    b <- unname(coef(fit)[term]); se <- sqrt(vcov(fit)[term, term])
    df <- if (weighted) fit$df.residual else Inf
    critical <- if (weighted) qt(0.975, df) else qnorm(0.975)
    p <- if (weighted) 2 * pt(abs(b / se), df, lower.tail = FALSE) else 2 * pnorm(abs(b / se), lower.tail = FALSE)
    data.frame(segment = label, log_OR = b, SE = se, OR = exp(b),
      OR_low = exp(b - critical * se), OR_high = exp(b + critical * se), P = p, residual_df = df)
  }
  slope_all <- extract(linear, "threshold_x", "Overall linear")
  slope_below <- extract(piece, "threshold_below", "At or below threshold")
  slope_above <- extract(piece, "threshold_above", "Above threshold")
  # For separate left/right slope terms, the correct null is beta_right=beta_left.
  # Testing beta_right=0 (as in original sections 4.2.3 and 6.1.7) is a different null.
  contrast <- survey::svycontrast(piece, c(threshold_below = -1, threshold_above = 1))
  difference <- as.numeric(coef(contrast)); se_difference <- as.numeric(survey::SE(contrast))
  statistic <- (difference / se_difference)^2
  ddf <- if (weighted) piece$df.residual else Inf
  p_diff <- if (weighted) pf(statistic, 1, ddf, lower.tail = FALSE) else pchisq(statistic, 1, lower.tail = FALSE)
  hinge_test <- survey::regTermTest(hinge, ~threshold_above, df = ddf, method = "Wald")
  prediction_error <- max(abs(as.numeric(predict(piece, type = "link")) - as.numeric(predict(hinge, type = "link"))))
  contrast_error <- abs(difference - unname(coef(hinge)["threshold_above"]))
  variance_error <- abs(se_difference^2 - vcov(hinge)["threshold_above", "threshold_above"])
  p_error <- abs(p_diff - as.numeric(hinge_test$p))
  stopifnot(all(is.finite(c(difference, se_difference, p_diff))), se_difference > 0,
            prediction_error < 1e-7, contrast_error < 1e-7, variance_error < 1e-7, p_error < 1e-7)
  meta <- data.frame(exposure = exposure, constant = constant, method = method,
                       threshold = threshold, threshold_source = source, increment = increment, units = units)
  sm <- cbind(meta, data.frame(n = nrow(d), cases = sum(d$kidney_stone),
    n_below = sum(x <= threshold), cases_below = sum(d$kidney_stone[x <= threshold]),
    n_above = sum(x > threshold), cases_above = sum(d$kidney_stone[x > threshold]),
    OR_all = slope_all$OR, CI_all_low = slope_all$OR_low, CI_all_high = slope_all$OR_high, P_all = slope_all$P,
    OR_below = slope_below$OR, CI_below_low = slope_below$OR_low, CI_below_high = slope_below$OR_high, P_below = slope_below$P,
    OR_above = slope_above$OR, CI_above_low = slope_above$OR_low, CI_above_high = slope_above$OR_high, P_above = slope_above$P,
    slope_difference_above_minus_below = difference, SE_difference = se_difference,
    Wald_statistic = statistic, numerator_df = 1, denominator_df = ddf, P_difference = p_diff,
    test = if (weighted) "Survey-adjusted Wald F" else "Wald chi-square",
    inference = "Exploratory; conditional on a same-sample selected threshold; no selection-adjusted inference",
    piecewise_formula = paste(deparse(f_piece), collapse = " "), row.names = NULL))
  slopes <- cbind(meta[rep(1L, 3L), ], rbind(slope_all, slope_below, slope_above)); rownames(slopes) <- NULL
  cf <- do.call(rbind, lapply(c("linear", "piecewise"), function(label) {
    fit <- if (label == "linear") linear else piece
    cbind(meta[rep(1L, length(coef(fit))), ], data.frame(model = label,
      term = names(coef(fit)), beta = unname(coef(fit)), SE = unname(sqrt(diag(vcov(fit)))), row.names = NULL))
  })); rownames(cf) <- NULL
  ck <- cbind(meta, data.frame(n_linear = nobs(linear), n_piecewise = nobs(piece),
    piecewise_converged = piece$converged, piecewise_full_rank = piece$rank == length(coef(piece)),
    equivalent_hinge_prediction_error = prediction_error, slope_contrast_error = contrast_error,
    variance_error = variance_error, P_difference_error = p_error, passed = TRUE))
  list(summary = sm, slopes = slopes, coefficients = cf, validation = ck, models = list(linear = linear, piecewise = piece, hinge = hinge))
}

fit_main_rcs <- function(d, exposure, model = RCS_MODEL, constant = .1) {
  variable <- switch(exposure, PA = "PAQ", ST = "sedentary", ASBI = "ASBI")
  units <- switch(exposure, PA = "MET-min/day", ST = "min/day", ASBI = "ASBI score")
  if (exposure == "ASBI") d$ASBI <- asbi_value(d, constant)
  x <- d[[variable]]
  ref <- rcs_weighted_median(x, d)
  knots <- as.numeric(quantile(x, KNOT_PROBS, type = 7))
  xr <- as.numeric(quantile(x, DISPLAY_PROBS, type = 7))
  stopifnot(length(unique(knots)) == 4L, all(diff(knots) > 0), xr[1] < xr[2])
  basis <- function(xx) {
    b <- Hmisc::rcspline.eval(xx, knots = knots, inclx = TRUE, norm = 2)
    colnames(b) <- SPLINE_NAMES; b
  }
  d[SPLINE_NAMES] <- as.data.frame(basis(x))
  covars <- model_covars(model, exposure)
  fit <- fit_survey(reformulate(c(SPLINE_NAMES, covars), "kidney_stone"), make_design(d))
  beta <- coef(fit)[SPLINE_NAMES]; V <- vcov(fit)[SPLINE_NAMES, SPLINE_NAMES, drop = FALSE]
  critical <- qt(.975, fit$df.residual)
  # 仅显示1%-99%区间，拟合始终使用全样本；不裁剪观测值。
  grid <- sort(unique(c(seq(xr[1], xr[2], length.out = 800L), ref)))
  delta <- sweep(basis(grid), 2L, as.numeric(basis(ref)), FUN = "-")
  eta <- drop(delta %*% beta)
  variance <- rowSums((delta %*% V) * delta)
  stopifnot(all(variance > -1e-10))
  se <- sqrt(pmax(variance, 0))
  curve <- data.frame(value = grid, log_OR = eta, SE = se, OR = exp(eta),
    OR_low = exp(eta - critical * se), OR_high = exp(eta + critical * se))
  at_ref <- curve[curve$value == ref, ]
  stopifnot(nrow(at_ref) == 1L, max(abs(unlist(at_ref[c("OR", "OR_low", "OR_high")]) - 1)) < 1e-10,
    all(is.finite(as.matrix(curve))), all(curve$OR_low <= curve$OR & curve$OR <= curve$OR_high))
  # 用完整设计矩阵独立核对曲线的对比估计与协方差。
  nd <- d[rep(1L, length(grid)), , drop = FALSE]; nd[SPLINE_NAMES] <- as.data.frame(basis(grid))
  nr <- d[1L, , drop = FALSE]; nr[SPLINE_NAMES] <- as.data.frame(basis(ref))
  X <- model.matrix(delete.response(terms(fit)), nd, contrasts.arg = fit$contrasts)
  Xref <- model.matrix(delete.response(terms(fit)), nr, contrasts.arg = fit$contrasts)
  D <- sweep(X, 2L, as.numeric(Xref), FUN = "-")
  prediction_error <- max(abs(drop(D %*% coef(fit)) - eta))
  covariance_error <- max(abs(rowSums((D %*% vcov(fit)) * D) - variance))
  stopifnot(prediction_error < 1e-8, covariance_error < 1e-8)
  imin <- which.min(curve$OR)
  overall <- wald_test(fit, SPLINE_NAMES, model, exposure)
  nonlinear <- wald_test(fit, SPLINE_NAMES[-1], model, exposure, "Nonlinear")
  sm <- data.frame(exposure = exposure, model = model,
    PA_adjustment = if (exposure != "ASBI") "Not_applicable" else "Without_PA",
    constant = if (exposure == "ASBI") constant else NA_real_, n = nrow(d), cases = sum(d$kidney_stone),
    median_reference = ref, minimum_value = curve$value[imin], minimum_OR = curve$OR[imin],
    minimum_at_boundary = imin %in% c(1L, nrow(curve)),
    P_overall = overall$P, P_nonlinear = nonlinear$P,
    df_overall = overall$numerator_df, df_nonlinear = nonlinear$numerator_df,
    denominator_df = overall$denominator_df, knot_5 = knots[1], knot_35 = knots[2],
    knot_65 = knots[3], knot_95 = knots[4], display_lower = xr[1], display_upper = xr[2],
    units = units, formula = formula_text(formula(fit)),
    direct_prediction_error = prediction_error, full_covariance_error = covariance_error,
    reference_OR = at_ref$OR, reference_low = at_ref$OR_low, reference_high = at_ref$OR_high)
  list(fit = fit, summary = sm, curve = curve, covars = covars)
}

run_nhanes_main <- function(loaded) {
  check_packages()
  oldopt <- options(survey.lonely.psu = "adjust", contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(oldopt), add = TRUE)
  message("正在核验固定主分析样本...")
  stopifnot(is.list(loaded), !is.null(loaded$data))
  d <- prepare_main_data(loaded$data)
  design <- make_design(d)
  message("分析样本: ", nrow(d), " 人，结石病史 ", sum(d$kidney_stone), " 人。")
  models <- model_meta <- category_rows <- category_tests <- descriptions <- list()
  register_model <- function(key, fit, purpose) {
    models[[key]] <<- fit
    model_meta[[key]] <<- data.frame(id = key, purpose = purpose, n = nobs(fit),
      cases = sum(d$kidney_stone), model_rank = fit$rank, residual_df = fit$df.residual,
      converged = fit$converged, same_participants = identical(rownames(fit$model), rownames(d)),
      formula = formula_text(formula(fit)))
    invisible(NULL)
  }
  definitions <- list(
    PA = list(variable = "PA_group", categories = paste0("P", 1:4), reference = "P1"),
    ST = list(variable = "ST_group", categories = paste0("S", 1:4), reference = "S1"),
    ASBI = list(variable = "ASBI_group", categories = paste0("Q", 1:4), reference = "Q3"))
  message("正在拟合PA、ST、ASBI的分类回归...")
  for (exposure in names(definitions)) {
    z <- definitions[[exposure]]
    desc <- category_descriptives(d, design, z$variable, z$categories)
    descriptions[[exposure]] <- cbind(exposure = exposure, desc)
    for (model in c("M1", "M2", "M3")) {
      key <- paste(exposure, model, sep = "_")
      f <- reformulate(c(z$variable, model_covars(model, exposure)), "kidney_stone")
      fit <- fit_survey(f, design)
      register_model(key, fit, "Main categorical association")
      tab <- category_estimates(fit, z$variable, z$categories, z$reference, model, exposure, "Main")
      category_rows[[key]] <- dplyr::left_join(tab, desc, by = "category")
      category_tests[[key]] <- wald_test(fit, z$variable, model, exposure)
    }
  }
  categories <- bind_rows_base(category_rows)

  message("正在拟合16组联合模型及正式交互检验...")
  canonical <- unlist(lapply(paste0("P", 1:4), function(p) paste(p, paste0("S", 1:4), sep = "_")))
  joint_desc <- category_descriptives(d, design, "joint_group", canonical)
  joint_rows <- joint_tests <- interaction_tests <- list()
  for (model in c("M1", "M2", "M3")) {
    covars <- model_covars(model, "Joint")
    fit <- fit_survey(reformulate(c("joint_group", covars), "kidney_stone"), design)
    register_model(paste0("Joint_", model), fit, "16 joint categories of PA and ST; P4&S1 reference")
    tab <- category_estimates(fit, "joint_group", canonical, "P4_S1", model, "Joint", "Main")
    tab$q_BH <- NA_real_
    tab$q_BH[tab$category != "P4_S1"] <- p.adjust(tab$P[tab$category != "P4_S1"], "BH", n = 15)
    joint_rows[[model]] <- dplyr::left_join(tab, joint_desc, by = "category")
    joint_tests[[model]] <- wald_test(fit, "joint_group", model, "Joint", "16-group overall association")
    interaction_fit <- fit_survey(reformulate(c("PA_group * ST_group", covars), "kidney_stone"), design)
    register_model(paste0("Interaction_", model), interaction_fit, "Multiplicative odds-scale PA-by-ST interaction")
    test <- wald_test(interaction_fit, "PA_group:ST_group", model, "PA_by_ST", "Multiplicative interaction")
    test$equivalent_joint_prediction_error <- max(abs(fitted(fit) - fitted(interaction_fit)))
    stopifnot(test$numerator_df == 9L, joint_tests[[model]]$numerator_df == 15L,
              test$equivalent_joint_prediction_error < 1e-8)
    interaction_tests[[model]] <- test
  }
  joint <- bind_rows_base(joint_rows)
  # 对与主图相同的调整层级计算标准化患病比例及复杂抽样方差。
  margin_covars <- model_covars(RCS_MODEL, "Joint")
  adj_formula <- if (length(margin_covars)) reformulate(margin_covars, "kidney_stone") else kidney_stone ~ 1
  adjustment_fit <- fit_survey(adj_formula, design)
  margins <- survey::svypredmeans(adjustment_fit, ~joint_group,
    predictat = factor(canonical, levels = levels(d$joint_group)))
  mu <- coef(margins)[canonical]; Vmargin <- vcov(margins)[canonical, canonical, drop = FALSE]
  crit <- qt(.975, models[[paste0("Joint_", RCS_MODEL)]]$df.residual)
  se_logit <- sqrt(diag(Vmargin)) / (mu * (1 - mu))
  adjusted <- data.frame(model = RCS_MODEL, category = canonical,
    adjusted_prevalence_pct = 100 * unname(mu),
    adjusted_low_pct = 100 * plogis(qlogis(mu) - crit * se_logit),
    adjusted_high_pct = 100 * plogis(qlogis(mu) + crit * se_logit))
  direct_mu <- vapply(canonical, function(g) {
    nd <- d; nd$joint_group <- factor(rep(g, nrow(d)), levels = levels(d$joint_group))
    weighted.mean(as.numeric(predict(models[[paste0("Joint_", RCS_MODEL)]], newdata = nd,
      type = "response", se.fit = FALSE)), d$WTMEC12YR)
  }, numeric(1))
  margin_error <- max(abs(unname(mu) - direct_mu))
  stopifnot(margin_error < 1e-8, all(is.finite(as.matrix(adjusted[, 3:5]))))

  message("正在拟合加权RCS及阈值模型...")
  rcs <- list()
  for (exposure in c("PA", "ST", "ASBI")) rcs[[exposure]] <- fit_main_rcs(d, exposure)
  threshold_records <- threshold_skips <- list()
  for (id in names(rcs)) {
    z <- rcs[[id]]; sm <- z$summary
    register_model(paste0("RCS_", id), z$fit, if (id %in% c("PA", "ST", "ASBI")) "Main RCS" else "Sensitivity RCS")
    if (RUN_THRESHOLD) {
      x <- switch(sm$exposure, PA = d$PAQ, ST = d$sedentary, ASBI = asbi_value(d, sm$constant))
      cutpoint <- sm$minimum_value
      below <- x <= cutpoint
      counts <- c(sum(below), sum(!below), sum(d$kidney_stone[below]), sum(d$kidney_stone[!below]))
      if (sm$minimum_at_boundary || cutpoint <= min(x) || cutpoint >= max(x) || min(counts) < 10L) {
        threshold_skips[[id]] <- data.frame(id = id, threshold = cutpoint,
          reason = "Minimum at display boundary or insufficient support for two-segment fit")
      } else {
        tr <- fit_threshold_effect(d, x, z$covars, cutpoint, sm$exposure, "Survey-weighted",
          switch(sm$exposure, PA = 100, ST = 60, ASBI = 1), sm$units,
          "Own RCS minimum on the 1st-99th percentile grid; selected in the same sample", sm$constant)
        for (kind in names(tr$models))
          register_model(paste0("Threshold_", id, "_", kind), tr$models[[kind]], "Exploratory threshold model")
        tr$models <- NULL
        for (nm in names(tr)) tr[[nm]] <- cbind(id = id, model = RCS_MODEL, PA_adjustment = sm$PA_adjustment, tr[[nm]])
        threshold_records[[id]] <- tr
      }
    }
  }
  rcs_summary <- bind_rows_base(lapply(names(rcs), function(id) cbind(id = id, rcs[[id]]$summary)))
  rcs_curves <- bind_rows_base(lapply(names(rcs), function(id) cbind(id = id, rcs[[id]]$curve)))
  plots <- list(PA = make_rcs_plot(rcs$PA), ST = make_rcs_plot(rcs$ST),
                ASBI = make_rcs_plot(rcs$ASBI), Joint = make_joint_plot(joint))
  prevalence <- survey::svyciprop(~kidney_stone, design, method = "logit", df = survey::degf(design))
  sample_summary <- data.frame(n = nrow(d), cases = sum(d$kidney_stone),
    weighted_prevalence_pct = 100 * as.numeric(coef(prevalence)),
    weighted_prevalence_low_pct = 100 * confint(prevalence)[1],
    weighted_prevalence_high_pct = 100 * confint(prevalence)[2],
    ST_480_group = ST_480_GROUP, n_ST_480 = sum(d$sedentary == 480),
    n_strata = length(unique(d$SDMVSTRA)), design_df = survey::degf(design),
    joint_margin_verification_error = margin_error)
  model_info <- bind_rows_base(model_meta)
  coefficients <- bind_rows_base(lapply(names(models), function(id) {
    fit <- models[[id]]
    data.frame(id = id, term = names(coef(fit)), beta = unname(coef(fit)),
      SE = unname(sqrt(diag(vcov(fit)))), row.names = NULL)
  }))
  notes <- analysis_notes()
  tables <- list(Readme = notes, Sample = sample_summary,
    Baseline = baseline_table(d, design), Category_descriptives = bind_rows_base(descriptions),
    Categorical_OR = categories[categories$role == "Main", ],
    Category_overall_tests = bind_rows_base(category_tests),
    Joint16_OR = joint_display_table(joint), Joint16_overall_tests = bind_rows_base(joint_tests),
    Interaction_tests = bind_rows_base(interaction_tests), Joint_adjusted_prevalence = joint_display_table(adjusted),
    RCS_summary = rcs_summary, RCS_curves = rcs_curves,
    Threshold_results = bind_rows_base(lapply(threshold_records, `[[`, "summary")),
    Threshold_slopes = bind_rows_base(lapply(threshold_records, `[[`, "slopes")),
    Threshold_coefficients = bind_rows_base(lapply(threshold_records, `[[`, "coefficients")),
    Threshold_validation = bind_rows_base(lapply(threshold_records, `[[`, "validation")),
    Threshold_skipped = bind_rows_base(threshold_skips),
    Model_formulas = model_info, Coefficients = coefficients)
  result <- list(data = d, design = design, models = models, tables = tables,
    rcs = rcs, plots = plots,
    joint_margins = margins, excel_file = NULL)
  
  message("Main analyses complete.")
  result
}

run_final_interaction <- function(d) {
  ns<-splines::ns
  des<-make_design(d)
  f<-reformulate(c("ns(sedentary, df=4) * PA_group",FULL_COVARS),"kidney_stone")
  fit<-fit_survey(f,des)
  tt<-survey::regTermTest(fit,~ns(sedentary,df=4):PA_group,method="Wald")
  b<-coef(fit);V<-vcov(fit);X<-model.matrix(fit)
  vary<-grepl("^\\(Intercept\\)$|ns\\(sedentary|PA_group",names(b));fixed<-!vary
  offset<-drop(X[,fixed,drop=FALSE]%*%b[fixed]);w<-d$WTMEC12YR/sum(d$WTMEC12YR)
  rg<-quantile(d$sedentary[d$PA_group=="P4"],c(.01,.99),names=FALSE)
  cr<-expand.grid(sedentary=seq(rg[1],rg[2],length.out=161),PA_group=paste0("P",1:4))
  cr$PA_group<-factor(cr$PA_group,levels=levels(d$PA_group))
  nd<-d[rep(1,nrow(cr)),];nd$sedentary<-cr$sedentary;nd$PA_group<-cr$PA_group
  Xn<-model.matrix(delete.response(terms(fit)),nd,contrasts.arg=fit$contrasts)
  critical<-qt(.975,fit$df.residual)
  pred<-t(vapply(seq_len(nrow(cr)),function(i) {
    pi<-plogis(offset+sum(Xn[i,vary]*b[vary]));mu<-sum(w*pi);a<-w*pi*(1-pi)
    g<-numeric(length(b));g[vary]<-Xn[i,vary]*sum(a);g[fixed]<-drop(crossprod(X[,fixed,drop=FALSE],a))
    se<-sqrt(drop(t(g)%*%V%*%g))/(mu*(1-mu))
    c(prevalence_pct=100*mu,low_pct=100*plogis(qlogis(mu)-critical*se),high_pct=100*plogis(qlogis(mu)+critical*se))
  },numeric(3)))
  cr<-cbind(cr,pred)
  test<-data.frame(n=nrow(d),cases=sum(d$kidney_stone),spline_df=4,interaction_df=tt$df,
    denominator_df=tt$ddf,F=as.numeric(tt$Ftest),P=as.numeric(tt$p),formula=formula_text(f))
  stopifnot(test$interaction_df==12L)
  cols<-c(P1="#ADD4E5",P2="#C6D0D4",P3="#EBBCB9",P4="#C1E2D4")
  pl<-ggplot(cr,aes(sedentary,prevalence_pct,colour=PA_group,fill=PA_group))+
    geom_ribbon(aes(ymin=low_pct,ymax=high_pct),alpha=.18,colour=NA)+geom_line(linewidth=.9)+
    scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
    scale_y_continuous(breaks=c(0,10,20,30),labels=function(x)paste0(x,"%"))+
    annotate("text",x=Inf,y=26,label=sprintf('italic(P)~"for interaction"==%.3f',test$P),parse=TRUE,
      hjust=1.05,family=FONT_FAMILY,size=FIGURE_LABEL_SIZE/.pt)+coord_cartesian(ylim=c(0,30))+
    labs(x="ST (min/day)",y="Adjusted prevalence (%)",colour=NULL,fill=NULL)+sf_theme()+
    theme(legend.position="inside",legend.position.inside=c(.64,.96),legend.direction="horizontal",legend.background=element_blank())
  list(model=fit,tables=list(Test=test,Curves=cr,Knots=data.frame(type=c(rep("Interior",3),rep("Boundary",2)),
    value=c(attr(ns(d$sedentary,df=4),"knots"),attr(ns(d$sedentary,df=4),"Boundary.knots")))),plot=style_final_plot(pl))
}

final_baseline_tests <- function(d,design) {
  vars<-c("age_group","sex","race","education","marital","pir","bmi","waist","drink","smoke",
    "water_g","diabetes","hyperlipid","hypertension","uric_g","ST_group","PA_group","ASBI_group")
  bind_rows_base(lapply(vars,function(v) {
    if(v%in%c("bmi","waist")) {
      tt<-survey::svyranktest(reformulate("kidney_stone",response=v),design,test="wilcoxon")
      method<-"Survey-weighted Wilcoxon rank-sum"
    } else {
      tt<-survey::svychisq(reformulate(c(v,"kidney_stone")),design,statistic="F")
      method<-"Rao-Scott F"
    }
    data.frame(variable=v,method=method,P=unname(tt$p.value),display=p_text(tt$p.value))
  }))
}

final_baseline_wide <- function(base,tests) {
  z<-base[base$variable%in%tests$variable,];key<-unique(z[c("variable","level")])
  key<-key[order(match(key$variable,tests$variable),match(key$level,c("",0:9,paste0("S",1:4),paste0("P",1:4),paste0("Q",1:4)))),]
  key$id<-paste(key$variable,key$level)
  for(g in c("All","No_history","History")) {
    zz<-z[z$group==g,];key[[g]]<-zz$display[match(key$id,paste(zz$variable,zz$level))]
  }
  key$P<-tests$P[match(key$variable,tests$variable)];key$P[duplicated(key$variable)]<-NA_real_
  key$id<-NULL;rownames(key)<-NULL;key
}

final_asbi_distribution <- function(d,design) {
  q<-as.numeric(coef(survey::svyquantile(~ASBI,design,quantiles=c(.25,.5,.75),ci=FALSE,qrule="math")))
  cuts<-round(q,2)
  tab<-category_descriptives(d,design,"ASBI_group",paste0("Q",1:4))
  tab$median<-vapply(tab$category,function(g)rcs_weighted_median(d$ASBI,d,d$ASBI_group==g),numeric(1))
  br<-seq(floor(min(d$ASBI)*2)/2,ceiling(max(d$ASBI)*2)/2,by=.5)
  h<-data.frame(bin=cut(d$ASBI,br,right=FALSE,include.lowest=TRUE,labels=FALSE),w=d$WTMEC12YR)
  hist<-data.frame(left=head(br,-1),right=tail(br,-1),center=head(br,-1)+.25,pct=0)
  ws<-tapply(h$w,h$bin,sum);hist$pct[as.integer(names(ws))]<-100*ws/sum(h$w)
  list(cutpoints=data.frame(quantile=c(.25,.5,.75),raw_cutpoint=q,weighted_cutpoint=cuts),histogram=hist,categories=tab)
}

# 3. ASBI sensitivity analyses

fit_asbi_sensitivity <- function(base_design, constant, sample_name) {
  stopifnot(sample_name %in% c("Full", "Exclude_P1S1"))
  design_full <- base_design
  d <- design_full$variables
  x <- asbi_value(d, constant)
  knots <- as.numeric(quantile(x, KNOT_PROBS, type = 7))
  # 每个分析样本使用自身的加权中位数。结点和显示区间仍沿用全样本设定。
  analysis_keep <- if (sample_name == "Full") rep(TRUE, nrow(d)) else !d$exclude_P1S1
  ref <- rcs_weighted_median(x, d, keep = analysis_keep)
  xr <- as.numeric(quantile(x, DISPLAY_PROBS, type = 7))
  stopifnot(all(is.finite(x)), length(unique(knots)) == 4L, all(diff(knots) > 0))
  basis <- function(xx) {
    b <- Hmisc::rcspline.eval(xx, knots = knots, inclx = TRUE, norm = 2)
    colnames(b) <- SPLINE_NAMES; b
  }
  design_full$variables$ASBI <- x
  design_full$variables[SPLINE_NAMES] <- as.data.frame(basis(x))
  design_full$variables$analysis_keep <- analysis_keep
  # 从完整分析样本的设计对象中取子人群，不从删行后的数据重新建设计。
  design <- subset(design_full, analysis_keep)
  dd <- design$variables
  observed_range <- range(dd$ASBI)
  stopifnot(ref > observed_range[1], ref < observed_range[2],
    xr[1] >= observed_range[1], xr[2] <= observed_range[2])
  covars <- model_covars(RCS_MODEL, "ASBI")
  fit <- fit_survey(reformulate(c(SPLINE_NAMES, covars), "kidney_stone"), design)
  critical <- qt(.975, fit$df.residual)
  grid <- sort(unique(c(seq(xr[1], xr[2], length.out = 800L), ref)))
  delta <- sweep(basis(grid), 2L, as.numeric(basis(ref)), FUN = "-")
  beta <- coef(fit)[SPLINE_NAMES]
  V <- vcov(fit)[SPLINE_NAMES, SPLINE_NAMES, drop = FALSE]
  eta <- drop(delta %*% beta)
  variance <- rowSums((delta %*% V) * delta)
  stopifnot(all(variance > -1e-10))
  se <- sqrt(pmax(variance, 0))
  curve <- data.frame(value = grid, log_OR = eta, SE = se, OR = exp(eta),
    OR_low = exp(eta - critical * se), OR_high = exp(eta + critical * se))
  at_ref <- curve[curve$value == ref, ]
  stopifnot(nrow(at_ref) == 1L, max(abs(unlist(at_ref[c("OR", "OR_low", "OR_high")]) - 1)) < 1e-10,
    all(is.finite(as.matrix(curve))), all(curve$OR_low <= curve$OR & curve$OR <= curve$OR_high))
  # 独立核对完整设计矩阵与协方差的对比结果。
  nd <- dd[rep(1L, length(grid)), , drop = FALSE]; nd[SPLINE_NAMES] <- as.data.frame(basis(grid))
  nr <- dd[1L, , drop = FALSE]; nr[SPLINE_NAMES] <- as.data.frame(basis(ref))
  X <- model.matrix(delete.response(terms(fit)), nd, contrasts.arg = fit$contrasts)
  Xr <- model.matrix(delete.response(terms(fit)), nr, contrasts.arg = fit$contrasts)
  contrast <- sweep(X, 2L, as.numeric(Xr), FUN = "-")
  prediction_error <- max(abs(drop(contrast %*% coef(fit)) - eta))
  variance_error <- max(abs(rowSums((contrast %*% vcov(fit)) * contrast) - variance))
  stopifnot(prediction_error < 1e-8, variance_error < 1e-8)
  overall <- wald_test(fit, SPLINE_NAMES, RCS_MODEL, "ASBI")
  nonlinear <- wald_test(fit, SPLINE_NAMES[-1], RCS_MODEL, "ASBI", "Nonlinear")
  imin <- which.min(curve$OR)
  sm <- data.frame(exposure = "ASBI", sample = sample_name, constant = constant, model = RCS_MODEL,
    n = nrow(dd), cases = sum(dd$kidney_stone),
    median_reference = ref, reference_source = "Analysis-sample survey-weighted median for this constant; qrule=math",
    sample_median = ref, sample_unweighted_median = median(dd$ASBI), minimum_value = curve$value[imin], minimum_OR = curve$OR[imin],
    minimum_at_boundary = imin %in% c(1L, nrow(curve)),
    P_overall = overall$P, P_nonlinear = nonlinear$P,
    F_overall = overall$F, F_nonlinear = nonlinear$F, df_overall = overall$numerator_df,
    df_nonlinear = nonlinear$numerator_df, denominator_df = overall$denominator_df,
    knot_5 = knots[1], knot_35 = knots[2], knot_65 = knots[3], knot_95 = knots[4],
    display_lower = xr[1], display_upper = xr[2], units = "ASBI score",
    design_df = survey::degf(design), formula = formula_text(formula(fit)),
    direct_prediction_error = prediction_error, full_covariance_error = variance_error,
    reference_OR = at_ref$OR, reference_low = at_ref$OR_low, reference_high = at_ref$OR_high,
    domain_subset = sample_name == "Exclude_P1S1", PA_adjustment = "Without_PA")
  list(fit = fit, summary = sm, curve = curve, covars = covars)
}

fit_sensitivity_threshold <- function(base_design, record, id) {
  sm <- record$summary
  x <- asbi_value(base_design$variables, sm$constant)
  cutpoint <- sm$minimum_value
  design_full <- base_design
  design_full$variables$threshold_x <- x - cutpoint
  design_full$variables$threshold_below <- pmin(x - cutpoint, 0)
  design_full$variables$threshold_above <- pmax(x - cutpoint, 0)
  design_full$variables$analysis_keep <- if (sm$sample == "Full") TRUE else !design_full$variables$exclude_P1S1
  design <- subset(design_full, analysis_keep)
  dd <- design$variables
  below <- dd$threshold_x <= 0
  side_support <- c(sum(below), sum(!below), sum(dd$kidney_stone[below]),
    sum(dd$kidney_stone[!below]), sum(1L - dd$kidney_stone[below]), sum(1L - dd$kidney_stone[!below]))
  if (sm$minimum_at_boundary || min(side_support) < 10L)
    return(list(skipped = data.frame(id = id, sample = sm$sample, constant = sm$constant,
      threshold = cutpoint, reason = "Minimum at display boundary or insufficient support on either side")))
  covars <- record$covars
  fits <- list(
    Linear = fit_survey(reformulate(c("threshold_x", covars), "kidney_stone"), design),
    Piecewise = fit_survey(reformulate(c("threshold_below", "threshold_above", covars), "kidney_stone"), design),
    Hinge = fit_survey(reformulate(c("threshold_x", "threshold_above", covars), "kidney_stone"), design))
  extract <- function(fit, term, segment) {
    b <- unname(coef(fit)[term]); se <- sqrt(vcov(fit)[term, term])
    critical <- qt(.975, fit$df.residual)
    data.frame(id = id, sample = sm$sample, constant = sm$constant, threshold = cutpoint,
      segment = segment, increment = 1, log_OR = b, SE = se, OR = exp(b),
      OR_low = exp(b - critical * se), OR_high = exp(b + critical * se),
      P = 2 * pt(-abs(b / se), fit$df.residual))
  }
  slopes <- rbind(extract(fits$Linear, "threshold_x", "Overall linear"),
    extract(fits$Piecewise, "threshold_below", "At or below threshold"),
    extract(fits$Piecewise, "threshold_above", "Above threshold"))
  con <- survey::svycontrast(fits$Piecewise, c(threshold_below = -1, threshold_above = 1))
  difference <- as.numeric(coef(con)); se <- as.numeric(survey::SE(con))
  Fstat <- (difference / se)^2
  p <- pf(Fstat, 1, fits$Piecewise$df.residual, lower.tail = FALSE)
  tt <- survey::regTermTest(fits$Hinge, ~threshold_above, method = "Wald")
  pred_error <- max(abs(as.numeric(predict(fits$Piecewise, type = "link")) -
                        as.numeric(predict(fits$Hinge, type = "link"))))
  b_error <- abs(difference - unname(coef(fits$Hinge)["threshold_above"]))
  v_error <- abs(se^2 - vcov(fits$Hinge)["threshold_above", "threshold_above"])
  p_error <- abs(p - as.numeric(tt$p))
  stopifnot(pred_error < 1e-7, b_error < 1e-7, v_error < 1e-7, p_error < 1e-7,
    all(vapply(fits, nobs, integer(1)) == sm$n))
  summary <- data.frame(id = id, sample = sm$sample, constant = sm$constant, model = RCS_MODEL,
    n = sm$n, cases = sm$cases, threshold = cutpoint,
    n_below = sum(below), cases_below = sum(dd$kidney_stone[below]),
    n_above = sum(!below), cases_above = sum(dd$kidney_stone[!below]),
    slope_difference = difference, SE_difference = se, F = Fstat,
    numerator_df = 1, denominator_df = fits$Piecewise$df.residual, P_difference = p,
    inference = "Exploratory, conditional on a same-sample selected RCS minimum; no selection-adjusted inference")
  check <- data.frame(id = id, prediction_error = pred_error, slope_contrast_error = b_error,
    variance_error = v_error, P_difference_error = p_error, passed = TRUE)
  list(summary = summary, slopes = slopes, validation = check, models = fits)
}

run_final_sensitivity <- function(d,design) {
  records<-thresholds<-list()
  for(sample_name in c("Full","Exclude_P1S1")) for(cc in c(.01,.1,.5)) {
    id<-paste0(sample_name,"_c",gsub("\\.","p",format(cc,scientific=FALSE,trim=TRUE)))
    records[[id]]<-fit_asbi_sensitivity(design,cc,sample_name)
    thresholds[[id]]<-fit_sensitivity_threshold(design,records[[id]],id)
  }
  tables<-list(RCS_summary=bind_rows_base(lapply(names(records),function(id)cbind(id=id,records[[id]]$summary))),
    RCS_curves=bind_rows_base(lapply(names(records),function(id)cbind(id=id,records[[id]]$curve))),
    Threshold_results=bind_rows_base(lapply(thresholds,`[[`,"summary")),
    Threshold_slopes=bind_rows_base(lapply(thresholds,`[[`,"slopes")),
    Threshold_validation=bind_rows_base(lapply(thresholds,`[[`,"validation")),
    Threshold_skipped=bind_rows_base(lapply(thresholds,`[[`,"skipped")))
  list(records=records,thresholds=thresholds,tables=tables,plots=lapply(records,make_rcs_plot))
}

# 4. Mutually adjusted PA component models

PC_DEMO <- c("sex", "age_group", "race", "education", "marital", "pir")

PC_LIFESTYLE <- c("smoke", "drink", "water_g")

PC_METABOLIC <- c("diabetes", "hypertension", "hyperlipid", "uric_g")

PC_COVARS <- list(A5 = c(PC_DEMO, PC_LIFESTYLE, "ST_group", "bmi", "waist", PC_METABOLIC))

PC_FAMILIES <- list(Intensity=c("VPA","NonVPA"), Domain=c("Work","Transport","Leisure"),
  DomainIntensity=c("WorkVig","WorkMod","Transport","LeisureVig","LeisureMod"))

PC_LABELS <- c(VPA="Vigorous",NonVPA="Moderate incl. transport",Work="Work / household",Transport="Transport",
  Leisure="Leisure",WorkVig="Vigorous work / household",WorkMod="Moderate work / household",
  LeisureVig="Vigorous leisure",LeisureMod="Moderate leisure")

pc_bind <- function(x) {
  if (!length(x)) return(data.frame())
  out <- do.call(rbind, x); rownames(out) <- NULL; out
}

pc_design <- function(d) survey::svydesign(ids=~SDMVPSU, strata=~SDMVSTRA,
  weights=~WTMEC12YR, nest=TRUE, data=d)

pc_fit <- function(f, ds) {
  # Only named data columns occur in formulas; avoid retaining the complete calling workspace.
  environment(f) <- baseenv()
  fit <- survey::svyglm(f, design=ds, family=quasibinomial(), na.action=na.fail,
    control=glm.control(epsilon=1e-10, maxit=100), model=TRUE)
  stopifnot(isTRUE(fit$converged), fit$rank==length(coef(fit)), !anyNA(coef(fit)),
    nobs(fit)==nrow(ds$variables), identical(rownames(fit$model),rownames(ds$variables)),
    all(is.finite(vcov(fit))), fit$df.residual>0)
  fit
}

pc_wald <- function(fit,terms) {
  t <- survey::regTermTest(fit,reformulate(terms),method="Wald")
  data.frame(F=as.numeric(t$Ftest),numerator_df=as.numeric(t$df),denominator_df=as.numeric(t$ddf),P=as.numeric(t$p))
}

pc_contrast <- function(fit,L) {
  stopifnot(length(L)==length(coef(fit)))
  beta <- drop(L%*%coef(fit)); se <- sqrt(drop(L%*%vcov(fit)%*%L)); df <- fit$df.residual
  crit <- qt(.975,df)
  data.frame(log_OR=beta,SE=se,OR=exp(beta),OR_low=exp(beta-crit*se),OR_high=exp(beta+crit*se),
    P=if(se>0)2*pt(-abs(beta/se),df) else NA_real_,denominator_df=df,
    OR_95CI=if(se>0)sprintf("%.3f (%.3f-%.3f)",exp(beta),exp(beta-crit*se),exp(beta+crit*se)) else "1.000 (Reference)")
}

pc_L <- function(fit,positive,negative=character()) {
  L <- setNames(rep(0,length(coef(fit))),names(coef(fit)))
  stopifnot(all(c(positive,negative)%in%names(L)))
  L[positive] <- 1; L[negative] <- -1; L
}

pc_prepare <- function(d) {
  input <- c(WorkVig="PAQ_work_vig", WorkMod="PAQ_work_mod",
    Transport="PAQ_transport", LeisureVig="PAQ_leisure_vig", LeisureMod="PAQ_leisure_mod")
  stopifnot(all(input %in% names(d)))
  for(nm in names(input)) {
    x <- d[[input[[nm]]]]
    stopifnot(is.numeric(x), all(is.finite(x)), all(x >= 0))
    d[[paste0(nm,"_met")]] <- x / 7 # supplied MET-min/week -> MET-min/day
  }
  d$VPA_met <- d$WorkVig_met + d$LeisureVig_met
  d$NonVPA_met <- d$WorkMod_met + d$LeisureMod_met + d$Transport_met
  d$Work_met <- d$WorkVig_met + d$WorkMod_met
  d$Leisure_met <- d$LeisureVig_met + d$LeisureMod_met
  for(nm in unique(unlist(PC_FAMILIES))) {
    x <- d[[paste0(nm,"_met")]]
    d[[paste0(nm,"_per100")]] <- x / 100
    d[[paste0(nm,"_any_num")]] <- as.integer(x > 0)
  }
  for(fam in PC_FAMILIES) stopifnot(
    max(abs(d$PAQ-rowSums(d[paste0(fam,"_met")]))) < 1e-8)
  stopifnot(nrow(d)==EXPECTED_N, sum(d$kidney_stone)==EXPECTED_CASES)
  rownames(d) <- as.character(d$SEQN)
  list(data=d)
}

pc_describe <- function(ds) {
  d <- ds$variables; w <- d$WTMEC12YR; comp <- list(); distributions <- list(); participation <- list()
  for(fam in names(PC_FAMILIES)) for(g in levels(d$PA_group)) for(nm in PC_FAMILIES[[fam]]) {
    sub <- ds[d$PA_group==g,]; dd <- sub$variables; v <- paste0(nm,"_met")
    ratio <- survey::svyratio(reformulate(v),~PAQ,sub)
    est <- as.numeric(coef(ratio)); se <- as.numeric(survey::SE(ratio)); crit <- qt(.975,survey::degf(sub))
    comp[[paste(fam,g,nm)]] <- data.frame(family=fam,PA_group=g,component=nm,n=nrow(dd),
      cases=sum(dd$kidney_stone),volume_share_pct=100*est,SE_pct=100*se,
      low_pct=100*(est-crit*se),high_pct=100*(est+crit*se))
    pr <- survey::svyciprop(reformulate(paste0(nm,"_any_num")),sub,method="logit",df=survey::degf(sub))
    ci <- as.numeric(confint(pr))
    participation[[paste(fam,g,nm)]] <- data.frame(family=fam,PA_group=g,component=nm,
      active_n=sum(dd[[paste0(nm,"_any_num")]]),active_cases=sum(dd$kidney_stone[dd[[paste0(nm,"_any_num")]]==1]),
      weighted_active_pct=100*as.numeric(coef(pr)),low_pct=100*ci[1],high_pct=100*ci[2])
  }
  for(nm in unique(unlist(PC_FAMILIES))) {
    x <- d[[paste0(nm,"_met")]]
    distributions[[nm]] <- data.frame(component=nm,n=length(x),zero_n=sum(x==0),zero_cases=sum(d$kidney_stone[x==0]),
      weighted_zero_pct=100*sum(w[x==0])/sum(w),weighted_mean_MET=weighted.mean(x,w),
      positive_median_MET=rcs_weighted_median(x,d,x>0),
      positive_unweighted_median_MET=median(x[x>0]),positive_max_MET=max(x))
  }
  list(composition=pc_bind(comp),participation=pc_bind(participation),distributions=pc_bind(distributions))
}

pc_definitions <- function() data.frame(
  component=c("VPA","NonVPA","Work","Transport","Leisure","WorkVig","WorkMod","LeisureVig","LeisureMod"),
  definition=c("WorkVig + LeisureVig","WorkMod + Transport + LeisureMod","WorkVig + WorkMod",
    "Walking/cycling for transport, assigned 4 MET","LeisureVig + LeisureMod",
    "Vigorous work/household/yard, assigned 8 MET","Moderate work/household/yard, assigned 4 MET",
    "Vigorous leisure, assigned 8 MET","Moderate leisure, assigned 4 MET"),
  questionnaire=c("PAQ605 + PAQ650","PAQ620 + PAQ635 + PAQ665","PAQ605 + PAQ620",
    "PAQ635 / PAQ640 / PAD645","PAQ650 + PAQ665","PAQ605 / PAQ610 / PAD615",
    "PAQ620 / PAQ625 / PAD630","PAQ650 / PAQ655 / PAD660","PAQ665 / PAQ670 / PAD675"),
  unit="MET-min/day = days/week * minutes/day * MET / 7")

run_nhanes_three_schemes_final <- function(loaded) {
  d <- pc_prepare(loaded$data)$data
  ds <- pc_design(d)
  models <- linear <- tests <- audit <- list()
  for(fam in names(PC_FAMILIES)) {
    components <- PC_FAMILIES[[fam]]
    terms <- paste0(components,"_per100")
    fit <- pc_fit(reformulate(c(terms,PC_COVARS$A5),"kidney_stone"),ds)
    models[[fam]] <- fit
    for(nm in components) linear[[paste(fam,nm)]] <- cbind(
      data.frame(family=fam,model="A5",component=nm,unit="per 100 MET-min/day"),
      pc_contrast(fit,pc_L(fit,paste0(nm,"_per100"))))
    tests[[fam]] <- cbind(data.frame(family=fam,model="A5"),pc_wald(fit,terms))
    audit[[fam]] <- data.frame(family=fam,n=nobs(fit),cases=sum(d$kidney_stone),
      converged=fit$converged,residual_df=fit$df.residual,formula=formula_text(formula(fit)))
  }
  desc <- pc_describe(ds)
  list(data=d,design=ds,models=models,tables=list(
    Linear_OR=pc_bind(linear),Composition=desc$composition,
    Component_definitions=pc_definitions(),Association_tests=pc_bind(tests),Model_audit=pc_bind(audit)))
}

# 5. Continuous-PA mediation with survey jackknife variance

M6C_PA_INCREMENT <- 100 # MET-min/day; 固定比较，不依据P值选择

M6C_PA_SCALE <- 100     # exposure unit in both regression models

M6C_QUADRATURE <- 20L

M6C_EXPECTED_N <- 19977L

M6C_EXPECTED_CASES <- 1713L

M6C_COVARS <- c("sex","age_group","race","education","marital","pir",
                "smoke","drink","water_g","ST_group")

M6C_MEDIATORS <- c(BMI="bmi", Waist="waist", Hypertension="hypertension",
                   Diabetes="diabetes", Hyperlipidemia="hyperlipid", Uric_acid="uric_acid")

M6C_BINARY <- c("hypertension","diabetes","hyperlipid")

m6c_bind <- function(x) { z <- do.call(rbind,x); rownames(z) <- NULL; z }

m6c_prepare <- function(d) {
  d$PA_per100 <- d$PAQ / M6C_PA_SCALE
  d$kidney_stone <- as.numeric(as.character(d$kidney_stone))
  d$PA_group <- cut(d$PAQ,c(-Inf,20,200,600,Inf),right=FALSE,labels=paste0("P",1:4))
  d$ST_group <- cut(d$sedentary,c(-Inf,240,360,480,Inf),right=FALSE,labels=paste0("S",1:4))
  for(v in M6C_BINARY) d[[v]] <- as.numeric(as.character(d[[v]]))
  needed <- c("SEQN","kidney_stone","PAQ","PA_per100",M6C_COVARS,unname(M6C_MEDIATORS),
              "WTMEC12YR","SDMVSTRA","SDMVPSU")
  stopifnot(nrow(d)==M6C_EXPECTED_N,sum(d$kidney_stone)==M6C_EXPECTED_CASES,
    !anyDuplicated(d$SEQN),all(complete.cases(d[needed])),
    all(d$WTMEC12YR>0),all(d$kidney_stone%in%0:1),all(is.finite(d$PAQ) & d$PAQ>=0))
  for(v in M6C_BINARY) stopifnot(all(d[[v]]%in%0:1))
  for(v in setdiff(unname(M6C_MEDIATORS),M6C_BINARY)) stopifnot(all(is.finite(d[[v]]) & d[[v]]>0))
  d
}

m6c_gh <- function(k) {
  # Golub-Welsch quadrature for a standard normal distribution, not exp(-x^2).
  J <- matrix(0,k,k); J[cbind(seq_len(k-1L),2:k)] <- sqrt(seq_len(k-1L))
  e <- eigen(J+t(J),symmetric=TRUE)
  list(nodes=e$values,weights=e$vectors[1,]^2)
}

m6c_cache <- function(d,variable,contrasts,exposure_values) {
  binary <- variable%in%M6C_BINARY
  center <- if(binary) 0 else median(log(d[[variable]]))
  dd <- d; dd$M <- if(binary) d[[variable]] else log(d[[variable]])-center
  fm <- reformulate(c("PA_per100",M6C_COVARS),response="M")
  fy <- reformulate(c("PA_per100*M",M6C_COVARS),response="kidney_stone")
  environment(fm) <- baseenv(); environment(fy) <- baseenv()
  Xm <- model.matrix(fm,dd); Xy <- model.matrix(fy,dd)
  mx <- y0 <- ydiff <- list()
  for(a in names(exposure_values)) {
    new <- dd; new$PA_per100 <- unname(exposure_values[a]) / M6C_PA_SCALE
    mx[[a]] <- model.matrix(fm,new)
    new$M <- 0; y0[[a]] <- model.matrix(fy,new)
    new$M <- 1; ydiff[[a]] <- model.matrix(fy,new)-y0[[a]]
  }
  pairs <- unique(do.call(rbind,lapply(seq_len(nrow(contrasts)),function(i) {
    a <- contrasts$active[i]; r <- contrasts$reference[i]
    data.frame(outcome_PA=c(a,a,r),mediator_PA=c(a,r,r))
  })))
  pairs$key <- paste0(pairs$outcome_PA,"_M",pairs$mediator_PA)
  list(data=dd,variable=variable,binary=binary,center=center,fm=fm,fy=fy,
       contrasts=contrasts,exposure_values=exposure_values,Xm=Xm,Xy=Xy,mx=mx,y0=y0,ydiff=ydiff,pairs=pairs,gh=m6c_gh(M6C_QUADRATURE))
}

m6c_fit <- function(X,y,w,family) {
  f <- glm.fit(X,y,weights=w/mean(w),family=family,
               control=glm.control(epsilon=1e-9,maxit=80))
  if(!f$converged || f$rank!=ncol(X) || any(!is.finite(f$coefficients)))
    stop("模型未收敛或出现秩缺失，停止计算，不忽略失败重复。")
  f
}

m6c_probabilities <- function(cache,bm,by,sigma,w,gh=cache$gh) {
  mu_m <- lapply(cache$mx,function(x) drop(x%*%bm))
  b0 <- lapply(cache$y0,function(x) drop(x%*%by))
  bM <- lapply(cache$ydiff,function(x) drop(x%*%by))
  out <- setNames(numeric(nrow(cache$pairs)),cache$pairs$key)
  for(i in seq_len(nrow(cache$pairs))) {
    a <- cache$pairs$outcome_PA[i]; b <- cache$pairs$mediator_PA[i]
    if(cache$binary) {
      pm <- plogis(mu_m[[b]])
      p <- (1-pm)*plogis(b0[[a]])+pm*plogis(b0[[a]]+bM[[a]])
    } else {
      # Integrate the outcome probability over the fitted log-mediator distribution.
      location <- b0[[a]]+bM[[a]]*mu_m[[b]]
      spread <- bM[[a]]*sigma
      p <- numeric(length(w))
      for(j in seq_along(gh$nodes)) p <- p+gh$weights[j]*plogis(location+spread*gh$nodes[j])
    }
    out[i] <- sum(w*p)/sum(w)
  }
  stopifnot(all(is.finite(out)),all(out>0 & out<1))
  out
}

m6c_decompose <- function(mu,contrasts) {
  out <- mu
  for(i in seq_len(nrow(contrasts))) {
    x <- contrasts[i,]; a <- x$active; r <- x$reference
    paa <- unname(mu[paste0(a,"_M",a)]); par <- unname(mu[paste0(a,"_M",r)])
    prr <- unname(mu[paste0(r,"_M",r)])
    # PNDE + TNIE = TE on the probability-difference scale, with A*M interaction.
    rd <- c(TE=paa-prr,PNDE=par-prr,TNIE=paa-par)
    lo <- c(TE=qlogis(paa)-qlogis(prr),PNDE=qlogis(par)-qlogis(prr),TNIE=qlogis(paa)-qlogis(par))
    names(rd) <- paste(x$contrast,names(rd),"RD",sep="__")
    names(lo) <- paste(x$contrast,names(lo),"logOR",sep="__")
    out <- c(out,rd,lo)
  }
  out
}

m6c_theta <- function(w,cache,details=FALSE) {
  fm <- m6c_fit(cache$Xm,cache$data$M,w,if(cache$binary) quasibinomial() else gaussian())
  fy <- m6c_fit(cache$Xy,cache$data$kidney_stone,w,quasibinomial())
  sigma <- if(cache$binary) NA_real_ else sqrt(sum(w*(cache$data$M-fm$fitted.values)^2)/sum(w))
  mu <- m6c_probabilities(cache,fm$coefficients,fy$coefficients,sigma,w)
  th <- m6c_decompose(mu,cache$contrasts)
  if(details) list(theta=th,mediator_fit=fm,outcome_fit=fy,sigma=sigma) else th
}

m6c_effect_table <- function(nm,theta,V,df,contrasts) {
  z <- qt(.975,df); rows <- list()
  for(i in seq_len(nrow(contrasts))) for(effect in c("TE","PNDE","TNIE")) {
    co <- contrasts[i,]; kr <- paste(co$contrast,effect,"RD",sep="__")
    ko <- paste(co$contrast,effect,"logOR",sep="__")
    est <- unname(theta[kr]); se <- sqrt(V[kr,kr]); logor <- unname(theta[ko]); sol <- sqrt(V[ko,ko])
    rows[[length(rows)+1L]] <- data.frame(mediator=nm,co,effect=effect,
      difference_pp=100*est,low_pp=100*(est-z*se),high_pp=100*(est+z*se),SE_pp=100*se,
      P=2*pt(-abs(est/se),df),OR=exp(logor),OR_low=exp(logor-z*sol),OR_high=exp(logor+z*sol),
      logOR=logor,SE_logOR=sol,P_OR=2*pt(-abs(logor/sol),df),df=df)
  }
  m6c_bind(rows)
}

run_nhanes_continuous_PA_mediation <- function(loaded,do_replicates=TRUE) {
  needed <- c("survey","ggplot2")
  missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("请先安装R包: ",paste(missing,collapse=", "))
  suppressPackageStartupMessages(library(ggplot2))
  old_opt <- options(survey.lonely.psu="adjust"); on.exit(options(old_opt),add=TRUE)
  d <- m6c_prepare(loaded$data)
  des <- survey::svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC12YR,nest=TRUE,data=d)
  pa_reference <- as.numeric(coef(survey::svyquantile(~PAQ,des,quantiles=.5,ci=FALSE,qrule="math")))
  stopifnot(length(pa_reference)==1L,is.finite(pa_reference),M6C_PA_INCREMENT>0)
  exposure_values <- c(a0=pa_reference,a1=pa_reference+M6C_PA_INCREMENT)
  stopifnot(all(exposure_values>=min(d$PAQ)),all(exposure_values<=max(d$PAQ)))
  contrasts <- data.frame(contrast="PA_median_plus100",active="a1",reference="a0",
    role="Continuous PA: weighted median plus 100 MET-min/day versus weighted median",
    reference_MET_min_day=pa_reference,active_MET_min_day=pa_reference+M6C_PA_INCREMENT,
    increment_MET_min_day=M6C_PA_INCREMENT,regression_unit_MET_min_day=M6C_PA_SCALE)
  # Exposure levels are fixed at the full-sample values across all JKn replicates.
  # Inference is conditional on this chosen contrast; median uncertainty is not propagated.
  cat("Continuous PA contrast:", exposure_values["a0"], "to", exposure_values["a1"], "MET-min/day\n")
  # Stratified delete-one-PSU jackknife: refit BOTH models and restandardize each replicate.
  repdes <- survey::as.svrepdesign(des,type="JKn",mse=TRUE)
  w <- weights(des); df <- survey::degf(repdes); nr <- ncol(repdes$repweights)
  cat("Common analytic sample:",nrow(d),"; cases:",sum(d$kidney_stone),"; JKn replicates:",nr,"; df:",df,"\n")
  fits <- reps <- covs <- points <- list(); effect_rows <- audit_rows <- formulas <- diagnostics <- list()
  for(nm in names(M6C_MEDIATORS)) {
    v <- M6C_MEDIATORS[[nm]]; cat("Mediator:",nm,"\n")
    ca <- m6c_cache(d,v,contrasts,exposure_values); full <- m6c_theta(w,ca,details=TRUE)
    # Public svyglm fits retained in R; naive glm standard errors are never reported.
    dm <- des; dm$variables$M <- ca$data$M
    sfm <- survey::svyglm(ca$fm,design=dm,family=if(ca$binary)quasibinomial() else gaussian(),
                         control=glm.control(epsilon=1e-9,maxit=80))
    sfy <- survey::svyglm(ca$fy,design=dm,family=quasibinomial(),control=glm.control(epsilon=1e-9,maxit=80))
    delta <- max(abs(coef(sfm)-full$mediator_fit$coefficients),abs(coef(sfy)-full$outcome_fit$coefficients))
    stopifnot(delta<1e-6)
    if(!ca$binary) {
      finer <- m6c_probabilities(ca,coef(sfm),coef(sfy),full$sigma,w,gh=m6c_gh(2L*M6C_QUADRATURE))
      integration_error <- max(abs(finer-full$theta[names(finer)]))
      stopifnot(integration_error<1e-7)
    } else integration_error <- 0
    if(do_replicates) {
      counter <- 0L
      estimator <- function(w,data) {
        counter <<- counter+1L
        if(counter%%40L==0L) cat(" ",nm,"replicate",counter-1L,"/",nr,"\n")
        m6c_theta(w,ca)
      }
      rr <- survey::withReplicates(repdes,estimator,return.replicates=TRUE)
      point <- coef(rr); V <- vcov(rr)
      stopifnot(max(abs(point-full$theta))<1e-10,all(is.finite(V)),all(is.finite(rr$replicates)))
      reps[[nm]] <- rr$replicates; covs[[nm]] <- V
      effect_rows[[nm]] <- m6c_effect_table(nm,point,V,df,contrasts)
    } else point <- full$theta
    points[[nm]] <- point
    fits[[nm]] <- list(mediator=sfm,outcome=sfy,sigma=full$sigma,center_logM=ca$center,
                       binary=ca$binary,variable=v)
    audit_rows[[nm]] <- data.frame(mediator=nm,n=nrow(d),cases=sum(d$kidney_stone),
      mediator_converged=sfm$converged,outcome_converged=sfy$converged,
      mediator_rank=sfm$rank,outcome_rank=sfy$rank,coefficient_check=delta,
      integration_error=integration_error,replicates=if(do_replicates)nr else 0,df=df,
      sigma=full$sigma,center_logM=ca$center)
    formulas[[nm]] <- data.frame(mediator=nm,part=c("Mediator","Outcome"),
      scale=c(if(ca$binary)"Binary 0/1" else "Log mediator centered at its sample median", "Binary stone history 0/1"),
      formula=c(paste(deparse(ca$fm),collapse=" "),paste(deparse(ca$fy),collapse=" ")))
    for(a in levels(d$PA_group)) {
      ii <- d$PA_group==a; xx <- d[[v]][ii]; ww <- w[ii]
      rr0 <- ca$data$M[ii]-predict(sfm,type="response")[ii]
      diagnostics[[paste(nm,a)]] <- data.frame(mediator=nm,PA_group=a,n=sum(ii),cases=sum(d$kidney_stone[ii]),
        weighted_mean=weighted.mean(xx,ww),minimum=min(xx),maximum=max(xx),
        zero_n=sum(xx==0),one_n=if(ca$binary)sum(xx==1)else NA_integer_,
        weighted_residual_SD=if(!ca$binary)sqrt(weighted.mean(rr0^2,ww))else NA_real_)
    }
  }
  if(!do_replicates) return(list(points=points,fits=fits,audit=m6c_bind(audit_rows),data=d,design=des,contrasts=contrasts,exposure_values=exposure_values))
  effects <- m6c_bind(effect_rows)
  effects$q_TNIE_across_6 <- NA_real_
  for(co in contrasts$contrast) {
    ix <- effects$contrast==co & effects$effect=="TNIE"
    effects$q_TNIE_across_6[ix] <- p.adjust(effects$P[ix],method="BH")
  }
  pm <- m6c_bind(lapply(names(M6C_MEDIATORS),function(nm) m6c_bind(lapply(contrasts$contrast,function(co) {
    ee <- effects[effects$mediator==nm & effects$contrast==co,]
    te <- ee[ee$effect=="TE",]; de <- ee[ee$effect=="PNDE",]; ie <- ee[ee$effect=="TNIE",]
    # Do not convert cancellation or imprecise total effects into misleading percentages.
    reason <- if(te$low_pp<=0 && te$high_pp>=0) "Total-effect CI includes zero" else if(de$difference_pp*ie$difference_pp<=0)
      "Direct and indirect estimates have opposite signs or include a zero point estimate" else "Not a primary estimand; only descriptive"
    eligible <- te$low_pp*te$high_pp>0 && de$difference_pp*ie$difference_pp>0
    data.frame(mediator=nm,contrast=co,proportion_mediated=if(eligible)ie$difference_pp/te$difference_pp else NA_real_,
      reportable=eligible,reason=reason)
  }))))
  probs <- m6c_bind(lapply(names(points),function(nm) {
    z <- points[[nm]]; kk <- !grepl("__",names(z),fixed=TRUE)
    data.frame(mediator=nm,scenario=names(z)[kk],probability=unname(z[kk]))
  }))
  note <- c(
    Purpose="Exploratory mediation analyses supporting the manuscript and the detailed supplementary mediation table.",
    Primary_contrast=sprintf("Total PA from the survey-weighted median %.10f to %.10f MET-min/day (+%.1f); reference levels held fixed across JKn replicates. Confidence intervals are conditional on these chosen exposure values.",pa_reference,pa_reference+M6C_PA_INCREMENT,M6C_PA_INCREMENT),
    Exposure="Raw continuous total PA (PAQ) divided by 100 MET-min/day, without categorization or logarithmic transformation. Same 19977 participants, including zero PA. Descriptive PA groups are retained only for checking the sample; they are not exposure terms in any mediation model.",
    Mediators="BMI, waist, hypertension, diabetes, hyperlipidemia and serum uric acid, fitted in six separate models. BMI/waist/urate are analyzed on a centered log scale with a Gaussian conditional mediator distribution; binary mediators use logistic models.",
    Outcome="Self-reported lifetime kidney-stone history. Its time of occurrence relative to typical-week PA and contemporaneously measured mediators is not established.",
    Adjustments=paste(M6C_COVARS,collapse=" + "),
    Adjustment_rationale="Retain demographic, smoking, alcohol, water-category and ST covariates for consistency with the pre-body/metabolic adjustment model. Other candidate mediators are not mutually adjusted. This relies on a working covariate-role assumption, not empirically established baseline temporality; water and ST could themselves be affected by activity or prior stones.",
    Models="Mediator ~ continuous PA_per100 + covariates. Outcome logistic model ~ PA_per100 * mediator + covariates. Exposure-mediator interaction retained for every mediator. PA is linear on the Gaussian mediator scale or logistic linear predictor, not necessarily on the outcome probability scale. No significance-based variable selection.",
    Integration="Parametric g-computation integrates the outcome probability over the fitted mediator distribution at each exposure level, then averages over the common survey-weighted covariate distribution. Binary mediators are integrated exactly; continuous log-mediators use Gaussian quadrature with a common residual SD. Twice-as-many-node results are checked.",
    Definitions="mu(a,b)=weighted mean of predicted Y under outcome exposure a and the modeled mediator distribution under b. TE=mu(a,a)-mu(r,r); PNDE=mu(a,r)-mu(r,r); TNIE=mu(a,a)-mu(a,r). PNDE+TNIE=TE. Natural-effect names describe the assumed causal model; in these cross-sectional data they are exploratory model-based decompositions.",
    Effect_scales="Primary: standardized probability difference in percentage points of lifetime stone history, not incident risk. Also report marginal OR decomposition from the same three probabilities; OR_TE=OR_PNDE*OR_TNIE. These are not conditional logistic coefficients or their differences.",
    Variance="NHANES WTMEC2YR/6, strata SDMVSTRA, PSUs SDMVPSU, nest=TRUE. JKn stratified delete-one-PSU replicate weights; both models and target-population averaging refitted in EVERY replicate. MSE variance around full-weight estimates; t intervals and P values use replicate-design degrees of freedom. No individual-level iid bootstrap or naive GLM SEs.",
    Multiplicity="All six candidate mediators reported. BH q values adjust the six TNIE probability-difference tests for the single continuous-PA contrast. TE and PNDE P values are nominal.",
    Proportion="Not used as the main result. Suppressed when the total-effect CI crosses zero or direct and indirect point estimates oppose one another. Opposing paths can occur without a significant total effect; this alone does not prove suppression or mediation.",
    Identification="A causal interpretation would require correct exposure-before-mediator-before-outcome ordering, consistency, positivity, correct model forms, no unmeasured exposure-outcome/exposure-mediator/mediator-outcome confounding, and no exposure-induced mediator-outcome confounding for these natural effects. These assumptions cannot be verified here.",
    Multiple_mediators="Separate mediator models do not estimate mutually independent pathways. Adiposity and metabolic variables may cause or confound each other; omission or adjustment of another candidate can change the estimand or induce bias. Do not sum individual indirect effects or mediation proportions, and do not interpret results as a joint pathway decomposition.",
    Limitations="Cross-sectional temporal ambiguity, reverse causation, residual confounding, a linear continuous-PA working model, Gaussian log-mediator assumptions, and exposure/mediator measurement error remain. This script does not test or prove linearity of every path. A significant indirect estimate does not establish that PA changed the mediator before stones.",
    Inherited="Original cohort selection, disease/lipid definitions, water coding and MEC weights are retained for comparability. This module does not resolve their measurement or weighting limitations.",
    References=paste("AGReMA: https://jamanetwork.com/journals/jama/fullarticle/2784353;",
      "Cross-sectional mediation: https://pubmed.ncbi.nlm.nih.gov/17402810/;",
      "Binary-outcome decomposition: https://pmc.ncbi.nlm.nih.gov/articles/PMC11991894/;",
      "Replicate inference: https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/withReplicates.html"))
  display <- effects
  display$Difference_95CI_pp <- sprintf("%.3f (%.3f to %.3f)",display$difference_pp,display$low_pp,display$high_pp)
  display$OR_95CI <- sprintf("%.4f (%.4f to %.4f)",display$OR,display$OR_low,display$OR_high)
  display <- display[c("mediator","contrast","effect","Difference_95CI_pp","OR_95CI","P","q_TNIE_across_6")]
  tables <- list(Notes=data.frame(item=names(note),details=unname(note)),
    Primary=display,All_effects=effects,All_display=display,
    Probabilities=probs,Proportion=pm,Descriptive=m6c_bind(diagnostics),
    Model_audit=m6c_bind(audit_rows),Formulas=m6c_bind(formulas),Contrasts=contrasts,
    Sample=data.frame(n=nrow(d),cases=sum(d$kidney_stone),PA_zero_n=sum(d$PAQ==0),
      n_strata=length(unique(d$SDMVSTRA)),n_PSU=length(unique(interaction(d$SDMVSTRA,d$SDMVPSU))),
      JKn_replicates=nr,df=df,ST_480_group="S4"),
    Versions=data.frame(package=c("R",needed),version=c(R.version.string,vapply(needed,function(p)as.character(utils::packageVersion(p)),character(1)))))
  excel <- NULL
  list(tables=tables,models=fits,replicates=reps,covariances=covs,points=points,data=d,
       design=des,replicate_design=repdes,contrasts=contrasts,exposure_values=exposure_values,
       plots=list(),excel_file=excel)
}

mediation_ratio_summary <- function(mediation_result) {
  z <- mediation_result
  te <- "PA_median_plus100__TE__RD"
  ie <- "PA_median_plus100__TNIE__RD"
  d <- z$replicate_design
  df <- survey::degf(d)
  out <- lapply(names(z$replicates), function(nm) {
    b <- z$points[[nm]]; r <- z$replicates[[nm]]
    stopifnot(all(is.finite(r[, c(te, ie)])), all(r[, te] != 0), b[te] != 0)
    estimate <- unname(100 * b[ie] / b[te])
    paired <- 100 * r[, ie] / r[, te]
    variance <- survey::svrVar(matrix(paired, ncol = 1), scale = d$scale,
                              rscales = d$rscales, mse = d$mse, coef = estimate)
    se <- sqrt(as.numeric(variance))
    data.frame(mediator = nm, estimate = estimate, low = estimate - qt(.975, df) * se,
               high = estimate + qt(.975, df) * se, SE = se, df = df,
               replicates = nrow(r), min_abs_replicate_TE = min(abs(r[, te])),
               replicate_TE_sign_changes = sum(sign(r[, te]) != sign(b[te])))
  })
  do.call(rbind, out)
}

# 6. Standalone figures

SF_FONT <- SF_NUM_FONT <- FONT_FAMILY

SF_SIZE <- 22

SF_COL <- c(ink="#555555",ribbon="#EFF4F6",stripe="#EBF1F3",
            blue="#DCE7EB",red="#EBBCB9",mint="#C7E1D7",median="#C6D0D4")

COMPONENT_COLORS <- c(Work="#DCE7EB",Transport="#C7E1D7",Leisure="#EBBCB9",
 VPA="#EBBCB9",NonVPA="#DCE7EB",WorkVig="#AFC6CF",WorkMod="#DCE7EB",
 LeisureVig="#D69E9A",LeisureMod="#EBBCB9")

MED_COLORS<-c(TE="#DCE7EB",PNDE="#EBBCB9",TNIE="#C7E1D7")

MED_LABELS<-c(BMI="BMI",Waist="Waist circumference",Hypertension="Hypertension",
 Diabetes="Diabetes",Hyperlipidemia="Hyperlipidemia",Uric_acid="Serum uric acid")

sf_p <- function(x) ifelse(is.na(x), "", ifelse(x<.001,"<0.001",sprintf("%.3f",x)))

sf_p_label <- function(x, label="P") paste0(label, ifelse(x<.001," < 0.001",paste0(" = ",sprintf("%.3f",x))))

sf_theme <- function(grid=TRUE) theme_classic(base_size=SF_SIZE, base_family=SF_FONT)+
  theme(axis.text=element_text(colour="black",size=FIGURE_LABEL_SIZE), axis.line=element_line(linewidth=.65),
    axis.ticks=element_line(linewidth=.6),
    panel.grid.major.y=if(grid) element_line(colour="grey94",linewidth=.4) else element_blank(),
    legend.position="bottom", legend.title=element_blank(),
    legend.text=element_text(size=18), plot.margin=margin(14,20,14,14),
    plot.background=element_rect(fill="white",colour=NA),
    plot.title=element_blank(),plot.subtitle=element_blank())

sf_rcs <- function(cr, ref, minimum, minimum_or, p1, pnl, xlabel,
                   p1_label="P overall", ref_label="Weighted median", units="", digits=2,
                   boundary=FALSE, limits=NULL) {
  stopifnot(all(is.finite(as.matrix(cr[c("value","OR","OR_low","OR_high")]))),
    all(cr$OR_low<=cr$OR+1e-10),all(cr$OR<=cr$OR_high+1e-10))
  if(is.null(limits)) limits <- c(min(c(cr$OR_low,1))*.85,max(c(cr$OR_high,1))*1.52)
  yr <- diff(limits)
  num <- function(z) paste0(formatC(z,format="f",digits=digits),if(nzchar(units))paste0(" ",units)else "")
  minlab <- if(boundary) "Minimum (boundary)" else "Minimum"
  p <- ggplot(cr,aes(value,OR))+
    geom_ribbon(aes(ymin=OR_low,ymax=OR_high),fill=SF_COL["ribbon"],colour=NA)+
    geom_hline(yintercept=1,linetype="dashed",colour="grey65",linewidth=.7)+
    geom_vline(xintercept=ref,linetype="dashed",colour=SF_COL["median"],linewidth=.8)+
    {if(!isTRUE(boundary)) geom_vline(xintercept=minimum,linetype="dashed",colour="grey35",linewidth=.8)}+
    geom_line(colour="black",linewidth=.9)+
    annotate("point",x=ref,y=1,colour=SF_COL["median"],size=2.5)+
    {if(!isTRUE(boundary)) annotate("point",x=minimum,y=minimum_or,colour="black",size=2.5)}+
    # 原图规则：文字左缘对应各自竖线，向右起排；两条说明上下错开。
    annotate("text",x=ref,y=limits[2]-.27*yr,label=paste0(ref_label,": ",num(ref)),
      hjust=0,vjust=1,family=SF_NUM_FONT,size=RCS_ANNOTATION_SIZE/.pt)+
    {if(!isTRUE(boundary)) annotate("text",x=minimum,y=limits[2]-.35*yr,label=paste0(minlab,": ",num(minimum)),
      hjust=0,vjust=1,family=SF_NUM_FONT,size=RCS_ANNOTATION_SIZE/.pt)}+
    annotate("text",x=Inf,y=Inf,label=paste(sf_p_label(p1,p1_label),sf_p_label(pnl,"P non-linear"),sep="\n"),
      hjust=1.05,vjust=1.25,family=SF_NUM_FONT,size=RCS_ANNOTATION_SIZE/.pt,lineheight=1.1)+
    scale_y_continuous(expand=expansion(mult=c(0,0)))+
    coord_cartesian(ylim=limits)+labs(x=xlabel,y="Odds ratio (95% CI)")+sf_theme()
  p
}

sf_forest <- function(d, labels, xlabel="Odds ratio", null=1, digits=2,
                      estimate="OR", low="OR_low", high="OR_high", effect_title="OR (95% CI)", extra=NULL, ink=SF_COL["ink"], stripe=SF_COL["stripe"],
                      p_column="P", p_header="P value", groups=NULL) {
  stopifnot(length(labels)==nrow(d),all(is.finite(d[[estimate]])))
  z<-d; z$.y<-rev(seq_len(nrow(z))); z$.label<-labels
  group_labels<-NULL
  if(!is.null(groups)) {
    stopifnot(length(groups)==nrow(z),!anyNA(groups))
    groups<-as.character(groups)
    starts<-which(c(TRUE,groups[-1L]!=head(groups,-1L)))
    group_id<-cumsum(seq_len(nrow(z))%in%starts)
    z$.y<-z$.y+length(starts)-group_id
    group_labels<-data.frame(.y=z$.y[starts]+.85,.label=groups[starts])
  }
  z$.est<-z[[estimate]];z$.lo<-z[[low]];z$.hi<-z[[high]]
  ref<-!is.finite(z$.lo)|!is.finite(z$.hi)
  if(all(c("category","reference")%in%names(z))) ref<-ref | as.character(z$category)==as.character(z$reference)
  if(any(ref & abs(z$.est-null)>1e-10)) stop("A non-reference CI is missing")
  z$.lo[ref]<-null;z$.hi[ref]<-null
  r<-range(c(z$.lo,z$.hi,null)); pad<-max(diff(r)*.12,if(null==1).025 else .02)
  r<-r+c(-pad,pad); if(null==1) r[1]<-max(0,r[1])
  tr<-function(x)(x-r[1])/diff(r)
  z$.x<-tr(z$.est);z$.xl<-tr(z$.lo);z$.xh<-tr(z$.hi)
  z$.value<-sprintf(paste0("%.",digits,"f (%.",digits,"f-%.",digits,"f)"),z$.est,z$.lo,z$.hi)
  if(null==0) z$.value<-sprintf(paste0("%.",digits,"f (%.",digits,"f, %.",digits,"f)"),z$.est,z$.lo,z$.hi)
  z$.value[ref]<-sprintf(paste0("%.",digits,"f (Reference)"),null)
  z$.p<-if(p_column%in%names(z))sf_p(z[[p_column]])else rep("",nrow(z));z$.p[ref]<-""
  ticks<-pretty(r,4);ticks<-ticks[ticks>=r[1]&ticks<=r[2]]
  td<-data.frame(x=tr(ticks),label=format(ticks,trim=TRUE,scientific=FALSE))
  n<-if(is.null(groups))nrow(z)else max(z$.y)+1
  bg<-z[seq_len(nrow(z))%%2==1,,drop=FALSE]
  label_x<-if(is.null(groups))-1.03 else -.90
  null_top<-if(is.null(groups))n+.6 else max(z$.y)+.46
  p<-ggplot()+geom_rect(data=bg,aes(ymin=.y-.46,ymax=.y+.46),xmin=-1.07,xmax=2.72,fill=stripe,colour=NA)+
    annotate("segment",x=tr(null),xend=tr(null),y=.53,yend=null_top,linetype="dashed",colour="grey60",linewidth=.55)+
    geom_segment(data=z[!ref,,drop=FALSE],aes(x=.xl,xend=.xh,y=.y,yend=.y),linewidth=.8,colour=ink)+
    geom_segment(data=z[!ref,,drop=FALSE],aes(x=.xl,xend=.xl,y=.y-.09,yend=.y+.09),linewidth=.65,colour=ink)+
    geom_segment(data=z[!ref,,drop=FALSE],aes(x=.xh,xend=.xh,y=.y-.09,yend=.y+.09),linewidth=.65,colour=ink)+
    geom_point(data=z,aes(.x,.y),size=2.7,colour=ink)+
    geom_text(data=z,aes(y=.y,label=.label),x=label_x,hjust=0,family=SF_FONT,size=FIGURE_LABEL_SIZE/.pt)+
    geom_text(data=z,aes(x=1.14,y=.y,label=.value),hjust=0,family=SF_NUM_FONT,size=FIGURE_LABEL_SIZE/.pt)+
    geom_text(data=z,aes(x=2.68,y=.y,label=.p),hjust=1,family=SF_NUM_FONT,size=FIGURE_LABEL_SIZE/.pt)+
    annotate("text",x=1.14,y=n+.85,label=effect_title,hjust=0,family=SF_NUM_FONT,size=FIGURE_LABEL_SIZE/.pt)+
    annotate("text",x=2.68,y=n+.85,label=if(p_column%in%names(z))p_header else "",hjust=1,family=SF_NUM_FONT,size=FIGURE_LABEL_SIZE/.pt)+
    annotate("segment",x=0,xend=1,y=.35,yend=.35,linewidth=.7)+
    geom_segment(data=td,aes(x=x,xend=x),y=.35,yend=.22,linewidth=.6)+
    geom_text(data=td,aes(x=x,label=label),y=-.05,family=SF_FONT,size=17/.pt)+
    annotate("text",x=.5,y=-.63,label=xlabel,family=SF_FONT,size=20/.pt)+
    coord_cartesian(xlim=c(-1.1,2.75),ylim=c(-1.12,n+1.15),expand=FALSE,clip="off")+
    theme_void(base_family=SF_FONT)+theme(plot.margin=margin(12,16,12,12),plot.background=element_rect(fill="white",colour=NA))
  if(!is.null(group_labels)) p<-p+geom_text(data=group_labels,
    aes(y=.y,label=.label),x=-1.03,hjust=0,fontface="bold",family=SF_FONT,size=FIGURE_LABEL_SIZE/.pt)
  if(!is.null(extra)) p<-p+labs(caption=extra)+theme(plot.caption=element_text(size=13,hjust=0,margin=margin(t=10)))
  p
}

plot_pa_st_categories <- function(tab) {
  # 共用OR坐标轴，将两个已拟合模型的分类结果排入同一张森林图。
  # PA以P1、ST以S1为参考；不重新估计、合并或平均效应量。
  needed<-c("model","exposure","category","reference","OR","OR_low","OR_high","P")
  stopifnot(is.data.frame(tab),all(needed%in%names(tab)))
  z<-tab[tab$model=="M3" & tab$exposure%in%c("PA","ST"),,drop=FALSE]
  wanted<-c(paste("PA",paste0("P",1:4)),paste("ST",paste0("S",1:4)))
  key<-paste(z$exposure,z$category)
  stopifnot(nrow(z)==8L,!anyDuplicated(key),setequal(key,wanted))
  z<-z[match(wanted,key),,drop=FALSE]
  stopifnot(all(as.character(z$reference)==rep(c("P1","S1"),each=4)))
  style_final_plot(sf_forest(z,as.character(z$category),groups=z$exposure))
}

sf_composition <- function(d, levels, labels, colors) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("请先安装ggplot2。")
  groups <- paste0("P", 1:4)
  required <- c("PA_group", "component", "volume_share_pct")
  stopifnot(all(required %in% names(d)), length(levels) == length(labels),
    all(levels %in% names(colors)))
  z <- d
  z$PA_group <- as.character(z$PA_group)
  z$component <- as.character(z$component)
  if (nrow(z) != length(groups) * length(levels) ||
      anyDuplicated(paste(z$PA_group, z$component)) ||
      !setequal(z$PA_group, groups) || !setequal(z$component, levels) ||
      any(!is.finite(z$volume_share_pct)) || any(z$volume_share_pct < 0) ||
      any(abs(tapply(z$volume_share_pct, z$PA_group, sum) - 100) > 1e-6)) {
    stop("组成数据应包含P1-P4的全部成分，各组未取整百分比之和应为100%。")
  }
  z <- z[order(match(z$PA_group, groups), match(z$component, levels)), , drop = FALSE]
  # 色块从左至右的顺序与数字区从上至下的成分顺序一致。
  z$.xmax <- ave(z$volume_share_pct, z$PA_group, FUN = cumsum)
  z$.xmin <- z$.xmax - z$volume_share_pct
  z$.bar_y <- 5 - match(z$PA_group, groups)
  z$.table_x <- c(40, 60, 80, 100)[match(z$PA_group, groups)]
  z$.table_y <- -2.15 - .88 * (match(z$component, levels) - 1)
  z$.percent <- sprintf("%.1f%%", z$volume_share_pct)
  z$component <- factor(z$component, levels = levels)
  key <- data.frame(component = factor(levels, levels = levels),
    label = vapply(as.character(labels), function(s)
      paste(strwrap(s, width = 23), collapse = "\n"), character(1)),
    y = -2.15 - .88 * (seq_along(levels) - 1))
  group_key <- data.frame(label = groups, y = 4:1, x = c(40, 60, 80, 100))
  ticks <- data.frame(x = seq(0, 100, 25))
  font <- getOption("nhanes.composition.font", "Arial")
  pt <- 72.27 / 25.4
  bottom <- min(key$y) - .53
  p <- ggplot2::ggplot(z) +
    ggplot2::geom_segment(data = ticks,
      ggplot2::aes(x = x, xend = x), y = .59, yend = 4.42,
      inherit.aes = FALSE, colour = "#EEEEEE", linewidth = .4) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .xmin, xmax = .xmax,
      ymin = .bar_y - .31, ymax = .bar_y + .31, fill = component),
      colour = "white", linewidth = .55, show.legend = FALSE) +
    ggplot2::geom_text(data = group_key, ggplot2::aes(y = y, label = label),
      x = -3, hjust = 1, inherit.aes = FALSE, family = font, size = 18 / pt) +
    ggplot2::annotate("segment", x = 0, xend = 100, y = .53, yend = .53,
      colour = "#555555", linewidth = .5) +
    ggplot2::geom_segment(data = ticks, ggplot2::aes(x = x, xend = x),
      y = .53, yend = .44, inherit.aes = FALSE, colour = "#555555", linewidth = .5) +
    ggplot2::geom_text(data = ticks, ggplot2::aes(x = x, label = paste0(x, "%")),
      y = .13, inherit.aes = FALSE, family = font, size = 14 / pt) +
    ggplot2::annotate("text", x = 50, y = -.39,
      label = "Contribution to total PA", family = font, size = 16 / pt) +
    ggplot2::annotate("text", x = -6, y = -1.29, label = "Activity component",
      hjust = 0, family = font, size = 15 / pt) +
    ggplot2::geom_text(data = group_key, ggplot2::aes(x = x, label = label),
      y = -1.29, inherit.aes = FALSE, family = font, size = 18 / pt) +
    ggplot2::annotate("segment", x = -10, xend = 106, y = -1.67, yend = -1.67,
      colour = "#B8B8B8", linewidth = .45) +
    ggplot2::geom_point(data = key, ggplot2::aes(y = y, fill = component),
      x = -10, shape = 22, size = 3.5, stroke = 0, inherit.aes = FALSE,
      show.legend = FALSE) +
    ggplot2::geom_text(data = key, ggplot2::aes(y = y, label = label),
      x = -6, hjust = 0, lineheight = .98, inherit.aes = FALSE,
      family = font, size = 15 / pt) +
    ggplot2::geom_text(ggplot2::aes(x = .table_x, y = .table_y, label = .percent),
      family = font, size = 16 / pt) +
    ggplot2::scale_fill_manual(values = colors, breaks = levels, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-13, 108), ylim = c(bottom, 4.55),
      expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_family = font) +
    ggplot2::theme(aspect.ratio = .55 + .10 * length(levels),
      legend.position = "none", plot.margin = ggplot2::margin(12, 10, 8, 10),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  attr(p, "composition_layout") <- "horizontal_100_percent_with_values"
  p
}

style_final_plot <- function(p) {
  p <- p + theme(text=element_text(family=FONT_FAMILY),plot.title=element_blank(),
    plot.subtitle=element_blank(),plot.caption=element_blank())
  for(i in seq_along(p$layers)) {
    if(inherits(p$layers[[i]]$geom,c("GeomText","GeomLabel")))
      p$layers[[i]]$aes_params$family <- FONT_FAMILY
    lab <- p$layers[[i]]$aes_params$label
    if(is.character(lab) && length(lab)==1L) {
      if(grepl("^P (overall|positive dose) ",lab)) {
        a <- strsplit(lab,"\n",fixed=TRUE)[[1]]
        a <- vapply(a,function(v) {
          v <- sub("^P ([A-Za-z -]+) ([<=]) ","italic(P)~\"\\1\"\\2",v)
          sub('"=', '"==',v,fixed=TRUE)
        },character(1))
        p$layers[[i]]$aes_params$label <- paste0("atop(",paste(a,collapse=","),")")
        p$layers[[i]]$geom_params$parse <- TRUE
      } else if(lab %in% c("P value", "q value")) {
        p$layers[[i]]$aes_params$label <- if(lab=="q value") 'italic(q)~value' else 'italic(P)~value'
        p$layers[[i]]$geom_params$parse <- TRUE
      }
    }
  }
  p
}

make_rcs_plot <- function(record) {
  sm<-record$summary
  style_final_plot(sf_rcs(record$curve,sm$median_reference,sm$minimum_value,sm$minimum_OR,
    sm$P_overall,sm$P_nonlinear,
    if(sm$exposure=="ASBI")"ASBI" else paste0(sm$exposure," (",sm$units,")"),
    units=if(sm$exposure=="ASBI")"" else sm$units,
    digits=if(sm$exposure=="ASBI")2 else 0,boundary=isTRUE(sm$minimum_at_boundary)))
}

make_joint_plot <- function(tab) {
  z<-joint_display_table(tab[tab$model=="M3",]);z$x<-as.integer(substr(z$category,2,2));z$y<-as.integer(substr(z$category,5,5))
  ref<-z$category=="P4&S1"
  stopifnot(sum(!ref)==15L, all(is.finite(z$P[!ref])))
  # Fig. 3c uses nominal P values; Fig. 3d uses BH q values for the same 15 comparisons.
  # ORs and pointwise CIs are the unchanged covariate-adjusted model estimates.
  z$or_label<-sprintf("%.2f",z$OR)
  z$ci_label<-ifelse(ref,"(Reference)",sprintf("(%.2f-%.2f)",z$OR_low,z$OR_high))
  z$p_label<-ifelse(z$P<.001,'italic(P) < 0.001',sprintf('italic(P) == "%.3f"',z$P))
  z$or_y<-z$y+ifelse(ref,.12,.21);z$ci_y<-z$y-ifelse(ref,.12,0)
  style_final_plot(ggplot(z,aes(x,y,fill=OR))+geom_tile(colour="white",linewidth=1.5)+
    geom_text(aes(y=or_y,label=or_label),family=FONT_FAMILY,size=18/.pt)+
    geom_text(aes(y=ci_y,label=ci_label),family=FONT_FAMILY,size=17/.pt)+
    geom_text(data=z[!ref,],aes(y=y-.21,label=p_label),parse=TRUE,family=FONT_FAMILY,size=17/.pt)+
    scale_fill_gradient(low="#DCE7EB",high="white",name="OR")+
    scale_x_continuous(breaks=1:4,labels=paste0("P",1:4))+
    scale_y_continuous(breaks=1:4,labels=paste0("S",1:4))+
    coord_cartesian(xlim=c(.5,4.5),ylim=c(.5,4.5),expand=FALSE)+labs(x="PA",y="ST")+
    sf_theme(FALSE)+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position="right"))
}

plot_asbi_distribution <- function(a) {
  h<-a$histogram;z<-a$categories;cuts<-a$cutpoints$weighted_cutpoint
  style_final_plot(ggplot(h,aes(center,pct))+
    geom_col(width=.495,fill="#EFF4F6",colour="white",linewidth=.25)+
    geom_vline(xintercept=cuts,linetype="dashed",colour="grey40",linewidth=.7)+
    geom_line(data=z,aes(median,crude_prevalence_pct),inherit.aes=FALSE,colour="#B7D4E1",linewidth=.9)+
    geom_point(data=z,aes(median,crude_prevalence_pct),inherit.aes=FALSE,colour="#B7D4E1",size=3)+
    geom_text(data=z,aes(median,crude_prevalence_pct,label=sprintf("%.1f%%",crude_prevalence_pct)),
      inherit.aes=FALSE,vjust=-.8,family=FONT_FAMILY,colour="#B7D4E1",size=RCS_ANNOTATION_SIZE/.pt)+
    scale_x_continuous(breaks=cuts,labels=sprintf("%.2f",cuts))+
    scale_y_continuous(breaks=seq(0,12,3),limits=c(0,max(12.5,h$pct*1.12)),
      sec.axis=sec_axis(~.,name="Weighted KSD prevalence (%)"))+
    labs(x="ASBI",y="Weighted ASBI distribution (%)")+sf_theme()+
    theme(axis.text.y.right=element_text(colour="#9CBCCB"),axis.title.y.right=element_text(colour="#9CBCCB"),
      axis.line.y.right=element_line(colour="#9CBCCB")))
}

plot_joint_share <- function(tab) {
  z<-joint_display_table(tab[tab$model=="M3",]);z$x<-as.integer(substr(z$category,2,2));z$y<-as.integer(substr(z$category,5,5))
  mx<-aggregate(weighted_population_pct~x,z,sum);my<-aggregate(weighted_population_pct~y,z,sum)
  # 边际柱从相同的零基线起绘，柱长按加权比例线性缩放，两个方向使用同一比例尺。
  bar_base<-4.58
  bar_scale<-.65/max(c(mx$weighted_population_pct,my$weighted_population_pct))
  mx$tip<-bar_base+mx$weighted_population_pct*bar_scale
  my$tip<-bar_base+my$weighted_population_pct*bar_scale
  style_final_plot(ggplot(z,aes(x,y,fill=weighted_population_pct))+geom_tile(colour="white",linewidth=1.5)+
    geom_text(aes(label=sprintf("%.1f%%",weighted_population_pct)),size=FIGURE_LABEL_SIZE/.pt,family=FONT_FAMILY)+
    geom_rect(data=mx,aes(xmin=x-.45,xmax=x+.45,ymin=bar_base,ymax=tip),
      inherit.aes=FALSE,fill="#DCE7EB",colour="white")+
    geom_text(data=mx,aes(x,y=tip+.06,label=sprintf("%.1f",weighted_population_pct)),
      inherit.aes=FALSE,vjust=0,size=FIGURE_LABEL_SIZE/.pt,family=FONT_FAMILY)+
    geom_rect(data=my,aes(xmin=bar_base,xmax=tip,ymin=y-.45,ymax=y+.45),
      inherit.aes=FALSE,fill="#DCE7EB",colour="white")+
    geom_text(data=my,aes(x=tip+.07,y=y,label=sprintf("%.1f",weighted_population_pct)),
      inherit.aes=FALSE,hjust=0,size=FIGURE_LABEL_SIZE/.pt,family=FONT_FAMILY)+
    scale_fill_gradient2(low="#DCE7EB",mid="white",high="#EBBCB9",midpoint=6.25,name=NULL,labels=function(x)paste0(x,"%"))+
    scale_x_continuous(breaks=1:4,labels=paste0("P",1:4))+scale_y_continuous(breaks=1:4,labels=paste0("S",1:4))+
    coord_cartesian(xlim=c(.5,5.85),ylim=c(.5,5.4),expand=FALSE)+labs(x="PA",y="ST")+sf_theme(FALSE)+
    theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position="right"))
}

plot_joint_prevalence <- function(z) {
  z<-joint_display_table(z)
  z$category<-factor(z$category,levels=z$category[order(z$adjusted_prevalence_pct,decreasing=TRUE)])
  # 右侧恢复每组调整患病率，保留两位小数；为文字列预留空间，避免与CI重叠。
  right_limit<-max(19,1.4*max(z$adjusted_high_pct,na.rm=TRUE))
  label_x<-right_limit*.985
  style_final_plot(ggplot(z,aes(adjusted_prevalence_pct,category,fill=adjusted_prevalence_pct))+
    geom_col(width=.72)+geom_errorbar(aes(xmin=adjusted_low_pct,xmax=adjusted_high_pct),orientation="y",width=.2,linewidth=.3,colour="grey35")+
    geom_text(aes(label=sprintf("%.2f%%",adjusted_prevalence_pct)),x=label_x,
      hjust=1,family=FONT_FAMILY,size=FIGURE_LABEL_SIZE/.pt)+
    scale_x_continuous(limits=c(0,right_limit),breaks=pretty(c(0,right_limit),4),expand=expansion(mult=c(0,0)))+
    scale_fill_gradient(low="#DCE7EB",high="#EBBCB9",guide="none")+labs(x="Adjusted KSD prevalence (%)",y=NULL)+sf_theme(FALSE))
}

plot_mediation_path <- function(z,nm) {
  effect_labels<-c(TE="Total effect (TE)",PNDE="Direct effect (PNDE)",TNIE="Indirect effect (TNIE)")
  ef<-function(e) {r<-z[z$effect==e,];sprintf("%s\n%.3f (%.3f, %.3f)",effect_labels[e],r$difference_pp,r$low_pp,r$high_pp)}
  nodes<-data.frame(x=c(1.5,5,8.5),y=c(2.3,4.4,2.3),label=c("PA",MED_LABELS[nm],"KSD"),fill=unname(MED_COLORS[c("TE","TNIE","PNDE")]))
  nodes$half_width<-c(1.2,1.7,1.2) # 中介名称较长，单独为中间框保留文字空间。
  style_final_plot(ggplot()+
    annotate("segment",x=2.25,xend=4,y=2.72,yend=4.02,colour="grey35",linewidth=.8,arrow=grid::arrow(length=grid::unit(2,"mm"),type="closed"))+
    annotate("segment",x=6,xend=7.75,y=4.02,yend=2.72,colour="grey35",linewidth=.8,arrow=grid::arrow(length=grid::unit(2,"mm"),type="closed"))+
    annotate("segment",x=2.7,xend=7.3,y=2.3,yend=2.3,colour="grey35",linewidth=.8,arrow=grid::arrow(length=grid::unit(2,"mm"),type="closed"))+
    lapply(seq_len(nrow(nodes)),function(i) {
      b<-nodes[i,]
      annotation_custom(grid::roundrectGrob(r=grid::unit(.12,"snpc"),
        gp=grid::gpar(fill=b$fill,col="grey45",lwd=.5*72.27/25.4)),
        xmin=b$x-b$half_width,xmax=b$x+b$half_width,ymin=b$y-.42,ymax=b$y+.42)
    })+
    geom_text(data=nodes,aes(x,y,label=label),family=FONT_FAMILY,size=18/.pt)+
    annotate("text",x=5,y=5.3,label=ef("TNIE"),size=FIGURE_LABEL_SIZE/.pt,lineheight=1.1,family=FONT_FAMILY)+
    annotate("text",x=5,y=1.5,label=ef("PNDE"),size=FIGURE_LABEL_SIZE/.pt,lineheight=1.1,family=FONT_FAMILY)+
    annotate("text",x=5,y=.5,label=ef("TE"),size=FIGURE_LABEL_SIZE/.pt,lineheight=1.1,family=FONT_FAMILY)+
    coord_cartesian(xlim=c(0,10),ylim=c(0,5.8),expand=FALSE)+theme_void(base_family=FONT_FAMILY))
}

plot_mediation_summary <- function(z) {
  z$effect<-factor(z$effect,levels=c("TE","PNDE","TNIE"))
  z$mediator<-factor(z$mediator,levels=rev(names(MED_LABELS)))
  style_final_plot(ggplot(z,aes(difference_pp,mediator,colour=effect))+
    geom_vline(xintercept=0,linetype="dashed",colour="grey60")+
    geom_errorbar(aes(xmin=low_pp,xmax=high_pp),orientation="y",width=.15,position=position_dodge(width=.65),linewidth=.8)+
    geom_point(position=position_dodge(width=.65),size=3)+
    scale_colour_manual(values=c(TE="#94B2BF",PNDE="#D39590",TNIE="#8AB5A2"))+
    scale_y_discrete(labels=MED_LABELS)+labs(x="Probability difference (percentage points)",y=NULL,colour=NULL)+sf_theme(FALSE))
}

# 7. Variable definitions and manuscript outputs

variable_definitions <- function() data.frame(Variable=c("Outcome","PA components","PA categories","ST","ASBI","Age and sex","Race and ethnicity","Education and marital status","Poverty-to-income ratio","BMI and waist","Diabetes","Hypertension","Hyperlipidemia","Smoking","Alcohol","Plain-water intake","Uric acid","Exclusions","Weights and missingness"), Definition=c("KIQ026: yes=1; no=0; other/unknown responses are missing. KSD is defined by a self-reported history of kidney stones.","Vigorous work: PAQ605/610/PAD615; moderate work: PAQ620/625/PAD630; transport: PAQ635/640/PAD645; vigorous leisure: PAQ650/655/PAD660; moderate leisure: PAQ665/670/PAD675. Weekly component volume = days × minutes × MET (8 vigorous; 4 moderate/transport); divide by 7 for daily volume. Explicit nonparticipation is zero. Unknown participation and duration codes 7777/9999 produce missing total PA.","P1 <20; P2 20 to <200; P3 200 to <600; P4 ≥600 MET-min/day. Non-vigorous PA = moderate work + transport + moderate leisure.","PAD680 typical sitting time. Values <2000 are retained; other values are missing. S1 <240; S2 240 to <360; S3 360 to <480; S4 ≥480 min/day.","ln[(PA/60+c)/(ST/60+c)], using PA in MET-min/day and ST in min/day. The transformed inputs are numerical values in MET-h/day and h/day, respectively. ASBI is reported without a physical unit. Main c=0.1, with c=0.01 and 0.5 in sensitivity analyses. At c=0.1, survey-weighted quartiles rounded to two decimal places defined Q1 <-2.14, Q2 -2.14 to <-0.47, Q3 -0.47 to <0.91, and Q4 ≥0.91. Q3 was the reference.","Age groups: 20-39, 40-59, and ≥60 years. Sex: female/male.","Mexican American; non-Hispanic Black; non-Hispanic White; other Hispanic; other race, using RIDRETH1.","Education: 9-11th grade, college or above, high school graduate, less than 9th grade, and some college. Marital status: divorced, living with partner, married, never married, separated, and widowed. Unmapped responses are missing.","<1; 1 to <3.5; ≥3.5.","Measured BMI (kg/m²) and waist circumference (cm), entered continuously in the main models.","Positive if DIQ010=1, DIQ050=1, DIQ070=1, fasting glucose ≥126 mg/dL, or HbA1c ≥6.5%; zero otherwise, including records without a positive available component.","Positive if BPQ020=1 or BPQ040A=1; zero otherwise, including records without a positive available component.","Positive if BPQ080=1, BPQ090D=1, total cholesterol ≥240 mg/dL, triglycerides ≥200 mg/dL, LDL cholesterol ≥160 mg/dL, or HDL cholesterol <40 mg/dL in men or <50 mg/dL in women; zero otherwise. Available fasting-subsample values contribute to this composite.","SMQ020=1 (≥100 lifetime cigarettes) coded 1; other responses coded 0.","2007-2016: ALQ101=1 coded 1, otherwise 0. 2017-2018: ALQ111=1 and ALQ121 in 1-7 coded 1; ALQ121=77/99 coded missing; other responses coded 0. This is the original binary project coding.","Positive DR1_320Z and DR2_320Z values are averaged; if only one is positive that value is retained; otherwise missing. Categories: <550, 550 to <1000, 1000 to <1650, ≥1650 g/day.","Positive serum uric acid values retained. Main-model categories: <4.5, 4.5 to <5.4, 5.4 to <6.3, ≥6.3 mg/dL. Continuous log uric acid is used in its mediation model.","KIQ022=1 indicates reported weak/failing kidneys; MCQ220=1 indicates cancer. Only participants with both responses coded no were retained. Unknown exclusion status was excluded.","WTMEC2YR/6 with SDMVSTRA and SDMVPSU. Complete-case selection applies to the final coded variables. Composite indicators use the positive-evidence/default-zero rules described above; missing component records are not separately imputed."),check.names=FALSE)

run_nhanes_manuscript <- function(analysis_data) {
  validate_analysis_input(analysis_data)
  check_packages()
  suppressPackageStartupMessages(library(ggplot2))
  old<-options(survey.lonely.psu="adjust",contrasts=c("contr.treatment","contr.poly"))
  on.exit(options(old),add=TRUE)
  started<-Sys.time()
  loaded <- list(data=analysis_data)
  message("[1/5] 正文：分类、联合分组、RCS、阈值")
  main<-run_nhanes_main(loaded);d<-main$data;ds<-main$design
  baseline_tests<-final_baseline_tests(d,ds)
  dist<-final_asbi_distribution(d,ds)
  message("[2/5] 加权 ST spline × PA 交互及边际标准化")
  interaction<-run_final_interaction(d)
  message("[3/5] ASBI 常数及排除 P1&S1 敏感性")
  sensitivity<-run_final_sensitivity(d,ds)
  message("[4/5] 三个完全调整活动成分模型")
  components<-run_nhanes_three_schemes_final(list(data=d))
  message("[5/5] 六个连续 PA 中介模型及 JKn 方差")
  mediation<-run_nhanes_continuous_PA_mediation(list(data=d))
  mediation$tables$Signed_indirect_total<-mediation_ratio_summary(mediation)
  stopifnot(identical(as.character(d$SEQN),as.character(components$data$SEQN)),
    identical(as.character(d$SEQN),as.character(mediation$data$SEQN)))
  main$tables$Baseline_tests<-baseline_tests
  main$tables$ASBI_cutpoints<-dist$cutpoints
  main$tables$ASBI_histogram<-dist$histogram
  main$tables$ASBI_group_prevalence<-dist$categories
  main$tables$Interaction_spline_test<-interaction$tables$Test
  main$tables$Interaction_spline_curves<-interaction$tables$Curves
  main$tables$Interaction_spline_knots<-interaction$tables$Knots
  # 发表表格与完整数值结果同时保留。
  tab2<-main$tables$Categorical_OR
  k<-unique(tab2[c("exposure","category","n","cases")])
  for(m in c("M1","M2","M3")) {
    a<-tab2[tab2$model==m,];idx<-match(paste(k$exposure,k$category),paste(a$exposure,a$category))
    k[[paste0(m,"_OR_95CI")]]<-a$OR_95CI[idx]
    k[[paste0(m,"_P")]]<-a$P[idx]
  }
  s5<-mediation$tables$All_effects
  ratio<-mediation$tables$Signed_indirect_total
  s5$Indirect_total_pct<-s5$Ratio_low<-s5$Ratio_high<-NA_real_
  ix<-s5$effect=="TNIE";j<-match(s5$mediator[ix],ratio$mediator)
  s5$Indirect_total_pct[ix]<-ratio$estimate[j];s5$Ratio_low[ix]<-ratio$low[j];s5$Ratio_high[ix]<-ratio$high[j]
  pub_main<-list(Baseline=final_baseline_wide(main$tables$Baseline,baseline_tests),Categorical_OR=k)
  pub_supp<-list(Variable_definitions=variable_definitions(),Threshold=main$tables$Threshold_results,
    Joint16_OR=main$tables$Joint16_OR[main$tables$Joint16_OR$model=="M3",],
    Components_OR=components$tables$Linear_OR,Mediation=s5)
  plots<-list(p_PA_RCS=main$plots$PA,p_ST_RCS=main$plots$ST,p_ASBI_RCS=main$plots$ASBI)
  for(ex in "ASBI") {
    a<-tab2[tab2$model=="M3" & tab2$exposure==ex,]
    plots[[paste0("p_",ex,"_categories")]]<-style_final_plot(sf_forest(a,a$category))
  }
  plots$p_PA_ST_categories<-plot_pa_st_categories(tab2)
  plots$p_Fig2d_interaction<-interaction$plot
  plots$p_Fig3a_joint_share<-plot_joint_share(main$tables$Joint16_OR)
  plots$p_Fig3b_joint_prevalence<-plot_joint_prevalence(main$tables$Joint_adjusted_prevalence)
  plots$p_Fig3c_joint_OR<-main$plots$Joint
  a<-main$tables$Joint16_OR;a<-a[a$model=="M3",]
  # 参照组置顶，其余15组按未取整OR从小到大排列；保留原CI和BH校正q值。
  a<-a[order(as.character(a$category)!=as.character(a$reference),a$OR,as.character(a$category)),,drop=FALSE]
  plots$p_Fig3d_joint_forest<-style_final_plot(sf_forest(a,joint_category_label(a$category),p_column="q_BH",p_header="q value"))
  plots$p_Fig4a_ASBI_distribution<-plot_asbi_distribution(dist)
  for(id in names(sensitivity$plots))plots[[paste0("p_S_ASBI_",id)]]<-sensitivity$plots[[id]]
  for(fam in names(PC_FAMILIES)) {
    z<-components$tables$Composition;z<-z[z$family==fam,];lv<-PC_FAMILIES[[fam]]
    plots[[paste0("p_S_composition_",fam)]]<-style_final_plot(sf_composition(z,lv,PC_LABELS[lv],COMPONENT_COLORS))
  }
  for(nm in names(MED_LABELS))plots[[paste0("p_S_mediation_",nm)]]<-plot_mediation_path(s5[s5$mediator==nm,],nm)
  plots$p_S_mediation_summary<-plot_mediation_summary(s5)
  index <- data.frame(object=names(plots),location=ifelse(grepl("^p_S_",names(plots)),
    "Supplementary","Main"),display=paste0("print(",names(plots),")"),font=FONT_FAMILY)
  stopifnot(all(vapply(plots,inherits,logical(1),what="ggplot")))
  result <- list(data=d,design=ds,main=main,sensitivity=sensitivity,components=components,
    mediation=mediation,interaction=interaction,
    publication=list(Main=pub_main,Supplementary=pub_supp),plots=plots,figure_index=index,
    elapsed_minutes=as.numeric(difftime(Sys.time(),started,units="mins")))
  message("全部完成，耗时 ",round(result$elapsed_minutes,1)," 分钟。单图数：",length(plots),"。未保存图片。")
  result
}

export_manuscript_tables <- function(result,directory="tables") {
  if(!requireNamespace("openxlsx",quietly=TRUE))stop("请先安装 openxlsx")
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  stamp<-format(Sys.time(),"%Y%m%d_%H%M%S")
  groups <- result$publication
  paths<-character()
  for(g in names(groups)) {
    tab<-groups[[g]];tab<-tab[vapply(tab,is.data.frame,logical(1))]
    wb<-openxlsx::createWorkbook();openxlsx::modifyBaseFont(wb,fontName=FONT_FAMILY,fontSize=11)
    # A sheet index retains the full name when Excel's 31-character limit applies.
    nm<-make.unique(substr(names(tab),1,27),sep="_")
    map<-data.frame(sheet=nm,table=names(tab))
    openxlsx::addWorksheet(wb,"Index");openxlsx::writeData(wb,"Index",map)
    wb$worksheets[[1]]$dimension<-sprintf('<dimension ref="A1:B%d"/>',nrow(map)+1L)
    hs<-openxlsx::createStyle(textDecoration="bold",halign="center",border="bottom")
    for(i in seq_along(tab)) {
      openxlsx::addWorksheet(wb,nm[i]);z<-joint_display_table(tab[[i]])
      if(!ncol(z))z<-data.frame(status="No records / not applicable")
      openxlsx::writeData(wb,nm[i],z,headerStyle=hs)
      wb$worksheets[[i+1L]]$dimension<-sprintf('<dimension ref="A1:%s%d"/>',openxlsx::int2col(ncol(z)),nrow(z)+1L)
      openxlsx::freezePane(wb,nm[i],firstRow=TRUE)
      openxlsx::setColWidths(wb,nm[i],seq_len(ncol(z)),widths=20)
    }
    path<-file.path(directory,paste0("NHANES_",g,"_",stamp,".xlsx"))
    if(file.exists(path))stop("目标文件已存在，未覆盖：",path)
    # 数值工作簿不含图片。移除部分openxlsx版本创建的空drawing引用，
    # 避免Excel/其他读取器寻找不存在的drawing.xml。
    for(j in seq_along(wb$worksheets)) if(!length(wb$drawings[[j]])) {
      wb$worksheets[[j]]$drawing<-character(0)
      wb$worksheets_rels[[j]]<-wb$worksheets_rels[[j]][!grepl('/drawing"',wb$worksheets_rels[[j]],fixed=TRUE)]
    }
    openxlsx::saveWorkbook(wb,path,overwrite=FALSE);paths[g]<-normalizePath(path,winslash="/")
  }
  print(paths);invisible(paths)
}
