# Analysis 1.1

## parsing Control1_1.R ##
cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("control file input: Control1_1\n")
cat("type(s):", paste(types, collapse = "_and_"), "\n")
cat("plddt_thresh:", plddt_thresh, "\n")

# output file names
(fileName = paste("a1_1_type_", paste(types, collapse = "_and_"), "_pthr_", plddt_thresh, "_", option, sep=""))

# libraries
library(readr)
library(dplyr)
library(broom)
library(detectseparation)

## empty data frame
results <- tibble()
fisher_results <- tibble()

## loop over types
for (type in types) {
  
  # process dfA, dfB using type, dfA, dfB
  source("Functions/dataProcessing_df.R")
  dfs <- dataProcessing_df(type, plddt_thresh) # function for data frame extraction
  final_dfA <- dfs$final_dfA # disease data
  final_dfB <- dfs$final_dfB # entanglement data
  rm(dfs)
  
  ######################################
  #### Disgenet Score Percentiles ####
  sub_finalA <- final_dfA[,c("Association_ID","score")] # subset those from disgenet
  sub_finalA <- sub_finalA[sub_finalA$Association_ID>0,]
  
  # Calculate score percentiles
  percentiles <- quantile(sub_finalA$score, c(0.95, 0.75, 0.50))
  rm(sub_finalA)
  
  #### Contingency tables ####
  # Subset necessary columns
  UniScoDis <- final_dfA[, c("uniprotids", "score","Entanglement", "Cov_Entanglement", "Knot")]
  
  # Remove duplicates
  UniScoDis <- distinct(UniScoDis)
  
  # Keep Uniprots with the highest score
  max_df <- UniScoDis %>%
    group_by(uniprotids) %>%
    filter(score == max(score)) 
  
  # Determine if the score falls in the percentile
  final_df <- max_df %>%
    mutate(
      `95th_percentile` = ifelse(is.na(score), "No", ifelse(score >= percentiles[1], "Yes", "No")),
      `75th_percentile` = ifelse(is.na(score), "No", ifelse(score >= percentiles[2], "Yes", "No")),
      `50th_percentile` = ifelse(is.na(score), "No", ifelse(score >= percentiles[3], "Yes", "No"))
    )
  
  # pull out unique lengths and essentiality per uniprot ID
  length_df <- final_dfB %>%
    select(gene, Length, Essential) %>%
    distinct()
  
  # join onto final_df
  final_df <- final_df %>%
    left_join(length_df, by = c("uniprotids" = "gene"))
  
  if (type == "af") {
    final_df_af <- final_df
    final_dfB_af <- final_dfB
  } else if (type == "crystal") {
    final_df_crystal <- final_df
    final_dfB_crystal <- final_dfB
  }
  
  # Add binary outcomes (0/1) and ensure Entanglement is a factor
  df2 <- final_df %>%
    mutate(
      Disease_95 = as.integer(`95th_percentile` == "Yes"),
      Disease_75 = as.integer(`75th_percentile` == "Yes"),
      Disease_50 = as.integer(`50th_percentile` == "Yes"),
      Essential = factor(Essential, levels = c("No", "NT", "Yes")),
      Knot = factor(Knot, levels = c("No", "Yes")),
      Entanglement = factor(Entanglement, levels = c("No","Yes")),
      Cov_Entanglement = factor(Cov_Entanglement, levels = c("No","Yes"))
    )
  df2$Dis_Category <- rowSums(df2[, c("Disease_50","Disease_75","Disease_95")], na.rm=TRUE)
  
  code_50 <- 1
  code_75 <- 2
  code_95 <- 3
  
  df2$remove_75 <- ifelse(!(df2$Dis_Category >= code_75 | df2$Dis_Category == 0), 1L, 0L)
  df2$remove_95 <- ifelse(!(df2$Dis_Category >= code_95 | df2$Dis_Category == 0), 1L, 0L)
  
  # function to capture warnings
  capture_warnings <- function(expr) {
    w <- character() # empty character vector to store warnings
    res <- withCallingHandlers( # allows to check for warnings while evaluating expr (glm)
      expr,
      warning = function(ww) {
        w <<- c(w, conditionMessage(ww)) # extract warning message and store it
        invokeRestart("muffleWarning") # suppress printing the warning
      }
    )
    list(result = res, warnings = unique(w)) # return result + warnings
  }
  
  # helpers to detect single-category (after removing NAs)
  n_levels_nonNA <- function(x) length(unique(na.omit(as.character(x)))) # counts the number of unique non-NA values in a vector
  single_level   <- function(x) n_levels_nonNA(x) < 2 # checks for less than two unique non-NA values
  
  # --- detect separation check ---
  check_sep <- function(formula, data, type, threshold) {
    glm(
      formula = formula,
      data    = data,
      family  = binomial("logit"),
      method  = detect_separation
    )
  }
  
 write_csv(df2, paste0("Data/DataRanInAnalysis/a1_1_", type, format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"))
  
  build_models <- function(df2, option = c("remove", "all")) {
    
    responses <- c(`95th` = "95th_percentile",
                   `75th` = "75th_percentile",
                   `50th` = "50th_percentile")
    
    fisher_list <- vector("list", length(responses))
    names(fisher_list) <- names(responses)
    
    models <- vector("list", length(responses)) # creates an empty list with n slots.
    names(models) <- names(responses) # name the list elements to match responses
    counts_list <- vector("list", length(responses)) # empty list of the same length
    names(counts_list) <- names(responses) # name the list elements to match responses
    warn_rows <- list() # empty list
    
    for (nm in names(responses)) {
      resp <- responses[[nm]]
      
      dat <- switch( # to choose what values to assign to dat
        nm, # percentile
        `95th` = if (option == "remove") filter(df2, remove_95 != 1L) else df2, # if remove, keep only rows where remove_95 != 1L
        `75th` = if (option == "remove") filter(df2, remove_75 != 1L) else df2, # if remove, keep only rows where remove_75 != 1L
        `50th` = df2
      )
      
      # ---- Make response, Entanglement, and Essentiality factors (No = reference) for check ----
      dat <- dat %>%
        mutate(
          !!resp := factor(.data[[resp]], levels = c("No", "Yes")),
          Entanglement = factor(Entanglement, levels = c("No", "Yes")),
          Cov_Entanglement = factor(Cov_Entanglement, levels = c("No", "Yes")),
          Knot = factor(Knot, levels = c("No", "Yes")),
          Essential    = factor(Essential,    levels = c("No","NT", "Yes"))
        )
      
      # Build the same formulas as string
      f_non  <- as.formula(sprintf("`%s` ~ Entanglement + scale(Length) + Essential", resp))
      f_cov  <- as.formula(sprintf("`%s` ~ Cov_Entanglement + scale(Length) + Essential", resp))
      f_both <- as.formula(sprintf("`%s` ~ Entanglement + Cov_Entanglement + scale(Length) + Essential", resp))
      f_knot <- as.formula(sprintf("`%s` ~ Knot + scale(Length) + Essential", resp))
      
      out_non <- capture.output(check_sep(f_non, dat, type, nm)) # capture all printed output from separation function
      sep_value_non <- sub(".*Separation:\\s*", "", grep("Separation:", out_non, value = TRUE)) # find the line containing "Separation:", extract separation status
      
      out_cov <- capture.output(check_sep(f_cov, dat, type, nm))
      sep_value_cov <- sub(".*Separation:\\s*", "", grep("Separation:", out_cov, value = TRUE))
      
      out_both <- capture.output(check_sep(f_both, dat, type, nm))
      sep_value_both <- sub(".*Separation:\\s*", "", grep("Separation:", out_both, value = TRUE))
      
      out_knot <- capture.output(check_sep(f_knot, dat, type, nm))
      sep_value_knot <- sub(".*Separation:\\s*", "", grep("Separation:", out_knot, value = TRUE))
      
      # choose disease column for counts ("95th_percentile" etc.)
      dcol <- paste0(nm, "_percentile")
      
      # ---- counts (Yes/No) ----
      is_yes <- function(x) {
        x_chr <- as.character(x)
        (x_chr == "Yes")
      }
      is_no <- function(x) {
        x_chr <- as.character(x)
        (x_chr == "No")
      }
      
      ent_yes_d_yes <- sum(is_yes(dat$Entanglement) & is_yes(dat[[dcol]]), na.rm = TRUE)
      ent_no_d_yes  <- sum(is_no(dat$Entanglement)  & is_yes(dat[[dcol]]), na.rm = TRUE)
      ent_yes_d_no  <- sum(is_yes(dat$Entanglement) & is_no(dat[[dcol]]),  na.rm = TRUE)
      ent_no_d_no   <- sum(is_no(dat$Entanglement)  & is_no(dat[[dcol]]),  na.rm = TRUE)
      
      cov_ent_yes_d_yes <- sum(is_yes(dat$Cov_Entanglement) & is_yes(dat[[dcol]]), na.rm = TRUE)
      cov_ent_yes_d_no  <- sum(is_yes(dat$Cov_Entanglement) & is_no(dat[[dcol]]),  na.rm = TRUE)
      cov_ent_no_d_yes  <- sum(is_no(dat$Cov_Entanglement)  & is_yes(dat[[dcol]]), na.rm = TRUE)
      cov_ent_no_d_no   <- sum(is_no(dat$Cov_Entanglement)  & is_no(dat[[dcol]]),  na.rm = TRUE)
      
      knot_ent_yes_d_yes <- sum(is_yes(dat$Knot) & is_yes(dat[[dcol]]), na.rm = TRUE)
      knot_ent_yes_d_no  <- sum(is_yes(dat$Knot) & is_no(dat[[dcol]]),  na.rm = TRUE)
      knot_ent_no_d_yes  <- sum(is_no(dat$Knot)  & is_yes(dat[[dcol]]), na.rm = TRUE)
      knot_ent_no_d_no   <- sum(is_no(dat$Knot)  & is_no(dat[[dcol]]),  na.rm = TRUE)
      
      counts_list[[nm]] <- tibble(
        Threshold = nm,
        ent_yes_d_yes = ent_yes_d_yes,
        ent_no_d_yes  = ent_no_d_yes,
        ent_yes_d_no  = ent_yes_d_no,
        ent_no_d_no   = ent_no_d_no,
        cov_ent_yes_d_yes = cov_ent_yes_d_yes,
        cov_ent_yes_d_no  = cov_ent_yes_d_no,
        cov_ent_no_d_yes  = cov_ent_no_d_yes,
        cov_ent_no_d_no   = cov_ent_no_d_no,
        knot_ent_yes_d_yes = knot_ent_yes_d_yes,
        knot_ent_yes_d_no  = knot_ent_yes_d_no,
        knot_ent_no_d_yes  = knot_ent_no_d_yes,
        knot_ent_no_d_no   = knot_ent_no_d_no
      )
      
      # ---- Fisher tests (Entanglement / Cov_Entanglement / Knot) ----
      # helper to run fisher safely and return a 1-row tibble
      run_fisher_2x2 <- function(yy, yn, ny, nn, predictor, threshold, type) {
        tab <- matrix(c(yy, yn, ny, nn), nrow = 2, byrow = TRUE,
                      dimnames = list(Predictor = c("Yes","No"), Disease = c("Yes","No")))
        ft <- fisher.test(tab)
        
        tibble(
          type = type,
          Threshold = threshold,
          Predictor = predictor,
          a_yes_d_yes = yy,
          b_yes_d_no  = yn,
          c_no_d_yes  = ny,
          d_no_d_no   = nn,
          OR_fisher = unname(ft$estimate),
          CI_lower_fisher = unname(ft$conf.int[1]),
          CI_upper_fisher = unname(ft$conf.int[2]),
          p_value_fisher  = unname(ft$p.value)
        )
      }
      
      # Create fisher rows (NOTE: counts are already computed above)
      fisher_ent  <- run_fisher_2x2(ent_yes_d_yes, ent_yes_d_no, ent_no_d_yes, ent_no_d_no,
                                    predictor = "Entanglement", threshold = nm, type = type)
      
      fisher_cov  <- run_fisher_2x2(cov_ent_yes_d_yes, cov_ent_yes_d_no, cov_ent_no_d_yes, cov_ent_no_d_no,
                                    predictor = "Cov_Entanglement", threshold = nm, type = type)
      
      fisher_knot <- run_fisher_2x2(knot_ent_yes_d_yes, knot_ent_yes_d_no, knot_ent_no_d_yes, knot_ent_no_d_no,
                                    predictor = "Knot", threshold = nm, type = type)
      
      # store in an attribute so we can pull it out alongside counts/warnings
      fisher_list[[nm]] <- bind_rows(fisher_ent, fisher_cov, fisher_knot)
      
      
      # ---- checking if response or used predictor has only one category ----
      resp_single <- single_level(dat[[resp]])
      ent_single  <- single_level(dat$Entanglement)
      cov_single  <- single_level(dat$Cov_Entanglement)
      knot_single <- single_level(dat$Knot)
      
      skip_non  <- resp_single || ent_single # and / or
      skip_cov  <- resp_single || cov_single
      skip_both <- resp_single || ent_single || cov_single
      skip_knot <- resp_single || knot_single
      
      # helpers to build a "skipped" placeholder
      skipped_model <- function() list(result = NULL, warnings = character(0))
      skipped_row_info <- function(model_name, reason) tibble(
        Threshold = nm, 
        Model = model_name,
        had_warning = FALSE,
        warning_text = paste0("skipped: ", reason),
        converged = as.logical(NA),
        separation_flag = as.logical(NA),
        nonconv_flag = as.logical(NA),
        sep_value = NA_character_ 
      )
      
      
      
      # fit or skip (use braces so else is unambiguous)
      fit_non <- if (!skip_non) {
        capture_warnings(
          glm(f_non,
              data = dat, family = binomial)
        )
      } else {
        skipped_model()
      }
      
      fit_cov <- if (!skip_cov) {
        capture_warnings(
          glm(f_cov,
              data = dat, family = binomial)
        )
      } else {
        skipped_model()
      }
      
      fit_both <- if (!skip_both) {
        capture_warnings(
          glm(f_both,
              data = dat, family = binomial)
        )
      } else {
        skipped_model()
      }
      
      fit_knot <- if (!skip_knot) {
        capture_warnings(
          glm(f_knot,
              data = dat, family = binomial)
        )
      } else {
        skipped_model()
      }
      
      models[[nm]] <- list(
        non_ent  = fit_non$result,
        cov_ent  = fit_cov$result,
        both_ent = fit_both$result,
        knot     = fit_knot$result
      )
      
      # warnings/flags table (mark skipped explicitly)
      wrs <- list(
        if (skip_non)  skipped_row_info("non_ent",  if (resp_single) "response has one level" else "Entanglement has one level") else
          tibble(Threshold = nm, Model = "non_ent",
                         had_warning = length(fit_non$warnings)>0L,
                         warning_text = ifelse(length(fit_non$warnings)>0L, paste(fit_non$warnings, collapse=" | "), NA_character_),
                         converged = isTRUE(fit_non$result$converged),
                         sep_value = sep_value_non),
        if (skip_cov)  skipped_row_info("cov_ent",  if (resp_single) "response has one level" else "Cov_Entanglement has one level") else
          tibble(Threshold = nm, Model = "cov_ent",
                         had_warning = length(fit_cov$warnings)>0L,
                         warning_text = ifelse(length(fit_cov$warnings)>0L, paste(fit_cov$warnings, collapse=" | "), NA_character_),
                         converged = isTRUE(fit_cov$result$converged),
                         sep_value = sep_value_cov),
        if (skip_both) skipped_row_info("both_ent", if (resp_single) "response has one level" else "Entanglement/Cov_Entanglement has one level") else
          tibble(Threshold = nm, Model = "both_ent",
                         had_warning = length(fit_both$warnings)>0L,
                         warning_text = ifelse(length(fit_both$warnings)>0L, paste(fit_both$warnings, collapse=" | "), NA_character_),
                         converged = isTRUE(fit_both$result$converged),
                         sep_value = sep_value_both),
        if (skip_knot) skipped_row_info("knot",     if (resp_single) "response has one level" else "Knot has one level") else
          tibble(Threshold = nm, Model = "knot",
                         had_warning = length(fit_knot$warnings)>0L,
                         warning_text = ifelse(length(fit_knot$warnings)>0L, paste(fit_knot$warnings, collapse=" | "), NA_character_),
                         converged = isTRUE(fit_knot$result$converged),
                         sep_value = sep_value_knot)
      )
      warn_rows[[nm]] <- bind_rows(wrs) %>%
        mutate(
          separation_flag = ifelse(!is.na(warning_text) & grepl("fitted probabilities numerically 0 or 1|separation", warning_text, TRUE), TRUE, FALSE),
          nonconv_flag    = ifelse(!is.na(warning_text) & grepl("did not converge", warning_text, TRUE), TRUE, FALSE)
        )
    }
    
    attr(models, "counts")   <- bind_rows(counts_list)
    attr(models, "warnings") <- bind_rows(warn_rows)
    attr(models, "fishers")  <- bind_rows(fisher_list)
    models
  }
  
  models <- build_models(df2, option = option)
  
  counts <- attr(models, "counts")
  
  fishers <- attr(models, "fishers")
  
  warn_fit <- attr(models, "warnings") %>%
    rename(
      had_warning_fit   = had_warning,
      warning_text_fit  = warning_text
    )
  
  # tidy helper that returns NA row for skipped models
  tidy_or_na <- function(m, mname, pct) {
    if (is.null(m)) { 
      return(tibble(
        Model = mname, term = "(skipped)", Threshold = pct,
        Estimate = NA_real_, Std_Error = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        OR = NA_real_, OR_Std_Error = NA_real_,
        CI_lower = NA_real_, CI_upper = NA_real_, p_value = NA_real_,
        had_warning_tidy = FALSE, warning_text_tidy = NA_character_
      ))
    }
    # if model was skipped above we stored list(result=NULL,...)
    if (is.null(m$coefficients)) {
      return(tibble(
        Model = mname, term = "(skipped)", Threshold = pct,
        Estimate = NA_real_, Std_Error = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        OR = NA_real_, OR_Std_Error = NA_real_,
        CI_lower = NA_real_, CI_upper = NA_real_, p_value = NA_real_,
        had_warning_tidy = FALSE, warning_text_tidy = NA_character_
      ))
    }
    
    tw <- capture_warnings(tidy(m, conf.int = TRUE))
    as_tibble(tw$result) %>%
      transmute(
        Model = mname,
        term,
        Threshold = pct,
        Estimate  = estimate,
        Std_Error = std.error,
        CI_low    = conf.low,
        CI_high   = conf.high,
        OR        = exp(estimate),
        OR_Std_Error = exp(std.error),
        CI_lower  = exp(conf.low),
        CI_upper  = exp(conf.high),
        p_value   = p.value,
        had_warning_tidy  = length(tw$warnings) > 0L,
        warning_text_tidy = ifelse(length(tw$warnings) > 0L,
                                   paste(unique(tw$warnings), collapse = " | "),
                                   NA_character_)
      )
  }
  
  # build results, preserving skipped models as NA rows
  results_list <- list()
  for (pct in names(models)) {
    mdl_set <- models[[pct]]
    for (mname in names(mdl_set)) {
      res <- tidy_or_na(mdl_set[[mname]], mname, pct)
      results_list[[length(results_list) + 1L]] <- res
    }
  }
  
  results_this <- bind_rows(results_list) %>%
    left_join(counts,  by = "Threshold") %>%
    left_join(warn_fit, by = c("Threshold","Model")) %>%
    mutate(
      warning = (had_warning_tidy | had_warning_fit) %>% replace_na(FALSE),
      separation_flag = coalesce(separation_flag, FALSE),
      nonconv_flag    = coalesce(nonconv_flag, FALSE)
    ) %>%
    select(
      -matches("warning"),  # move warning cols to end
      everything(),
      matches("warning")
    )
  
  ## append type and accumulate into `results`
  results <- bind_rows(
    results,
    mutate(results_this, type = type)
  )
  fisher_results <- bind_rows(fisher_results, fishers)
  
  # Pick data frames based on type
  final_df  <- get(paste0("final_df_", type))
  final_dfB <- get(paste0("final_dfB_", type))
  
  ## --- clean up: keep only results, types, plddt_thresh, option (and the loop var can be dropped) ---
  keep_vars <- c("results", "fisher_results", "types", "plddt_thresh", "option", "fileName", "final_df_af", "final_dfB_af", "final_df_crystal", "final_dfB_crystal")
  rm(list = setdiff(ls(), keep_vars))
  
  
}

results <- results %>%
  group_by(Model, term, Threshold) %>%
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  select(type, Threshold, everything())

fisher_results <- fisher_results %>%
  group_by(Predictor, Threshold) %>%
  mutate(p_value_fisher_BH = p.adjust(p_value_fisher, method = "BH")) %>%
  ungroup()

## Save results
write_csv(results, paste0("Results/Dataframes/", fileName, format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"),  ".csv"))
write_csv(fisher_results,paste0("Results/Dataframes/fisher_", fileName, format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"),  ".csv"))


