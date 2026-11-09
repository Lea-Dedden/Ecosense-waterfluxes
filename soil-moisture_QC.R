
# QUALITY CRITERIA #
# =====================================================================
# SOIL MOISTURE QC PIPELINE  (following ISMN; Dorigo et al. 2013)
# =====================================================================

library(dplyr)
library(zoo)
library(signal)


# flag criteria - flag added to qflag column if crtieria are met
# names = ['C01', 'C02', 'C03', 'D01', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09', 'D10', 'G']

# C01 soil moisture < threshold
# C02 soil moisture > threshold
# C03 soil moisture > saturation point
# D01 negative soil temperature( in situ)
# D04 rise in soil moisture without precipitation( in situ)
# D06 spikes
# D07 negative breaks(drops)
# D08 positive breaks(jumps)
# D09.1 constant values following negative break
# D09.2 saturated plateaus
# G good

# ---------------------------------------------------------------------
# 1) CONFIGURATION: site-specific
# ---------------------------------------------------------------------
cfg <- list(
  tz            = "CET",  # timezone for date conversion (e.g., for daily precip)
  deriv_p       = 2,     # SG polynomial degree
  deriv_n       = 3,     # SG window (n = p+1 -> pure derivative, no smoothing)
  c01_max       = 0.05,  # SM below this (m3/m3)
  c03_sat       = 0.6,   # SM above saturation point
  d01_temp_max  = 1,     # negative soil temperature 
  d04_rise15    = 0.1,   # 15-min rise threshold
  d04_precip_win= 96,    # 24h precip window (15-min steps)
  d04_sm_win    = 100,   # SM sd window
  d04_precip_max= 0.1,   # max precip in window for "no rain"
  d06_var_win   = 25,    # 25h variance window (center)
  d06_peak_win  = 4,     # peak detection window
  d07_mean_win  = 100,   # eq8 rolling mean window
  d09_12h       = 48,    # 12h in 15-min steps
  d09_24h       = 96,
  d09_48h       = 192
)



# ---------------------------------------------------------------------
# 2) HELPERS
# ---------------------------------------------------------------------

# flag appender:  add_flag(df, <logical vector>, "TAG")
add_flag <- function(df, condition, tag) {
  if (is.null(df)) return(df)
  if (!"qflag" %in% names(df)) df$qflag <- NA_character_

  idx <- which(condition)                      # treats NA as FALSE
  if (length(idx) == 0) return(df)

  existing <- df$qflag[idx]
  has_tag  <- !is.na(existing) &
              grepl(paste0("(^|,)", tag, "(,|$)"), existing)
  new_idx  <- idx[!has_tag]

  df$qflag[new_idx] <- ifelse(
    is.na(df$qflag[new_idx]),
    tag,
    paste(df$qflag[new_idx], tag, sep = ",")
  )
  df
}

# 1 & 2 derivatives of soil moisture (Savitzky-Golay filter following Dorigo et al. 2013)
add_savgol_deriv <- function(df, col = "SoilMoisture", m = 1,
                             p = cfg$deriv_p, n = cfg$deriv_n) {
  if (is.null(df) || !col %in% names(df)) return(df)

  out_col <- paste0("deriv", m)
  df[[out_col]] <- NA_real_

  valid_idx <- which(!is.na(df[[col]]))
  if (length(valid_idx) >= n) {
    df[[out_col]][valid_idx] <- signal::sgolayfilt(
      df[[col]][valid_idx], p = p, n = n, m = m
    )
  }
  df
}

# variance / mean over 25h window excluding current value
rolling_var_excl <- function(x) {
  if (length(x) <= 1) return(NA_real_)
  var(x[-13])
}
rolling_mean_excl <- function(x) {
  if (length(x) <= 1) return(NA_real_)
  mean(x[-13])
}

# peak detector: 0 = none, 1 = 1-step peak, 2 = 2-step (flat-top) peak
peak_D06 <- function(x) {
  if (anyNA(x)) return(0L)
  if (length(x) == 3) {
    if ((x[1] > x[2] & x[2] < x[3]) |
        (x[1] < x[2] & x[2] > x[3])) return(1L)
  } else if (length(x) > 3) {
    if ((x[1] > x[2] & x[2] == x[3] & x[3] < x[4]) |
        (x[1] < x[2] & x[2] == x[3] & x[3] > x[4])) return(2L)
  }
  0L
}

# anchored pattern for a tag (avoids partial matches like "D1" in "D10")
tag_pattern <- function(tag) paste0("(^|,)", tag, "(,|$)")

count_flag <- function(list_of_df, tag) {
  sum(vapply(list_of_df, function(df) {
    if (is.null(df) || !"qflag" %in% names(df)) return(0L)
    sum(grepl(tag_pattern(tag), df$qflag), na.rm = TRUE)
  }, integer(1)))
}


# ---------------------------------------------------------------------
# 3) QC FUNCTIONS
# ---------------------------------------------------------------------
qc_plot <- function(df, rain = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)

  # prep
  if ("TIMESTAMP" %in% names(df)) df$day <- as.Date(df$TIMESTAMP, tz = cfg$tz)
  if (!"qflag" %in% names(df))    df$qflag <- NA_character_

  # derivatives
  df <- add_savgol_deriv(df, m = 1)
  df <- add_savgol_deriv(df, m = 2)

  # C01: SM below threshold
  df <- add_flag(df, df$SoilMoisture < cfg$c01_max, "C01")

  # C02: SM above threshold
  # df <- add_flag(df, df$SoilMoisture > 60, "C02")

  # C03: SM above saturation point
  df <- add_flag(df, df$SoilMoisture > cfg$c03_sat, "C03")

  # D01: negative (or zero) soil temperature
  if ("SoilTemp" %in% names(df)) {
    df <- add_flag(df, df$SoilTemp < cfg$d01_temp_max, "D01")
  }

  # D04: SM rise without precipitation
  if (!is.null(rain) && "TIMESTAMP" %in% names(df)) {
    # NOTE: verify the precip column name in df_rain_47m_15min!
    #       and that TIMESTAMP values are unique in the rain table
    df <- dplyr::left_join(df, rain, by = "TIMESTAMP")
    if (!"precip" %in% names(df)) {
      stop("D04: rain table has no 'precip' column -- check the join.")
    }
    df$precip[is.na(df$precip)] <- 0

    df <- df %>% mutate(
      precip2    = rollapplyr(precip, width = cfg$d04_precip_win,
                              FUN = sum, fill = NA),
      std_x2     = rollapplyr(SoilMoisture, width = cfg$d04_sm_win,
                              FUN = function(x) sd(x, na.rm = FALSE) * 2,
                              fill = NA),
      rise24h    = SoilMoisture - dplyr::lag(SoilMoisture, 96),
      rise15min  = SoilMoisture - dplyr::lag(SoilMoisture, 1)
    )
    df <- add_flag(df,
      df$rise15min > cfg$d04_rise15 &
      df$rise24h   > df$std_x2 &
      !is.na(df$rise24h) & !is.na(df$std_x2) &
      !is.na(df$precip2) & df$precip2 < cfg$d04_precip_max,
      "D04")
  }

  # D06: spikes
  df <- df %>% mutate(
    eq4 = round(SoilMoisture / dplyr::lag(SoilMoisture), 3),
    eq5 = dplyr::lag(round(abs(savgol2 / dplyr::lead(savgol2, n = 2)), 3)),
    eq6_I   = abs(rollapply(SoilMoisture, width = cfg$d06_var_win,
                            FUN = rolling_var_excl,  align = "center", fill = NA) /
                  rollapply(SoilMoisture, width = cfg$d06_var_win,
                            FUN = rolling_mean_excl, align = "center", fill = NA)),
    eq6_II  = rollapply(SoilMoisture, width = cfg$d06_peak_win,
                        FUN = peak_D06, align = "center", fill = 0L, partial = TRUE),
    eq6_spike2h = as.numeric(dplyr::lag(eq6_II) > 1),
    eq6_spike = as.numeric(
      !is.na(eq4) & !is.na(eq5) & !is.na(eq6_I) & !is.na(eq6_spike2h) &
      ((eq4 > 1.04 | eq4 < 0.96) | eq6_spike2h > 0) &
      eq5 > 0.5 & eq5 < 1.5 &
      eq6_I < 1 & eq6_II > 0
    )
  )
  df <- add_flag(df,
    df$eq6_spike > 0 |
      (dplyr::lag(df$eq6_spike) > 0 & df$eq6_spike2h > 0),
    "D06")

  # D07 (drops) / D08 (jumps)
  df <- df %>% mutate(
    abs_change = c(NA, diff(SoilMoisture)),
    rel_change = abs(abs_change / SoilMoisture),
    eq8  = abs(rollapply(savgol1, width = cfg$d07_mean_win,
                         FUN = mean, na.rm = TRUE,          # FIXED: rollmean has no na.rm
                         align = "center", fill = NA) * 10),
    eq9a = round(abs(savgol2 / dplyr::lag(savgol2)), 1),
    eq9b = abs(savgol2 / dplyr::lead(savgol2))
  )
  df <- add_flag(df,
    df$rel_change > 0.1 &
    abs(df$abs_change) > 1 &
    abs(df$savgol1) > abs(df$eq8) &
    abs(df$savgol1 - df$eq8) > 1e-2 &
    abs(df$eq9a - 1) <= 1e-2 &
    df$eq9b > 10,
    "D07")

  # D09.1: invariant values for >= 12h
  df$abs_var_12h <- rollapplyr(df$SoilMoisture, width = cfg$d09_12h,
                               FUN = var, fill = NA)
  df$abs_var_24h <- rollapplyr(df$SoilMoisture, width = cfg$d09_24h,
                               FUN = var, fill = NA)
  df <- add_flag(df,
    df$abs_var_12h == 0 & df$rel_change == 0,
    "D09.1")

  # D09.2: exact same value for 48h
  if (nrow(df) >= cfg$d09_48h) {
    sm_chr <- ifelse(is.na(df$SoilMoisture), "__NA__", as.character(df$SoilMoisture))
    r <- rle(sm_chr)
    df$const48_run <- rep(r$lengths, r$lengths)   # run length back-expanded
    df <- add_flag(df,
      df$const48_run >= cfg$d09_48h & df$rel_change == 0,
      "D09.2")
  } else {
    df$const48_run <- NA_integer_
  }

  # G: good quality for unflagged values
  df$qflag[is.na(df$qflag)] <- "G"

  df
}

# ---------------------------------------------------------------------
# 4) RUN
# ---------------------------------------------------------------------
list_flagged <- lapply(list_SM_5, qc_plot, rain = df_rain_47m_15min)

# ---------------------------------------------------------------------
# 5) SANITY CHECK
# ---------------------------------------------------------------------

stopifnot(all(vapply(list_flagged, function(df)
  !any(grepl("(^|,)NA(,|$)", df$qflag), na.rm = TRUE), logical(1))))

# every frame has the essential columns
stopifnot(all(vapply(list_flagged, function(df)
  all(c("TIMESTAMP", "SoilMoisture", "qflag") %in% names(df)), logical(1))))

# ---------------------------------------------------------------------
# 6) SUMMARY
# ---------------------------------------------------------------------
tags <- c("C01", "C03", "D01", "D04", "D06", "D07", "D09.1", "D09.2")
flag_summary <- sapply(tags, count_flag, list_of_df = list_flagged)
print(flag_summary)

# flag-combination distribution across all plot
print(sort(table(unlist(lapply(list_flagged, function(df) df$qflag))),
           decreasing = TRUE))

# ---------------------------------------------------------------------
# 7) EXPORT
# ---------------------------------------------------------------------
list_flagged <- lapply(list_flagged, function(df) {
  df <- df %>%
    mutate(
      qflag_aqua = case_when(
        grepl("(^|,)(C01|C03|D01|D09\\.2)(,|$)", qflag) ~ "11",
        grepl("(^|,)(D04|D06|D07|D09\\.1)(,|$)", qflag) ~ "21",
        TRUE ~ "31"
      )
    ) %>%
    dplyr::select(dplyr::any_of(c("TIMESTAMP", "SoilMoisture", "ID",
                                  "qflag", "precip", "qflag_aqua")))
  df
})
