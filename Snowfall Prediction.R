
library(readxl)


xlsx_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Kaggle/Datathon/Copy of 2025 Allianz Datathon Dataset.xlsx"
sheet     <- "Climate Data"


wk_mmdd <- c("06-09","06-16","06-23","06-30",
             "07-07","07-14","07-21","07-28",
             "08-04","08-11","08-18","08-25",
             "09-01","09-08","09-15")

clim <- read_excel(xlsx_path, sheet = sheet)

# Normalize column names
names(clim) <- sub("\\s+", " ", names(clim))
station_col <- "Bureau of Meteorology station number"
rain_col    <- "Rainfall amount (millimetres)"


clim$Date <- as.Date(sprintf("%04d-%02d-%02d", clim$Year, clim$Month, clim$Day))
clim <- clim[!is.na(clim$Date) & clim$Year >= 2014, ]


clim[[station_col]] <- as.character(clim[[station_col]])

# Day-of-year & seasonal harmonics
clim$doy  <- as.integer(strftime(clim$Date, "%j"))
twopi     <- 2*pi
clim$sin1 <- sin(twopi * clim$doy / 365.25)
clim$cos1 <- cos(twopi * clim$doy / 365.25)
clim$sin2 <- sin(2*twopi * clim$doy / 365.25)
clim$cos2 <- cos(2*twopi * clim$doy / 365.25)

# Wet-day flag (>= 1 mm)
clim$wet <- as.integer(clim[[rain_col]] >= 1)


make_ski_weeks <- function(yr) {
  starts <- as.Date(paste0(yr, "-", wk_mmdd))
  data.frame(Year = yr, Week = seq_along(starts), start = starts, end = starts + 6)
}

years_for_weeks <- sort(unique(c(clim$Year, 2026)))
wk_windows <- do.call(rbind, lapply(years_for_weeks, make_ski_weeks))

assign_week <- function(d) {
  yr <- as.integer(format(d, "%Y"))
  ww <- wk_windows[wk_windows$Year == yr, ]
  i  <- which(d >= ww$start & d <= ww$end)
  if (length(i) == 1) return(ww$Week[i])
  return(NA_integer_)
}

# Keep only dates inside Ski Weeks (for model fitting)
clim$SkiWeek <- vapply(clim$Date, assign_week, FUN.VALUE = integer(1))
clim_sw <- clim[!is.na(clim$SkiWeek), ]


stations <- sort(unique(clim_sw[[station_col]]))

wet_models <- vector("list", length(stations)); names(wet_models) <- stations
amt_models <- vector("list", length(stations)); names(amt_models) <- stations
use_lognorm_fallback <- setNames(rep(FALSE, length(stations)), stations)

for (s in stations) {
  ds <- clim_sw[clim_sw[[station_col]] == s, ]
  
  # Wet-day probability (binomial logit)
  wet_models[[s]] <- glm(
    wet ~ sin1 + cos1 + sin2 + cos2 + Year,
    family = binomial(link = "logit"),
    data = ds
  )
  
  # Amount on wet days (Gamma log) – fallback to Gaussian on log(rain)
  ds_wet <- ds[ds$wet == 1 & !is.na(ds[[rain_col]]), ]
  if (nrow(ds_wet) >= 30) {
    amt_models[[s]] <- glm(
      ds_wet[[rain_col]] ~ sin1 + cos1 + sin2 + cos2 + Year,
      family = Gamma(link = "log"),
      data = ds_wet
    )
  } else {
    ds_wet$log_r <- log(pmax(ds_wet[[rain_col]], 1e-6))
    amt_models[[s]] <- glm(
      log_r ~ sin1 + cos1 + sin2 + cos2 + Year,
      family = gaussian(),
      data = ds_wet
    )
    use_lognorm_fallback[s] <- TRUE
  }
}

# BUILD 2026 DAILY NEWDATA (ski weeks only) 
dates_2026 <- as.Date("2026-01-01") + 0:366
dates_2026 <- dates_2026[format(dates_2026, "%Y") == "2026"]

doy_2026 <- as.integer(strftime(dates_2026, "%j"))
df_2026 <- data.frame(
  Date = dates_2026, Year = 2026, doy = doy_2026,
  sin1 = sin(twopi * doy_2026 / 365.25),
  cos1 = cos(twopi * doy_2026 / 365.25),
  sin2 = sin(2*twopi * doy_2026 / 365.25),
  cos2 = cos(2*twopi * doy_2026 / 365.25)
)
df_2026$SkiWeek <- vapply(df_2026$Date, assign_week, FUN.VALUE = integer(1))
df_2026 <- df_2026[!is.na(df_2026$SkiWeek), ]
stopifnot(nrow(df_2026) == 15*7)  # 105 days

# PREDICT DAILY 2026 BY STATION=
pred_rows <- list()
for (s in stations) {
  wm <- wet_models[[s]]
  am <- amt_models[[s]]
  
  nd <- df_2026
  p_wet <- as.numeric(predict(wm, newdata = nd, type = "response"))
  
  if (!use_lognorm_fallback[s]) {
    # Gamma(log): response is mean on original scale
    mu_wet <- as.numeric(predict(am, newdata = nd, type = "response"))
  } else {
    # Gaussian on log(rain): back-transform with lognormal correction
    mu_log <- as.numeric(predict(am, newdata = nd, type = "response"))
    sigma2 <- summary(am)$sigma^2
    mu_wet <- exp(mu_log + 0.5 * sigma2)
  }
  
  exp_rain <- p_wet * mu_wet
  
  pred_rows[[s]] <- data.frame(
    Station = s, Date = nd$Date, SkiWeek = nd$SkiWeek,
    p_wet = p_wet, E_rain_wet = mu_wet, E_rain = exp_rain
  )
}
pred_daily_2026 <- do.call(rbind, pred_rows)
row.names(pred_daily_2026) <- NULL

# Weekly Aggregate by station
weekly_station <- aggregate(E_rain ~ Station + SkiWeek, data = pred_daily_2026, sum)
names(weekly_station)[names(weekly_station) == "E_rain"] <- "weekly_total_mm"
weekly_station$weekly_mean_mmday <- weekly_station$weekly_total_mm / 7

# map stations to 9 resorts

# 71032 Thredbo, 71075 Perisher, 83024 Mt Buller, 83084 Falls Creek,
# 83085 Mt Hotham, 85291 Mt Baw Baw
map <- rbind(
  data.frame(Resort="Thredbo",        Station="71032", weight=1),
  data.frame(Resort="Perisher",       Station="71075", weight=1),
  data.frame(Resort="Mt Buller",      Station="83024", weight=1),
  data.frame(Resort="Falls Creek",    Station="83084", weight=1),
  data.frame(Resort="Mt Hotham",      Station="83085", weight=1),
  data.frame(Resort="Mt Baw Baw",     Station="85291", weight=1),
  data.frame(Resort="Selwyn",         Station="71075", weight=1),
  data.frame(Resort="Mt Stirling",    Station="83024", weight=1),
  data.frame(Resort="Charlotte Pass", Station="71032", weight=0.5),
  data.frame(Resort="Charlotte Pass", Station="71075", weight=0.5)
)

tmp <- merge(weekly_station, map, by = "Station", all.x = TRUE, all.y = FALSE)
tmp$weighted_mean <- tmp$weekly_mean_mmday * tmp$weight

weekly_resort <- aggregate(weighted_mean ~ Resort + SkiWeek, data = tmp, sum)
names(weekly_resort)[names(weekly_resort) == "weighted_mean"] <- "avg_mm_per_day"

# Ensure all Resort×Week present
all_pairs <- expand.grid(Resort = unique(map$Resort), SkiWeek = 1:15)
weekly_resort <- merge(all_pairs, weekly_resort, by = c("Resort","SkiWeek"), all.x = TRUE)
weekly_resort$avg_mm_per_day[is.na(weekly_resort$avg_mm_per_day)] <- 0

# final table (2025 week starts + 2026 values
wk2025 <- data.frame(SkiWeek = 1:15, WeekStart = as.Date(paste0("2025-", wk_mmdd)))

tmp2 <- merge(weekly_resort[, c("Resort","SkiWeek","avg_mm_per_day")],
              wk2025, by = "SkiWeek", all.x = TRUE)

tmp2$start_date_YY_mm_dd <- format(tmp2$WeekStart, "%y-%m-%d")

final_long <- tmp2[order(tmp2$WeekStart, tmp2$Resort),
                   c("start_date_YY_mm_dd", "Resort", "avg_mm_per_day")]
names(final_long) <- c("start_date_YY_mm_dd", "mountain", "pred_avg_mm_per_day_2026")

# ---- This is the dataframe you asked for ----
print(head(final_long, 20))
# Optional: save it
# dir.create("outputs", showWarnings = FALSE)
write.csv(final_long, "weekly_avg_mm_per_day_2026_long_with_2025_weekstarts.csv", row.names = FALSE)


