# function inputs to test easier

ms_stanford = read.csv("~/OSU/PhD/Mass_Shooting/Data/ms_stanford_1.csv")
ms_stanford = ms_stanford[which(ms_stanford$State != "AK"),]
ms_stanford = ms_stanford[which(ms_stanford$State != "HI"),]

k_bins = c(0, 5, 7, 9, 12, 100, 800) # standard
k_bins1 = c(0, 5, 7, 9, 12, 800) # truncated for time grouping
g_bins = c(0, 14, 93, 183, 365, 8000) # standard
h_bins = c(0, 0.5, 200)
k_bins2 = c(0,5,7,10,800)

dates = ms_stanford$Date
ref_date = "01/01/1999"
marks = ms_stanford$Total.Number.of.Victims
g_bins = g_bins
k_bins = k_bins1
k_quantile = TRUE
k_length = 5
lat = rep(0, length(dates))
lon = rep(0, length(dates))
h_bins = c(0,1)
g_quantile = FALSE
h_quantile = FALSE
k_length = 6
h_length = 6
g_length = 6
stopwhen = 1e-3
time_unit = "day"
dist_unit = "mile"
