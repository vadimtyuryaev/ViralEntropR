source("setup.R")
source("plot_results.R")
config <- build_config()
config$STUDY_DIR <- "."
generate_all_plots(config)

# without running bootstrapped CIs (if they are already cached)
# source("plot_results.R")
# config <- build_config()
# config$STUDY_DIR <- "."
# df <- load_all_summaries(config)
# 
# plots_dir <- file.path(output_dir(config), "plots")
# 
# plot_histograms(df,            out_path = file.path(plots_dir, "histograms.png"))
# plot_boxplots(df,              out_path = file.path(plots_dir, "boxplots.png"))
# plot_mean_ci(df,               out_path = file.path(plots_dir, "mean_ci.png"))
# plot_km_curves(df,             out_path = file.path(plots_dir, "km_curves.png"))
# plot_detection_rate_heatmaps(df, out_dir = plots_dir)
# plot_detection_curves(df,      out_path = file.path(plots_dir, "detection_curves.png"))