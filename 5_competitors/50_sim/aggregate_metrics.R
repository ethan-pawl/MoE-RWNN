metrics_dir <- file.path("5_competitors", "50_sim", "metrics")

metrics <- rbind(
  read.csv(file.path(metrics_dir, "gam.csv")), 
  read.csv(file.path(metrics_dir, "linear_flowmix.csv")), 
  read.csv(file.path(metrics_dir, "nonlinear_flowmix.csv")), 
  read.csv(file.path(metrics_dir, "mixdistreg.csv")), 
  read.csv(file.path(metrics_dir, "rf_wts.csv"))
)

library(ggplot2)

metrics$Model <- factor(metrics$Model, 
  levels = c("Nonlinear Model", "Linear Model", "GAM", "mixdistreg", "Random Forest"), 
  labels = c("Nonlinear Model", "Linear Model", "kMeans, GAM", "MoE-DR", "sidCluster, RF")
)

metrics$Metric <- factor(metrics$Metric,
  levels = c("RMSE, Mean", "RMSE, Probability", "Frobenius Error, Covariance"), 
  labels = c("RMSE(Mean)", "RMSE(Probability)", "FE(Covariance)")
)

# Mark RF covariance error as NA since using the residual covariance is not a reasonable comparison
metrics$Value[metrics$Model == "sidCluster, RF" & metrics$Metric == "FE(Covariance)"]

ggplot(metrics) + 
  geom_col(aes(Cluster, -log(Value), fill = Model), position = "dodge") + 
  facet_wrap(~ Metric)


pd <- position_dodge(width = 0.9)

set1 <- RColorBrewer::brewer.pal(9, "Set1")

# TODO: grid.arrange this so there's no weird hole in Probability
# TODO: add gray theme

base_size <- 18
p <- ggplot(metrics, aes(Cluster, -log(Value), color = Model)) +
  geom_segment(
    aes(x = Cluster, y = 0, yend = -log(Value)),
    position = pd,
    linewidth = 0.7
  ) +
  geom_point(
    position = pd,
    size = 2.5
  ) + 
  facet_wrap(~ Metric) + 
  scale_color_manual(
    values = c(
      "Nonlinear Model" = "black", 
      "Linear Model" = set1[5], 
      "GAM" = set1[4], 
      "MoE-DR" = set1[1], 
      "Random Forest" = set1[2]
    )
  ) + 
  labs(y = "-log(Error)") + 
  theme_gray() + 
  theme(
    plot.title    = ggplot2::element_text(size = base_size * 1.35),
    plot.subtitle = ggplot2::element_text(size = base_size * 1.10),
    plot.caption  = ggplot2::element_text(size = base_size * 0.85),
    axis.title    = ggplot2::element_text(size = base_size * 1.00),
    axis.text     = ggplot2::element_text(size = base_size * 0.80),
    legend.title  = ggplot2::element_text(size = base_size * 0.95),
    legend.text   = ggplot2::element_text(size = base_size * 0.90),
    strip.text    = ggplot2::element_text(size = base_size * 0.95)
  )


plist <- lapply(levels(metrics$Metric), function(cur_metric) {
  ggplot(subset(metrics, Metric == cur_metric), aes(Cluster, -log(Value), color = Model)) +
    geom_segment(
      aes(x = Cluster, y = 0, yend = -log(Value)),
      position = pd,
      linewidth = 0.7
    ) +
    geom_point(
      position = pd,
      size = 2.5
    ) + 
    facet_wrap(~ Metric) + 
    scale_color_manual(
      values = c(
        "Nonlinear Model" = "black", 
        "Linear Model" = set1[5], 
        "GAM" = set1[4], 
        "MoE-DR" = set1[1], 
        "Random Forest" = set1[2]
      )
    ) + 
    labs(x = "", y = "") + 
    theme_gray() + 
    theme(
      plot.title    = element_text(size = base_size * 1.35),
      plot.subtitle = element_text(size = base_size * 1.10),
      plot.caption  = element_text(size = base_size * 0.85),
      axis.title    = element_text(size = base_size * 1.00),
      axis.text     = element_text(size = base_size * 0.80),
      legend.title  = element_text(size = base_size * 0.95),
      legend.text   = element_text(size = base_size * 0.90),
      strip.text    = element_text(size = base_size * 0.95), 
      plot.margin = margin(0, 0, 0, 0)
    )
})

library(ggpubr)

leg <- get_legend(plist[[1]]) |> as_ggplot()
plist <- lapply(plist, function(p) p + theme(legend.position = "none"))
plist[[1]] <- plist[[1]] + labs(y = "-log(Error)")
plist[[2]] <- plist[[2]] + labs(x = "Cluster")
plist[[4]] <- leg

gridded_plot <- gridExtra::grid.arrange(grobs = plist, ncol = 4, widths = c(1, 0.4, 1, 0.35))
ggsave(file.path("5_competitors", "50_sim", "metrics_plot.pdf"), gridded_plot, width = 20.7, height = 7, units = "in")
