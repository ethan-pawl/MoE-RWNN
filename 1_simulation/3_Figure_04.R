library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggtext)

nlpl_df <- lapply(list.files(file.path("1_simulation", "NLPLs"), full.names = TRUE), function(fname) {
  readRDS(fname)
}) |> bind_rows()

load(file.path("1_simulation", "aux_simdata", "signals.Rdata"))

nlpl_df$`Signal Size Index` <- nlpl_df$`Signal Size`

nlpl_df <- nlpl_df |> 
  mutate(
    `Signal Size` = recode_values(
      Scenario,
      "Interaction in Mean" ~ signals[["inter_mu"]][`Signal Size Index`], 
      "Quadratic Mean" ~ signals[["quad_mu"]][`Signal Size Index`], 
      "Logistic Mean" ~ signals[["sig_mu"]][`Signal Size Index`], 
      default = signals[["linear"]][`Signal Size Index`]
    )
  )

scenarios <- c(
  "Linear", 
  "Interaction in Mean", "Quadratic Mean", "Logistic Mean", 
  "Interaction in Logit", "Quadratic Logit", "Logistic Logit"
)
p_list <- lapply(1:length(scenarios), function(i_scenario) {
  cur_scenario <- scenarios[i_scenario]
  p <- ggplot(filter(nlpl_df, Scenario == cur_scenario)) + 
    geom_line(
      aes(
        `Signal Size`, `NLPL (Above Oracle)`, group = Model, col = Model, 
        linetype = Model
      ), 
      linewidth = 1
    ) + 
    geom_point(aes(`Signal Size`, `NLPL (Above Oracle)`, group = Model, col = Model)) + 
    scale_linetype_manual(values = c("Linear" = "dotdash", "Nonlinear" = "solid")) + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.05))
  
  if(i_scenario == 1) {
    p <- p + 
      theme(
        legend.title = element_text(size = 18, hjust = 0.4), 
        legend.text = element_text(size = 16),
        legend.key.size = unit(2.5, "lines"), 
        legend.margin = margin(0, 0, 20, 60), 
        plot.title = element_textbox_simple(
          size = 20,
          margin = margin(t = 5, b = 10),
          face = "bold",
          fill = "grey85",
          box.color = NA,
          padding = margin(3, 6, 0.5, 6), 
          halign = 0.5
        ),
        axis.title.y.left = element_text(margin = margin(r = 10), size = 16), 
        axis.text = element_text(size = 14)
      ) + 
      labs(
        title = "Linear Means\nand Probabilities", y = "NLPL (Above Oracle)", x = ""
      ) + 
      guides(color = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2))

    leg <<- get_legend(p) |> as_ggplot()
  } else {
    p <- p + 
      theme(
        plot.title = element_text(size = 16, margin = margin(t = 5, b = 10)),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
      ) + 
      labs(title = strsplit(cur_scenario, " ")[[1]][1], y = "", x = "")
  }
  
  p <- p + theme(legend.position = "none")
  return(p)
})

p_list[[length(p_list) + 1]] <- leg

p_list[[length(p_list) + 1]] <- text_grob(expression("Signal Size (" * Delta * ")"), size = 14, x = 0.53, y = 0.9)

p_list[[length(p_list) + 1]] <- grobTree(
    rectGrob(gp = gpar(fill = "grey85", col = NA), height = unit(0.75, "cm"), 
             y = 0.275, width = unit(0.92, "npc"), x = unit(0.534, "npc")),
    text_grob("Nonlinear Means", size = 20, face = "bold", just = "top", 
              x = 0.5075)
)

p_list[[length(p_list) + 1]] <- grobTree(
    rectGrob(gp = gpar(fill = "grey85", col = NA), height = unit(0.75, "cm"), 
             y = 0.275, width = unit(0.91, "npc"), x = unit(0.5375, "npc")),
    text_grob("Nonlinear Probabilities", size = 20, face = "bold", just = "top", 
        x = 0.5325, y = 0.8)
)

p_list[[length(p_list) + 1]] <- nullGrob()
layout_mat <- matrix(c(12, 10, 10, 10,
                       12, 2, 3, 4, 
                       1, 2, 3, 4, 
                       1, 11, 11, 11,
                       1, 12, 12, 12,
                       1, 5, 6, 7,
                       8, 5, 6, 7, 
                       9, 9, 9, 9), ncol = 4, byrow = TRUE)

pdf(file.path("plots", "Figure04.pdf"), 16, 7.9)
grid.arrange(grobs = p_list, layout_matrix = layout_mat, 
             heights = c(0.25, 0.89, 1.11, 0.1, 0.05, 1.11, 0.89, 0.15), 
             widths = c(1.05, 1, 1, 1))
graphics.off()
