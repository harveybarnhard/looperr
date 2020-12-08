library(data.table)
library(ggplot2)
library(scales)

# Load results =================================================================
results_R = fread("data-raw/results_R.csv")
results_S = fread("data-raw/results_stata.csv")

# Merge results and calculate relative speed ===================================
results = rbind(results_R, results_S)
results = merge(results, results_S[, c("seconds", "nObs", "nGroups")],
                by=c("nObs", "nGroups"))
results[, relative_speed := seconds.y /seconds.x]


# Rename results ===============================================================
results[expr=="baseR", expr:="Base R loop with lm()"]
results[expr=="cpp_onecore", expr:="C++ One Thread"]
results[expr=="cpp_twocore", expr:="C++ Two Threads"]
results[expr=="cpp_threecore", expr:="C++ Three Threads"]
results[expr=="cpp_fourcore", expr:="C++ Four Threads"]
results[expr=="regressby", expr:="Stata regressby Four Threads"]

# Plot results =================================================================
g = ggplot(results, aes(x=as.character(nObs), y=as.character(nGroups),
                        fill=relative_speed)) +
  geom_tile(color="black") +
  geom_text(aes(label=round(relative_speed,2)),fontface = "bold") +
  scale_fill_gradientn(colours=c("white", "yellow", "orange","red"),
                       values=rescale(c(1,5,18,200), from=c(min(results$relative_speed), max(results$relative_speed)))) +
  facet_wrap(~expr) +
  theme_minimal() +
  theme(legend.position="none",
        strip.text = element_text(size = 20),
        axis.title = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x="Number of Observations Per Group",
       y="Number of Groups")
ggsave("examples/benchmark_regressby.svg", plot=g)
