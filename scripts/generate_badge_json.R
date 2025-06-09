# generate_badge_json.R
library(cranlogs)
library(jsonlite)

pkg <- "PartialNetwork"
downloads <- cran_downloads(pkg, from = "2023-01-01")
total <- sum(downloads$count)

# Format number
label <- if (total >= 1e6) {
  paste0(round(total / 1e6, 1), "M")
} else if (total >= 1e3) {
  paste0(round(total / 1e3), "K")
} else {
  as.character(total)
}

json <- list(
  schemaVersion = 1,
  label = "CRAN",
  message = paste0(label, " downloads"),
  color = "blue"
)

# Save JSON
# Remove everything in badges/
unlink("badges", recursive = TRUE)
dir.create("badges", showWarnings = FALSE)
write_json(json, "badges/cran_downloads.json", auto_unbox = TRUE)
