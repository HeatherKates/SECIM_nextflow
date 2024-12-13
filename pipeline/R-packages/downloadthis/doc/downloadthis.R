## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(downloadthis)

mtcars %>%
  download_this(
    output_name = "mtcars dataset",
    output_extension = ".csv",
    button_label = "Download data as csv",
    button_type = "default",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## -----------------------------------------------------------------------------
list(mtcars, iris) %>%
  download_this(
    output_name = "mtcars and iris datasets",
    output_extension = ".xlsx",
    button_label = "Download datasets as xlsx",
    button_type = "warning",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## -----------------------------------------------------------------------------
list("mtcars dataset" = mtcars, "iris dataset" = iris) %>%
  download_this(
    output_name = "mtcars and iris datasets",
    output_extension = ".xlsx",
    button_label = "Download datasets as xlsx",
    button_type = "warning",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## -----------------------------------------------------------------------------
vector_example <- 1:10
linear_model <- lm(mpg ~ gear, data = mtcars)

list(mtcars, iris, vector_example, linear_model) %>%
  download_this(
    output_name = "datasets, vector, and linear model",
    output_extension = ".rds",
    button_label = "Download as rds",
    button_type = "warning",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## -----------------------------------------------------------------------------
## Link in Github repo
download_link(
  link = "https://github.com/fmmattioni/downloadthis/raw/master/inst/example/file_1.pdf",
  button_label = "Download pdf file",
  button_type = "danger",
  has_icon = TRUE,
  icon = "fa fa-save",
  self_contained = FALSE
)

## -----------------------------------------------------------------------------
download_file(
  path = system.file("assets/css/all.min.css", package = "downloadthis"),
  output_name = "CSS file from downloadthis",
  button_label = "Download css file",
  button_type = "danger",
  has_icon = TRUE,
  icon = "fa fa-save",
  self_contained = FALSE
)

## -----------------------------------------------------------------------------
path_files <- list.files(path = system.file("assets/css", package = "downloadthis"), full.names = TRUE)

download_file(
  path = path_files,
  output_name = "Files from downloadthis",
  button_label = "Download files",
  button_type = "danger",
  has_icon = TRUE,
  icon = "fa fa-save",
  self_contained = FALSE
)

## -----------------------------------------------------------------------------
download_dir(
  path = system.file("assets", package = "downloadthis"),
  output_name = "example dir",
  button_label = "Download directory",
  button_type = "success",
  has_icon = TRUE,
  icon = "fa fa-save",
  self_contained = FALSE
)

## -----------------------------------------------------------------------------
library(ggplot2)

plot <- ggplot(cars, aes(speed, dist)) +
  geom_point()

plot

## -----------------------------------------------------------------------------
download_this(plot)

