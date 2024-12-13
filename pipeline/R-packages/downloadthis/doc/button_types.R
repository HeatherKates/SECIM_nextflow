## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(downloadthis)

## ----example_1----------------------------------------------------------------
mtcars %>%
  download_this(
    output_name = "mtcars dataset",
    output_extension = ".csv",
    button_label = "Download data as csv",
    button_type = "default",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## ----example_2----------------------------------------------------------------
mtcars %>%
  download_this(
    output_name = "mtcars dataset",
    output_extension = ".xlsx",
    button_label = "Download data as xlsx",
    button_type = "primary",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## ----example_3----------------------------------------------------------------
mtcars %>%
  download_this(
    output_name = "mtcars dataset",
    output_extension = ".csv",
    button_label = "Download data as csv",
    button_type = "success",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## ----example_4----------------------------------------------------------------
mtcars %>%
  download_this(
    output_name = "mtcars dataset",
    output_extension = ".xlsx",
    button_label = "Download data as xlsx",
    button_type = "info",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## ----example_5----------------------------------------------------------------
mtcars %>%
  download_this(
    output_name = "mtcars dataset",
    output_extension = ".csv",
    button_label = "Download data as csv",
    button_type = "warning",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

## ----example_6----------------------------------------------------------------
mtcars %>%
  download_this(
    output_name = "mtcars dataset",
    output_extension = ".xlsx",
    button_label = "Download data as xlsx",
    button_type = "danger",
    has_icon = TRUE,
    icon = "fa fa-save"
  )

