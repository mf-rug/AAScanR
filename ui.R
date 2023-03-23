library(tidyverse)
library(shiny)
library(shinyFiles)
library(DT)
library(shinyjs)
library(tools)
library(Biostrings)

ui <- fluidPage(
  useShinyjs(),
  titlePanel('AAScanR'), 
  HTML('<strong>Input</strong>'),hr(),
  fluidRow(
    column(2, fileInput('aafiles', '20 Files from AAscan', multiple = TRUE)),
    column(2, fileInput('mutfile', 'File with a list of mutations', multiple = FALSE))
  ),
  HTML('<strong>Output</strong>'),hr(),
  fluidRow(
    column(12,
           disabled(downloadButton('download', 'Download tsv file', icon = icon('download'))),
           align = 'right')
  ), br(),
  DTOutput('table'),
  br(),hr(),br(),
  htmlOutput('info')
)

