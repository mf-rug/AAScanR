library(shiny)
library(shinyjs)
library(DT)
library(ggplot2)
library(magrittr)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(tools)
library(reticulate)
library(zip)

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(".butt{background-color:#1f002e;} .butt{color: #ff0000;}"),
    tags$style(".butt2{background-color:#ffffff;} .butt2{color: #616161;}"),
    tags$style(HTML("
      .header {
          display: flex;
          justify-content: space-between;
          margin: 0 auto;
        }
      "))
  ),
  titlePanel('AAScanR'), 
   title = 'Generate Primers with pyAAScan',
   br(),
  div(
    style = "display: inline-block;vertical-align:top;line-height: 2.5;",
    HTML(
      '<strong>List of mutations:&nbsp&nbsp</strong>'
    )
  ),   div(style="display: inline-block;vertical-align:top", fileInput('mutfile', NULL, multiple = FALSE)), 
  div(style="display: inline-block;vertical-align:top", textAreaInput('mutinput', NULL)), 
  br(),
   div(style="display: inline-block;vertical-align:top;line-height: 2;", HTML('<Strong>Target Sequence (including flanking region)&nbsp&nbsp&nbsp</strong>')),
   div(style="display: inline-block;vertical-align:top", actionButton('example', 'load example', style = "line-height: 0.8;", class = 'butt2')),
   div(style="display: inline-block;vertical-align:top;line-height: 2;", HTML('&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp')),
   div(style="display: inline-block;vertical-align:top", disabled(actionButton('orf', 'select longest ORF', style = "line-height: 0.8;", class = 'butt2'))),
   textAreaInput('seq', NULL, width = '100%', height = '25vh', 
   ),
   fluidRow(
     column(3,
            HTML('<strong>position of first nucleotide of first codon</strong>'),
            textInput('start', NULL, width = '10ch')
     ),
     column(3,
            HTML('<strong>position of last nucleotide of last codon</strong>'),
            textInput('end', NULL, width = '10ch')
     ),
     column(3,
            HTML('<strong>number of selected codons</strong>'),
            htmlOutput('codnum')
     ),
     br()
   ),
  htmlOutput('info'),br(),
   div(style="display: inline-block;vertical-align:top", HTML('<strong>Output&nbsp</strong>')),
   div(style="display: inline-block;vertical-align:top", 
       tags$a(href = 'https://github.com/matteoferla/DirEvo_tools',
              target = "_blank",
              tags$p(
                tags$i(
                  class = "glyphicon glyphicon-info-sign",
                  style = "color:#0072B2;",
                  title = 'Primers are calculated with AAScan\'s python implementation: https://github.com/kt-korbeld/pyAAscan'
                )
              )
       )),
   HTML('<hr style="margin: 3px 0 5px" />'),
   fluidRow(
     div(
       style = 'line-height:0.3',
       br()),
       div(style = 'height:0.3vh', class='header', 
           div(style="display: inline-block;vertical-align:top", HTML('<strong>&nbsp&nbspSelected&nbspsequence:</strong>')),
           div(style="display: inline-block;vertical-align:top", plotOutput('plot', width = '70vw', height = '3em')),
           div(style="display: inline-block;vertical-align:top", disabled(downloadButton('download_prims', 'Download all primers', icon = icon('download'), class = "butt"))),
   )), br(),
   HTML('<hr style="margin: 19px 0 0px" />'),
   br(), column(12, DTOutput('prim_table')), br()
)

