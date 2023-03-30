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
  tabsetPanel(
   tabPanel(
      br(),
      title = 'Select Primers',
      div(style="display: inline-block;vertical-align:top", HTML('<strong>Input&nbsp</strong>')),
      div(style="display: inline-block;vertical-align:top", 
          tags$a(href = 'https://github.com/dmitryveprintsev/AAScan',
                 target = "_blank",
                 tags$p(
                   tags$i(
                     class = "glyphicon glyphicon-info-sign",
                     style = "color:#0072B2;",
                     title = 'AAScan is available at https://github.com/dmitryveprintsev/AAScan (Windows) and https://github.com/dbv123w/AAScan (MacOS)'
                   )
                 )
          )),
      HTML('<hr style="margin: 3px 0 25px" />'),
      div(style="display: inline-block;vertical-align:top;line-height: 2.5;", HTML('<strong>20 Files from AAscan / deep_mut_scan:&nbsp&nbsp</strong>')),
      div(style="display: inline-block;vertical-align:top", fileInput('aafiles', NULL, multiple = TRUE)),
      div(style="display: inline-block;vertical-align:top", 
          tags$p(
            tags$i(
              style = "line-height: 2.5;",
              class = "glyphicon glyphicon-info-sign",
              style = "color:#0072B2;",
              title = 'Required input: 20 plain text files named A.txt, C.txt, D.txt, ... for each amino acid. Each file must contain two columns, in the format <residueNum>_<F/R> <primerSequence>, e.g.:\n\n1_F aaaagatgcatgaagatcggaca\n1_R agcctgtaatagtagcatcagtatata\n2_F agtctgatctagtttattattaaa\n...'
            )
          )
      ),      
      div(style="display: inline-block;vertical-align:top;line-height: 2.5;", HTML('<strong>&nbsp&nbsp&nbsp&nbsp&nbspFile with a list of mutations:&nbsp&nbsp</strong>')),
      div(style="display: inline-block;vertical-align:top", fileInput('mutfile', NULL, multiple = FALSE)), 
      div(style="display: inline-block;vertical-align:top", 
          tags$p(
            tags$i(
              style = "line-height: 2.5;",
              class = "glyphicon glyphicon-info-sign",
              style = "color:#0072B2;",
              title = 'Required input: 1 plain text file with the desired mutations. Must contain one column, in the format <originalAA><residueNum><mutatedAA>, e.g.:\n\nD3A\nG6F\nH20K\n...'
            )
          )
      ), br(),  
      HTML('<strong>Output</strong>'),
      HTML('<hr style="margin: 8px 0 25px" />'),
      fluidRow(
        column(12,
               disabled(downloadButton('download', 'Download tsv file', icon = icon('download'), class = "butt")),
               align = 'right')
      ), br(),
      DTOutput('table'),
      br(),hr(),br(),
      htmlOutput('info')
    ),
   tabPanel(
     title = 'Generate Primers with deep_mut_scan',
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
     div(style="display: inline-block;vertical-align:top", HTML('<strong>Output&nbsp</strong>')),
     div(style="display: inline-block;vertical-align:top", 
         tags$a(href = 'https://github.com/matteoferla/DirEvo_tools',
                target = "_blank",
                tags$p(
                  tags$i(
                    class = "glyphicon glyphicon-info-sign",
                    style = "color:#0072B2;",
                    title = 'Primers are calculated with DirEvo tools\' deep_mut_scanning function: https://github.com/matteoferla/DirEvo_tools'
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
  )
)

