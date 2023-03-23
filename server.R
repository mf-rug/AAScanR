server <- function(input, output, session) {
  
  primers <- eventReactive(input$aafiles,{
    msgp <- ''
    filecount <- nrow(input$aafiles)
    filenames <- tools::file_path_sans_ext(input$aafiles$name) %>% sort()
    if (filecount == 20) {
      msgp <- paste(msgp, '\n', "✓ Correctly selected 20 files","<br>")
      if (identical(filenames, sort(AA_STANDARD))) {
        msgp <- paste(msgp, '\n', paste("✓ Filenames match", paste(sort(AA_STANDARD), collapse = ' ')),"<br>")
        prim.l <- list()
        for (i in 1:filecount) {
          aa <- file_path_sans_ext(input$aafiles[i, 'name'])
          df <- read_delim(input$aafiles[i, 'datapath'], delim = ' ', col_names = c('name', 'sequence'), show_col_types = F, skip_empty_rows = T) %>% as.data.frame()
          # get rid of rows with only NA
          prim.l[[aa]] <- df[!!rowSums(!is.na(df)),]
        }
        if (lapply(prim.l, nrow) %>% unlist() %>% unique() %>% length() == 1) {
          msgp <- paste(msgp, '\n', '✓ All input files contain the same number of primers',"<br>")
        } else {
           msgp <- paste(msgp, '\n', '× ERROR: Input files do not contain the same number of primers',"<br>")
           msgp <- paste(msgp, '\n', lapply(prim.l, nrow) %>% unlist() %>% paste(., collapse = ' '),"<br>")
        }
      } else {
         msgp <- paste(msgp, '\n', paste("× ERROR: Filenames do not match", paste(sort(AA_STANDARD), collapse = ',')),"<br>")
         msgp <- paste(msgp, '\n', paste(filenames, collapse = ','),"<br>")
      }
    } else {
      msgp <- paste(msgp, '\n', "× ERROR: Did not select 20 files, instead got:")
      msgp <- paste(msgp, '\n', filecount,"<br>")
    }
    if (exists('prim.l')) {
      list(prim.l, msgp)
    } else {
      list(NA, msgp)
    }
  })
  
  muts <- eventReactive(input$mutfile,{
    msg <- ''
    filename <- file_path_sans_ext(input$mutfile$name)
    df <- read_csv(input$mutfile$datapath, col_names = 'Mutation', show_col_types = F, skip_empty_rows = T)
    if (all(sapply(df$Mutation, function(x) {
      str_detect(x, paste0(
        '^[',
        paste0(AA_STANDARD, collapse = ''),
        '][0-9]+[',
        paste0(AA_STANDARD, collapse = ''),
        ']$'
      ))
    }))) {
      msg <- paste(msg, '\n', '✓ Mutant input file correct')
      df$aa <- str_extract(df$Mutation, '^[A-Z]')
      df$num <- str_extract(df$Mutation, '[0-9]+')
      df$mut <- str_extract(df$Mutation, '[A-Z]$')
      list(df, msg)
    } else {
      msg <- paste(msg, '\n', '× ERROR: Mutant input file incorrect')
      msg <- paste(msg, '\n', sapply(df$Mutation, function(x)
       str_detect(
         x, paste0(
           '^[',
           paste0(AA_STANDARD, collapse = ''),
           '][0-9]+[',
           paste0(AA_STANDARD, collapse = ''),
           ']$'
         )
       )) %>% head(10) %>% paste(., collapse = ' '))
      list(NA, msg)
    }
})
  
  out <- reactive({
    if (is.list(primers()[[1]]) && !all(is.na(muts()[[1]]))) {
      enable('download')
      df <- muts()[[1]]
      prim.l <- primers()[[1]]
      for (i in 1:nrow(df)) {
        aa <- df[i, 'mut'] %>% as.character()
        prim.df <- prim.l[[aa]]
        df[i, 'fw'] <- prim.df[prim.df$name == paste0(df[i, 'num'],'_F'), 'sequence']
        df[i, 'rv'] <- prim.df[prim.df$name == paste0(df[i, 'num'],'_R'), 'sequence']
      }
      df
    }
  })
  
  output$download <- downloadHandler(
      filename = function() {
        paste('data-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write_tsv(out(), con)
      }
    )
  
  output$info <- renderUI({
    HTML(
      primers()[[2]], '<br>',
      muts()[[2]]
    )
  })
  output$table <- renderDT({
    out()
  })
}

