server <- function(input, output, session) {
  AA_STANDARD <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  dms <- import("deep_mut_scanning_mf")
  cods <- data.frame('AA' = c("A","D","E","F","C","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"),
             'cod' = c("GCG","CGC","AAC","GAT","TGC","GAA","CAG","GGC","CAT","ATT","CTG","AAA","ATG","TTT","CCG","AGC","ACC","TGG","TAT","GTG"))
  full_data <- reactive({
    if (str_detect(input$seq, '^[ATGCatgc]+$') && str_detect(input$start, '^[0-9]+$') && str_detect(input$end, '^[0-9]+$')) {
      if (as.numeric(input$start) < as.numeric(input$end) && as.numeric(input$start) > 25 && as.numeric(input$end) < (str_length(input$seq) -25)) {
        full <- data.frame()
        withProgress(
          message = 'Calculating primers', value = 0, {
            for (i in 1:nrow(cods)) {
              out <- dms$deep_mutation_scan(input$seq, c(as.integer(input$start),as.integer(input$end)), mutation = cods[i, 'cod'], overlap_len=as.integer(15))
              out.df <-
                lapply(out, function(x) lapply(x, as.character) %>% unlist() %>% as.data.frame)
              out.df <- do.call('cbind', out.df) %>% t() %>% as.data.frame()
              out.df$mut <- cods[i, 'AA']
              out.df$mutation <- paste0(out.df$AA, out.df$mut)
              full <- rbind(full, out.df)
              incProgress(1/nrow(cods), detail = HTML("<br>current amino acid is", cods[i, 'AA']))
            }
          }
        )
        enable('download_prims')
      } else { disable('download_prims') }
    } else { disable('download_prims') }
    if (exists('full')) {
      full
    } else {
      NA
    }
  })
  
  output$codnum <- renderUI({
    if (str_detect(input$seq, '^[ATGCatgc]+$')) {
      if (str_detect(input$start, '^[0-9]+$') && str_detect(input$end, '^[0-9]+$')) {
        if (as.numeric(input$start) < as.numeric(input$end)) {
          if (as.numeric(input$start) > 25) {
            if (as.numeric(input$end) < (str_length(input$seq) -25)) {
              selected_nts <- as.numeric(input$end) - as.numeric(input$start) + 1
              if ((selected_nts /3) %% 1 == 0) {
                out <- paste0(selected_nts /3, ' (✓ selected nucleotides a multiple of 3)')
              } else {
                out <- paste0(round(selected_nts /3, 0), ' (× selected nucleotides NOT a multiple of 3)')
              } 
            } else {
              out <- '× Error: need at least 25 nucleotides downstream of end of mutated region'
            }
          } else {
            out <- '× Error: need at least 25 nucleotides upstream of start of mutated region'
          }
        } else {
          out <- '× Error: end number smaller than start number'
        }
      } else {
        out <- 'Enter numbers for start and end'
        } 
      } else {
        out <- '× Error: input sequence can only contain any of: A T G C a t g c'
      }
    
    full <- full_data()
    if (!all(is.na(full))) {
      warn <- !sum(!str_detect(full$AA, '^\\*')) == length(full$AA)
    } else {
      warn <- NA
    }
    if (is.na(warn) || !warn) {
      HTML('<i>',out,'</i>')
    } else {
      HTML('<i>',out,'<br>× WARNING: Stop codons deteceted! Selected the wrong frame? (Check start codon)</i>')
    }
  })
  
  output$prim_table <- renderDT({
    full <- full_data()
    if (!all(is.na(full))) {
      datatable(
        full[,c('base', 'mutation', 'homology_Tm', 'fw_primer', 'fw_len_primer', 'fw_anneal_Tm','rv_primer', 'rv_len_primer', 'rv_anneal_Tm' )],
        options = list('pageLength' = 50, lengthMenu = c(10,25,50,100,200)),
        rownames= FALSE
      )
    }
  })
  
  output$plot <- renderPlot({
    full <- full_data()
    if (!all(is.na(full))) {
      df <- data.frame('num' = 1:str_length(input$seq))
      df$nt <- str_split(input$seq, '')[[1]]
      df$target <- 'no'
      df[input$start:input$end, 'target'] <- 'yes'
      p <- ggplot(df, aes(x=num, y=1)) +
        # geom_point(show.legend = F) +
        geom_line(aes(x = num - 0.5, group = target, color = target), linewidth = 10, show.legend = F) + 
        geom_line(aes(x = num + 0.5, group = target, color = target), linewidth = 10, show.legend = F) + 
        scale_color_manual(values = c('yes' = 'grey', 'no' = 'grey96')) +
        labs(x = NULL, y = NULL) +
        annotate('text', x=1, y=0, label = '1', size = 2.5) +
        annotate('text', x=str_length(input$seq), y=0, label = str_length(input$seq), size = 2.5, hjust=0.85) +
        annotate('text', x=as.numeric(input$start), y=0, label = input$start, size = 2.5) +
        annotate('text', x=as.numeric(input$end), y=0, label = input$end, size = 2.5) +
        scale_y_continuous(limits = c(-1,2), expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        theme_void() 

      if (str_length(input$seq) < 2000) {
        p <- p +
          geom_text(aes(label = nt), size = 500 / str_length(input$seq) , show.legend = F)
      }
      p
    }
  })
  
  observeEvent(input$example, {
    updateTextAreaInput('seq', session = getDefaultReactiveDomain(), value = 'ATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCCTGGTGCCGCGCGGCAGCCATTGACTTGGGCCCGAACAAAAACTCATCTCAGAAGAGGA')
    updateTextInput('start', session = getDefaultReactiveDomain(), value = 112)
    updateTextInput('end', session = getDefaultReactiveDomain(), value = 171)
  })
  
  observe({
    full <- full_data()
    if (!all(is.na(full))) {
      for (aa in cods$AA) {
        writeout <- full[full$mut == cods$AA[1],c('mutation', 'fw_primer', 'rv_primer')] %>% pivot_longer(cols=-1) %>% select(2:3)
        writeout$name <- paste0(rep(1:(nrow(writeout) /2), each=2), c('_F', '_R'))
        write_delim(writeout, file.path(tempdir(),paste0(aa, '.txt')), col_names = FALSE, delim = ' ')
      }
    }
  })
  
  output$download_prims <- downloadHandler(
    filename = function() {
      paste("output", "zip", sep = ".")
      paste0('AAScanR_all_primers_', Sys.Date(), '.zip')
    },
    content = function(con) {
      fs <- c()
      for (i in 1:length(cods$AA)) {
        path <- file.path(tempdir(), paste0(cods[i, 'AA'], '.txt'))
        fs <- c(fs, path)
      }
      zip(zipfile = con, files = fs, mode="cherry-pick")
    },
    contentType = "application/zip"
  )

  
  primers <- eventReactive(input$aafiles,{
    msgp <- ''
    filecount <- nrow(input$aafiles)
    filenames <- file_path_sans_ext(input$aafiles$name) %>% sort()
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
      if (any(lapply(prim.l, function(x) str_remove(x$name, '_.$') %>% as.numeric() %>% max()) < max(as.numeric(df$num)))) {
        df <- data.frame('Error' = 'Problem with input: mutation files with residue numbers not provided in the primer files.')
      } else {
        for (i in 1:nrow(df)) {
          aa <- df[i, 'mut'] %>% as.character()
          prim.df <- prim.l[[aa]]
          df[i, 'fw'] <- prim.df[prim.df$name == paste0(df[i, 'num'],'_F'), 'sequence']
          df[i, 'rv'] <- prim.df[prim.df$name == paste0(df[i, 'num'],'_R'), 'sequence']
        }
      }
      df
    }
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste('AAScanR_select_primers_', Sys.Date(), '.tsv', sep='')
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

