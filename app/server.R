# create venv if on shinyapps.io
# code from https://stackoverflow.com/questions/74637485/shinyapps-says-python-package-not-installed-after-just-installing-it
# if (Sys.info()[['user']] == 'shiny'){
  # When running on shinyapps.io, create a virtualenv
  envs<-reticulate::virtualenv_list()
  if(!'venv_aascanr_app' %in% envs)
  {
    reticulate::virtualenv_create(envname = 'venv_aascanr_app',
                                  python = '/opt/homebrew/bin/python3')
    reticulate::virtualenv_install('venv_aascanr_app',
                                   packages = c('Bio'))
  }
  # https://github.com/ranikay/shiny-reticulate-app
  # Set environment BEFORE this
  reticulate::use_virtualenv('venv_aascanr_app', required = TRUE)
  
# }

aas <- c("A","A","A","A","R","R","R","R","R","R","N","N","D","D","C","C","Q","Q","E","E","G","G","G","G","H","H","*","*","*","I","I","I","L","L","L","L","L","L","K","K","M","F","F","P","P","P","P","S","S","S","S","S","S","T","T","T","T","W","Y","Y","V","V","V","V")
cods <- c("GCT","GCC","GCA","GCG","CGT","CGC","CGA","CGG","AGA","AGG","AAT","AAC","GAT","GAC","TGT","TGC","CAA","CAG","GAA","GAG","GGT","GGC","GGA","GGG","CAT","CAC","TAA","TAG","TGA","ATT","ATC","ATA","CTT","CTC","CTA","CTG","TTA","TTG","AAA","AAG","ATG","TTT","TTC","CCT","CCC","CCA","CCG","TCT","TCC","TCA","TCG","AGT","AGC","ACT","ACC","ACA","ACG","TGG","TAT","TAC","GTT","GTC","GTA","GTG")

new_translate <- function(x) {
  p <- substring(x, 1, str_count(x) -(str_count(x) %% 3))
  p <- strsplit(p, "")[[1]]
  p <- paste0(p[c(TRUE, FALSE, FALSE)], p[c(FALSE, TRUE, FALSE)], p[c(FALSE, FALSE, TRUE)])
  r <- grep('[^AGTC]', p)
  p <- aas[match(p, cods)]
  p[r] <- "X"
  paste(p, collapse = "")
}

translate_frames <- function(dna_sequence) {
  # Store results for each frame
  translations <- list()
  
  # Length of the DNA sequence
  dna_length <- nchar(dna_sequence)
  
  # Process each frame
  for (frame in 1:3) {
    start <- frame
    # Initialize amino acid sequence for this frame
    aa_sequence <- c()
    
    # Translate each codon in this frame
    while (start <= dna_length - 2) {
      codon <- substr(dna_sequence, start, start + 2)
      if (nchar(codon) == 3) { # Ensure full codon
        aa_sequence <- c(aa_sequence, new_translate(codon))
      }
      start <- start + 3
    }
    
    # Store the translated sequence for this frame
    translations[[frame]] <- aa_sequence
  }
  
  return(translations)
}


# Function to find the longest ORF in a list of amino acid sequences
find_longest_orf <- function(translations) {
  # Function to find longest ORF in a single amino acid sequence
  longest_orf_in_frame <- function(aa_seq) {
    start_positions <- which(aa_seq == "M")
    end_positions <- which(aa_seq == "*")
    
    # Initialize variables to track the longest ORF
    longest_length <- 0
    longest_orf <- NULL
    
    # Iterate over all start positions
    for (start in start_positions) {
      # Find the first end position that is greater than the start position
      valid_ends <- end_positions[end_positions > start]
      if (length(valid_ends) > 0) {
        end <- valid_ends[1]
        orf_length <- end - start + 1
        
        # Check if this ORF is the longest found so far
        if (orf_length > longest_length) {
          longest_length <- orf_length
          longest_orf <- list(start = start, end = end)
        }
      }
    }
    
    return(longest_orf)
  }
  
  # Apply the function to each frame and store results
  longest_orfs <- lapply(translations, longest_orf_in_frame)
  
  return(longest_orfs)
}

# Function to find the overall longest ORF across all frames and convert to nucleotide positions
find_overall_longest_orf <- function(dna_seq) {
  translations <- translate_frames(dna_seq)
  longest_orfs <- find_longest_orf(translations)
  
  # Initialize a variable to store the longest ORF info
  overall_longest <- list(length = 0, frame = NULL, start = NULL, end = NULL, nt_start = NULL, nt_end = NULL)
  
  # Check each frame's longest ORF and find the overall longest
  for (i in seq_along(longest_orfs)) {
    orf <- longest_orfs[[i]]
    if (!is.null(orf)) {
      orf_length <- orf$end - orf$start + 1
      if (orf_length > overall_longest$length) {
        overall_longest$length <- orf_length
        overall_longest$frame <- i
        overall_longest$start <- orf$start
        overall_longest$end <- orf$end
      }
    }
  }
  
  # Convert amino acid positions back to nucleotide positions
  if (!is.null(overall_longest$frame)) {
    overall_longest$nt_start <- (overall_longest$start - 1) * 3 + overall_longest$frame
    overall_longest$nt_end <- (overall_longest$end - 1) * 3 + overall_longest$frame + 2
    overall_longest$seq <- str_sub(dna_seq, overall_longest$nt_start, overall_longest$nt_end)
  }
  
  return(overall_longest)
}



parse_data_to_df <- function(text) {
  # Extract key-value pairs and sequences
  key_vals <- str_extract_all(text, "\\b(?:[a-zA-Z]+=[^ ]+|\\b[acgtACGTnNxX]+\\b)")
  key_vals <- lapply(key_vals, function(x) x[-1])
  
  # add 'f' and 'r' to distinguish data for fw and rv
  key_vals[[1]][str_detect(key_vals[[1]], '=')] <- paste0('f', key_vals[[1]][str_detect(key_vals[[1]], '=')])
  key_vals[[2]][str_detect(key_vals[[2]], '=')] <- paste0('r', key_vals[[2]][str_detect(key_vals[[2]], '=')])

  # add a key to the primer seq, which AAScan spits out without key
  key_vals[[1]][!str_detect(key_vals[[1]], '=')] <- paste0('fw_primer=', key_vals[[1]][!str_detect(key_vals[[1]], '=')])
  fkv <- str_split(key_vals[[1]], ' ') %>% unlist() %>% str_split('=')
  fwdat <- sapply(fkv, function(x) x[2])
  names(fwdat) <- sapply(fkv, function(x) x[1])
  
  # same for rv
  key_vals[[2]][!str_detect(key_vals[[2]], '=')] <- paste0('rv_primer=', key_vals[[2]][!str_detect(key_vals[[2]], '=')])
  rkv <- str_split(key_vals[[2]], ' ') %>% unlist() %>% str_split('=')
  revdat <- sapply(rkv, function(x) x[2])
  names(revdat) <- sapply(rkv, function(x) x[1])
  
  df <- cbind(fwdat %>% as.matrix() %>% t() %>% as.data.frame(),
              revdat %>% as.matrix() %>% t() %>% as.data.frame())
  return(df)
}


# server
server <- function(input, output, session) {
  muts <- reactiveVal(NULL)
  muts_good <- reactiveVal(NULL)
  muts_msg <- reactiveVal(NULL)
  allsubstr <- function(x, n) unique(substring(x, 1:(nchar(x) - n + 1), n:nchar(x)))
  AA_STANDARD <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  dms <- reticulate::import("deep_mut_scanning_mf")
  aas <- reticulate::import("pyAAscan")
  cods <- data.frame('AA' = c("A","D","E","F",  "C",  "G",  "H",  "I",  "K",  "L",  "M",  "N",  "P",  "Q",  "R",  "S",  "T",  "V",  "W",  "Y"),
             'cod1' = c("GCG","GAT","GAA","TTT","TGC","GGC","CAT","ATT","AAA","CTG","ATG","AAC","CCG","CAG","CGT","AGC","ACC","GTG","TGG","TAT"),
             'cod2' = c("GCC","GAC","GAG","TTC","TGT","GGT","CAG","ATC","AAG","TTA","ATG","AAT","CCA","CAA","CGC","TCT","ACG","GTT","TGG","TAC"))

  full_data <- reactive({
    if (str_detect(input$seq, '^[ATGCatgc]+$') && str_detect(input$start, '^[0-9]+$') && str_detect(input$end, '^[0-9]+$')) {
      if (as.numeric(input$start) < as.numeric(input$end) && as.numeric(input$start) > 25 && as.numeric(input$end) < (str_length(input$seq) -25)) {
        full <- data.frame()
        full_old <- data.frame()
        withProgress(
          message = 'Calculating primers', value = 0, {
            # loop over i 20 AAs
            for (i in 1:nrow(cods)) {

              out.df <- data.frame()
              # loop over j codons in the selected sequence
              for (j in 1:(floor((as.integer(input$end) - as.integer(input$start) -1) /3))) {
                out <- aas$Mutate(seq_in = input$seq,
                         cod1pos = as.integer(as.numeric(input$start) -1),
                         mutpos = as.integer(j),
                         codon1 = cods[i, 'cod1'],
                         codon2 = cods[i, 'cod2'])
                add_df <- parse_data_to_df(out)
                cod_end <- as.numeric(input$start) -1 + (j*3)
                codon <- str_sub(input$seq,cod_end -2, cod_end)
                add_df$codon <- codon
                # add_df$AA <- Biostrings::translate(DNAString(codon)) %>% as.character()
                add_df$AA <- new_translate(codon)
                add_df$pos <- j
                out.df <- rbind(out.df, add_df)
              }
              out.df$mut <- cods[i, 'AA']
              out.df$mutation <- paste0(out.df$AA, out.df$pos, out.df$mut)
              full <- rbind(full, out.df)
              incProgress(1/nrow(cods), detail = HTML("current amino acid is", cods[i, 'AA']))
            }
          }
        )
        enable('download_prims')
      } else { disable('download_prims') }
    } else { disable('download_prims') }
    if (exists('full')) {
      if (!is.null(muts_good()) && muts_good()) {
        full_subset <- full[full$mutation %in% muts()$Mutation,]
        if (nrow(full_subset) == 0) {
          muts_msg('× Warning: none of the supplied mutations matched the selected sequence!')
        } else if (nrow(full_subset) < nrow(muts())) {
          muts_msg('× Warning: not all the supplied mutations matched the selected sequence!')
        }
        full_subset
      } else {
        full
      }
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
              selected_nts <- as.numeric(input$end) - as.numeric(input$start) +1
              if ((selected_nts /3) %% 1 == 0) {
                out <- paste0(selected_nts /3, ' (✓ ', selected_nts, ' selected nucleotides are a multiple of 3)')
              } else {
                out <- paste0(round(selected_nts /3, 0), ' (× ', selected_nts, ' selected nucleotides are NOT a multiple of 3)')
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
  observeEvent(input$seq,{
    if (str_detect(input$seq, '^[ATGCatgc]+$')) {
      enable('orf')
    } else {
      disable('orf')
    }
  })
  observeEvent(input$orf,{
    if (str_detect(input$seq, '^[ATGCatgc]+$')) {
      # u <- str_extract(input$seq, regex('atg(?:[atgc][atgc][atcg])*?(?:taa|tga|tag)', ignore_case = TRUE))
      # updateTextInput('start', session = getDefaultReactiveDomain(), value = str_locate(input$seq, u)[1,]['start'] %>% as.numeric() -1)
      # updateTextInput('end', session = getDefaultReactiveDomain(), value = str_locate(input$seq, u)[1,]['end'] %>% as.numeric() -3)
      u <- find_overall_longest_orf(input$seq)
      updateTextInput('start', session = getDefaultReactiveDomain(), value = u$nt_start)
      updateTextInput('end', session = getDefaultReactiveDomain(), value = u$nt_end -3)
    }
  }, ignoreInit = TRUE)
  
  output$prim_table <- renderDT({
    full <- full_data()
    if (!all(is.na(full))) {
      datatable(
        full[,c('mutation', 'codon', 'fw_primer', 'flen', 'fTm','fTmfull', 'fGCcontent', 'rv_primer', 'rlen', 'rTm','rTmfull', 'rGCcontent')],
        # full[,c('base', 'mutation', 'homology_Tm', 'fw_primer', 'fw_len_primer', 'fw_anneal_Tm','rv_primer', 'rv_len_primer', 'rv_anneal_Tm' )],
        options = list('pageLength' = 50, lengthMenu = c(10,25,50,100,200)),
        rownames= TRUE
      )
    }
  })
  
  output$plot <- renderPlot({
    full <- full_data()
    if (!all(is.na(full))) {
      df <- data.frame('num' = 1:str_length(input$seq))
      df$nt <- str_split(input$seq, '')[[1]]
      df$target <- 'no'
      df[(as.numeric(input$start) -1):as.numeric(input$end), 'target'] <- 'yes'
      p <- ggplot(df, aes(x=num, y=1)) +
        # geom_point(show.legend = F) +
        geom_line(aes(x = num - 0.5, group = target, color = target), linewidth = 10, show.legend = F) + 
        geom_line(aes(x = num + 0.5, group = target, color = target), linewidth = 10, show.legend = F) + 
        scale_color_manual(values = c('yes' = 'grey', 'no' = 'grey96')) +
        labs(x = NULL, y = NULL) +
        annotate('text', x=1, y=0, label = '1', size = 2.5) +
        annotate('text', x=str_length(input$seq), y=0, label = str_length(input$seq), size = 2.5, hjust=0.85) +
        annotate('text', x=as.numeric(input$start) -1, y=0, label = input$start, size = 2.5) +
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
    updateTextAreaInput('seq', session = getDefaultReactiveDomain(), value = 'TTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGAACACGATTAACATCGCTAAGAACGACTTCGAGATGGGTGAAGCACGCTTCGAGGTTGCGGATAACGCTGCCGCCTAACTTGGGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATAGC')
    updateTextInput('start', session = getDefaultReactiveDomain(), value = 43)
    updateTextInput('end', session = getDefaultReactiveDomain(), value = 68)
  })

  
  output$download_prims <- downloadHandler(
    filename = function() {
      paste0('AAScanR_primers_', Sys.Date(), '.tsv')
    },
    content = function(con) {
      df <- full_data()
      muts <- muts()$Mutation
      if (!is.null(muts)) {
        df <- df[df$mutation %in% muts, ]
      }
      df <- df[, c("mutation", "fw_primer", "rv_primer")]
      write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE)
    },
    contentType = "text/tab-separated-values"
  )
  
  
  observeEvent(input$mutinput,{
    msg <- ''
    df <- data.frame('Mutation' = str_split(input$mutinput, '\n')[[1]])
    if (all(sapply(df$Mutation, function(x) {
      str_detect(x, paste0(
        '^[',
        paste0(AA_STANDARD, collapse = ''),
        '][0-9]+[',
        paste0(AA_STANDARD, collapse = ''),
        ']$'
      ))
    }))) {
      msg <- paste(msg, '\n', '✓ Mutant input correct')
      df$aa <- str_extract(df$Mutation, '^[A-Z]')
      df$num <- str_extract(df$Mutation, '[0-9]+')
      df$mut <- str_extract(df$Mutation, '[A-Z]$')
      muts_good(TRUE)
    } else {
      msg <- paste(msg, '\n', '× ERROR: Mutant input incorrect. Showing all Primers.')
      muts_good(FALSE)
    }
    muts(df)
    muts_msg(msg)
  },ignoreInit = TRUE)
  
  
  observeEvent(input$mutfile,{
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
      muts_good(TRUE)
    } else {
      msg <- paste(msg, '\n', '× ERROR: Mutant input file incorrect. Showing all Primers.')
      muts_good(FALSE)
    }
    muts(df)
    muts_msg(msg)
  },ignoreInit = TRUE)
  
  
  output$info <- renderUI({
    HTML(
      muts_msg()
    )
  })
  
}

