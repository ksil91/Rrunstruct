


#### Reading the stdout files  ####

#' slurp out the traces from a file and return them as a data frame with headers
#' @param file the name of the stdout file to slurp
#' @param K the value of K that this was run at
#' @param rep Which rep number was this structure run.
slurp_traces <- function(file, K, rep) {
  x <- readLines(file)
  x <- x[1:which(str_detect(x, "MCMC completed"))]  # ditch part of file after traces

  # pick out lines that have a sweep number followed by a colon
  # and then remove the colon and the leading and trailing spaces
  # and also toss the LnPD "0", if it on there.
  x2 <- x[str_detect(x, "^ *[0-9]+:")]  %>%
    str_replace_all(":", "") %>%
    str_trim %>%
    str_replace_all(" +0$", "")


  # get the header and fix spaces in names, etc
  header <- x[str_detect(x, "Rep#:")][1] %>%
    str_replace_all("Ln Like", "LnLike") %>%
    str_replace_all("Rep#:", "Sweep") %>%
    str_trim

  # then catenate those and read.table them through a textConnection, add the K to
  # it and tbl_df it.
  tbl_df(cbind(K = K, Rep = rep, read.table(textConnection(c(header, x2)), header = TRUE, na.strings = "--")))
}


#' make the trace data frame tidy and long format
gather_traces <- function(df) {
  gather(df, "variable", "value", -(K:Sweep))
}



#' read all the traces in from all the files in the specified output directory
#'
#' this function also puts everything into long (tidy) format
#'
#' @param dir  the directory holding all the text output from the structure_run() function
#' @export
read_traces <- function(D) {
  files <- dir(D, pattern = "*stdout*", full.names = TRUE)  # files to slurp



  # cycle over files, slurp em in, and rbind it all at the end
  lapply(files, function(x) {
    # some stuff to get the K and rep index
    fn <- str_replace(x, "^.*stdout_K", "")
    K <- str_split(fn, "_")[[1]][2]
    rep <- str_split(fn, "_")[[1]][4]

    slurp_traces(x, K = K, rep = rep) %>%
      gather_traces
  }) %>%
    do.call(rbind, .)

}


#### Functions to read the results files ####

slurp_results <- function(file, K, Rep) {
  # i/o checks
  if (!is.na(file) & !is.character(file)) stop("filename must be a string.")

  s_f <- readr::read_lines(file)
  s_f <- s_f[s_f!=""]

  # Check for predefined populations
  p <- grep("^Proportion of membership.*", s_f)
  pre_pop <- TRUE
  if (length(p) == 0) {
    p = grep("^Overall proportion of membership*", s_f)
    pre_pop <- FALSE
  }
  # Check for USEPOPINFO
  if (length(grep("^USEPOPINFO", s_f)) ==0) {
    upi <- FALSE
  } else { upi <- TRUE}

  # Get cluster memberships and number of predefined pops
  mem_lines <- s_f[(p+3):(grep("^Allele-freq", s_f)-2)]
  mem_lines <- str_trim(mem_lines)
  pops=NULL
  if (!pre_pop) {
    mem_lines <- str_split_fixed(mem_lines, "\\s+", n=as.integer(K))
    mem_lines <- t(mem_lines)
    class(mem_lines) <- "numeric"
    mem_df <- data.frame("Cluster" = mem_lines[,1],
                         "Proportion" = mem_lines[,2],
                         stringsAsFactors = FALSE)
  } else {
    mem_lines <- str_split_fixed(mem_lines, "\\s+", n=as.integer(K)+2)
    mem_lines[,1] <- str_replace(mem_lines[,1], ":", "")
    mem_df <- data.matrix(data.frame(mem_lines[-1, ], stringsAsFactors = FALSE))
    colnames(mem_df) <- mem_lines[1,]
    pops = nrow(mem_df)
  }

  # Get inferred ancestry of individuals

  if (!upi){
  ances_lines <- s_f[(grep("^Inferred ancestry of.*", s_f)+1):(grep("^Estimated Allele Frequencies .*", s_f)-1)]
  ances_lines <- str_trim(ances_lines)
  } else {
    ances_lines <- s_f[(grep("^Inferred ancestry of.*", s_f)+2):(grep("^Estimated Allele Frequencies .*", s_f)-1)]
    ances_lines <- str_trim(ances_lines)
  }
  # extract genotype missing proportions for individuals
  header <- gsub("\\(|\\)|:","", ances_lines[1], " ")
  header <- str_split(header, " ")[[1]] %>% str_replace_all("[^a-zA-Z ]","")
  header <- c("Index",header)

  ances_lines <- ances_lines[-1]
  missing_proportions <- as.numeric(gsub("[\\(\\)]", "",
                                         regmatches(ances_lines, regexpr("\\(.*?\\)", ances_lines))))
  sample_label <- str_trim(gsub("\\(", "",
                                regmatches(ances_lines, regexpr(".*\\(", ances_lines))))
  sample_label <- str_split_fixed(sample_label, "\\s+", n = 2)

  if (!pre_pop) {
    sample_summary <- data.frame(sample_label, missing_proportions)
    split_n <- 3
  } else {
    population_assignment <- as.integer(str_trim(gsub("\\)|:", "",
                                                      regmatches(ances_lines, regexpr("\\).*?:", ances_lines)))))
    sample_summary <- data.frame(sample_label, missing_proportions, population_assignment)
    split_n <- 4
  }

  if (upi) {
    ancest_matrix <- gsub(":  ", "", regmatches(ances_lines, regexpr("..:(.*)", ances_lines)))
    ancest_matrix <- str_split_fixed(ancest_matrix, "\\s+", n=as.integer(K)+1)
    suppressWarnings(class(ancest_matrix) <- "numeric")

    x <- ancest_matrix[!complete.cases(ancest_matrix),]
    for(i in 1:ncol(x)-1){
      v = rep(0,ncol(x))
      v[1] <-i
      v[i+1] <- 1
      x[which(x[,1] == i),] <-  rep(v,each = nrow(x[which(x[,1] == i),]))
    }
    ancest_matrix[!complete.cases(ancest_matrix),] <- x
    ancest_matrix <- ancest_matrix[,-1]
  } else {
    ancest_matrix <- gsub(":  ", "", regmatches(ances_lines, regexpr(":(.*)", ances_lines)))
    ancest_matrix <- str_split_fixed(ancest_matrix, "\\s+", n=as.integer(K))
    suppressWarnings(class(ancest_matrix) <- "numeric")

  }

  ancest_df <- data.frame(sample_summary,
                          ancest_matrix, stringsAsFactors = FALSE)
  colnames(ancest_df)[1:split_n] <- header[1:split_n]
  colnames(ancest_df)[(split_n+1):ncol(ancest_df)] <- seq(1,ncol(ancest_matrix))

  tbl_df(cbind(K = K, Rep = Rep, ancest_df))
}


#' makes a long format data frame of the results
gather_results <- function(df) {
 gather_(df, "cluster", "probability", gather_cols = names(df)[str_detect(names(df), "^[1-9]")])
}




#' read all the results in from all the files in the specified output directory
#'
#' this function also puts everything into long (tidy) format
#'
#' @param D  the directory holding all the text output from the structure_run() function
#' @export
read_results <- function(D) {
  files <- dir(D, pattern = "*results*", full.names = TRUE)  # files to slurp



  # cycle over files, slurp em in, and rbind it all at the end
  lapply(files, function(x) {
    # some stuff to get the K and rep index
    fn <- str_replace(x, "^.*results_K", "")
    K <- str_split(fn, "_")[[1]][2]
    Rep <- str_split(fn, "_")[[1]][4]

    slurp_results(x, K = K, Rep = Rep) %>%
      gather_results
  }) %>%
    do.call(rbind, .)

}


