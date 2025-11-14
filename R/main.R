#' Encode DNA Sequences into Feature Matrix
#'
#' This function converts DNA 5-mer sequences into a position-wise categorical feature matrix.
#'
#' @param dna_strings A character vector of DNA 5-mer strings.
#' @return A data frame with nucleotide positions (nt_pos1...nt_posN) as columns and factor levels (A, T, C, G).
#' @examples
#' dna_encoding(c("ATCGA", "TGCAT"))
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A Sites for Multiple Samples
#'
#' Predict m6A modification probability and status for multiple samples using a trained ML model.
#'
#' @param ml_fit A trained random forest model (e.g. rf_fit.rds).
#' @param feature_df A data frame containing features including gc_content, RNA_type, RNA_region,
#' exon_length, distance_to_junction, evolutionary_conservation, and DNA_5mer.
#' @param positive_threshold Numeric. Threshold above which samples are labeled as "Positive" (default = 0.5).
#' @return A data frame with predicted probabilities and labels added.
#' @examples
#' \dontrun{
#' mdl <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' df  <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' prediction_multiple(mdl, df, positive_threshold = 0.6)
#' }
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  feature_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  prob <- predict(ml_fit, feature_df, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_prob <- prob
  feature_df$predicted_m6A_status <- ifelse(prob > positive_threshold, "Positive", "Negative")
  return(feature_df)
}


#' Predict m6A Sites for a Single Sample
#'
#' This function predicts m6A probability and classification for a single RNA feature input.
#'
#' @param ml_fit A trained random forest model (e.g. rf_fit.rds).
#' @param gc_content GC content (numeric).
#' @param RNA_type RNA type (character: "mRNA", "lincRNA", "lncRNA", or "pseudogene").
#' @param RNA_region RNA region (character: "CDS", "intron", "3'UTR", or "5'UTR").
#' @param exon_length Exon length (numeric).
#' @param distance_to_junction Distance to exon–intron junction (numeric).
#' @param evolutionary_conservation Evolutionary conservation score (numeric).
#' @param DNA_5mer DNA 5-mer sequence (character).
#' @param positive_threshold Threshold for classifying "Positive"/"Negative" (default = 0.5).
#' @return A named vector with predicted probability and status.
#' @examples
#' \dontrun{
#' mdl <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' prediction_single(
#'   mdl,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "ATCGT",
#'   positive_threshold = 0.5
#' )
#' }
#' @export
#' @import randomForest
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5){
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
    RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  pred_df <- prediction_multiple(ml_fit, feature_df, positive_threshold)
  returned_vector <- c(
    predicted_m6A_prob = pred_df$predicted_m6A_prob[1],
    predicted_m6A_status = pred_df$predicted_m6A_status[1]
  )
  return(returned_vector)
}
