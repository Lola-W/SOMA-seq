#' Load Required Libraries
library(Seurat)
library(tidyverse)
library(parallel)
library(biomaRt)
library(dplyr)

#' Read and Process a Single File for checking
#'
#' This function reads a single file (either RDS or CSV) and processes it to extract
#' donor metadata, which includes donor ID, cell subclass, and UMI count.
#' It also generates a barcode from the metadata row names.
#'
#' @param filePath A string containing the path to the file.
#' @return A data frame with columns for donor ID, subclass, UMI count, and barcode.
#' @export
#' @examples
#' readAndProcessFile("path/to/your/file.rds")
readAndProcessFile <- function(filePath) {
  fileExtension <- tools::file_ext(filePath)

  if (fileExtension == "rds") {
    print(filePath)
    seuratObj <- readRDS(filePath)
    metadata <- seuratObj@meta.data
  } else if (fileExtension == "csv") {
    metadata <- read.csv(filePath, row.names = 1)
  } else {
    stop("Unsupported file type: ", fileExtension)
  }

  info_meta <- metadata %>%
    select(donor_id, Subclass, `Number of UMIs`) %>%
    mutate(barcode = sub("-.*", "", rownames(metadata)))

  return(info_meta)
}

#' Process Multiple RDS Files
#'
#' This function processes multiple RDS files within a given directory to extract
#' and combine single cell metadata from each file.
#' 
#' @param directoryPath The directory path containing the RDS files of Seurat objects.
#' @return A combined data frame of metadata from all processed RDS files.
readAndProcessFiles <- function(directoryPath) {
  # List all RDS files in the directory
  rdsFiles <- list.files(directoryPath, pattern = "\\.rds$", full.names = TRUE)
  
  # Initialize an empty list to store metadata from each file
  metadataList <- list()
  
  # Loop over RDS files and process each file
  for (filePath in rdsFiles) {
    # Print the file path being processed
    print(paste("Processing file:", filePath))
    
    # Read the RDS file
    seuratObj <- readRDS(filePath)
    
    # Extract metadata
    metadata <- seuratObj@meta.data
    
    # Process metadata
    info_meta <- metadata %>%
      dplyr::mutate(barcode = sub("-.*", "", rownames(metadata)), read_count = `Number of UMIs`) %>% 
      dplyr::select(barcode, donor_id, Subclass, read_count)
    
    # Add the processed metadata to the list
    metadataList[[length(metadataList) + 1]] <- info_meta
  }
  
  # Combine all metadata into one data.frame
  combinedMetadata <- bind_rows(metadataList)
  
  return(combinedMetadata)
}

#' Process Barcode Information from Monopogen Directories
#'
#' This function aggregates barcode information from specified directories
#' under the Monopogen base directory, aligning it with sample names extracted
#' from the directory structure.
#'
#' @param monopogenBaseDir A string representing the base directory containing Monopogen output.
#' @param barcodeBaseDir A string representing the base directory where barcode count files are stored.
#' @param pattern A string pattern to identify the relevant barcode count files.
#' @return A data frame with barcode, read count, and sample information.
#' @examples
#' processBarcodeInfo("/path/to/monopogen/output/")
processBarcodeInfo <- function(monopogenBaseDir,
                               barcodeBaseDir = "Soma_seq/STARsolo/STAR_results/",
                               pattern = "/soloOut/GeneFull/filtered/barcode_counts.csv"){
  # Read a directory of barcodes
  subDirs <- list.dirs(monopogenBaseDir, full.names = FALSE, recursive = FALSE)
  print(subDirs)
  result <- c()
  for (dir in subDirs) {
    barcode_info <- data.table::fread(paste0(barcodeBaseDir, dir, pattern))
    barcode_info$sample <- dir
    result <- rbind(result,barcode_info)
  }
  result<- result %>%
    rename(barcode = cell, read_count = id)
  return(result)
}

#' Process Sample Information, Helper Function
#'
#' Processes and cleans sample information, typically involving modification
#' of sample file names and splitting on specified delimiters.
#'
#' @param sample_info A data frame containing sample information including file names.
#' @return A data frame with cleaned and separated sample information.
#' @examples
#' processSampleInfo(sample_info)
processSampleInfo <- function(sample_info) {
  # Remove the "_I1_001.fastq.gz" part and split by "|"
  processed <- sample_info %>%
    mutate(sample = gsub("_I1_001\\.fastq\\.gz", "", RNAseq_file_names)) %>%
    separate_rows(sample, sep = "\\|")

  return(processed)
}

#' Read, Combine RDS Files with Donor IDs, and Filter Mutations
#'
#' This function reads multiple RDS files from given directories, combines them,
#' and filters the data based on somatic mutation criteria and a minimum
#' B-allele frequency (BAF).
#'
#' @param baseDir The base directory containing RDS files grouped by sample.
#' @param sample_info A data frame containing sample and donor information.
#' @param BAF The B-allele frequency threshold used for filtering, default 0.3.
#' @return A list of data frames, each corresponding to combined data from one directory.
#' @examples
#' readAndCombineRDSWithDonorId("/path/to/data", sample_info)
readAndCombineRDSWithDonorId <- function(baseDir, sample_info, BAF = 0.3) {
  # Get barcode and sample info, filter out somatic mutation.
  # Process the sample_info data.frame
  processed_sample_info <- processSampleInfo(sample_info)

  # Find all subdirectories in the base directory
  subDirs <- list.dirs(baseDir, full.names = TRUE, recursive = FALSE)
  
  # Initialize a list to temporarily store the data frames
  dataList <- list()

  for (dir in subDirs) {
    # List all files matching the pattern in the directory
    filePaths <- list.files(path = file.path(dir, "somatic"),
                            pattern = "chr.*\\.SNV_mat\\.RDS",
                            full.names = TRUE)
    # Check if there are any files
    if (length(filePaths) > 0) {
      folderName <- basename(dir)
      if (folderName %in% processed_sample_info$sample) {
        # Read all RDS files and store them in a list
        dataListTmp <- lapply(filePaths, function(filePath) {
          # print(filePath)
          data <- readRDS(filePath)
          if (!is.data.frame(data)) {
            warning(paste("File does not contain a data frame:", filePath))
            return(NULL)
          }

          # get somatic mutation
          chrPattern <- gsub("SNV_mat\\.RDS", "putativeSNVs\\.csv", basename(filePath))
          putativeFilePath <- file.path(dir, "somatic", chrPattern)
          putative <- read.csv(putativeFilePath)
          somatic <- putative %>%
            # threshold of 0.3 used, can be changed
            filter(SVM_pos_score > 0.5 & LDrefine_merged_score > 0.25 & BAF_alt < BAF) %>%
            dplyr::rename("ref" = "Ref_allele",
                   "alt" = "Alt_allele")

          common_cols <- c("chr", "pos", "ref", "alt")
          data <- semi_join(data, somatic, by = common_cols)

          # Map folderName to donor_id
          donor_id <- processed_sample_info$external_donor_name[processed_sample_info$sample == folderName]

          print(paste0(filePath, ": ", nrow(data)))
          # Add donor_id column to the data
          data %>% mutate(donor_id = donor_id, sample = folderName) %>%
            dplyr::select(donor_id, sample, everything())
        })
        # print(head(dataListTmp))

        # Combine all data frames in the list into one
        combinedData <- bind_rows(dataListTmp)

        # Add the combined data frame to the dataList
        dataList[[length(dataList) + 1]] <- combinedData
      } else {
        print(paste0("Not found ", folderName))
      }
    } else {
      print(paste0("No matching files in ", dir))
    }
  }

  return(dataList)
}

#' Preprocess Monopogen Data
#'
#' This function preprocesses a list of data tables containing Monopogen output, merges them with
#' barcode and metadata information, and groups the data according to specified factors.
#'
#' @param dataList A list of data frames or data tables; each contains Monopogen output for a sample.
#' @param info_meta A data frame containing metadata for samples such as donor ID and barcode.
#' @param barcode_info A data frame containing read counts and barcode identifiers.
#' @param groupByFactors A character vector specifying the columns to group by after processing.
#' @return A data frame containing the grouped and processed data.
#' @examples
#' preprocessMonopogen(dataList, info_meta, barcode_info, c("donor_id", "sample"))
preprocessMonopogen <- function(dataList, info_meta, barcode_info, groupByFactors) {
  # combine with Mopnopogen data
  # Ensure groupByFactors is a character vector
  if (!is.character(groupByFactors)) {
    stop("groupByFactors must be a character vector.")
  }

  # Initialize an empty data frame for the summary
  groupedData <- data.frame()

  # Process each data frame in the list
  for (i in seq(length(dataList))) {
    # Convert data.table to data.frame if necessary
    combinedDataTable <- as.data.table(dataList[[i]])

    # Melt the combined data table, each row is one variant for one barcode,
    combinedDataTable <- combinedDataTable %>%
      rename(org.genotype = genotype)

    meltedData <- melt(combinedDataTable, id.vars = 1:20, variable.name = "barcode", value.name = "genotype") %>%
      filter(genotype != "0/0" & genotype != "1/1")

    # Join with info_meta
    joinedData <- inner_join(meltedData, info_meta, by = c("barcode","donor_id"))
    # get read count
    joinedData <- left_join(joinedData, barcode_info, by = c("barcode", "sample"))

    groupedData <- rbind(groupedData, joinedData)
  }

  return(groupedData)
}

#' Merge Processed Data with Gene Information
#'
#' This function merges processed variant data with gene information, extracting chromosomal positions
#' and merging on genetic coordinates to attach gene annotations to each variant.
#'
#' @param processed A data frame containing processed variant data with chromosomal identifiers.
#' @param gene_info A data frame containing gene information including locations and alleles.
#' @return A data frame with variants annotated with gene information and calculated gene sizes.
#' @examples
#' mergeProcessedWithGeneInfo(processed, gene_info)
mergeProcessedWithGeneInfo <- function(processed, gene_info) {
  # Parse the Location to extract chr and pos
  gene_info <- gene_info %>%
    mutate(chr = as.integer(gsub(":.+$", "", Location)),
           pos = as.integer(gsub("^.*-", "", Location))) %>%
    filter(Gene != "-") %>%
    distinct(chr, pos, Allele, .keep_all = TRUE)

  # Adjust the chr column to be numeric for comparison in processed
  processed <- processed %>%
    mutate(chr_num = as.integer(gsub("chr", "", chr)))

  # Get variant count
  summarized <- processed %>%
    dplyr::group_by(chr_num, pos, ref, alt, Subclass, donor_id) %>%
    #dplyr::group_by(chr_num, pos, ref, alt) %>%
    dplyr::summarize(freq = n(), .groups = 'drop')

  # Merge the two datasets
  merged_data <- summarized %>%
    left_join(gene_info, by = c("chr_num" = "chr", "pos", "alt" = "Allele")) %>%
    #dplyr::select(chr = chr_num, pos, ref, alt, freq, Gene) %>%
    dplyr::select(chr = chr_num, pos, ref, alt, Subclass, donor_id, freq, Gene) %>%
    filter(Gene != "-")

  # Get gene size
  message("Start getting gene size")
  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="ensembl_gene_id", values=merged_data$Gene, mart=human)
  gene_coords$size=gene_coords$end_position - gene_coords$start_position
  gene_coords <- gene_coords[, c("hgnc_symbol","ensembl_gene_id", "size" )]

  # Merge gene details back into merged_data
  merged_data <- merged_data %>%
    left_join(gene_coords, by = c("Gene" = "ensembl_gene_id"))

  return(merged_data)
}

#' Process RNA editing data and merge with whole genome data
#'
#' This function reads RNA editing data, renames the columns to a standard format,
#' and performs an anti-join operation with a dataset of processed whole genome data.
#' It then merges the result with sample information and returns the selected columns 
#' along with disease information.
#'
#' @param rna_editing_file Path to the RNA editing data file.
#' @param processed_whole A data frame of processed whole genome data.
#' @param sample_info A data frame containing sample information.
#' @return A data frame containing the filtered genome data merged with sample information.
#' @importFrom data.table fread
#' @importFrom dplyr anti_join left_join rename select
#' @export
#' @examples
#' process_rna_editing("path/to/rna_editing_data.txt", processed_whole, sample_info)
process_rna_editing <- function(rna_editing_file, processed_whole, sample_info) {
  RNA_editing <- fread(rna_editing_file)
  colnames(RNA_editing) <- trimws(colnames(RNA_editing))
  RNA_editing <- RNA_editing %>%
    dplyr::rename(chr = Region, pos = Position, ref = Ref, alt = Ed)
  
  result <- anti_join(processed_whole, RNA_editing, by = c("chr", "pos", "ref", "alt"))
  rows_wanted <- colnames(result)
  
  processed_whole_disease <- result %>%
    left_join(sample_info, by = c("donor_id" = "external_donor_name")) %>%
    dplyr::rename(disease = Disease) %>%
    dplyr::select(c(rows_wanted, "disease"))
  
  return(processed_whole_disease)
}

#' Prepare data for VEP input
#'
#' This function takes a data frame with variant information and 
#' transforms it to a format suitable for input to the Variant Effect Predictor (VEP).
#' The resulting data is written to a specified file.
#'
#' @param processed_whole A data frame containing columns chr, pos, ref, and alt with chromosome,
#'        position, reference base, and alternative base information.
#' @param output_file Path to the output file where the transformed data will be saved.
#' @importFrom dplyr select
#' @importFrom data.table fread
#' @export
#' @examples
#' prepare_input_for_vep(processed_whole, "/path/to/output/transformed.csv")
prepare_input_for_vep <- function(processed_whole, output_file) {
  transformedData <- unique(processed_whole[, c("chr", "pos", "ref", "alt")])
  transformedData$chr <- as.numeric(gsub("chr", "", transformedData$chr))
  transformedData <- transformedData[order(transformedData$chr, transformedData$pos)]
  transformedData <- unique(transformedData[, .(paste(chr, pos, ".", ref, alt, ".", ".", "."), collapse="\t")])
  
  write.table(transformedData, file = output_file, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
}

#' Process VEP results and merge with whole genome and sample information
#'
#' This function reads the VEP result file, merges it with the variant data,
#' and then further combines it with sample information.
#'
#' @param processed_whole A data frame containing columns chr, pos, ref, and alt with chromosome,
#'        position, reference base, and alternative base information.
#' @param gene_info_file Path to the VEP output file.
#' @param sample_info A data frame containing sample information.
#' @param output_file Path to save the processed combined data as an RDS file.
#' @importFrom dplyr left_join rename select
#' @importFrom data.table fread
#' @return A summary data frame containing the genotype information for its input
#'         processed_whole with ensembl id, hgnc symbol, size, and frequency in cell
#'         subclass.
#' @export
#' @examples
#' process_vep_results(processed_whole, "/path/to/ensembl_vep_output.txt", sample_info, "/path/to/output/processed_data.rds")
process_vep_results <- function(processed_whole, gene_info_file, sample_info, output_file) {
  gene_info <- fread(gene_info_file)
  processed_gene_whole <- mergeProcessedWithGeneInfo(processed_whole, gene_info)
  
  processed_gene_whole_genome_disease <- processed_gene_whole %>%
    left_join(sample_info, by = c("donor_id" = "external_donor_name")) %>%
    rename(disease = Disease) %>%
    select(chr, pos, ref, alt, Subclass, donor_id, freq, Gene, hgnc_symbol, size, disease)
  
  return(processed_gene_whole_genome_disease)
  saveRDS(processed_gene_whole_genome_disease, file = output_file)
}
