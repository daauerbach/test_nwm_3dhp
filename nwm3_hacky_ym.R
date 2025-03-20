library(tidyverse)
library(sf)
dir_data_common <- "~/T/DFW-Team WDFW Watershed Synthesis - data_common"

#for COMIDs all of WA
sf_nhdp_wa <- readRDS("~/T/DFW-Team WDFW Watershed Synthesis - flow_trees_heat/sf_nhdp_wa.rds")
sf_nhdp_wa_comid <- sf_nhdp_wa |> 
  as_tibble() |> 
  dplyr::filter(COMID > 0, areasqkm < 1500) |> 
  select(comid = COMID)
rm(sf_nhdp_wa)

library(Rarr)


my_format_chunk <- function (decompressed_chunk, metadata, alt_chunk_dim) 
{
  datatype <- metadata$datatype
  actual_chunk_size <- length(decompressed_chunk)/datatype$nbytes
  if ((datatype$base_type == "py_object") || (actual_chunk_size == 
                                              prod(unlist(metadata$chunks)))) {
    chunk_dim <- unlist(metadata$chunks)
  }
  else if (actual_chunk_size == prod(alt_chunk_dim)) {
    chunk_dim <- alt_chunk_dim
  }
  else {
    stop("Decompresed data doesn't match expected chunk size.")
  }
  if (metadata$order == "C") {
    chunk_dim <- rev(chunk_dim)
  }
  if (datatype$base_type == "string") {
    converted_chunk <- .format_string(decompressed_chunk, 
                                      datatype)
    dim(converted_chunk[[1]]) <- chunk_dim
  }
  else if (datatype$base_type == "unicode") {
    converted_chunk <- .format_unicode(decompressed_chunk, 
                                       datatype)
    dim(converted_chunk[[1]]) <- chunk_dim
  }
  else if (datatype$base_type == "py_object") {
    converted_chunk <- .format_object(decompressed_chunk, 
                                      metadata, datatype)
    dim(converted_chunk[[1]]) <- chunk_dim
  }
  else {
    output_type <- switch(datatype$base_type, boolean = 0L, 
                          #            int = 1L, uint = 1L, float = 2L)
                          int = 1L, uint = 1L, float = 1L) #force float to 1L, not a general fix
    converted_chunk <- .Call("type_convert_chunk", decompressed_chunk, 
                             output_type, datatype$nbytes, datatype$is_signed, 
                             chunk_dim, PACKAGE = "Rarr")
  }
  if (metadata$order == "C") {
    converted_chunk[[1]] <- aperm(converted_chunk[[1]])
  }
  names(converted_chunk) <- c("chunk_data", "warning")
  return(converted_chunk)
}

#changes: add Rarr::: throughout; update chunk_file declaration to include '/' 
my_read_chunk <- function (zarr_array_path, chunk_id, metadata, s3_client = NULL, 
                           alt_chunk_dim = NULL) 
{
  if (missing(metadata)) {
    metadata <- Rarr:::read_array_metadata(zarr_array_path, s3_client = s3_client)
  }
  dim_separator <- ifelse(is.null(metadata$dimension_separator), 
                          yes = ".", no = metadata$dimension_separator)
  chunk_id <- paste(chunk_id, collapse = dim_separator)
  datatype <- metadata$datatype
  #    chunk_file <- paste0(zarr_array_path, chunk_id)
  chunk_file <- paste0(zarr_array_path,"/", chunk_id)
  
  if (nzchar(Sys.getenv("RARR_DEBUG"))) {
    message(chunk_file)
  }
  if (is.null(s3_client)) {
    size <- file.size(chunk_file)
    if (file.exists(chunk_file)) {
      compressed_chunk <- readBin(con = chunk_file, what = "raw", 
                                  n = size)
    }
    else {
      compressed_chunk <- NULL
    }
  }
  else {
    parsed_url <- Rarr:::parse_s3_path(chunk_file)
    if (Rarr:::.s3_object_exists(s3_client, parsed_url$bucket, parsed_url$object)) {
      compressed_chunk <- s3_client$get_object(Bucket = parsed_url$bucket, 
                                               Key = parsed_url$object)$Body
    }
    else {
      compressed_chunk <- NULL
    }
  }
  if (!is.null(compressed_chunk)) {
    decompressed_chunk <- Rarr:::.decompress_chunk(compressed_chunk, 
                                                   metadata)
    #        converted_chunk <- Rarr:::.format_chunk(decompressed_chunk, 
    converted_chunk <- my_format_chunk(decompressed_chunk, 
                                       metadata, alt_chunk_dim)
  }
  else {
    converted_chunk <- list(chunk_data = array(metadata$fill_value, 
                                               dim = unlist(metadata$chunks)), warning = 0L)
  }
  return(converted_chunk)
}

#need to declare here as not even ':::' available
#also changing to:   chunk <- my_read_chunk(
my_extract_elements <- function(i, metadata, index, required_chunks, zarr_array_path, s3_client, chunk_idx) {
  ## find elements to select from the chunk and what in the output we replace
  index_in_result <- index_in_chunk <- list()
  alt_chunk_dim <- unlist(metadata$chunks)
  
  for (j in seq_len(ncol(required_chunks))) {
    index_in_result[[j]] <- which(chunk_idx[[j]] == required_chunks[i, j])
    ## are we requesting values outside the array due to overhanging chunks?
    outside_extent <- index_in_result[[j]] > metadata$shape[[j]]
    if (any(outside_extent))
      index_in_result[[j]] <- index_in_result[[j]][-outside_extent]
    if (any(index_in_result[[j]] == metadata$shape[[j]])) 
      alt_chunk_dim[j] <- length(index_in_result[[j]])
    
    index_in_chunk[[j]] <- ((index[[j]][index_in_result[[j]]] - 1) %% metadata$chunks[[j]]) + 1
  }
  
  ## read this chunk
  #chunk <- Rarr:::read_chunk(
  chunk <- my_read_chunk(
    zarr_array_path,
    chunk_id = required_chunks[i, ],
    metadata = metadata,
    s3_client = s3_client,
    alt_chunk_dim = alt_chunk_dim
  )
  
  warn <- chunk$warning[1]
  chunk_data <- chunk$chunk_data
  
  ## extract the required elements from the chunk
  selection <- R.utils::extract(chunk_data, indices = index_in_chunk, drop = FALSE)
  
  return(list(selection, index_in_result, warning = warn))
}

#change to FUN = .extract_elements
my_read_data <- function(required_chunks, zarr_array_path, s3_client, 
                         index, metadata) {
  warn <- 0L
  
  ## determine which chunk each of the requests indices belongs to
  chunk_idx <- mapply(\(x,y) { (x-1) %/% y }, index, metadata$chunks, SIMPLIFY = FALSE)
  
  ## hopefully we can eventually do this in parallel
  chunk_selections <- lapply(seq_len(nrow(required_chunks)), 
                             #FUN = .extract_elements,
                             FUN = my_extract_elements,
                             metadata = metadata, index = index,
                             required_chunks = required_chunks,
                             zarr_array_path = zarr_array_path,
                             s3_client = s3_client,
                             chunk_idx = chunk_idx)
  
  ## predefine our array to be populated from the read chunks
  output <- array(metadata$fill_value, dim = vapply(index, length, integer(1)))
  
  ## proceed in serial and update the output with each chunk selection in turn
  for (i in seq_along(chunk_selections)) {
    index_in_result <- chunk_selections[[i]][[2]]
    cmd <- Rarr:::.create_replace_call(x_name = "output", idx_name = "index_in_result",
                                       idx_length = length(index_in_result), 
                                       y_name = "chunk_selections[[i]][[1]]")
    eval(parse(text = cmd))
    warn <- max(warn, chunk_selections[[i]]$warning[1])
  }
  return(list(output = output, warn = warn))
}

my_read_zarr_array <- function (zarr_array_path, index, s3_client) 
{
  zarr_array_path <- Rarr:::.normalize_array_path(zarr_array_path)
  if (missing(s3_client)) 
    s3_client <- Rarr:::.create_s3_client(path = zarr_array_path)
  metadata <- Rarr:::read_array_metadata(zarr_array_path, s3_client = s3_client)
  #yikes
  metadata$datatype$base_type <- "float"
  
  if (missing(index)) {
    index <- vector(mode = "list", length = length(metadata$shape))
  }
  index <- Rarr:::check_index(index = index, metadata = metadata)
  required_chunks <- as.matrix(Rarr:::find_chunks_needed(metadata, index))
  res <- my_read_data(required_chunks, zarr_array_path, s3_client, index, metadata)
  if (isTRUE(res$warn > 0)) {
    warning("Integer overflow detected in at least one chunk.\n", 
            "Overflowing values have been replaced with NA", 
            call. = FALSE)
  }
  return(res$output)
}

#Now runs...to 56K tibble
feat_id <- tibble(
  # #v2.1: comid = read_zarr_array("https://noaa-nwm-retrospective-2-1-zarr-pds.s3.amazonaws.com/chrtout.zarr/feature_id/")
  comid = my_read_zarr_array("https://noaa-nwm-retrospective-3-0-pds.s3.amazonaws.com/CONUS/zarr/chrtout.zarr/feature_id")
) |> 
  tibble::rowid_to_column(var = "i") |> 
  inner_join(sf_nhdp_wa_comid, by = "comid")


pull_nwm_zarr <- function(
    #s3_address = "https://noaa-nwm-retrospective-2-1-zarr-pds.s3.amazonaws.com/chrtout.zarr/streamflow/",
  s3_address =   "https://noaa-nwm-retrospective-3-0-pds.s3.amazonaws.com/CONUS/zarr/chrtout.zarr/streamflow/",
  y = 2022,
  m = 8,
  h = 13,
  feat = feat_id,
  dir_write = file.path(dir_data_common, "nwm_retro/v3_0")
) {
  sclfctr <- 0.009999999776482582
  cms2cfs <- 35.31468492
  
  #vector of all time available, v2.1 length was 367439
  time_dim <- seq(as.POSIXct("1979-02-01 01:00:00", tz = 'UTC'), length.out = 385704, by = "1 hour")
  #indices of desired subset
  time_i <- which(
    #year(time_dim) %in% 2023 & month(time_dim) %in% 8 & hour(time_dim) %in% 13
    year(time_dim) %in% y &
      month(time_dim) %in% m &
      hour(time_dim) %in% h
  )
  
  print(Sys.time())
  z <- my_read_zarr_array(s3_address, index = list(time_i, feat$i)) |> 
    as.data.frame() |> 
    set_names(feat$comid) |> 
    as_tibble() |> 
    mutate(
      #correct 'scale factor' brings to cms, then to cfs
      across(everything(), ~.* sclfctr * cms2cfs),
      across(everything(), ~if_else(.<0, NA_real_, .)),
      i = time_dim[time_i], #col of POSIXct
      m = month(i),
      yday = yday(i)
    )
  print(Sys.time())
  
  if(!is.null(dir_write)){
    if(identical(m, 1:12)) {
      f = file.path(dir_write, paste0("nwm_wa_",y, ".rds"))
    } else {
      f = file.path(dir_write, paste0("nwm_wa_",y,"_",paste0(str_pad(m,width = 2, pad = "0"),collapse = ""), ".rds"))
    }
    print(paste("data written to ", f))
    saveRDS(z, f)
  } else {
    return(z)
  }
}


args <- commandArgs(trailingOnly = TRUE)

#x <- expand_grid(y = args[1], m = args[2])

pull_nwm_zarr(y = args[1], m = args[2])

