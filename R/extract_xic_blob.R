#' Extract XIC data from an SQLite blob
#'
#' Certain mass spec data alignment formats store feature extracted ion
#' chromatograms as SQLite blobs containing a compressed binary file.
#' This function extracts the stored data from such blobs.
#'
#' The precise format of the binary file has not yet been completely elucidated.
#' Chunk 1 is unknown, but directly related to retention time. Chunk 2 is the
#' retention time in minutes. Chunk 3 is the XIC itself. Chunk 4 looks like
#' probably the instantaneous noise level, but might be something else.
#'
#'
#' @param blb A raw vector containing an SQLite blob of XIC data
#'
#' @return A tibble containing the data
#'
#' @export
extract_xic_blob <- function(blb) {
  # XIC blobs are gzipped binary files
  # still need to figure out what exactly it means and whether lengths can be
  # different. Probably header bytes encode that info somehow.

  # decompress the blob
  df <- memDecompress(from = blb, type = "gzip")

  # First blob has header:
  # BC 01 95 39 4C EA 01 48 B5 27 74 4B BD 23 4C E3 02 [18] 04 00 00 01 17 53
  # last 4 bytes might be actually first u32 of chunk 1?
  # Second blob has header:
  # BC 01 95 39 4C EA 01 48 B5 27 74 4B BD 23 4C E3 02 [52] 04 00 00 01 17 53
  # note single byte is different
  # 0x1804 = 1048, the chunk length
  # 0x5204 = 1106, presumably chunk length for second blob

  # obvious data chunks have fixed length in at least one example
  seg_len <- readBin(df[18:19], "integer", n = 1, size = 2)
  # 4 byte (32 bit) integers, probably unsigned
  dsize <- 4

  # FIXME: need to figure out how chunk offsets are calculated!
  # Not clear whether the gaps between chunks are constant or forcing positions
  # to some multiple of some value or something else
  # Limited testing suggests constant gaps for some weird reason

  # chunk1 is some odd measure of time
  start_chunk1 <- 21
  end_chunk1 <- start_chunk1 + seg_len*dsize - 1
  chunk1 <- readBin(df[start_chunk1:end_chunk1], "integer",
                    n = seg_len, size = dsize)

  # chunk2 is definitely time in minutes
  start_chunk2 <- end_chunk1 + 4
  end_chunk2 <- start_chunk2 + seg_len*dsize - 1
  chunk2 <- readBin(df[start_chunk2:end_chunk2], "double",
                    n = seg_len, size = dsize)
  # chunk3 is the actual XIC
  start_chunk3 <- end_chunk2 + 2
  end_chunk3 <- start_chunk3 + seg_len*dsize - 1
  chunk3 <- readBin(df[start_chunk3:end_chunk3], "double",
                    n = seg_len, size = dsize)

  # might be the noise level? always less than chunk3, slightly correlated
  start_chunk4 <- end_chunk3 + 2
  end_chunk4 <- start_chunk4 + seg_len*dsize - 1
  chunk4 <- readBin(df[start_chunk4:end_chunk4], "double",
                    n = seg_len, size = dsize)

  xic_data <- dplyr::tibble(Chunk1 = chunk1, Time = chunk2,
                             XICData = chunk3, NoiseThresh = chunk4)

  return(xic_data)
}
