library(nexonsR)

local({
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  for (type in c('gene', 'partial', 'unique')) {
    header <- 'Gene_ID\tGene_Name\tChr\tStart\tEnd\tStrand\ta.bam\tb.bam'
    writeLines(if (type == 'gene') header else paste0('Transcript_ID\t', header),
               file.path(folder, paste0(type, '.txt')))
  }
  write_flex <- function(sample, rows) {
    con <- gzfile(file.path(folder, paste0(sample, '_flexout.txt.gz')), 'wt')
    on.exit(close(con))
    writeLines(c('Transcript_ID\tStart_Flex\tEnd_Flex', rows), con)
  }
  error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))
  write_flex('a', c('t1\t-20000\t3', 't2\t0\t-1', 't1\t-20000\t4',
                    't1\t5\t3', 't1\t-20000\t3'))
  write_flex('b', character())
  x <- read_nexons(folder, '')
  messages <- character()
  withCallingHandlers(x$parse_flexout(chunk_size = 2L), message = function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart('muffleMessage')
  })
  stopifnot(length(messages) == 2L, identical(x$flexout_loaded, c('a.bam', 'b.bam')),
    identical(x$get_flexout('t1')[[1]], data.frame(value = c(-20000L, 5L), count = c(3, 1))),
    identical(x$get_flexout('t1', 'End_Flex')[[1]], data.frame(value = c(3L, 4L), count = c(3, 1))),
    nrow(x$get_flexout('t1')[[2]]) == 0L,
    nrow(x$get_flexout('absent')[[1]]) == 0L, length(x$loaded) == 0L)
  error(x$get_flexout('t1', samples = 'unknown'))
  error(x$parse_flexout(chunk_size = 0))
  clone <- x$clone(deep = TRUE)
  write_flex('a', 't1\t1.5\t2')
  x$parse_flexout() # Cached samples are not reopened.
  x$clear_cache()
  stopifnot(length(x$flexout_loaded) == 0L, length(clone$flexout_loaded) == 2L)
  for (row in c('t1\t1.5\t2', 't1\t1\t', 't1\t1\t2\t3', '\t1\t2',
                't1\t2147483648\t2', 't1\tNA\t2')) {
    write_flex('a', row)
    suppressMessages(error(x$parse_flexout('a.bam', chunk_size = 1L)))
    stopifnot(length(x$flexout_loaded) == 0L)
  }
  # Chunk boundaries must not change aggregates.
  write_flex('a', c('t1\t-20000\t3', 't2\t0\t-1', 't1\t-20000\t4',
                    't1\t5\t3', 't1\t-20000\t3'))
  suppressMessages(x$parse_flexout(chunk_size = 100L))
  stopifnot(identical(x$get_flexout('t1'), clone$get_flexout('t1')))
})
