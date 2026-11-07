.libPaths(c("library", .libPaths()))
library(data.table)
options(timeout = 600)

for (folder in c("data/raw", "data/processed", "results/tables", "results/models", "results/figures")) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

# The series matrix contains expression values and sample information together.
sources <- read.csv("data/source_files.csv")
for (i in seq_len(nrow(sources))) {
  path <- file.path("data/raw", sources$file[i])
  if (!file.exists(path)) {
    download.file(sources$url[i], path, mode = "wb", method = "libcurl")
  }
  if (unname(tools::md5sum(path)) != sources$md5[i]) {
    stop("The GEO file differs from the version used here: ", path,
         ". Check it before updating data/source_files.csv.")
  }
}

connection <- gzfile("data/raw/GSE56035_series_matrix.txt.gz", "rt")
header <- character()
repeat {
  line <- readLines(connection, n = 1)
  if (!length(line)) stop("Could not find the expression table.")
  if (line == "!series_matrix_table_begin") break
  header <- c(header, line)
}
close(connection)

read_field <- function(name) {
  line <- header[startsWith(header, paste0(name, "\t"))]
  stopifnot(length(line) == 1)
  scan(text = line, what = "", sep = "\t", quote = '"', quiet = TRUE)[-1]
}
characteristics <- header[startsWith(header, "!Sample_characteristics_ch1\t")]
characteristics <- lapply(characteristics, function(line) {
  scan(text = line, what = "", sep = "\t", quote = '"', quiet = TRUE)[-1]
})
read_characteristic <- function(pattern) {
  found <- vapply(characteristics, function(x) any(grepl(pattern, x, ignore.case = TRUE)), logical(1))
  stopifnot(sum(found) == 1)
  sub("^[^:]+: *", "", characteristics[[which(found)]])
}

metadata <- data.frame(
  sample_id = read_field("!Sample_geo_accession"),
  title = read_field("!Sample_title"),
  Age = as.numeric(read_characteristic("^age ")),
  Sex = tools::toTitleCase(tolower(read_characteristic("^(sex|gender):"))),
  CellType = ifelse(grepl("monocyte", read_characteristic("^cell type:"), ignore.case = TRUE),
                    "Monocyte", "CD4_T"),
  Batch = read_characteristic("^batch:")
)
metadata$Individual <- sub("\\..*$", "", metadata$title)
stopifnot(!anyNA(metadata), !anyDuplicated(metadata$sample_id),
          all(grepl("^IGTB[0-9]+$", metadata$Individual)))

# Donor prefixes link cell types. Keep repeat suffixes in the original title.
donors <- as.data.table(metadata)[, .(
  n_arrays = .N, n_cell_types = uniqueN(CellType),
  n_ages = uniqueN(Age), n_sexes = uniqueN(Sex),
  recorded_ages = paste(sort(unique(Age)), collapse = ";")
), by = Individual]
donors$paired <- donors$n_cell_types == 2
donors$included <- donors$paired & donors$n_ages == 1 & donors$n_sexes == 1
metadata$paired <- donors$paired[match(metadata$Individual, donors$Individual)]
metadata$included <- donors$included[match(metadata$Individual, donors$Individual)]
for (name in c("Sex", "CellType", "Batch", "Individual")) {
  metadata[[name]] <- factor(metadata[[name]])
}
rownames(metadata) <- metadata$sample_id
fwrite(donors, "results/tables/donor_qc.csv")
fwrite(metadata, "results/tables/sample_metadata.csv")

# Unzip to a temporary file so fread does not need an extra gzip package.
connection <- gzfile("data/raw/GSE56035_series_matrix.txt.gz", "rb")
temporary <- tempfile(fileext = ".txt")
output <- file(temporary, "wb")
repeat {
  block <- readBin(connection, "raw", n = 1024 * 1024)
  if (!length(block)) break
  writeBin(block, output)
}
close(connection)
close(output)
table <- fread(temporary, skip = length(header) + 1,
               nrows = as.integer(read_field("!Sample_data_row_count")[1]))
unlink(temporary)
probe_ids <- as.character(table[[1]])
expression <- as.matrix(table[, -1])
storage.mode(expression) <- "double"
rownames(expression) <- probe_ids
stopifnot(identical(colnames(expression), metadata$sample_id),
          !anyDuplicated(probe_ids), all(is.finite(expression)))

# Check batch structure before changing the supplied scale or selecting samples.
batch_r2 <- function(values, batch) {
  y <- t(values)
  fit <- lm.fit(model.matrix(~ batch), y)
  total <- colSums(sweep(y, 2, colMeans(y))^2)
  1 - colSums(fit$residuals^2) / total
}
set.seed(20261104)
batch_checks <- list()
for (cell in unique(as.character(metadata$CellType))) {
  samples <- metadata$CellType == cell
  batch <- droplevels(metadata$Batch[samples])
  batch_checks[[cell]] <- data.frame(
    probe_id = probe_ids, cell_type = cell,
    observed_r2 = batch_r2(expression[, samples], batch),
    shuffled_r2 = batch_r2(expression[, samples], sample(batch))
  )
}
fwrite(do.call(rbind, batch_checks), "results/tables/batch_structure.csv")

# GEO describes these as linear-scale RMA values. Some are nevertheless negative.
# Exclude the affected probes instead of clipping values or adding a pseudocount.
stopifnot(all(grepl("Output on 2\\^ scale", read_field("!Sample_data_processing"))))
probe_qc <- data.frame(probe_id = probe_ids, n_nonpositive = rowSums(expression <= 0))
keep_probe <- probe_qc$n_nonpositive == 0
probe_qc$included <- keep_probe
expression <- log2(expression[keep_probe, ])
stopifnot(all(apply(expression[, metadata$included], 1, var) > 1e-8))
fwrite(probe_qc, "results/tables/probe_qc.csv")

# Keep probe IDs as row names; annotation is for interpretation, not aggregation.
connection <- gzfile("data/raw/GPL6244.annot.gz", "rt")
lines <- readLines(connection)
close(connection)
start <- which(lines == "!platform_table_begin") + 1
end <- which(lines == "!platform_table_end") - 1
platform <- fread(text = paste(lines[start:end], collapse = "\n"), sep = "\t", quote = "", fill = TRUE)
platform <- platform[match(probe_ids, as.character(platform$ID)), ]
annotation <- data.frame(probe_id = probe_ids, symbol = platform[["Gene symbol"]],
                         gene_title = platform[["Gene title"]], gene_id = platform[["Gene ID"]])
annotation$symbol[is.na(annotation$symbol) | annotation$symbol == ""] <- NA_character_
annotation$unique_symbol <- !is.na(annotation$symbol) & !grepl("///", annotation$symbol, fixed = TRUE)
fwrite(annotation, "results/tables/probe_annotation.csv")

saveRDS(list(expression = expression, metadata = metadata, annotation = annotation),
        "data/processed/immvar.rds", compress = FALSE)
counts <- data.frame(
  item = c("Input arrays", "Input probes", "Donor prefixes", "Paired donors", "Primary donors", "Primary arrays", "Retained probes"),
  n = c(nrow(metadata), length(probe_ids), nrow(donors), sum(donors$paired),
        sum(donors$included), sum(metadata$included), nrow(expression))
)
fwrite(counts, "results/tables/sample_counts.csv")
print(counts, row.names = FALSE)
