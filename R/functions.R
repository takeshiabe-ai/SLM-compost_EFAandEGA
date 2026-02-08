numify <- function(x){
  if(is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

numify_df <- function(df){
  as.data.frame(lapply(df, numify))
}

clean_syntax <- function(s){
  s <- as.character(s)
  s <- gsub("_x000D_", "\n", s, fixed=TRUE)
  s <- gsub("\r\n", "\n", s)
  s <- gsub("\r", "\n", s)
  s <- gsub("[ \t]+", " ", s)
  trimws(s)
}

extract_items_from_syntax <- function(model_syntax){
  lines <- unlist(strsplit(model_syntax, "\n", fixed=TRUE))
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#", lines)]
  meas  <- lines[grepl("=~", lines, fixed=TRUE)]
  if(length(meas) == 0) return(character(0))

  rhs <- trimws(gsub("^.*=~", "", meas))
  parts <- trimws(unlist(strsplit(rhs, "\\+")))
  parts <- parts[parts != ""]
  items <- trimws(sub(".*\\*", "", parts))
  unique(items)
}

collapse_sparse_ordinal <- function(x, min_n=5){
  x <- suppressWarnings(as.numeric(x))
  if(all(is.na(x))) return(x)
  tab <- table(x, useNA="no")
  if(length(tab) < 2) return(x)

  repeat{
    tab <- table(x, useNA="no")
    if(length(tab) < 2) break
    low <- which(tab < min_n)
    if(length(low)==0) break

    k <- as.numeric(names(tab))[low[1]]
    levs <- sort(as.numeric(names(tab)))
    pos <- match(k, levs)

    if(pos == 1) target <- levs[2]
    else if(pos == length(levs)) target <- levs[length(levs)-1]
    else {
      left <- levs[pos-1]; right <- levs[pos+1]
      target <- ifelse(abs(k-left) <= abs(k-right), left, right)
    }
    x[x == k] <- target
  }
  x
}

make_ordered_observed_levels <- function(x){
  x <- suppressWarnings(as.numeric(x))
  lev <- sort(unique(x[!is.na(x)]))
  ordered(x, levels=lev)
}
