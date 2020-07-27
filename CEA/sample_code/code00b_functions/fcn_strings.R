#'Various functions to work with strings

#' Ask for the last $n$ characters of a string
substrRight = function(x, n){substr(x, nchar(x)-n+1, nchar(x))}
#' Trim white spaces in the left or right of strings, but not in between.
trim = function (x) gsub("^\\s+|\\s+$", "", x)
trim_between = function (x) gsub("\\s+", " ", x)
#' To split a text into two 
wrapper = function(x, w) {
  paste(strwrap(x, width=w), collapse = "\n")
}

#' example:
wrapper("To be or not to be, that is the question", 22)

wrapper_latex = function(x, w, ch) {
  paste(strwrap(x, width=w), collapse = ch)
}

delete_prefix = function(x, ch){return(strsplit(x, ch)[[1]][2])}
#' example:
delete_prefix("hat.mepp", "[.]")
sapply(c("hat.mepp", "break.fast"), "delete_prefix", "[.]")