# R reader for the central run configuration (framework_config.json).

# Walk up the directory tree until framework_config.json is found, so scripts
# resolve the repo root regardless of their working directory.
get_framework_root <- function(start_path = getwd()) {
  path <- normalizePath(start_path, mustWork = TRUE)
  if (file.exists(path) && !dir.exists(path)) {
    path <- dirname(path)
  }

  repeat {
    if (file.exists(file.path(path, "framework_config.json"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate framework_config.json from: ", start_path, call. = FALSE)
    }
    path <- parent
  }
}

read_framework_config <- function(start_path = getwd()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read framework_config.json.", call. = FALSE)
  }

  root <- get_framework_root(start_path)
  cfg_path <- file.path(root, "framework_config.json")

  cfg <- jsonlite::fromJSON(cfg_path, simplifyVector = FALSE)
  cfg$framework_root <- root  # remember where the config lives for later path resolution
  cfg$framework_config_path <- cfg_path
  cfg
}

cfg_get <- function(cfg, keys, default = NULL, required = TRUE) {
  value <- cfg
  for (key in keys) {
    if (is.null(value[[key]])) {
      if (required) {
        stop("Missing config entry: ", paste(keys, collapse = "."), call. = FALSE)
      }
      return(default)
    }
    value <- value[[key]]
  }
  value
}

path_from_framework_root <- function(cfg, relative_path) {
  if (is.null(relative_path) || !nzchar(relative_path)) {
    return(relative_path)
  }

  is_absolute <- grepl("^/", relative_path) || grepl("^[A-Za-z]:[/\\\\]", relative_path)  # POSIX or Windows drive paths
  if (is_absolute) {
    return(normalizePath(relative_path, mustWork = FALSE))
  }

  normalizePath(file.path(cfg$framework_root, relative_path), mustWork = FALSE)
}

current_script_path <- function(default = getwd()) {
  match <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(match)) {
    return(normalizePath(default, mustWork = FALSE))
  }
  normalizePath(sub("--file=", "", match[1]), mustWork = TRUE)
}
