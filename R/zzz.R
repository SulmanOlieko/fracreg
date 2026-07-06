fracregStartupMessage <- function() {
  art_lines <- c(
    "     ,...                                                ",
    "  .d' \"\"                                                ",
    "  dM`                                                   ",
    " mMMmm`7Mb,od8 ,6\"Yb.  ,p6\"bo `7Mb,od8 .gP\"Ya   .P\"Ybmmm",
    "  MM    MM' \"'8)   MM 6M'  OO   MM' \"',M'   Yb :MI  I8  ",
    "  MM    MM     ,pm9MM 8M        MM    8M\"\"\"\"\"\"  WmmmP\"  ",
    "  MM    MM    8M   MM YM.    ,  MM    YM.    , 8M       ",
    ".JMML..JMML.  `Moo9^Yo.YMbmd' .JMML.   `Mbmmd'  YMMMMMb ",
    "                                               6'     dP",
    "                                               Ybmmmd'  "
  )
  base_msg <- paste0(
    "\n\n* Please cite the 'fracreg' package as:\n",
    "Owili, SO. (2026). fracreg: Fractional Response Regressions. R package version 1.0.0.\n\n",
    "See also: citation(\"fracreg\")\n\n",
    "* For any questions, suggestions, or comments on the 'fracreg' package, you can contact the authors directly or visit:\n",
    "  https://github.com/SulmanOlieko/fracreg/issues\n"
  )

  console_width <- getOption("width")
  art_width <- max(nchar(art_lines))

  if (console_width < art_width) {
    # 1. Truncate all lines to fit the console width
    trimmed_art <- substr(art_lines, 1, console_width - 3)

    # 2. Add the ASCII ellipsis strictly to the last row
    ellipsis <- "..."
    e_len <- nchar(ellipsis)
    cut_point <- max(0, console_width - e_len)

    last_line_index <- length(trimmed_art)
    trimmed_art[last_line_index] <- paste0(
      substr(art_lines[last_line_index], 1, cut_point),
      ellipsis
    )

    # 3. Set the version alignment target
    align_width <- console_width
  } else {
    trimmed_art <- art_lines
    align_width <- art_width
  }

  art_msg <- paste(trimmed_art, collapse = "\n")

  # 4. Pad the version text so it hits the exact right edge
  version_msg <- paste0(
    "\n",
    sprintf(paste0("%", align_width, "s"), "version 1.0.0")
  )

  return(paste0(art_msg, version_msg, base_msg))
}

.onAttach <- function(libname, pkgname) {
  # We still skip printing during the CRAN check to keep the test logs clean
  in_chk <- Sys.getenv("_R_CHECK_PACKAGE_NAME_") != ""

  if (!in_chk) {
    msg <- fracregStartupMessage()
  } else {
    msg <- paste0(
      "\n* Please cite the 'fracreg' package as:\n",
      "Owili, SO. (2026). fracreg: Fractional Response Regressions. R package version 1.0.0.\n\n",
      "See also: citation(\"fracreg\")\n\n",
      "* For any questions, suggestions, or comments on the 'fracreg' package, you can contact the authors directly or visit:\n",
      "  https://github.com/SulmanOlieko/fracreg/issues\n"
    )
  }

  packageStartupMessage(msg)
  invisible()
}
