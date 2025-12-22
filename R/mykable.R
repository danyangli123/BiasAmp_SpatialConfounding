mykable <- function(x, format = "html", caption = NULL, ...) {
    kable(x, format = format, caption = caption) %>%
        kable_styling(full_width = F, ...)
}

mykable2 <- function(x, ...) {
    mykable(x, bootstrap_options = "condensed", font_size = 12, ...)
}

#https://rdrr.io/github/kaz-yos/tableone/src/R/kableone.R
#https://stackoverflow.com/questions/31430140/print-html-or-word-table-in-knitr-so-that-whitespaces-in-strings-are-respected
mykableone <- function(x, missing = TRUE, nonnormal = NULL, smd = FALSE, return_obj = FALSE, format = "fp", catDigits = 1, contDigits = 2) {
    capture.output(x <- print(x, missing = missing, nonnormal = nonnormal, smd = smd, format = format, catDigits = catDigits, contDigits = contDigits, printToggle = FALSE))
    nms <- rownames(x)
    x %<>% as_tibble() %>%
        mutate(Variable = gsub(" ", "&nbsp;", nms, fixed = TRUE)) %>%
        #mutate(Variable = nms) %>%
        select(Variable, everything())
    #knitr::kable(x, ...)
    #x2 <-  mutate_all(x, function(var) gsub(" ", "&nbsp;", var, fixed = TRUE))
    if (return_obj) return(x)
    knitr::kable(x, format = "html", escape = FALSE, align = c("l", rep("r", ncol(x) - 1))) %>%
        kableExtra::kable_styling(full_width = F, bootstrap_options = c("condensed"), font_size = 12)
}
# #https://rdrr.io/github/kaz-yos/tableone/src/R/kableone.R
# #https://stackoverflow.com/questions/31430140/print-html-or-word-table-in-knitr-so-that-whitespaces-in-strings-are-respected
# mykableone <- function(x, missing = TRUE, return_obj = FALSE) {
#     capture.output(x <- print(x, missing = missing, printToggle = FALSE))
#     nms <- rownames(x)
#     x %<>% as_tibble() %>%
#         mutate(Variable = gsub(" ", "&nbsp;", nms, fixed = TRUE)) %>%
#         #mutate(Variable = nms) %>%
#         select(Variable, everything())
#     #knitr::kable(x, ...)
#     #x2 <-  mutate_all(x, function(var) gsub(" ", "&nbsp;", var, fixed = TRUE))
#     if (return_obj) return(x)
#     kable(x, format = "html", escape = FALSE, align = c("l", rep("r", ncol(x) - 1))) %>%
#         kable_styling(full_width = F, bootstrap_options = c("condensed"), font_size = 12)
# }
