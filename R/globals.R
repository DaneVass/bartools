# fixes issues with no visible binding for global variables. See here
# https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887
utils::globalVariables(c("."))
