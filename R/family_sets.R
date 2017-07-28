family_set_archimedean <- c(
  "clayton", "gumbel", "frank", "joe"
)

family_set_elliptical <- c(
  "gaussian", "student"
)

family_set_bb <- c(
  "bb1", "bb6", "bb7", "bb8"
)

family_set_onepar <- c(
  "gaussian", family_set_archimedean
)

family_set_twopar <- c(
  "t", family_set_bb
)

family_set_parametric <- c(
  "indep", family_set_onepar, family_set_twopar
)

family_set_all <- family_set_parametric
family_set_rotations <- setdiff(c(family_set_archimedean,
                                  family_set_bb), "frank")

family_set_defs <- c(
  "archimedean", "elliptical", "bbs", "oneparametric", "twoparametric", 
  "parametric", "all"
)

family_set_all_defs <- c(
  family_set_all, family_set_defs
)

check_family_set <- function(family_set) {
  i_wrong <- which(!(family_set %in% family_set_all_defs))
  if (length(i_wrong) > 0) {
    stop(
      "unknown families in family_set: ",
      paste0('"', family_set[i_wrong], '"', collapse = ", ")
    )
  }
}

get_bounds <- function(family) {
  switch(family,
         indep = c(0, 0, 0),
         gaussian = c(0, 0.5, 1),
         t = rbind(c(0, 0.5, 1),
                   c(2, 5, 50)),
         clayton = c(0, 1, 200),
         gumbel = c(1, 2, 50),
         frank = c(0, 1, 50),
         joe = c(1, 2, 50),
         bb1 = rbind(c(0, 1, 200),
                     c(1, 2, 200)),
         bb6 = rbind(c(1, 2, 200),
                     c(1, 2, 200)),
         bb7 = rbind(c(1, 2, 200),
                     c(0, 1, 200)),
         bb8 = rbind(c(1, 2, 200),
                     c(0, 0.5, 1)))
}