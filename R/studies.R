setClass("studySize")

setClass("oneWayStudySize",
  representation(size = "numeric"),
  contains = "studySize",
  validity = function(object) {
    if (!is.natural(object@size)) {
      "size must be a natural number"
    } else if (length(object@size) != 1) {
      "size must be a scalar"
    } else {
      TRUE
    }
  }
)

setClass("twoWayStudySize",
  representation(
    grp1 = "numeric",
    grp2 = "numeric"
  ),
  contains = "studySize",
  validity = function(object) {
    grp1.size <- object@grp1
    grp2.size <- object@grp2
    if (!is.natural(grp1.size)) {
      "grp1 must be a natural number"
    } else if (!is.natural(grp2.size)) {
      "grp2 must be a natural number"
    } else if (length(grp1.size) != 1) {
      "grp1 must be a scalar"
    } else if (length(grp2.size) != 1) {
      "grp2 must be a scalar"
    } else {
      TRUE
    }
  }
)

##' @export
study.size <- function(total.size = NULL, grp1.size = NULL, grp2.size = NULL) {
  # get the study sizing
  # options are
  # total.size given, even positive integer, other not given
  # grp1.size given positive integer, others not given
  # grp1.size and grp2.size distinct positive integers, others not given
  if (!is.null(total.size) && (is.null(grp1.size) && is.null(grp2.size))) {
    sizing <- total.size / 2
    new("twoWayStudySize", grp1 = sizing, grp2 = sizing)
  } else if ((is.null(total.size) && is.null(grp2.size)) &&
    (!is.null(grp1.size))) {
    new("twoWayStudySize", grp1 = grp1.size, grp2 = grp1.size)
  } else if (is.null(total.size) && !(is.null(grp1.size) || is.null(grp2.size))) {
    new("twoWayStudySize", grp1 = grp1.size, grp2 = grp2.size)
  } else {
    stop("bad parameters")
  }
}

setMethod(
  "show",
  signature("oneWayStudySize"),
  function(object) {
    cat("size : ", object@size, "\n", sep = "")
  }
)

setMethod(
  "show",
  signature("twoWayStudySize"),
  function(object) {
    cat("group 1 size : ", object@grp1, "\n",
      "group 2 size : ", object@grp2, "\n",
      sep = ""
    )
  }
)
