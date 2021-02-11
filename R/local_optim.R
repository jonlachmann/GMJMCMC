# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-02-11

simulated.annealing <- function (indices) {

}

greedy.optim <- function () {

}

local.optim <- function (model, features, type) {
  if (type == 1) {
    return(simulated.annealing(indices))
  }
  if (type == 2) {
    return(greedy.optim(indices))
  }
  if (type == 3) {
    return("not implemented")
  }
  stop("Invalid local optimizer chosen")
}