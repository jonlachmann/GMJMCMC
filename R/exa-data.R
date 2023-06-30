#' Excerpt from the Open Exoplanet Catalogue data set
#'
#' Data fields include planet and host star attributes.
#'
#' The variables are as follows:
#'
#' \itemize{
#' \item TypeFlag: Flag indicating the type of data
#' \item PlanetaryMassJpt: Mass of the planetary object in Jupiter masses
#' \item RadiusJpt: Radius of the planetary object in Jupiter radii
#' \item PeriodDays: Orbital period of the planetary object in days
#' \item SemiMajorAxisAU: Semi-major axis of the planetary object's orbit in astronomical units
#' \item Eccentricity: Eccentricity of the planetary object's orbit
#' \item HostStarMassSlrMass: Mass of the host star in solar masses
#' \item HostStarRadiusSlrRad: Radius of the host star in solar radii
#' \item HostStarMetallicity: Metallicity of the host star
#' \item HostStarTempK: Effective temperature of the host star in Kelvin
#' \item PlanetaryDensJpt: Density of the planetary object up to a constant
#' }
#'
#' @docType data
#' @keywords datasets
#' @name exoplanet
#' @usage data(exoplanet)
#' @format A data frame with 223 rows and 11 variables
#' @source Dataset downloaded from the Open Exoplanet Catalogue Repository.
#' \url{https://github.com/OpenExoplanetCatalogue/oec_tables/}
#'
#' Creators:
#'
#' 1. Prof. Hanno Rein, Department for Physical and Environmental Sciences.
#' University of Toronto at Scarborough
#' Toronto, Ontario M1C 1A4
#' hanno.rein 'at' utoronto.ca
#'
"exoplanet"