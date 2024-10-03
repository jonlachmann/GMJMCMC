#' Gene expression data lymphoblastoid cell lines of all 210 unrelated HapMap 
#' individuals from four populations
#'
#' A 210 times 47293 matrix with indviduals along the rows and expression data along the columns
#'
#' The 
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
#' @name SangerData
#' @usage data(SangerData)
#' @format A numerical matrix  with 210 rows and 47293 variables
#' @source Dataset downloaded from 
#' \url{https://ftp.sanger.ac.uk/pub/genevar/}
#'
#' Reference:
#'
#' Stranger, BE et al (2007): Relative impact of nucleotide and copy number variation on gene expression phenotypes
#' Science, 2007•science.org
#'
"SangerData"