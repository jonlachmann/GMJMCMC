#' Breast Cancer Wisconsin (Diagnostic) Data Set
#'
#' Features are computed from a digitized image of a fine needle aspirate (FNA) of a breast mass.
#' They describe characteristics of the cell nuclei present in the image.
#'
#' Separating plane described above was obtained using Multisurface Method-Tree (MSM-T)
#' (K. P. Bennett, "Decision Tree Construction Via Linear Programming." Proceedings of the 4th Midwest Artificial
#' Intelligence and Cognitive Science Society, pp. 97-101, 1992), a classification method which uses linear programming to
#' construct a decision tree. Relevant features were selected using an exhaustive search in the space of 1-4 features
#' and 1-3 separating planes.
#'
#' The actual linear program used to obtain the separating plane in the 3-dimensional space is that described in:
#' (K. P. Bennett and O. L. Mangasarian: "Robust Linear Programming Discrimination of Two Linearly Inseparable Sets",
#' Optimization Methods and Software 1, 1992, 23-34).
#'
#' The variables are as follows:
#'
#' \itemize{
#' \item ID number
#' \item Diagnosis (1 = malignant, 0 = benign)
#' \item Ten real-valued features are computed for each cell nucleus
#' }
#'
#' @docType data
#' @keywords datasets
#' @name breastcancer
#' @usage data(breastcancer)
#' @format A data frame with 569 rows and 32 variables
#' @source Dataset downloaded from the UCI Machine Learning Repository.
#' \url{http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)}
#'
#' Creators:
#'
#' 1. Dr. William H. Wolberg, General Surgery Dept.
#' University of Wisconsin, Clinical Sciences Center
#' Madison, WI 53792
#' wolberg 'at' eagle.surgery.wisc.edu
#'
#' 2. W. Nick Street, Computer Sciences Dept.
#' University of Wisconsin, 1210 West Dayton St., Madison, WI 53706
#' street 'at' cs.wisc.edu 608-262-6619
#'
#' 3. Olvi L. Mangasarian, Computer Sciences Dept.
#' University of Wisconsin, 1210 West Dayton St., Madison, WI 53706
#' olvi 'at' cs.wisc.edu
#'
#' Donor: Nick Street
#' @references W.N. Street, W.H. Wolberg and O.L. Mangasarian. Nuclear feature extraction for breast tumor diagnosis.
#' IS&T/SPIE 1993 International Symposium on Electronic Imaging: Science and Technology, volume 1905, pages 861-870, San Jose, CA, 1993.
#' @references Lichman, M. (2013). UCI Machine Learning Repository \url{http://archive.ics.uci.edu/ml}.
#' Irvine, CA: University of California, School of Information and Computer Science.
"breastcancer"