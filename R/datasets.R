#' @name MS
#' @docType data
#' @title Multiple Sclerosis data
#'
#' @description A \code{data.frame} with a design and proteomic data.
#'
#' @usage
#' data(MS)
#'
#' @details
#' Data from biobank are analysed a study population of 101 patients, 37 were diagnosed with
#' multiple sclerosis, and 64 without multiple sclerosis. Of the patients without multiple
#' sclerosis, 50 were diagnosed with other neurological disorders and 14 were neurologically
#' healthy patients who had undergone spinal anaesthesia for orthopaedic surgery on the knee
#' or ankle, i.e. neurologically healthy controls. Unless otherwise stated, all the patients
#' without multiple sclerosis were considered as controls for this study. All patients with
#' multiple sclerosis had relapsing remitting multiple sclerosis. The proteome were obtained
#' on cerebrospinal fluid samples from all patients prior medical treatment for multiple
#' sclerosis. It was discovered the patients separated into two clusters, called cluster 1 and
#' cluster 2. This is utilised in the data analysis by considering the data as 2-way factorial
#' design with the two factors: MS and clusters both on two levels.
#'
#' @author Ellen Færgestad Mosleth
#' @references
#' * Opsahl, J.A. et al. Label-free analysis of human cerebrospinal fluid addressing various normalization strategies and revealing protein groups affected by multiple sclerosis. Proteomics 16, 1154-1165 (2016).
#'
#' * Ellen Færgestad Mosleth, Christian Alexander Vedeler, Kristian Hovde Liland, Anette McLeod, Gerd Haga Bringland, Liesbeth Kroondijk, Frode Berven, Artem Lysenko, Christopher J. Rawlings, Karim El-Hajj Eid, Jill Anette Opsahl, Bjørn Tore Gjertsen, Kjell-Morten Myhr and Sonia Gavasso, Cerebrospinal fluid proteome shows disrupted neuronal development in multiple sclerosis. Scientific Reports – Nature 11(4087), (2021).
#'
#' @examples
#' data(MS)
#' str(MS)
#'
NULL

#' @name Lactobacillus
#' @docType data
#' @title Lactobacillus data
#'
#' @description A \code{data.frame} with a design and proteomic data, transcriptomic data and phenotypic data.
#'
#' @usage
#' data(Lactobacillus)
#'
#' @details
#' Experiment on Lactobacillus sakei was performed as a 2-way factorial design with two
#' factors both on two levels: strain (L. sakei strains LS25 and 23K) (factor A) and growth
#' condition (high and low glucose availability) (factor B) both on two levels, and their
#' interaction term (factor AB). There were three biological replicates within each group.
#' Transcriptome, proteome and end product profile (lactate, formate, acetate and ethanol)
#' were observed.
#'
#' @author Ellen Færgestad Mosleth
#' @references McLeod et al. 2017. Effects of glucose availability in Lactobacillus sakei; metabolic change and regulation of the proteome and transcriptome. Plos One 12, e0187542.
#'
#' @examples
#' data(Lactobacillus)
#' str(Lactobacillus)
#'
NULL

#' @name Diabetes
#' @docType data
#' @title Diabetes data
#'
#' @description A \code{data.frame} with a design and transcriptomic data.
#'
#' @usage
#' data(Diabetes)
#'
#' @details
#' Clinical study on humans was performed as a 2-way factorial design with two factors
#' both on two levels: bariatric surgery on two levels (before and after the bariatric
#' surgery) and type 2 diabetes (T2D) on two levels (with and without T2D). There were
#' 8 patients without T2D and 7 with T2D. It was discovered that the patients with T2D
#' would be separated in two groups: 3 patients in the group called T2D1 and 4 patients
#' in the group called T2D2. The experiment can therefore also analysed as 2 way factorial
#' design where the disease factor is on three levels. All patients were obese before
#' bariatric surgery (BMI >45). Transcriptome in the subcutaneous adipose tissue were
#' obtained before and one year after bariatric surgery.
#'
#' @author Ellen Færgestad Mosleth
#' @references Dankel et al. 2010. Switch from Stress Response to Homeobox Transcription Factors in Adipose Tissue After Profound Fat Loss. Plos One 5.
#'
#' @examples
#' data(Diabetes)
#' str(Diabetes)
#'
NULL
