utils::globalVariables(
  c(
    "A1_ADT_normal",
    "A1_coord",
    "A1_gene_pca",
    "A1_true_label",
    "MMTV_ADT_normal",
    "MMTV_gene_pca",
    "MMTV_coord",
    "sim1_1_data_1",
    "sim1_1_data_2",
    "sim1_1_data_coord",
    "sim1_1_true_label",
    "sim1_2_data_1",
    "sim1_2_data_2",
    "sim1_2_data_coord",
    "sim1_2_true_label",
    "sim2_1_data_1",
    "sim2_1_data_2",
    "sim2_1_data_coord",
    "sim2_1_true_label",
    "sim2_2_data_1",
    "sim2_2_data_2",
    "sim2_2_data_coord",
    "sim2_2_true_label",
    "sim2_3_data_1",
    "sim2_3_data_2",
    "sim2_3_data_coord",
    "sim2_3_true_label",
    "sim2_4_data_1",
    "sim2_4_data_2",
    "sim2_4_data_coord",
    "sim2_4_true_label"
  )
)

#' Load A1 real application dataset
#' @export
load_A1_data <- function() {
  list(
    ADT_normal = A1_ADT_normal,
    coord     = A1_coord,
    gene_pca   = A1_gene_pca,
    true_label = A1_true_label
  )
}

#' Load MMTV real application dataset
#' @export
load_MMTV_data <- function() {
  list(ADT_normal = MMTV_ADT_normal,
       coord = MMTV_coord,
       gene_pca = MMTV_gene_pca)
}

#' Load sim1_1 simulation dataset
#' @export
load_sim1_1_data <- function() {
  list(
    data1 = sim1_1_data_1,
    data2 = sim1_1_data_2,
    coord = sim1_1_data_coord,
    true_label = sim1_1_true_label
  )
}

#' Load sim1_2 simulation dataset
#' @export
load_sim1_2_data <- function() {
  list(
    data1 = sim1_2_data_1,
    data2 = sim1_2_data_2,
    coord = sim1_2_data_coord,
    true_label = sim1_2_true_label
  )
}

#' Load sim2_1 simulation dataset
#' @export
load_sim2_1_data <- function() {
  list(
    data1 = sim2_1_data_1,
    data2 = sim2_1_data_2,
    coord = sim2_1_data_coord,
    true_label = sim2_1_true_label
  )
}

#' Load sim2_2 simulation dataset
#' @export
load_sim2_2_data <- function() {
  list(
    data1 = sim2_2_data_1,
    data2 = sim2_2_data_2,
    coord = sim2_2_data_coord,
    true_label = sim2_2_true_label
  )
}

#' Load sim2_3 simulation dataset
#' @export
load_sim2_3_data <- function() {
  list(
    data1 = sim2_3_data_1,
    data2 = sim2_3_data_2,
    coord = sim2_3_data_coord,
    true_label = sim2_3_true_label
  )
}

#' Load sim2_4 simulation dataset
#' @export
load_sim2_4_data <- function() {
  list(
    data1 = sim2_4_data_1,
    data2 = sim2_4_data_2,
    coord = sim2_4_data_coord,
    true_label = sim2_4_true_label
  )
}
