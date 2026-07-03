#' Generate Sim2_2 simulation data (normal distribution, mixed L1 and L2)
#'
#' This function generates a synthetic spatial multi-omics dataset for Simulation 2_2.
#'  Feature values are drawn from normal distributions with domain-specific means.
#'

#' @param seed Random seed for reproducibility. Default is NULL (no seed set).
#' @return A list with four components:
#'   \item{data_1}{Data frame of first modality (n spots x 20 genes).}
#'   \item{data_2}{Data frame of second modality (n spots x 20 genes).}
#'   \item{coord}{Data frame of spatial coordinates (columns: x_spot, y_spot).}
#'   \item{true_label}{Data frame with column 'label' (integer domain labels 0-5).}
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rt
#' @export
#'
#' @examples
#' \dontrun{
#' sim_data <- generate_sim2_2_data(seed = 2011)
#' str(sim_data)
#' }

generate_sim2_2_data <- function(seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  
  # =========================================================================================
  # 22x22 spots
  # =========================================================================================
  x <- rep(1:22, each = 22)
  y <- rep(1:22, times = 22)
  coord <- data.frame(x, y)
  
  # Generate C_label
  coord$C <- ifelse(
    (coord$x %in% 5:8 & coord$y == 22) |
      (coord$x %in% 3:10 & coord$y == 21) |
      (coord$x %in% 4:11 & coord$y == 20) |
      (coord$x %in% 5:12 & coord$y == 19) |
      (coord$x %in% 6:12 & coord$y == 18) |
      (coord$x %in% 7:12 & coord$y == 17) |
      (coord$x %in% 8:11 & coord$y == 16) |
      (coord$x %in% 9:11 & coord$y == 15) |
      (coord$x %in% 10:11 & coord$y %in% 10:14) |
      (coord$x %in% 3:4 & coord$y == 10) |
      (coord$x %in% 8:11 & coord$y == 10) |
      (coord$x %in% 3:9 & coord$y == 9) |
      (coord$x %in% 2:8 & coord$y == 8) |
      (coord$x %in% 2:6 & coord$y %in% 6:7) |
      (coord$x %in% 2:5 & coord$y %in% 4:5) |
      (coord$x %in% 2:4 & coord$y == 3) |
      (coord$x %in% 3:4 & coord$y %in% 1:2),
    1,
    NA
  )
  
  coord$C <- ifelse(
    (coord$x == 2 & coord$y %in% 13:18) |
      (coord$x == 3 & coord$y %in% 11:20) |
      (coord$x == 4 & coord$y %in% 11:19) |
      (coord$x == 5 & coord$y %in% 10:18) |
      (coord$x == 6 & coord$y %in% 10:17) |
      (coord$x == 7 & coord$y %in% 10:16) |
      (coord$x == 8 & coord$y %in% 11:15) |
      (coord$x == 9 & coord$y %in% 11:14),
    2,
    coord$C
  )
  
  coord$C <- ifelse(
    (coord$x == 5 & coord$y %in% 2:3) |
      (coord$x == 6 & coord$y %in% 2:5) |
      (coord$x %in% 7:8 & coord$y %in% 2:7) |
      (coord$x == 9 & coord$y %in% 3:8) |
      (coord$x == 10 & coord$y %in% 3:9) |
      (coord$x == 11 & coord$y %in% 4:9) |
      (coord$x == 12 & coord$y %in% 4:11) |
      (coord$x == 13 & coord$y %in% 3:10) |
      (coord$x == 14 & coord$y %in% 3:9) |
      (coord$x == 15 & coord$y %in% 2:7) |
      (coord$x == 16 & coord$y %in% 2:5) |
      (coord$x == 17 & coord$y %in% 2:3) |
      (coord$x == 18 & coord$y == 2),
    3,
    coord$C
  )
  
  coord$C <- ifelse(
    (coord$x %in% 12:20 & coord$y == 12) |
      (coord$x %in% 13:20 & coord$y == 11) |
      (coord$x %in% 14:21 & coord$y == 10) |
      (coord$x %in% 15:21 & coord$y %in% 8:9) |
      (coord$x %in% 16:21 & coord$y %in% 6:7) |
      (coord$x %in% 17:22 & coord$y %in% 4:5) |
      (coord$x %in% 18:21 & coord$y == 3) |
      (coord$x %in% 19:20 & coord$y == 2),
    4,
    coord$C
  )
  
  coord$C <- ifelse(
    (coord$x == 12 & coord$y %in% 13:16) |
      (coord$x %in% 13:20 & coord$y %in% 13:19) |
      (coord$x %in% 17:20 & coord$y == 20) |
      (coord$x %in% 15:21 & coord$y == 20) |
      (coord$x == 21 & coord$y %in% 15:20) |
      (coord$x == 22 & coord$y %in% 17:19),
    5,
    coord$C
  )
  
  coord$C[is.na(coord$C)] <- 0
  coord <- coord[coord$C != 0, ]
  true_label <- data.frame(label = coord$C)
  
  # Generate label_L1
  data_1 <- coord[, c("x", "y")]
  
  data_1$L1 <- ifelse(
    (data_1$x %in% 5:8 & data_1$y == 22) |
      (data_1$x %in% 3:10 & data_1$y == 21) |
      (data_1$x %in% 4:11 & data_1$y == 20) |
      (data_1$x %in% 5:12 & data_1$y == 19) |
      (data_1$x %in% 6:12 & data_1$y == 18) |
      (data_1$x %in% 7:12 & data_1$y == 17) |
      (data_1$x %in% 8:11 & data_1$y == 16) |
      (data_1$x %in% 9:11 & data_1$y == 15) |
      (data_1$x %in% 10:11 & data_1$y %in% 10:14) |
      (data_1$x %in% 3:4 & data_1$y == 10) |
      (data_1$x %in% 8:11 & data_1$y == 10) |
      (data_1$x %in% 3:9 & data_1$y == 9) |
      (data_1$x %in% 2:8 & data_1$y == 8) |
      (data_1$x %in% 2:6 & data_1$y %in% 6:7) |
      (data_1$x %in% 2:5 & data_1$y %in% 4:5) |
      (data_1$x %in% 2:4 & data_1$y == 3) |
      (data_1$x %in% 3:4 & data_1$y %in% 1:2),
    1,
    NA
  )
  
  data_1$L1 <- ifelse(
    (data_1$x == 2 & data_1$y %in% 13:18) |
      (data_1$x == 3 & data_1$y %in% 11:20) |
      (data_1$x == 4 & data_1$y %in% 11:19) |
      (data_1$x == 5 & data_1$y %in% 10:18) |
      (data_1$x == 6 & data_1$y %in% 10:17) |
      (data_1$x == 7 & data_1$y %in% 10:16) |
      (data_1$x == 8 & data_1$y %in% 11:15) |
      (data_1$x == 9 & data_1$y %in% 11:14),
    2,
    data_1$L1
  )
  
  data_1$L1 <- ifelse(
    (data_1$x == 5 & data_1$y %in% 2:3) |
      (data_1$x == 6 & data_1$y %in% 2:5) |
      (data_1$x %in% 7:8 & data_1$y %in% 2:7) |
      (data_1$x == 9 & data_1$y %in% 3:8) |
      (data_1$x == 10 & data_1$y %in% 3:9) |
      (data_1$x == 11 & data_1$y %in% 4:9) |
      (data_1$x == 12 & data_1$y %in% 4:11) |
      (data_1$x == 13 & data_1$y %in% 3:10) |
      (data_1$x == 14 & data_1$y %in% 3:9) |
      (data_1$x == 15 & data_1$y %in% 2:7) |
      (data_1$x == 16 & data_1$y %in% 2:5) |
      (data_1$x == 17 & data_1$y %in% 2:3) |
      (data_1$x == 18 & data_1$y == 2),
    3,
    data_1$L1
  )
  
  data_1$L1 <- ifelse(
    (data_1$x %in% 12:20 & data_1$y == 12) |
      (data_1$x %in% 13:20 & data_1$y == 11) |
      (data_1$x %in% 14:21 & data_1$y == 10) |
      (data_1$x %in% 15:21 & data_1$y %in% 8:9) |
      (data_1$x %in% 16:21 & data_1$y %in% 6:7) |
      (data_1$x %in% 17:22 & data_1$y %in% 4:5) |
      (data_1$x %in% 18:21 & data_1$y == 3) |
      (data_1$x %in% 19:20 & data_1$y == 2),
    3,
    data_1$L1
  )
  
  data_1$L1 <- ifelse(
    (data_1$x == 12 & data_1$y %in% 13:16) |
      (data_1$x %in% 13:20 & data_1$y %in% 13:19) |
      (data_1$x %in% 17:20 & data_1$y == 20) |
      (data_1$x %in% 15:21 & data_1$y == 20) |
      (data_1$x == 21 & data_1$y %in% 15:20) |
      (data_1$x == 22 & data_1$y %in% 17:19),
    5,
    data_1$L1
  )
  
  data_1$L1[which(data_1$L1 == 3)] <- sample(c(3, 4), length(which(data_1$L1 == 3)), replace = TRUE)
  random_idx_1 <- sample(1:nrow(data_1), size = 50, replace = FALSE)
  data_1$L1[random_idx_1] <- sample(c(3, 4), size = length(random_idx_1), replace = TRUE)
  
  # Generate label_L2
  data_2 <- coord[, c("x", "y")]
  
  data_2$L2 <- ifelse(
    (data_2$x %in% 5:8 & data_2$y == 22) |
      (data_2$x %in% 3:10 & data_2$y == 21) |
      (data_2$x %in% 4:11 & data_2$y == 20) |
      (data_2$x %in% 5:12 & data_2$y == 19) |
      (data_2$x %in% 6:12 & data_2$y == 18) |
      (data_2$x %in% 7:12 & data_2$y == 17) |
      (data_2$x %in% 8:11 & data_2$y == 16) |
      (data_2$x %in% 9:11 & data_2$y == 15) |
      (data_2$x %in% 10:11 & data_2$y %in% 10:14) |
      (data_2$x %in% 3:4 & data_2$y == 10) |
      (data_2$x %in% 8:11 & data_2$y == 10) |
      (data_2$x %in% 3:9 & data_2$y == 9) |
      (data_2$x %in% 2:8 & data_2$y == 8) |
      (data_2$x %in% 2:6 & data_2$y %in% 6:7) |
      (data_2$x %in% 2:5 & data_2$y %in% 4:5) |
      (data_2$x %in% 2:4 & data_2$y == 3) |
      (data_2$x %in% 3:4 & data_2$y %in% 1:2),
    1,
    NA
  )
  
  data_2$L2 <- ifelse(
    (data_2$x == 2 & data_2$y %in% 13:18) |
      (data_2$x == 3 & data_2$y %in% 11:20) |
      (data_2$x == 4 & data_2$y %in% 11:19) |
      (data_2$x == 5 & data_2$y %in% 10:18) |
      (data_2$x == 6 & data_2$y %in% 10:17) |
      (data_2$x == 7 & data_2$y %in% 10:16) |
      (data_2$x == 8 & data_2$y %in% 11:15) |
      (data_2$x == 9 & data_2$y %in% 11:14),
    1,
    data_2$L2
  )
  
  data_2$L2 <- ifelse(
    (data_2$x == 5 & data_2$y %in% 2:3) |
      (data_2$x == 6 & data_2$y %in% 2:5) |
      (data_2$x %in% 7:8 & data_2$y %in% 2:7) |
      (data_2$x == 9 & data_2$y %in% 3:8) |
      (data_2$x == 10 & data_2$y %in% 3:9) |
      (data_2$x == 11 & data_2$y %in% 4:9) |
      (data_2$x == 12 & data_2$y %in% 4:11) |
      (data_2$x == 13 & data_2$y %in% 3:10) |
      (data_2$x == 14 & data_2$y %in% 3:9) |
      (data_2$x == 15 & data_2$y %in% 2:7) |
      (data_2$x == 16 & data_2$y %in% 2:5) |
      (data_2$x == 17 & data_2$y %in% 2:3) |
      (data_2$x == 18 & data_2$y == 2),
    3,
    data_2$L2
  )
  
  data_2$L2 <- ifelse(
    (data_2$x %in% 12:20 & data_2$y == 12) |
      (data_2$x %in% 13:20 & data_2$y == 11) |
      (data_2$x %in% 14:21 & data_2$y == 10) |
      (data_2$x %in% 15:21 & data_2$y %in% 8:9) |
      (data_2$x %in% 16:21 & data_2$y %in% 6:7) |
      (data_2$x %in% 17:22 & data_2$y %in% 4:5) |
      (data_2$x %in% 18:21 & data_2$y == 3) |
      (data_2$x %in% 19:20 & data_2$y == 2),
    1,
    data_2$L2
  )
  
  data_2$L2 <- ifelse(
    (data_2$x == 12 & data_2$y %in% 13:16) |
      (data_2$x %in% 13:20 & data_2$y %in% 13:19) |
      (data_2$x %in% 17:20 & data_2$y == 20) |
      (data_2$x %in% 15:21 & data_2$y == 20) |
      (data_2$x == 21 & data_2$y %in% 15:20) |
      (data_2$x == 22 & data_2$y %in% 17:19),
    1,
    data_2$L2
  )
  
  random_idx_2 <- sample(1:nrow(data_2), size = 20, replace = FALSE)
  data_2$L2[random_idx_2] <- 3
  
  # =========================================================================================
  # Generate feature values using normal distributions
  # =========================================================================================
  n_features <- 20
  params_L1 <- data.frame(
    group = c(1, 2, 3, 4, 5),
    mean = c(-2.5, -1.2, 0.3, 1.8, 2.9),
    sd   = c(0.1, 0.1, 0.1, 0.1, 0.1)
  )
  params_L2 <- data.frame(
    group = c(1, 3),
    mean = c(-0.9, 1.1),
    sd   = c(0.1, 0.1)
  )
  
  for (i in 1:n_features) {
    new_col <- numeric(nrow(data_1))
    for (j in 1:nrow(params_L1)) {
      group_val <- params_L1$group[j]
      Factor <- data_1$L1 == group_val
      new_col[Factor] <- rnorm(sum(Factor), params_L1$mean[j], params_L1$sd[j])
    }
    data_1[[paste0("gene_", i)]] <- new_col
  }
  
  for (i in 1:n_features) {
    new_col <- numeric(nrow(data_2))
    for (j in 1:nrow(params_L2)) {
      group_val <- params_L2$group[j]
      Factor <- data_2$L2 == group_val
      new_col[Factor] <- rnorm(sum(Factor), params_L2$mean[j], params_L2$sd[j])
    }
    data_2[[paste0("gene_", i)]] <- new_col
  }
  
  
  coord <- coord[, -3]
  colnames(coord) <- c("x_spot", "y_spot")
  data_1 <- data_1[, !names(data_1) %in% c("x", "y", "L1")]
  data_2 <- data_2[, !names(data_2) %in% c("x", "y", "L2")]
  
  invisible(list(
    data_1 = data_1,
    data_2 = data_2,
    coord = coord,
    true_label = true_label
  ))
}
