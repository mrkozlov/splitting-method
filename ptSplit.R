# Подключение библиотек
library(FNN)      # Для KNN
library(ggplot2)  # Для визуализации
library(svglite)  # Для сохранения графиков
library(dplyr)    # Для манипуляций с данными
library(terra) # Для работы с пространственными данными
library(sf)
#pts<-train_reduced
#pts<- st_read("dra70_2000.gpkg")

#pts<-train_balanced
 #pts$Duplicates<-NULL
#pts$geom<-NULL

#pts$y<-pts$decimalLati
 #pts$x<-pts$decimalLong
#pts$decimalLong<-NULL
 #pts$decimalLati<-NULL
 
# pts<-as.data.frame(pts.sp1)
 #pts$layer<-NULL
 
params <- list(
  test_size = 0.2,            # Устанавливаем 20% для тестового набора (33-34 точки из 167)
  max_attempts = 90000,       # Увеличиваем до 90000 для более тщательного поиска оптимального разбиения
  n_clusters = 12,            # Устанавливаем 17 кластеров, чтобы распределить 33-34 тестовые точки (по 2 точки на кластер, с 1 кластером на 1-2 точки)
  w_cross_dist = 7.0,         # Увеличиваем до 7.0, чтобы значительно увеличить межгрупповое расстояние и снизить автокорреляцию
  w_coverage = 9.0,           # Увеличиваем до 9.0, чтобы гарантировать полное покрытие всех кластеров
  w_balance = 60.0,           # Увеличиваем до 60.0, чтобы добиться равномерного распределения тестовых точек
  w_dist_diff = 0.1,          # Уменьшаем до 0.1, чтобы сохранить малую разницу в расстояниях (текущая 0.079 уже хороша)
  max_points_per_cluster = 2  # Оставляем 2, чтобы строго ограничить число точек в кластере и избежать переполнения (как в кластерах 1, 8, 9 с 12 точками)
)

# 2. Функции для вычисления пространственных метрик

# Вычисление внутригрупповых расстояний
check_spatial_metrics <- function(data) {
  if (nrow(data) < 2) {
    return(list(mean_dist = 0, min_dist = 0, max_dist = 0))
  }
  
  coords <- as.matrix(data[, c("x", "y")])
  dist_matrix <- as.matrix(dist(coords))
  
  return(list(
    mean_dist = mean(dist_matrix[upper.tri(dist_matrix)]),
    min_dist = min(dist_matrix[upper.tri(dist_matrix)]),
    max_dist = max(dist_matrix[upper.tri(dist_matrix)])
  ))
}

# Вычисление межгруппового расстояния
check_cross_distance <- function(train_data, test_data) {
  if (nrow(train_data) == 0 || nrow(test_data) == 0) {
    return(0)
  }
  
  knn <- get.knnx(train_data[, c("x", "y")], 
                  test_data[, c("x", "y")], 
                  k = 1)
  return(mean(knn$nn.dist))
}

# 3. Улучшенная функция разбиения на основе кластеров (без нормализации)
balanced_cluster_split <- function(pts, params) {
  # Проверка наличия столбцов x и y
  if (!all(c("x", "y") %in% colnames(pts))) {
    stop("Данные должны содержать колонки 'x' и 'y'")
  }
  
  # Проверка числового типа данных
  if (!is.numeric(pts$x) || !is.numeric(pts$y)) {
    stop("Столбцы 'x' и 'y' должны быть числовыми")
  }
  
  # Преобразование в числовую матрицу для kmeans
  coords <- as.matrix(pts[, c("x", "y")])
  
  # Кластеризация для учета плотности
  set.seed(42)
  clusters <- kmeans(coords, centers = params$n_clusters)$cluster
  pts$cluster <- as.factor(clusters)
  
  # Вывод распределения всех точек по кластерам
  cat("\nРаспределение всех точек по кластерам:\n")
  print(table(pts$cluster))
  
  # Целевое количество точек в тестовом наборе
  target_test <- round(table(pts$cluster) * params$test_size)
  target_test[target_test == 0] <- 1  # Минимум 1 точка на кластер
  
  # Проверяем доступные кластеры
  cluster_counts <- table(pts$cluster)
  available_clusters <- names(cluster_counts[cluster_counts > 0])
  
  best_split <- NULL
  best_score <- -Inf
  
  for (attempt in 1:params$max_attempts) {
    # Отбор точек из каждого кластера
    test_indices <- unlist(lapply(names(target_test), function(cl) {
      cluster_pts <- which(pts$cluster == cl)
      sample(cluster_pts, min(target_test[cl], length(cluster_pts)), replace = FALSE)
    }))
    
    # Проверка покрытия кластеров
    temp_test_data <- pts[test_indices, ]
    cluster_coverage <- table(temp_test_data$cluster)
    if (!all(available_clusters %in% names(cluster_coverage))) {
      missing_clusters <- available_clusters[!available_clusters %in% names(cluster_coverage)]
      for (cl in missing_clusters) {
        cluster_pts <- which(pts$cluster == cl & !(1:nrow(pts) %in% test_indices))
        if (length(cluster_pts) > 0) {
          test_indices <- c(test_indices, sample(cluster_pts, 1))
        }
      }
    }
    
    # Ограничение максимального числа точек в кластере
    temp_test_data <- pts[test_indices, ]
    cluster_counts_test <- table(temp_test_data$cluster)
    for (cl in names(cluster_counts_test)) {
      if (cluster_counts_test[cl] > params$max_points_per_cluster) {
        # Удаляем лишние точки из кластера
        cluster_indices <- which(temp_test_data$cluster == cl)
        excess <- cluster_counts_test[cl] - params$max_points_per_cluster
        remove_indices <- sample(cluster_indices, excess)
        test_indices <- test_indices[!test_indices %in% temp_test_data$index[remove_indices]]
      }
    }
    
    # Добавляем индекс для отслеживания строк
    pts$index <- 1:nrow(pts)
    temp_test_data <- pts[test_indices, ]
    
    # Проверка покрытия после ограничения
    cluster_coverage <- table(temp_test_data$cluster)
    if (!all(available_clusters %in% names(cluster_coverage))) {
      missing_clusters <- available_clusters[!available_clusters %in% names(cluster_coverage)]
      for (cl in missing_clusters) {
        cluster_pts <- which(pts$cluster == cl & !(1:nrow(pts) %in% test_indices))
        if (length(cluster_pts) > 0) {
          test_indices <- c(test_indices, sample(cluster_pts, 1))
        }
      }
    }
    
    # Создание тренировочного и тестового наборов
    train_data <- pts[-test_indices, ]
    test_data <- pts[test_indices, ]
    
    # Удаляем временный столбец index
    train_data$index <- NULL
    test_data$index <- NULL
    pts$index <- NULL
    
    # Вычисление метрик
    train_metrics <- check_spatial_metrics(train_data)
    test_metrics <- check_spatial_metrics(test_data)
    cross_dist <- check_cross_distance(train_data, test_data)
    coverage <- mean(table(test_data$cluster) > 0)  # Доля заполненных кластеров
    balance <- 1 - sd(table(test_data$cluster)) / mean(table(test_data$cluster))
    dist_diff <- abs(train_metrics$mean_dist - test_metrics$mean_dist)
    
    # Комбинированная оценка
    score <- (cross_dist * params$w_cross_dist + 
                coverage * params$w_coverage + 
                balance * params$w_balance - 
                dist_diff * params$w_dist_diff)
    
    if (score > best_score) {
      best_split <- list(train = train_data, test = test_data)
      best_score <- score
    }
  }
  
  return(best_split)
}

# ... (предыдущий код остаётся без изменений до этого места)

# 4. Определение экстента точек
# Вычисляем минимальные и максимальные значения координат x и y
#e <- extent(min(pts$x), max(pts$x), min(pts$y), max(pts$y))
cat("\nЭкстент точек:\n")
print(e)
e <- extent(floor(min(pts$x)), ceiling(max(pts$x)), floor(min(pts$y)), ceiling(max(pts$y)))
print(e)

# 4. Применение функции разбиения
set.seed(42)
split_result <- balanced_cluster_split(pts, params)

# ... (дальнейший код остаётся без изменений)
train_data <- split_result$train
test_data <- split_result$test

# 5. Вычисление метрик
train_metrics <- check_spatial_metrics(train_data)
test_metrics <- check_spatial_metrics(test_data)
cross_dist <- check_cross_distance(train_data, test_data)
coverage <- length(unique(test_data$cluster)) / params$n_clusters  # Покрытие кластеров

# 6. Визуализация
# Вычисление ближайших соседей для линий
knn <- get.knnx(train_data[, c("x", "y")], test_data[, c("x", "y")], k = 1)
segment_data <- data.frame(
  x = test_data$x,
  y = test_data$y,
  xend = train_data$x[knn$nn.index],
  yend = train_data$y[knn$nn.index]
)

# Построение графика
final_plot <- ggplot() +
  geom_point(data = pts, aes(x, y), color = "grey90", size = 2) +
  geom_point(data = train_data, aes(x, y, color = "Training"), size = 3, alpha = 0.7) +
  geom_point(data = test_data, aes(x, y, color = "Testing"), shape = 17, size = 4) +
  geom_segment(
    data = segment_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "darkgrey", linetype = "dashed", alpha = 0.5
  ) +
  scale_color_manual(
    name = "Dataset",
    values = c("Training" = "#3575b5", "Testing" = "#e74c3c")
  ) +
  labs(
    title = "Balanced Cluster-Based Data Split",
    subtitle = sprintf(
      "Training: %d points | Testing: %d points | Cross-dist: %.2f | Cluster coverage: %.0f%%",
      nrow(train_data), nrow(test_data), cross_dist, coverage * 100
    ),
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Сохранение графика в PNG, TIFF и PDF
ggsave("balanced_cluster_split_plot.png", final_plot, width = 8, height = 6, dpi = 300)
ggsave("balanced_cluster_split_plot.tiff", final_plot, width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave("balanced_cluster_split_plot.pdf", final_plot, width = 8, height = 6, dpi = 300)

# Сохранение данных
write.table(test_data, file = "testAria3.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(train_data, file = "trainAria3.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 7. Создание текстового отчета
report_text <- paste(
  "=== SPLITTING RESULTS ===",
  "\nParameters:",
  sprintf("- Test size: %.1f", params$test_size),
  sprintf("- Max attempts: %d", params$max_attempts),
  sprintf("- Number of clusters: %d", params$n_clusters),
  sprintf("- Cross-distance weight: %.1f", params$w_cross_dist),
  sprintf("- Coverage weight: %.1f", params$w_coverage),
  sprintf("- Balance weight: %.1f", params$w_balance),
  sprintf("- Distance difference weight: %.1f", params$w_dist_diff),
  sprintf("- Max points per cluster: %d", params$max_points_per_cluster),
  "\nDataset sizes:",
  sprintf("- Training set: %d points (%.1f%%)", nrow(train_data), nrow(train_data)/nrow(pts)*100),
  sprintf("- Testing set: %d points (%.1f%%)", nrow(test_data), nrow(test_data)/nrow(pts)*100),
  "\nSpatial metrics:",
  sprintf("- Mean distance in train: %.3f", train_metrics$mean_dist),
  sprintf("- Mean distance in test: %.3f", test_metrics$mean_dist),
  sprintf("- Cross-group distance: %.3f", cross_dist),
  sprintf("- Cluster coverage: %.1f%%", coverage * 100),
  "\nTesting points per cluster:",
  capture.output(print(table(test_data$cluster))),
  sep = "\n"
)

# Сохранение отчета в файл
writeLines(report_text, "splitting_report.txt")

# Вывод графика и отчета в консоль
print(final_plot)
cat(report_text)
make_extent_from_pts <- function(pts, digits = 1) {
  if (!all(c("x", "y") %in% names(pts))) {
    stop("Датафрейм должен содержать колонки 'x' и 'y'")
  }
  
  round_down <- function(x, digits) floor(x * 10^digits) / 10^digits
  round_up   <- function(x, digits) ceiling(x * 10^digits) / 10^digits
  
  min_lon <- round_down(min(pts$x, na.rm = TRUE), digits)
  max_lon <- round_up(max(pts$x, na.rm = TRUE), digits)
  min_lat <- round_down(min(pts$y, na.rm = TRUE), digits)
  max_lat <- round_up(max(pts$y, na.rm = TRUE), digits)
  
  extent(min_lon, max_lon, min_lat, max_lat)
}


e_30s <- make_extent_from_pts(pts)
print(e_30s)

