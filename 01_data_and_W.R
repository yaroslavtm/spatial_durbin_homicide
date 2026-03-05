# 0. Библиотеки
library(readxl)
library(sf)
library(spdep)
library(plm)
library(stargazer)

# 1. Загрузка панельных данных
df <- read_excel("ваш/путь/Данные.xlsx")
# Переименование столбцов
colnames(df) <- c("year", "fips", "abbr", "state", "homicide", "population",
                  "income_nominal", "homicide_rate", "ethanol", "unemployment",
                  "income_real", "opioid_deaths", "opioid_rate", "fs_ratio",
                  "gini", "incarc_rate",
                  "u16", "u17", "u18", "u19", "u20", "u21")

# Оставляем только нужные переменные
df <- df[, c("year", "fips", "abbr", "state", "homicide_rate", "income_real",
             "ethanol", "unemployment", "opioid_rate", "fs_ratio", "gini",
             "incarc_rate")]
cat("Размерность панели:", nrow(df), "наблюдений\n")
cat("Штаты:", length(unique(df$state)), "\n")
cat("Годы:", min(df$year), "-", max(df$year), "\n")

# 2. Загрузка шейпфайла и построение матрицы Queen
states_sf <- st_read("C:/Users/Downloads/cb_2018_us_state_500k/cb_2018_us_state_500k.shp", quiet = TRUE)

# Фильтрация: 48 смежных штатов
exclude_fips <- c("02", "15", "11", "60", "66", "69", "72", "78")
states48 <- states_sf[!states_sf$STATEFP %in% exclude_fips, ]
states48 <- states48[order(states48$STATEFP), ]

cat("Штатов в шейпфайле:", nrow(states48), "\n")

# Построение матрицы Queen-contiguity
w_queen <- poly2nb(states48, queen = TRUE)

# Проверка: число соседей
cat("\nСтруктура матрицы смежности:\n")
summary(w_queen)

# Row-standardized listw объект (нужен для всех пространственных моделей)
W_listw <- nb2listw(w_queen, style = "W")

# Извлечение полной матрицы 48×48 (для визуализации)
W_matrix <- listw2mat(W_listw)
cat("\nРазмерность W:", dim(W_matrix), "\n")
cat("Суммы строк (должны быть 1):", range(rowSums(W_matrix)), "\n")

# 3. важно: порядок штатов в W и в данных должен совпадать
# иначе будут некорректные выводы
# Порядок штатов в матрице W
w_order <- states48$NAME
cat("\nПорядок штатов в W:\n")
print(w_order)

# Проверяем что все штаты из данных есть в W и наоборот
cat("\nВсе штаты совпадают:", setequal(unique(df$state), w_order), "\n")

# Сортируем данные: сначала по штату (в порядке W), потом по году
df$state_f <- factor(df$state, levels = w_order)
df <- df[order(df$state_f, df$year), ]
df$state_f <- NULL

cat("Первые строки отсортированных данных:\n")
print(head(df[, c("state", "year", "homicide_rate")], 22))

# 4. Описательная статистика
vars <- c("homicide_rate", "income_real", "ethanol", "unemployment",
          "opioid_rate", "fs_ratio", "gini", "incarc_rate")

# Приводим к numeric на всякий случай
for (v in vars) df[[v]] <- as.numeric(df[[v]])

desc <- data.frame(
  Variable     = vars,
  Mean         = sapply(df[vars], mean, na.rm = TRUE),
  SD           = sapply(df[vars], sd, na.rm = TRUE),
  Min          = sapply(df[vars], min, na.rm = TRUE),
  Max          = sapply(df[vars], max, na.rm = TRUE),
  row.names    = NULL
)

cat("\n=== Описательная статистика ===\n")
print(data.frame(desc[1], round(desc[-1], 3)))

# Таблица для LaTeX / Word (stargazer)
stargazer(as.data.frame(df[vars]),
          type = "text",
          title = "Описательная статистика панельных данных, 2008-2018",
          digits = 3)

# 5. Корреляционная матрица
cor_matrix <- cor(df[vars], use = "complete.obs")
cat("\n=== Корреляционная матрица ===\n")
print(round(cor_matrix, 3))

# 6. Сохранение
# Сохраняем W и отсортированные данные для следующих скриптов
save(df, W_listw, w_queen, states48, file = "ваш/путь/workspace_01.RData")
cat("\nРабочее пространство сохранено в workspace_01.RData\n")

