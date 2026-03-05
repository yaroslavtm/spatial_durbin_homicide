# 0. Библиотеки и загрузка данных
load("ваш/путь/workspace_01.RData")
library(spdep)
library(spatialreg)
library(splm)
library(plm)
library(SDPDmod)
library(car)
N <- 48
years <- 2008:2018

# 1. Pooled OLS
ols <- lm(
  homicide_rate ~ income_real + ethanol + unemployment +
    opioid_rate + fs_ratio + gini + incarc_rate,
  data = df
)
summary(ols)

# 2. VIF
print(vif(ols))

# 3. Panel Two-Way FE
pdf <- pdata.frame(df, index = c("state", "year"))

fe_twoway <- plm(
  homicide_rate ~ income_real + ethanol + unemployment +
    opioid_rate + fs_ratio + gini + incarc_rate,
  data = pdf,
  model = "within",
  effect = "twoways"
)
summary(fe_twoway)

# 4. Hausman Test (FE vs RE)
fe_ind <- plm(
  homicide_rate ~ income_real + ethanol + unemployment +
    opioid_rate + fs_ratio + gini + incarc_rate,
  data = pdf,
  model = "within",
  effect = "individual"
)
re_ind <- plm(
  homicide_rate ~ income_real + ethanol + unemployment +
    opioid_rate + fs_ratio + gini + incarc_rate,
  data = pdf,
  model = "random",
  effect = "individual"
)
print(phtest(fe_ind, re_ind))

# 5. Moran's I по годам (остатки cross-section OLS)
moran_results <- data.frame(year = integer(),
                            I = numeric(),
                            p_value = numeric())

for (yr in years) {
  sub <- df[df$year == yr, ]
  ols_yr <- lm(
    homicide_rate ~ income_real + ethanol + unemployment +
      opioid_rate + fs_ratio + gini + incarc_rate,
    data = sub
  )
  mt <- moran.test(residuals(ols_yr), W_listw)
  moran_results <- rbind(moran_results,
                         data.frame(
                           year = yr,
                           I = as.numeric(mt$estimate[1]),
                           p_value = mt$p.value
                         ))
}

row.names(moran_results) <- NULL
print(round(moran_results, 4))


# 6. LM-тесты (кросс-секционные средние)
vars <- c(
  "homicide_rate",
  "income_real",
  "ethanol",
  "unemployment",
  "opioid_rate",
  "fs_ratio",
  "gini",
  "incarc_rate"
)
df_mean <- aggregate(df[vars], by = list(state = df$state), FUN = mean)
df_mean <- df_mean[match(states48$NAME, df_mean$state), ]

ols_cs <- lm(
  homicide_rate ~ income_real + ethanol + unemployment +
    opioid_rate + fs_ratio + gini + incarc_rate,
  data = df_mean
)

print(lm.RStests(ols_cs, W_listw, test = "all"))

# 7. Панельные LM-тесты (splm)
print(
  slmtest(
    homicide_rate ~ income_real + ethanol + unemployment +
      opioid_rate + fs_ratio + gini + incarc_rate,
    data = pdf,
    listw = W_listw,
    test = "lml"
  )
)

print(
  slmtest(
    homicide_rate ~ income_real + ethanol + unemployment +
      opioid_rate + fs_ratio + gini + incarc_rate,
    data = pdf,
    listw = W_listw,
    test = "lme"
  )
)


# 8. SDM (Spatial Durbin Model, Two-Way FE)
df_sorted <- df[order(df$year, match(df$state, states48$NAME)), ]
sdm_data <- data.frame(
  state = df_sorted$state,
  year = df_sorted$year,
  y = df_sorted$homicide_rate,
  df_sorted[, c(
    "income_real",
    "ethanol",
    "unemployment",
    "opioid_rate",
    "fs_ratio",
    "gini",
    "incarc_rate"
  )]
)
W_bin <- listw2mat(nb2listw(w_queen, style = "B"))
W_rs <- W_bin / rowSums(W_bin)
sdm <- SDPDm(
  formula = y ~ income_real + ethanol + unemployment +
    opioid_rate + fs_ratio + gini + incarc_rate,
  data = sdm_data,
  W = W_rs,
  index = c("state", "year"),
  model = "sdm",
  effect = "twoways",
  LYtrans = TRUE
)
summary(sdm)

# 9. SAR (для сравнения с SDM)
sar <- SDPDm(
  formula = y ~ income_real + ethanol + unemployment +
    opioid_rate + fs_ratio + gini + incarc_rate,
  data = sdm_data,
  W = W_rs,
  index = c("state", "year"),
  model = "sar",
  effect = "twoways",
  LYtrans = TRUE
)
summary(sar)

# 10. Декомпозиция эффектов (с p значениями через метод Монте-Карло
library(MASS)
rho <- as.numeric(sdm$rho)
betas <- as.numeric(sdm$coefficients[1:7])
thetas <- as.numeric(sdm$coefficients[8:14])
VC <- sdm$varcov[1:15, 1:15]

var_names <- c(
  "income_real",
  "ethanol",
  "unemployment",
  "opioid_rate",
  "fs_ratio",
  "gini",
  "incarc_rate"
)
I_mat <- diag(N)
n_sim <- 10000
set.seed(42)

mean_vec <- c(rho, betas, thetas)
draws <- mvrnorm(n_sim, mean_vec, VC)

direct_sim   <- matrix(0, n_sim, 7)
indirect_sim <- matrix(0, n_sim, 7)
total_sim    <- matrix(0, n_sim, 7)

for (s in 1:n_sim) {
  rho_s   <- draws[s, 1]
  beta_s  <- draws[s, 2:8]
  theta_s <- draws[s, 9:15]
  A_inv   <- solve(I_mat - rho_s * W_rs)
  
  for (k in 1:7) {
    S_k <- A_inv %*% (beta_s[k] * I_mat + theta_s[k] * W_rs)
    direct_sim[s, k]   <- sum(diag(S_k)) / N
    total_sim[s, k]    <- sum(S_k) / N
    indirect_sim[s, k] <- total_sim[s, k] - direct_sim[s, k]
  }
}

# Точечные оценки, SE и p-values
cat("Direct effects:\n")
for (k in 1:7) {
  m <- mean(direct_sim[, k])
  s <- sd(direct_sim[, k])
  p <- 2 * pnorm(-abs(m / s))
  sig <- ifelse(p <0.001, "***", ifelse(p <0.01, "**", ifelse(p <0.05, "*", ifelse(p <0.1, ".", ""))))
  cat(sprintf("  %-15s  %10.6f  SE=%9.6f  p=%.4f %s\n", var_names[k], m, s, p, sig))
}

cat("\nIndirect effects:\n")
for (k in 1:7) {
  m <- mean(indirect_sim[, k])
  s <- sd(indirect_sim[, k])
  p <- 2 * pnorm(-abs(m / s))
  sig <- ifelse(p <0.001, "***", ifelse(p <0.01, "**", ifelse(p <0.05, "*", ifelse(p <0.1, ".", ""))))
  cat(sprintf("  %-15s  %10.6f  SE=%9.6f  p=%.4f %s\n", var_names[k], m, s, p, sig))
}

cat("\nTotal effects:\n")
for (k in 1:7) {
  m <- mean(total_sim[, k])
  s <- sd(total_sim[, k])
  p <- 2 * pnorm(-abs(m / s))
  sig <- ifelse(p <0.001, "***", ifelse(p <0.01, "**", ifelse(p <0.05, "*", ifelse(p <0.1, ".", ""))))
  cat(sprintf("  %-15s  %10.6f  SE=%9.6f  p=%.4f %s\n", var_names[k], m, s, p, sig))
}

# Таблица точечных оценок
A_inv_pt <- solve(I_mat - rho * W_rs)
Direct <- Indirect <- Total <- numeric(7)
for (k in 1:7) {
  S <- A_inv_pt %*% (betas[k] * I_mat + thetas[k] * W_rs)
  Direct[k]   <- sum(diag(S)) / N
  Total[k]    <- sum(S) / N
  Indirect[k] <- Total[k] - Direct[k]
}
effects_table <- data.frame(
  Variable = var_names,
  Direct   = round(Direct, 6),
  Indirect = round(Indirect, 6),
  Total    = round(Total, 6)
)
cat("\nТочечные оценки:\n")
print(effects_table)
# 11. LR-тест SDM vs SAR
ll_sdm <- sdm$likl
ll_sar <- sar$likl
lr_stat <- 2 * (ll_sdm - ll_sar)
p_val <- 1 - pchisq(lr_stat, 7)
cat("LR(SDM vs SAR) =",
    round(lr_stat, 3),
    " df = 7  p =",
    format(p_val, digits = 4),
    "\n")

# 12. Сохранение
save(
  effects_table,
  moran_results,
  rho,
  betas,
  thetas,
  W_rs,
  sdm_data,
  direct_sim,
  indirect_sim,
  total_sim,
  file = "ваш/путь/workspace_02.RData"
)


