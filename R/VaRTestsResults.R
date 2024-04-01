# Create vector of actual portfolio returns

setwd('~/UCL MSc Statistics/STAT0034 MSc Project/Code')

fig_dir <- 'C:/Users/jonat/OneDrive/Documents/UCL MSc Statistics/STAT0034 MSc Project/Report figures'

library(quantmod)
library(RColorBrewer)
library(stats)

colours <- brewer.pal(5, 'Set1')

start_date <- as.Date('01-01-2010', format = "%d-%m-%Y")
end_date <- as.Date('31-12-2019', format = "%d-%m-%Y")

getSymbols(Symbols = c('^FTSE', '^BVSP', '^GSPC'),
           src = 'yahoo',
           from = start_date,
           to = end_date,
           periodicity = 'daily')

# Check if data contains any NA values

if (any(is.na(GSPC))) {
  GSPC <- GSPC[complete.cases(GSPC), ]
}

if (any(is.na(FTSE))) {
  FTSE <- FTSE[complete.cases(FTSE), ]
}

if (length(setdiff(zoo::index(FTSE), zoo::index(GSPC))) != 0) {
  common_dates <- as.Date(intersect(zoo::index(FTSE), zoo::index(GSPC)), origin = "1970-01-01")
  GSPC <- GSPC[common_dates, ]
  FTSE <- FTSE[common_dates, ]
}

# Create vectors of returns
FTSENetRet <- diff(FTSE$FTSE.Close)/xts::lag.xts(FTSE$FTSE.Close, k = 1)
GSPCNetRet <- diff(GSPC$GSPC.Close)/xts::lag.xts(GSPC$GSPC.Close, k = 1)

# Create vector of portfolio net returns
weights <- c(0.5, 0.5)
portNetRet <- cbind(FTSENetRet, GSPCNetRet) %*% weights
portNetRet <- portNetRet[2:length(portNetRet)]

# Exprort net returns to file
portRet_filename <- paste0('~/UCL MSc Statistics/STAT0034 MSc Project/Code/VaRBacktest/portfolio_returns.csv')
write.table(portNetRet, portRet_filename, row.names = F, col.names = F, sep = ',')

plot(portNetRet[1001:2469], type = 'h')

portNetRetTest <- portNetRet[1001:2469]

HF_Var <- read.csv('./VaRBacktest/Results/HF_VaR_VaR.csv')
EM_Var <- read.csv('./VaRBacktest/Results/EM_VaR_VaR.csv')

HF_Var <- HF_Var[1:1469, ]
EM_Var <- EM_Var[1:1469, ]

# Check that EM and HF have the same missing data
all(which(complete.cases(HF_Var)) == which(complete.cases(EM_Var)))

HF_Var_Clean <- HF_Var[complete.cases(HF_Var), ]
EM_Var_Clean <- EM_Var[complete.cases(EM_Var), ]
portNetRetClean <- portNetRetTest[complete.cases(HF_Var)]

Var_dates <- index(FTSE)[1001:2469]
dates <- Var_dates[complete.cases(HF_Var)]

# Compute statistics for the Basel traffic light test

HF_exceedances <- matrix(nrow = nrow(HF_Var), ncol = ncol(HF_Var))
EM_exceedances <- matrix(nrow = nrow(EM_Var), ncol = ncol(EM_Var))
HF_exceedances_clean <- matrix(nrow = nrow(HF_Var_Clean), ncol = ncol(HF_Var_Clean))
EM_exceedances_clean <- matrix(nrow = nrow(EM_Var_Clean), ncol = ncol(EM_Var_Clean))
for (i in 1:6) {
  HF_exceedances[ , i] <- (portNetRetTest < HF_Var[ , i])
  EM_exceedances[ , i] <- (portNetRetTest < EM_Var[ , i])
  HF_exceedances_clean[ , i] <- (portNetRetClean < HF_Var_Clean[ , i])
  EM_exceedances_clean[ , i] <- (portNetRetClean < EM_Var_Clean[ , i])
}

HF_TLTest <- rollapply(HF_exceedances, FUN = sum, width = 250, by.column = T)
EM_TLTest <- rollapply(EM_exceedances, FUN = sum, width = 250, by.column = T)
HF_TLTest_clean <- rollapply(HF_exceedances_clean, FUN = sum, width = 250, by.column = T)
EM_TLTest_clean <- rollapply(EM_exceedances_clean, FUN = sum, width = 250, by.column = T)

families_nice <- c('Clayton', 'Gumbel', 'Gaussian')
legend_position <- c('topleft', 'topleft', 'top')

tl_fig_filename <- paste0(fig_dir, '/tl_test_fig_99.png')
png(tl_fig_filename, height = 1400, width = 1200, pointsize = 16)
par(mfrow = c(3, 1), cex = 1.1)
for (i in 1:3) {
  main <- paste0(families_nice[i], ' copula')
  plot(HF_TLTest_clean[ , 2*i - 1] ~ dates, type = 'l', col = 'red', ylim = c(0, max(c(HF_TLTest_clean[ , 2*i - 1]))),
       xlab = 'Date', ylab = 'Failures', main = main)
  lines(EM_TLTest_clean[ , 2*i - 1] ~ dates, col = 'blue')
  legend(legend_position[i], fill = c('red', 'blue'), legend = c('IFM', 'EM'), box.lty = 'blank', inset = 0.05)
}
mtext(text = 'Failures of 99% copula-GARCH Value-at-Risk models over a rolling window of 250 days', side = 3, outer = TRUE, at = 0.5, line = -1.2, cex = 1.6)
dev.off()

# Create table detailing how often each model is in each zone

HF_model_zones <- matrix(nrow = nrow(HF_TLTest_clean), ncol = ncol(HF_TLTest_clean))
EM_model_zones <- matrix(nrow = nrow(EM_TLTest_clean), ncol = ncol(EM_TLTest_clean))

# Specify the zone for 99% VAR
HF_model_zones[ , c(1, 3, 5)][HF_TLTest_clean[ , c(1, 3, 5)] <= 4] <- 1
EM_model_zones[ , c(1, 3, 5)][EM_TLTest_clean[ , c(1, 3, 5)] <= 4] <- 1
HF_model_zones[ , c(1, 3, 5)][(HF_TLTest_clean[ , c(1, 3, 5)] <= 9) & (HF_TLTest_clean[ , c(1, 3, 5)] > 4)] <- 2
EM_model_zones[ , c(1, 3, 5)][(EM_TLTest_clean[ , c(1, 3, 5)] <= 9) & (EM_TLTest_clean[ , c(1, 3, 5)] > 4)] <- 2
HF_model_zones[ , c(1, 3, 5)][HF_TLTest_clean[ , c(1, 3, 5)] > 9] <- 3
EM_model_zones[ , c(1, 3, 5)][EM_TLTest_clean[ , c(1, 3, 5)] > 9] <- 3

# Specify zones for 95% VAR
HF_model_zones[ , c(2, 4, 6)][HF_TLTest_clean[ , c(2, 4, 6)] <= 17] <- 1
EM_model_zones[ , c(2, 4, 6)][EM_TLTest_clean[ , c(2, 4, 6)] <= 17] <- 1
HF_model_zones[ , c(2, 4, 6)][(HF_TLTest_clean[ , c(2, 4, 6)] <= 26) & (HF_TLTest_clean[ , c(2, 4, 6)] > 17)] <- 2
EM_model_zones[ , c(2, 4, 6)][(EM_TLTest_clean[ , c(2, 4, 6)] <= 26) & (EM_TLTest_clean[ , c(2, 4, 6)] > 17)] <- 2
HF_model_zones[ , c(2, 4, 6)][HF_TLTest_clean[ , c(2, 4, 6)] > 26] <- 3
EM_model_zones[ , c(2, 4, 6)][EM_TLTest_clean[ , c(2, 4, 6)] > 26] <- 3

HF_zone_props <- list()
EM_zone_props <- list()
for(i in 1:6) {
  HF_zone_props[[i]] <- prop.table(table(HF_model_zones[ , i]))
  EM_zone_props[[i]] <- prop.table(table(EM_model_zones[ , i]))
}
zones <- c('Green', 'Yellow', 'Red')
indices <- c(2, 4, 6, 1, 3, 5)
cat(paste0('\\multirow{3}{*}{Zone} & \\multicolumn{6}{c}{95\\% Value-at-Risk} & \\multicolumn{6}{c}{99\\% Value-at-Risk} \\\\\n'))
cat(paste0(' & ', paste0(rep(paste0('\\multicolumn{2}{c}{', families_nice, '}'), 2), collapse = ' & '), '\\\\\n'))
cat(paste0(' & ', paste0(rep(c('\\textit{IFM}', '\\textit{EM}'), 6), collapse = ' & '), '\\\\\n'))
cat('\\hline\n')
for (i in 1:3) {
  cat(paste0(zones[i], ' & '))
  for (j in 1:6) {
    cat(paste0(formatC(HF_zone_props[[indices[j]]][i] * 100, digits = 3), ' & '))
    cat(paste0(formatC(EM_zone_props[[indices[j]]][i] * 100, digits = 3), ' & '))
  }
  cat('\\\\\n')
}

HF_failures <- colSums(HF_exceedances_clean)
EM_failures <- colSums(EM_exceedances_clean)
n <- 1469
alpha <- rep(c(0.99, 0.95), 3)

HF_zstat <- (HF_failures - n*(1 - alpha))/(sqrt(n * alpha * (1 - alpha)))
EM_zstat <- (EM_failures - n*(1 - alpha))/(sqrt(n * alpha * (1 - alpha)))
HF_pvalues <- 2 * (1 - pnorm(HF_zstat))
EM_pvalues <- 2 * (1 - pnorm(EM_zstat))

cat(paste0(' & ', paste0(paste0('\\multicolumn{2}{c|}{', families_nice, '}'), collapse = ' & '), '\\\\\n'))
cat(paste0(' & ', paste0(rep(c('\\textit{IFM}', '\\textit{EM}'), 3), collapse = ' & '), '\\\\\n'))

cat('No. of failures & ')
for(i in 1:3) {
  cat(paste0(formatC(HF_failures[indices[i]]), ' & '))
  cat(paste0(formatC(EM_failures[indices[i]]), ' & '))
}
cat('\\\\\n')
cat(paste0('Z-statistic & '))
for (i in 1:3) {
  cat(paste0(formatC(HF_zstat[indices[i]], digits = 4), ' & '))
  cat(paste0(formatC(EM_zstat[indices[i]], digits = 4), ' & '))
}
cat('\\\\\n')
cat('$p$-values & ')
for (i in 1:3) {
  cat(paste0(formatC(HF_pvalues[indices[i]], width = 6, digits = 4, format = 'f'), ' & '))
  cat(paste0(formatC(EM_pvalues[indices[i]], width = 6, digits = 4, format ='f'), ' & '))
}
cat('\n')
cat('\n')
cat(paste0(' & ', paste0(paste0('\\multicolumn{2}{c|}{', families_nice, '}'), collapse = ' & '), '\\\\\n'))
cat(paste0(' & ', paste0(rep(c('\\textit{IFM}', '\\textit{EM}'), 3), collapse = ' & '), '\\\\\n'))
cat('No. of failures & ')
for(i in 1:3) {
  cat(paste0(formatC(HF_failures[indices[i + 3]]), ' & '))
  cat(paste0(formatC(EM_failures[indices[i + 3]]), ' & '))
}
cat('\\\\\n')
cat(paste0('Z-statistic & '))
for (i in 1:3) {
  cat(paste0(formatC(HF_zstat[indices[i + 3]], digits = 4), ' & '))
  cat(paste0(formatC(EM_zstat[indices[i + 3]], digits = 4), ' & '))
}
cat('\\\\\n')
cat('$p$-values & ')
for (i in 1:3) {
  cat(paste0(formatC(HF_pvalues[indices[i + 3]], width = 6, digits = 4, format = 'f'), ' & '))
  cat(paste0(formatC(EM_pvalues[indices[i + 3]], width = 6, digits = 4, format ='f'), ' & '))
}

# Create time series plot of exceedances

exceedance_sample_filename <- paste0(fig_dir, '/exceedance_sample.png')
png(exceedance_sample_filename, height = 800, width = 1200, pointsize = 16)
par(cex = 1.15, cex.lab = 1.2)
plot(portNetRetClean[200:500] ~ dates[200:500], type = 'h',
     ylim = c(min(HF_Var_Clean[200:500, 1]), max(portNetRetClean[200:500])),
     xlab = 'Date', ylab = 'Log-Returns',
     main = '99% Value-at-Risk using the Clayton copula with exceedance points',
     bty = 'l')
lines(HF_Var_Clean[200:500, 1] ~ dates[200:500], col = colours[1], lwd = 2)
abline(h = 0)
exceedance_points <- which(HF_Var_Clean[200:500, 1] > portNetRetClean[200:500])
exceedance_dates <- dates[200:500][exceedance_points]
exceedance_vars <- HF_Var_Clean[200:500, 1][exceedance_points]
points(as.numeric(exceedance_dates), exceedance_vars, pch = 4, cex = 1.8, col = colours[2], lwd = 2)
dev.off()

# Conduct the test for Haas time between failures independence test
HF_violation_gaps <- list()
EM_violation_gaps <- list()
for (i in 1:6) {
  HF_violation_gaps[[i]] <- c(which(HF_exceedances_clean[ , i])[1], diff(which(HF_exceedances_clean[ , i])))
  EM_violation_gaps[[i]] <- c(which(EM_exceedances_clean[ , i])[1], diff(which(EM_exceedances_clean[ , i])))
}

HaasLRStat <- function(data, p) {
  powers <- data - 1
  num <- p * ((1 - p)^powers)
  denom <- (1/data) * ((1 - 1/data)^powers)
  LRStat <- -2 * sum(log(num/denom))
  return(LRStat)
}

level <- rep(c(0.99, 0.95), 3)

HF_Haas_stat <- numeric(6)
EM_Haas_stat <- numeric(6)
HF_Haas_pvalue <- numeric(6)
EM_Haas_pvalue <- numeric(6)
for (i in 1:6) {
  HF_Haas_stat[i] <- HaasLRStat(HF_violation_gaps[[i]], p = (1 - level[i]))
  EM_Haas_stat[i] <- HaasLRStat(EM_violation_gaps[[i]], p = (1 - level[i]))
  HF_Haas_pvalue[i] <- 1 - pchisq(HF_Haas_stat[i], df = length(HF_violation_gaps[[i]]))
  EM_Haas_pvalue[i] <- 1 - pchisq(EM_Haas_stat[i], df = length(EM_violation_gaps[[i]]))
}

cat(paste0(' & ', paste0(paste0('\\multicolumn{2}{c|}{', families_nice, '}'), collapse = ' & '), '\\\\\n'))
cat(paste0(' & ', paste0(rep(c('\\textit{IFM}', '\\textit{EM}'), 3), collapse = ' & '), '\\\\\n'))
cat('\\hline\n')
cat('LR statistic')
for (i in 4:6) {
  cat(paste0(' & ', formatC(HF_Haas_stat[indices[i]], width = 7, digits = 2, format = 'f')))
  cat(paste0(' & ', formatC(EM_Haas_stat[indices[i]], width = 7, digits = 2, format = 'f')))
}
cat('\\\\\n')
cat('$p$-value')
for (i in 4:6) {
  cat(paste0(' & ', formatC(HF_Haas_pvalue[indices[i]], width = 7, digits = 4, format = 'f')))
  cat(paste0(' & ', formatC(EM_Haas_pvalue[indices[i]], width = 7, digits = 4, format = 'f')))
}
