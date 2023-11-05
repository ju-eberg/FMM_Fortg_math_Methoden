x0 = 0
j = 2

f <- sin

t = \(x) x # äquiv. zu t= function (x) {return(x)}

xs = 10^seq(-8, 0, 0.1)

xs
#Liste, die bei 10^-8 startet, logarithmisch

errs = abs(f(xs) - t(xs))
errs # Liste der Fehler

length(xs)
length(errs) #jew. 81

plot(xs, errs, log="xy")
lines(xs, (xs-x0)^2)

# Frage mit k: Wie schnell wird mein Fehler kleiner?
# durchgezogene Linie hat Steigung 2
# schneller, weil Steigung steiler

# abs(f(x)-t(x)) = c (x-x0)^k
# log abs(( f(x)-t(x) ) = log(c) + k*log(x-x0)


#########################################################
# Blatt 1)
# 26.10.2023

# 1.1: ipad, rechnen

#1.2: 

x0 <- 1
f = function(x) x^2

# Funktion zur Berechnung der rechten Differenzenquotienten
ef_prime <- function(f, x0, h) {
  return((f(x0 + h) - f(x0)) / h)
}

# Wahre Ableitung von f(x) = x^2 an der Stelle x0 = 1
true_derivative <- 2

# Vektor von h-Werten (h wird kleiner)
h_values <- 10^seq(-1, -17, by = -0.1)

# Vektor zur Speicherung der Fehler
errors <- abs(ef_prime(function(x) x^2, 1, h_values) - true_derivative)

# Plot des Fehlers auf logarithmischen Achsen
plot(log10(h_values), log10(errors), type = "l", xlab = "log10(h)", ylab = "log10(Error)")

#### Ende 1.2


# 1.3:

# Mitternachtsformel
mitternacht <- function(p, q) {
  Diskri <- p^2 - q
  if (Diskri < 0) {
    return(c(NA, NA))
  }
  zwei_summ <- sqrt(Diskri)
  x1 <- -p + zwei_summ
  x2 <- -p - zwei_summ
  return(c(x1, x2))
}

# Satz von Vieta
Vieta <- function(p, q) {
  x2 <- -p
  x1 <- q / x2
  return(c(x1, x2))
}

# Testwerte für p und q
Werte_p <- c(1, 0.1, 0.01, 0.001, 0.0001)
Werte_q <- c(1, 0.1, 0.01, 0.001, 0.0001)

# Vergleiche die Ergebnisse der beiden Ansätze für verschiedene p und q
for (p in Werte_p) {
  for (q in Werte_q) {
    cat("p =", p, "q =", q, "\n")
    cat("Mitternachtsformel: ", mitternacht(p, q), "\n")
    cat("Satz von Vieta: ", mitternacht(p, q), "\n")
    cat("\n")
  }
}



####### Ende 1.3

######## BLATT 2

#2.1 und 2.2 nur in Lösung auf GoodNotest

##2.3:

#' @param X Matrix mit Zeilen X_1, ..., X_n.
#' @param beta Parametervektor

sigma_2 <- function(X, beta) {
  n <- nrow(X)
  d <- ncol(X)
  Sig <- matrix(0, d, d)
  for (i in 1:n) {
    Sig <- Sig + (X[i, ] - colMeans(X)) %*% t(X[i, ] - colMeans(X))
  }
  as.numeric(t(beta) %*% Sig %*% beta / n)
}

# 
# (a) Gebe die Laufzeit- und Speicher-Komplexität des Algorithmus in Landau-Notation an.
# 
# Lösung
# 
# Speicherplatz für X: O(nd)
# 
# Speicherplatz für beta: O(d)
# 
# Speicherplatz für Sig: O(d ^2)
# 
# Insgesamt haben wir also O(nd + d^2)
# 
# 
# Laufzeit für colMeans(X): O(nd + d) = O(nd) (n − 1 Additionen für jede der d Spalten
#                                              und d Multiplikationen für die Durchschnittsbildung)
# 
#  
# Laufzeit für (X[i, ] - colMeans(X)) %*% t(X[i, ] - colMeans(X)): O(nd + d + d^2) = O(nd + d^2)
#  (O(nd) Operationen für colMeans(X), 2d Subtraktionen, d^2 Multiplikationen für das äußere Produkt$)
# 
# • n Iterationen über i.
# 
# Insgesamt haben wir also O(n(nd + d^2)) = =(nd^2 + n^2d)
# 

#  # b)  
# Überprüfe die Laufzeiteigenschaften numerisch. Erzeuge dazu mehrere Datensätze
# mit n ∈ [10, 1000] und d ∈ [10, 1000] und benchmarke den Algorithmus (z.B. mit
#                                         microbenchmark::microbenchmark())

library(microbenchmark)
library(ggplot2)

grid <- expand.grid(
  n = seq.int(10, 1000, length = 4),
  d = seq.int(10, 1000, length = 4)
)

times <- numeric(nrow(grid))

for (k in 1:nrow(grid)) {
  n <- grid[k, 1]
  d <- grid[k, 2]
  X <- replicate(d, rnorm(n))
  cat("n = ", n, ", d = ", d, ", time = ", sep = "")
  times[k] <- median(microbenchmark(sigma_2(X, rep(1, d)), times = 10)$time)
  cat(times[k] / 1e9, "\n")
}

cbind(grid, time = times) |>
  ggplot(aes(n, time, color = as.factor(d))) +
  geom_line()

cbind(grid, time = times) |>
  ggplot(aes(d, time, color = as.factor(n))) +
  geom_line()


#######
# (c) Finde und benchmarke einen alternativen Algorithmus mit Laufzeit- und Speicherkomplexität O(nd)

# Kernideen in Lsg nachlesen

#' @param X Matrix mit Zeilen X_1, ..., X_n.
#' @param beta Parametervektor

#' @param X Matrix mit Zeilen X_1, ..., X_n.
#' @param beta Parametervektor
sigma_2_besser <- function(X, beta) {
  Xbeta <- X %*% beta
  mu <- mean(Xbeta)
  mean((Xbeta - mu)^2)
}

grid <- expand.grid(
  n = seq.int(10, 1000, length = 4),
  d = seq.int(10, 1000, length = 4)
)

times <- numeric(nrow(grid))

for (k in 1:nrow(grid)) {
  n <- grid[k, 1]
  d <- grid[k, 2]
  X <- replicate(d, rnorm(n))
  cat("n = ", n, ", d = ", d, ", time = ", sep = "")
  times[k] <- median(microbenchmark(sigma_2_besser(X, rep(1, d)), times = 1000)$time)
  cat(times[k] / 1e9, "\n")
}

cbind(grid, time = times) |>
  ggplot(aes(n, time, color = as.factor(d))) +
  geom_line()

cbind(grid, time = times) |>
  ggplot(aes(d, time, color = as.factor(n))) +
  geom_line()

# Wir sehen also, dass die Laufzeit von d und n hier tatsächlich linear sind.
