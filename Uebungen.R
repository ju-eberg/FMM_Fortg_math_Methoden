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



##################### BLATT 3

#3.1: goodnotest

# #3.2:
# # 
# # Um zu zeigen, dass die gegebene Matrix
# # 
# # A =
# #   | 1 1 |
# #   | 0 ε |
# #   
# #   numerisch singulär bezüglich δ ist, müssen wir zeigen, dass die Kondition der Matrix A, κ(A),
# multipliziert mit δ größer oder gleich 1 ist, wobei δ ≤ ε.
# # 
# # Die Kondition einer Matrix kann als Verhältnis der größten und kleinsten Singulärwerte 
# der Matrix ausgedrückt werden. Die Kondition κ(A) ist definiert als:
# #   
# #   κ(A) = σ_max / σ_min,
# # 
# # wobei σ_max der größte Singulärwert und σ_min der kleinste Singulärwert von A ist.
# # 
# # Wir können die Singulärwerte von A berechnen. Zuerst berechnen wir die Eigenwerte von AA, 
#   da die Singulärwerte von A die Quadratwurzeln der Eigenwerte von AA sind:
# #   
# #   A * A =
# #   | 1 1 | * | 1 0 | = | 2 1 |
# #   | 0 ε | | 1 ε | | ε ε² |
# #   
# #   Die Eigenwerte von AA sind die Lösungen der charakteristischen Gleichung det(AA - λI) = 0:
# #   
# #   det(A*A - λI) = | 2-λ 1 |
# #   | ε ε²-λ |
# #   
# #   Die Eigenwerte sind die Lösungen der Gleichung (2-λ)(ε²-λ) - ε = 0:
# #   
# #   (2-λ)(ε²-λ) - ε = 0
# # (2-λ)(λ² - ε²) - ε = 0
# # 2λ² - 2ε² - λ³ + λε² - ε = 0
# # 
# # Da ε > 0 und λ ist eine Lösung der Gleichung, und die anderen beiden Lösungen sind komplex, 
# können wir die größte und kleinste Eigenwerte bestimmen:
# #   
# #   Der größte Eigenwert (σ_max) ist der positive reale Wert von λ.
# # 
# # Der kleinste Eigenwert (σ_min) ist der nicht negative reale Wert von λ.
# # 
# # Die kleinste Eigenwert ist 0, und der größte Eigenwert ist der positive reale Wert. 
#   Daher ist σ_max / σ_min unendlich, da der kleinste Eigenwert 0 ist.
# # 
# # Nun, um die numerische Kondition von A zu bestimmen, teilen wir den größten 
# Singulärwert durch den kleinsten Singulärwert:
# #   
# #   κ(A) = σ_max / σ_min = ∞ / 0 = ∞.
# # 
# # Schließlich, da ε > 0 und δ ≤ ε, können wir sagen, dass κ(A) * δ ≥ ∞ * δ ≥ 1,
#   was bedeutet, dass A numerisch singulär bezüglich δ ist.



#3.3:

set.seed(5)
M <- matrix(rnorm(500 * 500), 500, 500)
A <- M %*% t(M)
B <- rnorm(500)
c(
  all(A == t(A)), # A ist symmetrisch
  min(eigen(A)$values) > 0 # A ist positiv definit
)

# a. 
# Schreibe eine Funktion spd_solve(A, B), welche die Cholesky-Faktorisierung von A
# verwendet um A−1B zu berechnen. Die Funktionen chol(), forwardsolve() und
# backsolve() sollten dabei nützlich sein.

spd_solve <- function(A, B) {
  # Cholesky-Faktorisierung von A
  L <- chol(A)
  
  # Löse das Gleichungssystem L^T * x = B
  x <- forwardsolve(t(L), B)
  
  # Löse das Gleichungssystem L * y = x
  y <- backsolve(L, x)
  
  return(y)
}

# In dieser Funktion wird zuerst die Cholesky-Faktorisierung von A durchgeführt, 
# indem chol(A) verwendet wird. Dann werden die Gleichungen L^T * x = B und L * y = x gelöst, 
# um den Vektor y zu erhalten, 
# der das Ergebnis von A⁻¹B darstellt.

# b. Vergleiche solve und spd_solve

# Erzeugen Sie eine zufällige SPD-Matrix A
set.seed(5)
M <- matrix(rnorm(500 * 500), 500, 500)
A <- M %*% t(M)

# Erzeugen Sie einen zufälligen Vektor B
B <- rnorm(500)

# Laufzeitmessung für solve(A, B)
time_solve <- system.time(solve(A, B))

# Laufzeitmessung für spd_solve(A, B)
time_spd_solve <- system.time(spd_solve(A, B))

# Ergebnisse ausgeben
print("Laufzeit solve(A, B):")
print(time_solve)

print("Laufzeit spd_solve(A, B):")
print(time_spd_solve)



### c.  . 
# Wie skalieren die Laufzeiten für k = 1, 10, 100? Wie lässt sich das Beobachtete erklären?
# 
# Die Laufzeiten für solve(A, B) und spd_solve(A, B) werden wahrscheinlich mit steigendem Wert von k skalieren, da die Größe der Matrix B sich ändert und sich dadurch die Anzahl der Operationen für die Berechnungen erhöht.
# 
# Hier ist eine allgemeine Erklärung:
#   
#   Für k = 1 (einen Vektor B):
#   
#   solve(A, B) löst das Gleichungssystem A⁻¹B für einen Vektor B.
# spd_solve(A, B) verwendet die Cholesky-Faktorisierung von A und führt zwei Vorwärts- und Rückwärtssubstitutionen durch. Die Cholesky-Faktorisierung hat eine geringere Laufzeitkomplexität als die allgemeine Matrixinversion.
# Für k = 10 (eine Matrix B mit 10 Spalten):
#   
#   Sowohl solve(A, B) als auch spd_solve(A, B) müssen eine Matrix B mit 10 Spalten behandeln. Die Laufzeiten beider Funktionen werden länger sein als bei k = 1, aber spd_solve wird wahrscheinlich schneller sein, da die Cholesky-Faktorisierung für A einmal durchgeführt wird, und dann wird die Vorwärts- und Rückwärtssubstitution für jede Spalte von B durchgeführt.
# Für k = 100 (eine Matrix B mit 100 Spalten):
#   
#   Mit steigendem k werden die Laufzeiten für beide Funktionen noch länger. solve(A, B) wird aufgrund der Matrixinversion insbesondere bei großen Matrizen erheblich länger dauern. spd_solve(A, B) wird ebenfalls länger dauern, aber die Cholesky-Faktorisierung wird immer noch schneller sein als die Matrixinversion, insbesondere bei großen k.
# Zusätzlich zur Größe von k hängt die Laufzeit auch von der Größe der Matrix A ab. Wenn A eine größere Matrix ist, werden die Berechnungen generell länger dauern. Die spezifischen Laufzeiten hängen auch von der Rechenleistung des verwendeten Computers und der Implementierung der verwendeten Funktionen ab.
# 
# Um genaue Aussagen über die Laufzeitkomplexität für verschiedene Werte von k und A zu treffen, wäre eine detaillierte Analyse der verwendeten Algorithmen und eine experimentelle Messung erforderlich. Es ist jedoch zu erwarten, dass die Cholesky-Faktorisierung insbesondere für große k und große A schneller ist als die Matrixinversion.
# 
# 


################################
# 4.1: in Goodnotes, hier Zusatz:

# Schreibe eine Funktion my_lu(A), welche die LU-Zerlegung von A berechnet.
A <- matrix(c(2,4,3,-4,-7,-5,6,8,2), nrow = 3, byrow = TRUE)
lu(A)

#4.2 in Goodnotes

# 4.3


# Lade die pracma-Bibliothek
library(pracma)

# Definiere die Matrix A mit einem kleinen Wert für ε
epsilon <- 1e-10
A <- matrix(c(epsilon, 1, 1, 1), nrow = 2, byrow = TRUE)

# Berechne die LU-Zerlegung ohne Pivoting
lu_result <- lu(A)

# Gib die LU-Zerlegung aus
print("Matrix L:")
print(lu_result$L)
print("Matrix U:")
print(lu_result$U)



