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

####

