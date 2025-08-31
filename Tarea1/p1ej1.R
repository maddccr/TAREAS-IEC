# Función que genera densidades de una normal
dist_normal <- function(mu = 0,sigma = 1,delta){
  x <- seq(-5, 5, by = delta) # valores de x de -5 a 5
  probs <- (1/(sqrt(2*pi)*sigma)) * exp(-((x - mu)^2)/(2*sigma^2)) 
  return(list(x = x, probs = probs)) # devuelve lista con x y densidades
}

# Probando diferentes valores delta para ver el efecto en la curva
delta <- c(0.5, 0.1, 0.01)
for(d in delta){
  dt <- dist_normal(0,1,d)
  plot(dt$x, dt$probs, col = "red", xlab = "x", ylab = "Densidad",
       cex = d, main = "Densidad Distribución Normal típica",
       sub = paste0("(con distancia entre valores = ", d, ")"))
}

# Preparar ventana para 3 gráficos juntos
par(mfrow = c(1,3)) 

# Diferentes combinaciones de media y desviación
parametros <- list(
  c(mu = 0, sigma = 1),
  c(mu = 5, sigma = 1),
  c(mu = 5, sigma = 4)
)

for(p in parametros){
  mu <- p["mu"]
  sigma <- p["sigma"]
  dt <- dist_normal(mu, sigma, delta = 0.01)
  plot(dt$x, dt$probs, col = "red", xlab = "x", ylab = "Densidad", type = "l",
       lwd = 2, ylim = c(0, 0.4), main = paste0("Normal (μ=", mu, ",σ=", sigma, ")"))
}

