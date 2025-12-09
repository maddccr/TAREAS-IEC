
#Primero suavizamos con metodo LOESS
t <- ABM_1$df_N$t
ajuste <- loess(ABM_1$N ~ t, span = 0.2)
ABM_1_suavizado <- predict(ajuste)

#Plot pasico para registrar resultados
plot(ABM_1_suavizado, type="p", 
     xlab="Tiempo", ylab="Población N(t)",
     main="ABM Base")





mean_interp <- apply(sim_splines, 2, mean)
var_interp  <- apply(sim_splines, 2, var)
se_interp   <- sqrt(var_interp / nSim)
mean(se_interp)

#Banda de confianza 95%
ribbon_df <- data.frame(
  t    = t200,
  ymin = mean_interp - 1.96 * se_interp,
  ymax = mean_interp + 1.96 * se_interp
)













(plot) {
  
  plot_grid <- rep(0, L)
  plot_grid[which(habitat == 1 & bio_sex == 0)] <- 1
  plot_grid[which(habitat == 1 & bio_sex == 1)] <- 2
  
  
  if(t %in% seq(1, iter - iter/2, by = 5)){
    
    habCuadrado <- matrix(plot_grid, nrow=habSize, byrow=TRUE)
    image(habCuadrado, col = c("white", "pink", "purple"), axes = FALSE, asp = 1)
    #title(main = paste("Iteración", t, " - N =", N[length(N)]))
  }
  
}


sensibilidad <- function(K = 500, iter = 60, p_birth = 0.9, reps = 50){
  t_inicial <- Sys.time()
  p_death <- seq(0, 0.5, by = 0.05)
  p_len <- length(p_death) #cauntos valores distintos de p
  
  #Monitorear
  extincion <- numeric(p_len)
  t_hasta_ext_prom <- numeric(p_len)
  
  #PARA LUEGO GRAFICAR NECESITAMOS PROMEDIAR LAS REPS
  mean_N <- matrix(0, nrow = iter, ncol = p_len)
  var_N <- matrix(0, nrow = iter, ncol = p_len)
  
  
  for(i in 1:p_len) {
    
    pd <- p_death[i]
    #guardar la serie de cada corrida
    runs <- matrix(NA, nrow = reps, ncol = iter)
    
    t_ext <- rep(Inf, reps) #si no extincion t = Inf
    
    for(r in 1:reps){
      N_t <- ABM_ext3(K, iter, pd, p_birth)
      runs[r, ] <- N_t
      
      if(any(N_t == 0)) { #primer instante donde la poblacion = 0
        t_ext[r] <- which(N_t == 0)[1]
      }
    }
    
    extincion[i] <- mean(t_ext < Inf)
    t_hasta_ext_prom[i] <- ifelse(all(extincion[i] == 0), 
                                  NA, 
                                  mean(t_ext[t_ext < Inf]))
    
    mean_N[ ,i] <- colMeans(runs)
    var_N[ ,i] <- apply(runs, 2, var)
    
  }
  
  list(p_death = p_death,
       mean_N = mean_N,
       var_N = var_N,
       extincion_prob = extincion,
       t_hasta_ext_prom = t_hasta_ext_prom,
       t_ej = (t_inicial - Sys.time()))
  
}

res <- sensibilidad()


res$t_ej




pd <- res$p_death #cols x
t <- 1:nrow(res$mean_N) #filas y
z <- res$mean_N

persp(y = pd, x = t, z = z, theta = 0, phi = -9, ticktype = "detailed", expand = 0.5)

image(t, pd, z, xlab="p_death", ylab="t", main="Mean N(t) heatmap")
contour(t, pd, z, add = TRUE, drawlabels = FALSE)

library(plotly)


pd <- res$p_death                 # x = columnas
t  <- 1:nrow(res$mean_N)          # y = filas
z  <- res$mean_N                  # matriz filas × columnas

plot_ly() |>
  add_surface(
    x = pd,     # eje X
    y = t,      # eje Y
    z = z       # matriz Z
  ) |>
  layout(
    scene = list(
      xaxis = list(title = "p_death"),
      yaxis = list(title = "t"),
      zaxis = list(title = "N promedio")
    )
  )
























dN/dt = rN[(K − N)/K]

r = b-d o r = (b − d) + (i − e)

#arranca optima // abundancia de recursos 

N <- ABM(100, 60, 0.3, 0.8)

plot(N$N, type="p", 
     
     xlab="Tiempo", ylab="Población N(t)",
     
     main="ABM - Crecimiento Logístico")

points((N$p_birth_eff)*100, col= "red")

N$p_birth_eff

plot(ABM_1$N, type="p", 
     xlab="Tiempo", ylab="Población N(t)",
     main="ABM - Crecimiento Logístico")
#points((N$p_birth_eff)*100, col= "red")



points(sim[2,], col="red")
points(sim[3,], col="yellow")
points(sim[4,], col="green")
points(sim[5,], col="violet")
points(sim[6,], col="blue")
points(sim[7,], col="red")
points(sim[8,], col="yellow")
points(sim[9,], col="green")
points(sim[10,], col="violet")

#grafico de interpolados
plot(sim_splines[1,], col="blue", pch=18)
points(sim_splines[2,], col="red", pch=18)
points(sim_splines[3,], col="yellow",pch=18)
points(sim_splines[4,], col="green", pch=18)
points(sim_splines[5,], col="violet", pch=18)
points(sim_splines[6,], col="blue", pch=18)
points(sim_splines[7,], col="red",pch=18)
points(sim_splines[8,], col="yellow", pch=18)
points(sim_splines[9,], col="green", pch=18)
points(sim_splines[10,], col="violet", pch=18)


#Ejemplo grafico de una simulacion 
plot(sim[1,], col="blue")
t200 <- seq(1, 60, length.out = 200)
sim_splines <- interpolar_sims(sim)
points(t200, sim_splines[1,], col="blue", pch=18, size = 0.3)


















#EXTENDER EL MODELO AGREGANDO SEXO BIOLOGICO
ABM_ext1 <- function(K, iter, p_death, p_birth) {
  
  #Inicializar el ABM
  #Creacion del ambiente
  
  habSize <- ceiling(sqrt(K)) #quemos una grilla LxL que pueda ocupar K individuos
  L <- habSize^2
  habitat <- rep(0, L) #habitad vacia
  
  #Cantidad inicial de agentes (1%)
  N0 <- round(0.05 * K) #incorporar como parametro? 
  
  #PARA EVITAR EXTINCION SUBITA > 5%
  
  #Posicionar de forma aleatoria
  habitat[sample(1:L, N0)] <- 1 #1 = lleno
  
  #Lo que se busca monitonear es la cantidad de agentes en el tiempo
  N <- c() #crecieminto dinamico
  
  #INICIALIZAR LA EXTENCION
  bio_sex <- rep(NA, L)
  bio_sex[which(habitat ==1)] <- rbinom(N0, 1, 0.5)
  #a cada agente inicial le asigno un sexo bio (1= fem 0 = masc) con 50% de prob
  
  #ahora si
  
  for(t in 1:iter){
    
    vivos <- which(habitat == 1) #recuperar los ind donde el hab esta ocupado
    
    #MUERTES: a cada agente vivo se le asigna una probabilidad a partir de una unif(0,1). Si p < p_death --> r.i.p
    
    for(agente in vivos){
      if( runif(1) < p_death){
        habitat[agente] <- 0
        bio_sex[agente] <- NA
      }
    }
    
    #NACIMIENTOS 
    #para simular el compartamiento logistico es necesario que la p_birth sea controlada por la carga del ambiente
    
    n <- sum(habitat == 1) #sumamos los vivos
    p_birth_K <- p_birth*(1 - n / K) #por la limitaciones de la carga, a medida que aumenta N, disminuye la prob de nacer
    
    
    for(agent in vivos){
      if(!is.na(bio_sex[agent]) & bio_sex[agent] == 1) {
        if (runif(1) < p_birth_K){
          libres <- which(habitat == 0) #ind de espacios libres en la grilla
          if(length(libres) > 0){ 
            i <- sample(libres, 1)
            habitat[i]  <- 1 #elijo un ind aleatorio para colocar
            bio_sex[i] <- rbinom(1, 1, 0.5)
          }
        }
      }
    }
    
    #Actualizar la cantidad de agentes
    
    N <- c(N, sum(habitat == 1))
    
    plot_grid <- rep(0, L)
    plot_grid[which(habitat == 1 & bio_sex == 0)] <- 1
    plot_grid[which(habitat == 1 & bio_sex == 1)] <- 2
    
    
    if(t %in% seq(1, iter - iter/2, by = 5)){
      
      habCuadrado <- matrix(plot_grid, nrow=habSize, byrow=TRUE)
      image(habCuadrado, col = c("white", "pink", "purple"), axes = FALSE, asp = 1)
      #title(main = paste("Iteración", t, " - N =", N[length(N)]))
    }
    
  }
  
  #intento de vizualizacion
  
  
  
  return(N)
  
}


N <- ABM_ext(1000, 60, 0.3, 0.9)


plot(N, type="l", 
     xlab="Tiempo", ylab="Población N(t)",
     main="ABM - Crecimiento Logístico")


#1000 03 - 09 extincion


#curvas por curva// parametros --> sacar los coeficientes 
#en promedio coso // intervalito // no correr sino usar muchas reps 