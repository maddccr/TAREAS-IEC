
#----
modelo logistico:
  dN/dt = rN[(K − N)/K]
r = b-d o r = (b − d) + (i − e)
N=individuos
K= capacidad de carga (en numero de individuos)
segun chat, r no iria en la ecuacion, se representaria con probs de nacimiento
y muerte por agente y puede variar segun la densidad de poblaicon

ABMLogistico <- function(K, years, p_birth, p_death) {
  #k capacidad de carga en individuos, anios iteraciones, p_birth, p_death prob de nacimiento y de muerte
  
  #tamanio del habitat (cada celda 0=vacía, 1=agente)
  habSize <- ceiling(sqrt(K))
  T <- habSize^2 #tamanio de cuadricula
  
  # Estado inicial: pocos agentes (1%)
  N0 <- round(0.01 * K) #numero inicial de individuos
  habitat <- rep(0, T) #habitat vacio
  habitat[sample(1:T, N0)] <- 1 #coloca los individuos en lugares aleatorios
  
  Nmonitor <- c() #vector que monitorea el numero de individuos
  
  for (t in 1:years) {
    
    living <- which(habitat == 1) #vector con indices de celdas con individuos vivos
    
    # MUERTES: a cada agente vivo le asigna una probabilidad a partir de una unif, si es menor a la p_death lo mata
    for (agent in living) {
      if (runif(1) < p_death) {
        habitat[agent] <- 0
      }
    }
    
    # NACIMIENTOS
    N <- sum(habitat == 1)  #N es el numero de individuos
    p_birth_eff <- p_birth * (1 - N / K)   # actualiza la prob de nacer segun numero de individuos, a medida que aumenta N, disminuye la prob de nacer
    
    for (agent in living) {
      if (runif(1) < p_birth_eff) { #genero num aleatorio, si es menor a la prob de nacer, nace un individuo
        empties <- which(habitat == 0) #vector con indices de celdas vacias
        if (length(empties) > 0) { #si hay lugar
          newPos <- sample(empties, 1) #creo nueva posicion tomando un valor aleatorio de los indices
          habitat[newPos] <- 1 #coloco al individuo en esa posicion
        }
      }
    }
    
    Nmonitor <- c(Nmonitor, sum(habitat == 1)) #actualizo el monitoreo de individuos
  }
  
  plot(Nmonitor, type="l", 
       xlab="Tiempo", ylab="Población N(t)",
       main="ABM - Crecimiento Logístico")
  
  return(Nmonitor)
}


  