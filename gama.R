# -------------------------------------------------------------------------------------------

# Animação da otimização IWLS distribução Gama

# -------------------------------------------------------------------------------------------

# libraries:
library(tibble)
library(ggplot2)
library(gganimate)
library(hrbrthemes)
library(viridis)
library(gifski)

# -------------------------------------------------------------------------------------------

# Declara funções

# -------------------------------------------------------------------------------------------

# Declara função calcula W de W
calcula_W <- function(eta, theta, alfa){
  
  # Calcula W
  #(1 / theta) * exp(2 * eta)
  (alfa/theta^2) * (-1/eta^2)^2
  
}

# Declara função Z
calcula_Z <- function(eta, Y){
  
  # Calcula Z
  #eta + ( Y * exp(-eta) ) - 1 
  2*eta - (Y*(eta^2))
  
}

# Funcao estimativa de betas1
EMV <- function(Y, Xs, beta, iter, tol, alfa){
  
  # Declara betas
  betas <- tibble(iter = 0:iter,
                  beta0 = as.numeric(beta[1]),
                  beta1 = as.numeric(beta[2]))
  
  for(r in 1:iter){
    
    beta <- t(betas[r, c(2,3) ] %>% as.matrix())
    
    # Calcula eta0
    eta <- Xs %*% beta
    
    # Calcula media inicial
    theta <- (1 / eta)
    
    # Calcula diagonal da matriz W
    W_i <- calcula_W(eta, theta, alfa)
    
    # Declara matriz W
    W <- diag(length(Y))
    diag(W) <- W_i
    
    # Calcula z
    z_i <- calcula_Z(eta, Y)
    
    # Calcula beta_0
    betas[r+1, "beta0" ] <- (solve( t(Xs) %*% W %*% Xs ) %*% t(Xs) %*% W %*% z_i  )[1]
    
    # Calcula beta_1
    betas[r+1, "beta1" ] <- (solve( t(Xs) %*% W %*% Xs ) %*% t(Xs) %*% W %*% z_i  )[2]
    
    # Calcula norma euclideana
    norma <- sqrt(sum((abs(betas[r+1, c(2,3)] - betas[r, c(2,3)]))^2))
    
    # Condição de parada
    
    if(norma < 10^(-tol)){
      
      # Resultado
      resultado <- betas[1:r,]
      
      # Iterrompe laço
      break
      
    }
    
  }
  
  # Retorno da função
  print(betas)
  return(resultado)
  
}

# -------------------------------------------------------------------------------------------

# Otimização Gamma

# -------------------------------------------------------------------------------------------

# Declara Y
#Y = sort(rgamma(100, shape = 2, rate = 2/5), decreasing = T)

# Declara X e Y
alpha = 5
b = 2
X = sort(rgamma(100, shape = alpha, rate = b), decreasing = T)
Y = -alpha/(0.5 - X*5)

# Declara matriz de Xs
Xs = cbind(rep(1, length(X)), X)

# Declara beta 0
beta0 = c(00000.1, 00000.1)

# Executa EMV
sol = EMV(Y, Xs, beta0, 100, 10, alpha)

# Define estados
estados = paste0(rep(sol$iter, each = length(Y)),"    ",
                 "Beta 0:   ", rep(round(sol$beta0,4), each = length(Y)),"    ",
                 "Beta 1:   ", rep(round(sol$beta1,4), each = length(Y)))

# Gera data frame
df_anime <- tibble(iter = rep(sol$iter, each = length(Y)),
                   estados = factor(estados, levels = unique(estados)),
                   beta0 = rep(sol$beta0, each = length(Y)),
                   beta1 = rep(sol$beta1, each = length(Y)),
                   X = rep(X, times = nrow(sol)),
                   Observado = rep(Y, times = nrow(sol)),
                   pred = 0.1,
                   Y = "Estimado")

# Laco para popular betas
for(i in 1:nrow(df_anime)){
  df_anime[i, "pred"] = 1/(df_anime$beta0[i] + df_anime$beta1[i] * df_anime$X[i])
}

# Data frame original
original <- tibble(iter = rep(df_anime$iter[nrow(df_anime)], length(Y)),
                   estados = rep(df_anime$estados[nrow(df_anime)], length(Y)),
                   beta0 = rep(df_anime$beta0[nrow(df_anime)], length(Y)),
                   beta1 = rep(df_anime$beta1[nrow(df_anime)], length(Y)),
                   X = X,
                   Observado = Y,
                   pred = Y,
                   Y = "Observado")

# Declara Data
dados <- df_anime %>% 
  bind_rows(original)

# Plot
p <- dados %>%
  ggplot( aes(x=X, y=pred, group = Y, color = Y)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "Optimização EMV IWLS Sigma",
       subtitle = paste("Iteração:  ","{closest_state}"),
       y = "Y")+
  transition_states(estados,
                    transition_length = 2,
                    state_length = 4,
                    wrap = FALSE) +
  view_follow(fixed_x = T)

# Salva animacao
animate(p, fps = 10, renderer = gifski_renderer(loop = F))


# -------------------------------------------------------------------------------------------