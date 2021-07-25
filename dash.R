# -------------------------------------------------------------------------------------------

# Animação da otimização IWLS  distribução Binomial

# -------------------------------------------------------------------------------------------

# Carrega Pacotes
library(tidyverse)
library(ggplot2)
library(gganimate)
library(hrbrthemes)
library(viridis)
library(gifski)
library(plotly)
library(transformr)
library(metR)

# -------------------------------------------------------------------------------------------

# Declara funções

# -------------------------------------------------------------------------------------------

# Declara função calcula W de W
calcula_W <- function(eta, theta){
  
  # Calcula W
  (1 / theta) * exp(2 * eta)
  
}

# Declara função Z
calcula_Z <- function(eta, Y){
  
  # Calcula Z
  eta + ( Y * exp(-eta) ) - 1 
  
}

# Funcao estimativa de betas
EMV <- function(Y, Xs, beta, iter, tol){
  
  # Declara betas
  betas <- tibble(iter = 0:iter,
                  beta0 = as.numeric(beta[1]),
                  beta1 = as.numeric(beta[2]))
  
  # Laco para popular betas calculados nas iterações
  for(r in 1:iter){
    
    # Recupera beta da vez
    beta <- t(betas[r, c(2,3) ] %>% as.matrix())
    
    # Calcula eta0
    eta <- Xs %*% beta
    
    # Calcula media inicial
    theta <- exp(eta)
    
    # Calcula diagonal da matriz Wi
    W_i <- calcula_W(eta, theta)
    
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
    if(norma < tol){
      
      # Resultado
      resultado <- betas[1:r,]
      
      # Iterrompe laço
      break
      
    }
    
  }
  
  # Retorno da função
  return(resultado)
  
}


# -------------------------------------------------------------------------------------------

# Otimização 

# -------------------------------------------------------------------------------------------

# Declara betas reais
betas_reais <- c(1, -0.5) 

# Declara X
X <- sort(rnorm(50, 0, 0.5))

# Declara matriz de Xs
Xs <- cbind(rep(1, length(X)), X)

# Calcula eta0
eta <- Xs %*% betas_reais

# Declara vetor Y
Y = exp(eta)

# Declara beta 0
beta0 <- c(0, 0)

# Executa EMV
sol <- EMV(Y, Xs, beta0, 200, 10e-5) 

# Visuaaliza beta estimado
sol %>%  tail()

# --------------------------------------------------------------------------------------------

# Espaço de busca

# --------------------------------------------------------------------------------------------

# Define data frame espaço de busca
espaco <- tibble(beta0 = rep(seq(-5,5,0.1),times = 101),
                 beta1 = rep(seq(-5,5,0.1), each = 101)) %>% 
  mutate(dif_betas = sqrt((beta0 - betas_reais[1])^2 + (beta1 - betas_reais[2])^2))

# Declara figura
plot_ly() %>% 
  add_trace(data = data,  x=espaco$beta0, y=espaco$beta1, z=espaco$dif_betas, type="mesh3d" ) 

# Declara grafico ggplot
espaco %>% 
  ggplot(mapping = aes(x = beta0, y = beta1, fill = dif_betas))+
  geom_tile()

# Declara grafico ggplot
espaco %>% 
  ggplot(mapping = aes(x = beta0, y = beta1, z = dif_betas))+
  geom_contour_filled()+
  geom_point(aes(x = betas_reais[1], y = betas_reais[2]), color = "red", size = 3)

# Declara data frame de animacao
anime <- tibble()

# Laco para construcao de dfs
for(i in 1:(nrow(sol)-1)){
  
  # Gera data frame de animacao
  df_anime <- tibble(iter = rep(i,10202),
                     beta0 = c(espaco$beta0,sol$beta0[i]),
                     beta1 = c(espaco$beta1,sol$beta1[i]),
                     dif_betas = c(espaco$dif_betas, 
                                   sqrt((sol$beta0[i] - betas_reais[1])^2 +
                                          (sol$beta1[i] - betas_reais[2])^2)),
                     Tipo = c(rep(10, 10201), 30))
  
  # Gera segunda iteracao
  anime <- anime %>% 
    bind_rows(df_anime)
  
}

# Declara grafico ggplot
p <- anime %>% 
  ggplot(mapping = aes(x = beta0,
                       y = beta1,
                       z = dif_betas,
                       color = dif_betas,
                       size = Tipo))+
  geom_point()+
  scale_color_gradient(low="blue", high="red")+
  transition_states(iter,
                    transition_length = 2,
                    state_length = 4,
                    wrap = FALSE) +
  view_follow(fixed_x = T, fixed_y = T)

# Salva animacao
animate(p, fps = 5, renderer = gifski_renderer(loop = F))

# --------------------------------------------------------------------------------------------

# Animação

# --------------------------------------------------------------------------------------------

# Declara data frame de animacao
anime <- tibble()

# Laco para construcao de dfs
for(i in 1:(nrow(sol)-1)){
  
  # Gera data frame de animacao
  df_anime <- tibble(iter = rep(i,10202),
                     beta0 = c(espaco$beta0,sol$beta0[i]),
                     beta1 = c(espaco$beta1,sol$beta1[i]),
                     dif_betas = c(espaco$dif_betas, 
                                   sqrt((sol$beta0[i] - betas_reais[1])^2 +
                                          (sol$beta1[i] - betas_reais[2])^2)),
                     Tipo = c(rep(0, 10201), 1))
  
  # Gera segunda iteracao
  anime <- anime %>% 
    bind_rows(df_anime)
  
}

sol

# Declara grafico ggplot
p <- ggplot(aes(x = beta0, y = beta1), data = sol)+
  geom_point(color = "red", size = 3)+
  geom_contour(data = espaco, mapping = aes(x = beta0, y = beta1, z = dif_betas) )+
  transition_time(iter) +
  view_follow(fixed_x = T, fixed_y = T)

# Salva animacao
animate(p, fps = 5, renderer = gifski_renderer(loop = F))

# --------------------------------------------------------------------------------------------

# Declara grafico ggplot
p <- ggplot(aes(x = beta0, y = beta1), data = sol)+
  geom_tile(show.legend = T,data = espaco, mapping = aes(x = beta0, y = beta1, fill = dif_betas) )+
  geom_contour(data = espaco, mapping = aes(x = beta0, y = beta1, z = dif_betas), color = "White" )+
  geom_text_contour(data = espaco,
                    aes(x = beta0, y = beta1, z = dif_betas),
                    color = "white",
                    size = 6,
                    min.size = 1,
                    skip = 0)+
  geom_point(color = "red", size = 4)+
  transition_time(iter) +
  view_follow(fixed_x = T, fixed_y = T)+
  labs(title = "Diferença de Betas Reais e Estimados",
       x = "Beta 0",
       y = "Beta 1",
       fill = "Diferença")

# Salva animacao
animate(p, fps = 5, renderer = gifski_renderer(loop = F))

# --------------------------------------------------------------------------------------------


# <div>Icons made by <a href="https://www.flaticon.com/authors/ultimatearm" title="ultimatearm">ultimatearm</a> from <a href="https://www.flaticon.com/" title="Flaticon">www.flaticon.com</a></div>