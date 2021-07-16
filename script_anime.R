# -------------------------------------------------------------------------------------------

# App Script

# -------------------------------------------------------------------------------------------

# libraries:
library(tibble)
library(ggplot2)
library(gganimate)
library(babynames)
library(hrbrthemes)
library(viridis)
library(gifski)

# -------------------------------------------------------------------------------------------


# Declara Y
Y = sort(rpois(25, 15), decreasing = T)

# Declara X
X = sort(rnorm(25, 0, 0.5))

# Declara matriz de Xs
Xs = cbind(rep(1, length(X)), X)

# Declara beta 0
beta0 = c(1,1)

# Executa EMV
sol = EMV(Y, Xs, beta0, 40, 10e-5) 

# Define estados
estados = paste0(rep(sol$iter, each = length(Y))," ",
                   "Beta0: ", rep(round(sol$beta0,4), each = length(Y))," ",
                   "Beta1: ", rep(round(sol$beta1,4), each = length(Y)))

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
  df_anime[i, "pred"] = exp(df_anime$beta0[i] + df_anime$beta1[i] * df_anime$X[i])  
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
anime <- dados %>%
  ggplot( aes(x=X, y=pred, group = Y, color = Y)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text(aes(label = paste0("Betas = (",round(beta0,3)," , ", round(beta1,3), ")")),
            y = min(original$Observado),
            x = -0.45,
            check_overlap = TRUE)+
  labs(title = "Optimização EMV IWLS Poisson",
       subtitle = "Iterações: {closest_state}",
       y = "Y")+
  transition_states(estados,
                    transition_length = 4,
                    state_length = 2,
                    wrap = FALSE) +
  view_follow()
  
  

animate(anime, renderer = gifski_renderer(loop = F))


# -------------------------------------------------------------------------------------------


