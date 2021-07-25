
# --------------------------------------------------------------------------

# Aplicativo Shiny Ajuste da Curva

# --------------------------------------------------------------------------

# Carrega Bibliotecas
library(shiny)
library(tidyverse)
library(tibble)
library(ggplot2)
library(gganimate)
library(hrbrthemes)
library(viridis)
library(gifski)

# --------------------------------------------------------------------------

# Define funções

# --------------------------------------------------------------------------

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
    if(norma < 10^(-tol)){
      
      # Resultado
      resultado <- betas[1:r,]
      
      # Iterrompe laço
      break
      
    }
    
  }
  
  # Retorno da função
  return(resultado)
  
}

# --------------------------------------------------------------------------

# Interface do Usuário

# --------------------------------------------------------------------------

# Declara interface do usuario
ui <- fluidPage(
  
  # Titulo do aplicativo
  titlePanel("Estimador IWLS Poisson"),
  
  # Barra lateral de entradas do usuário
  sidebarLayout(
    sidebarPanel(
      h3("Amostra"),
      sliderInput("obs",
                  "Número de Observações:",
                  min = 10,
                  max = 100,
                  value = 25),
      sliderInput("beta0_real",
                  "Beta 0 Real:",
                  min = 0,
                  max = 5,
                  value = 2),
      sliderInput("beta1_real",
                  "Beta 1 Real",
                  min = -3,
                  max = 3,
                  value = -2),
      actionButton(inputId = "amostre", 
                   label = "Amostrar"),
      h3("Optimização"),
      sliderInput("beta0",
                  "Beta 0 Inicial",
                  min = -5,
                  max = 5,
                  value = 1),
      sliderInput("beta1",
                  "Beta 1 Inicial",
                  min = -5,
                  max = 5,
                  value = 1),
      sliderInput("ints",
                  "Máximo de Iterações",
                  min = 1,
                  max = 200,
                  value = 100),
      sliderInput("tol",
                  "Tolerância (10^-x)",
                  min = 1,
                  max = 10,
                  value = 5),
      actionButton(inputId = "optimize", 
                   label = "Optimizar")
    ),
    
    # Desenha grafico
    mainPanel(
      imageOutput("grafico1"),
      imageOutput("grafico2")
    )
  )
)

# --------------------------------------------------------------------------

# Servidor

# --------------------------------------------------------------------------

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Observa clique do botao amostrar
  dados <- eventReactive(input$amostre, {
    
    # Declara betas reais
    betas_reais <- c(input$beta0_real, input$beta1_real) 
    
    # Declara X
    X <- sort(rnorm(input$obs, 0, 0.5))
    
    # Declara matriz de Xs
    Xs <- cbind(rep(1, length(X)), X)
    
    # Calcula eta0
    eta <- Xs %*% betas_reais
    
    # Declara vetor Y
    Y = exp(eta)
    
    # Declara data frame
    amostra <- tibble(Y, X)
    
  })
  
  # Observa clique do botao optimizar
  animacao1 <- eventReactive(input$optimize, {
    
    # Declara Y
    Y = dados()$Y
    
    # Declara X
    X = dados()$X
    
    # Declara matriz de Xs
    Xs = cbind(rep(1, length(X)), X)
    
    # Declara beta 0
    beta0 = c(input$beta0, input$beta1)
    
    # Executa EMV
    sol = EMV(Y, Xs, beta0, input$ints, input$tol) 
    
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
    
    # Retorno
    return(dados)
    
  })
  
  # Renderiza Grafico
  output$grafico1 <- renderImage({
    
    # Declara objeto gráfico
    p <- animacao1() %>%
      ggplot( aes(x=X, y=pred, group = Y, color = Y, shape = Y)) +
      geom_point(size = 5) +
      scale_shape_manual(values=c(1, 20))+
      scale_colour_manual(values = c("red", "#376db2"))+
      labs(title = "Optimização EMV IWLS Poisson",
           subtitle = paste("Iteração:  ","{closest_state}"),
           y = "Y")+
      transition_states(estados,
                        transition_length = 2,
                        state_length = 4,
                        wrap = FALSE) +
      view_follow(fixed_x = T)+
      theme(panel.grid = element_blank(),
            panel.background = element_blank())
    
    # Salva animacao
    anim_save("outfile.gif", animate(p, fps = 5, renderer = gifski_renderer(loop = F))) 
    
    # Retorna lista contendo o arquivo
    list(src = "outfile.gif",
         contentType = 'image/gif',
         width = 600,
         height = 400)
    
  }, deleteFile = TRUE)
  
  # Observa clique do botao optimizar
  animacao2 <- eventReactive(input$optimize, {
    
    # Declara Y
    Y = dados()$Y
    
    # Declara X
    X = dados()$X
    
    # Declara matriz de Xs
    Xs = cbind(rep(1, length(X)), X)
    
    # Declara betas reais
    betas_reais <- c(input$beta0_real, input$beta1_real) 
    
    # Declara beta 0
    beta0 = c(input$beta0, input$beta1)
    
    # Executa EMV
    sol = EMV(Y, Xs, beta0, input$ints, input$tol) 
    
    # Declara betas 0 e betas 1
    betas0 <- seq(round(min(sol$beta0) -10),round(max(sol$beta0)+10),1)
    betas1 <- seq(round(min(sol$beta1) -10),round(max(sol$beta1)+10),1)
    
    # Define data frame espaço de busca
    espaco <- tibble(beta0 = rep(betas0, times = length(betas1)),
                     beta1 = rep(betas1, each = length(betas0))) %>% 
      mutate(dif_betas = sqrt((beta0 - betas_reais[1])^2 +
                                (beta1 - betas_reais[2])^2))
    
    # Retorno
    return(list(sol,espaco))
    
  })
  
  # Renderiza Grafico
  output$grafico2 <- renderImage({
    
    # Recupera objetos
    sol <- animacao2()[[1]]
    espaco <- animacao2()[[2]]
    
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
      labs(title = "Distância Euclidiana de Betas Reais e Estimados",
           x = "Beta 0",
           y = "Beta 1",
           fill = "Distância",
           subtitle = paste("Iteração:  ","{frame_time}"))+
      theme(panel.grid = element_blank(),
            panel.background = element_blank())
    
    # Salva animacao
    anim_save("outfile.gif", animate(p, fps = 5, renderer = gifski_renderer(loop = F))) 
    
    # Retorna lista contendo o arquivo
    list(src = "outfile.gif",
         contentType = 'image/gif',
         width = 600,
         height = 400)
    
  }, deleteFile = TRUE)
  
}

# --------------------------------------------------------------------------

# Executa o aplicativo 
shinyApp(ui = ui, server = server)
