

# -----------------------------------------------------------------------------------------------

# Aplicativo IWLS

# -----------------------------------------------------------------------------------------------

# Carrega Bibliotecas
library(shiny)
library(shinydashboard)
library(tidyverse)
library(tibble)
library(ggplot2)
library(gganimate)
library(hrbrthemes)
library(viridis)
library(gifski)
library(transformr)
library(metR)

# --------------------------------------------------------------------------

# Define funções

# --------------------------------------------------------------------------

# Calcula Y
calculaY <- function(eta){
  
  # Calcula Y Poisson
  exp(eta)
  
}

# Declara função calcula W de W
calcula_W <- function(eta, theta){
  
  # Calcula W Poisson
  (1 / theta) * exp(2 * eta)
  
}

# Declara função Z
calcula_Z <- function(eta, Y){
  
  # Calcula z Poisson
  eta + ( Y * exp(-eta) ) - 1 
  
}

#Declara funcao theta
calcula_theta <- function(eta, dist){
  
  # Conficao de distribuicao
  if(dist == "poisson"){
    
    # Calcula theta poisson
    exp(eta)
    
  } else if(dist == "exp"){
    
    # Calcula theta gama
    (-1 / eta)
    
  } else if(dist == "bin"){
    
    # Calcula theta binomial
    exp(eta) / (1 + exp(eta))
    
  } else if(dist == "nbin"){
    
    # Calcula theta da normal
    exp(eta)
    
  }
  
}

# Funcao estimativa de betas
EMV <- function(Y, Xs, beta, tol, dist){
  
  # Declara betas
  betas <- tibble(iter = 1,
                  beta0 = as.numeric(beta[1]),
                  beta1 = as.numeric(beta[2]))
  
  # Declara norma
  norma <- 1
  
  # Declara iteracao
  r = 1
  
  # Laco para popular betas calculados nas iterações
  while(norma > tol){
    
    # Recupera beta da vez
    beta <- t(betas[r, c(2,3) ] %>% as.matrix())
    
    # Calcula eta0
    eta <- Xs %*% beta
    
    # Calcula media inicial
    theta <- calcula_theta(eta, dist)
    
    # Calcula diagonal da matriz Wi
    W_i <- calcula_W(eta, theta)
    
    # Declara matriz W
    W <- diag(length(Y))
    
    # Popula diagonal da matriz W
    diag(W) <- W_i
    
    # Calcula z
    z_i <- calcula_Z(eta, Y)
      
    # Calcula beta_0
    betas[r+1, "iter" ] <- r + 1
    
    # Calcula beta_0
    betas[r+1, "beta0" ] <- (solve( t(Xs) %*% W %*% Xs ) %*% t(Xs) %*% W %*% z_i  )[1]
    
    # Calcula beta_1
    betas[r+1, "beta1" ] <- (solve( t(Xs) %*% W %*% Xs ) %*% t(Xs) %*% W %*% z_i  )[2]
    
    # Calcula norma euclideana
    norma <- sqrt(sum((abs(betas[r+1, c(2,3)] - betas[r, c(2,3)]))^2))
    
    # Incrementa r
    r = r + 1
    
  }
  
  # Retorno da função
  return(betas)
  
}

# Declara dicionario de distribuicoes
distribuicoes <- c("Poisson" = "poisson",
                   "Exponêncial" = "exp",
                   "Binomial" = "bin",
                   "Binomial Negativa" = "nbin")


# --------------------------------------------------------------------------

# Interface do Usuário

# --------------------------------------------------------------------------

# Declara Interface do usuario
ui <- dashboardPage(
  
  # Cabecalho do aplicativo
  dashboardHeader(
    title = span(
      tags$img(src="imagens/optimizador_iwls.png", width = '100%'))
  ),
  
  # Declara barra lateral
  dashboardSidebar(
    
    # Seletores de Entrada
    selectInput("dist",
                "Distribuição",
                choices = c("Poisson" = "poisson",
                            "Exponêncial" = "exp",
                            "Binomial" = "bin",
                            "Binomial Negativa" = "nbin")),
    sliderInput("obs",
                "Número de Observações:",
                min = 10,
                max = 100,
                value = 25),
    sliderInput("beta0_real",
                "Beta 0 Real:",
                min = 0,
                max = 10,
                value = 2,
                step = 0.1),
    sliderInput("beta1_real",
                "Beta 1 Real",
                min = -10,
                max = 10,
                value = -2,
                step = 0.1),
    actionButton(inputId = "amostre", 
                 label = "Amostrar",
                 style="color: #fff; background-color: #008521; border-color: #2e6da4"),
    sliderInput("beta0",
                "Beta 0 Inicial",
                min = -10,
                max = 10,
                value = 1,
                step = 0.1),
    sliderInput("beta1",
                "Beta 1 Inicial",
                min = -10,
                max = 10,
                value = 1,
                step = 0.1),
    sliderInput("tol",
                "Tolerância (10^-)",
                min = 1,
                max = 10,
                value = 5),
    actionButton(inputId = "optimize", 
                 label = "Optimizar",
                 style="color: #fff; background-color: #008521; border-color: #2e6da4")
    
  ),
  
  # Declara corpo do aplicativo
  dashboardBody(
    # Define html style
    tags$style(".skin-blue .main-header .logo { background-color: #3e8fce; }
                .skin-blue .main-header .logo:hover { background-color: #3e8fce;}
                .skin-blue .main-header .navbar { background-color: #3e8fce;}"),
    fluidRow(
      box(width = 12, imageOutput("grafico1"), align = "center")
    ),
    fluidRow(
      box(width = 12, imageOutput("grafico2"), align = "center")
    )
  )
  
)

# -----------------------------------------------------------------------------------------------

# Declara Servidor

# -----------------------------------------------------------------------------------------------

# Declara função servidor
server <- function(input, output, session){
  
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
    Y = calculaY(eta)
    
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
    sol = EMV(Y, Xs, beta0, 10^(-input$tol), dist = input$dist) 
    
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
    resultados <- df_anime %>% 
                    bind_rows(original)
    
    # Retorno
    return(list(resultados))
    
  })
  
  # Renderiza Grafico
  output$grafico1 <- renderImage({
    
    # Declara animacao
    anime <- animacao1()[[1]]
    
    # Observa evento
    distribuicao <- isolate(input$dist)
    
    # Recupera nome da distribuicao
    nome <- names(distribuicoes[distribuicoes==distribuicao])
    
    # Declara objeto gráfico
    p <- anime %>%
      ggplot( aes(x=X, y=pred, group = Y, color = Y, shape = Y)) +
      geom_point(size = 5) +
      scale_shape_manual(values=c(1, 20))+
      scale_colour_manual(values = c("red", "blue"))+
      labs(title = paste0("Curva de Ajuste da Distribuição ", nome),
           subtitle = paste("Iteração:  ","{closest_state}"),
           y = "Y")+
      transition_states(estados,
                        transition_length = 2,
                        state_length = 4,
                        wrap = FALSE) +
      view_follow(fixed_x = T)+
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            text = element_text(size=13))
    
    # Salva animacao
    anim_save("outfile.gif", animate(p, fps = 5, renderer = gifski_renderer(loop = F))) 
    
    # Retorna lista contendo o arquivo
    list(src = "outfile.gif",
         contentType = 'image/gif',
         width = 550,
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
    sol = EMV(Y, Xs, beta0, 10^(-input$tol), dist = input$dist) 
    
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
           subtitle = paste("Tempo:  ","{frame_time}"))+
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            text = element_text(size=13))
    
    # Salva animacao
    anim_save("outfile.gif", animate(p, fps = 5, renderer = gifski_renderer(loop = F))) 
    
    # Retorna lista contendo o arquivo
    list(src = "outfile.gif",
         contentType = 'image/gif',
         width = 550,
         height = 400)
    
  }, deleteFile = TRUE)
  
}

# -----------------------------------------------------------------------------------------------

# Executa função

# -----------------------------------------------------------------------------------------------


# Executa aplicativo
shinyApp(ui, server)

# -----------------------------------------------------------------------------------------------