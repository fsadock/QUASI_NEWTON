#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
# --------------------------------------------------------------------------

# Aplicativo Shiny

# --------------------------------------------------------------------------
 
# Carrega Bibliotecas
library(shiny)
library(tidyverse)

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
    
    for(r in 1:iter){
        
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

# Declara interface do usuario
ui <- fluidPage(

    # Titulo do aplicativo
    titlePanel("Estimador IWLS Poisson"),

    # Barra lateral de entradas 
    sidebarLayout(
        sidebarPanel(
            h3("Amostra"),
            sliderInput("obs",
                        "Número de Observações:",
                        min = 10,
                        max = 100,
                        value = 25),
            sliderInput("theta",
                        "Media (Theta)",
                        min = 1,
                        max = 50,
                        value = 15),
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
                        max = 50,
                        value = 30),
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
           imageOutput("grafico")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    # Observa clique do botao optimizar
    dados <- eventReactive(input$amostre, {
        
        # Declara Y
        Y = sort(rpois(input$obs, input$theta), decreasing = T)
        
        # Declara X
        X = sort(rnorm(input$obs, 0, 0.5))
        
        # Declara data frame
        amostra <- tibble(Y, X)
        
        # Retorno
        return(amostra)
        
    })

    # Observa clique do botao optimizar
    animacao <- eventReactive(input$optimize, {
            
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
    output$grafico <- renderImage({
        
        # Plot
        p <- animacao() %>%
            ggplot( aes(x=X, y=pred, group = Y, color = Y)) +
            geom_point(size = 4, alpha = 0.7) +
            labs(title = "Optimização EMV IWLS Poisson",
                 subtitle = paste("Iteração:  ","{closest_state}"),
                 y = "Y")+
            transition_states(estados,
                              transition_length = 2,
                              state_length = 4,
                              wrap = FALSE) +
            view_follow(fixed_x = T)
        
        # Salva animacao
        anim_save("outfile.gif", animate(p, fps = 5, renderer = gifski_renderer(loop = F))) 
        
        # Retorna lista contendo o arquivo
        list(src = "outfile.gif",
             contentType = 'image/gif',
             width = 800,
             height = 600)
        
    }, deleteFile = TRUE)
        
}

# Run the application 
shinyApp(ui = ui, server = server)
