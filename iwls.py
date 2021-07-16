
# -------------------------------------------------------------------------------------------------------

# Implementação do Algoritmo IWLS Optimização MLG Estimador de Máxima Verossimilhança

# -------------------------------------------------------------------------------------------------------

# Importa Pacotes
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# -------------------------------------------------------------------------------------------------------

# Seja Y1, Y2, . . . , Yn uma amostra aleatória tal que Yi ∼ Poisson(θi) com θi > 0. Desejamos analisar
# o impacto do regressor X1i sobre a resposta Yi. Considere o preditor linear ηi = β0 + β1X1i.
# Denotando β = (β0, β1)' e Xi = (1, X1i)'
# podemos tambem escrever ηi = Xi β.
# No caso Poisson, lembre que a função de ligação é estabelecida por θi = exp{ηi}.
# A seguinte amostra devera ser utilizada para resolver esta questão:

# Declara resposta e covariável
Y = np.array([7, 2, 24, 33, 2, 5, 2, 5, 8, 8, 16, 20, 11, 4, 37, 1, 8, 8, 5, 7, 5, 4, 4, 2, 26])
X = np.array([-0.27, 0.57, -0.74, -0.94, 0.64, 0.86, 0.5, 0.12, -0.36, 0.12, -0.55, -0.81, -0.21,
       0.53, -1, 0.4, -0.07, 0.15, -0.13, 0.22, 0.96, 0.45, -0.06, 0.97, -0.71])

# Dados da questão 2
df = pd.DataFrame({"Y" : Y, "X" : X})

# Fecha janela de graficos
plt.close()

# Imprime grafico
sns.scatterplot(data = df, x = "X", y = "Y")
plt.title("Dispersão da Amostra Poisson")

# Exibe grafico
plt.show()

# -------------------------------------------------------------------------------------------------------

# a) Estabeleça o chute inicial βˆ(r=0) = (1, 1) a ser utilizado no IWLS.
# Usando este chute e o resultado amostral fornecido no enunciado.
# Construa no R (exibir o script) a matriz W(0) correspondente a este problema.
# Mostre a submatriz 5 × 5 formada pelas linhas 1 a 5 e colunas 1 a 5 da matriz W(0).

# -------------------------------------------------------------------------------------------------------

# Declara função calcula W de W
def calcula_W(eta, theta):
  
  # Calcula W
  return (1 / theta) * np.exp(2 * eta)
  
# Declara beta0 chute incial
beta = np.array([1, 1])

# Declara matriz de Xs
Xs = np.stack([np.transpose(np.repeat(1, X.shape)), np.transpose(X)], axis = 1)

# Calcula eta0
eta = np.matmul(Xs, beta)

# Calcula media inicial
theta = np.exp(eta)

# Calcula diagonal da matriz Wi
W_i = calcula_W(eta, theta)

# Declara matriz W
W = np.eye(Y.shape[0])
np.fill_diagonal(W, W_i)

# Imprime primeiros 5 elementos da matriz diagonal
W[0:5, 0:5]

# -------------------------------------------------------------------------------------------------------

# b) Ainda levando em conta o chute inicial βˆ(r=0) = (1, 1)>, calcule o vetor z (0)
# definido no IWLS deste problema relacionado a Poisson.
# Mostre os valores do seu vetor z(0).

# -------------------------------------------------------------------------------------------------------

# Declara função Z
def calcula_Z(eta, Y):
  
  # Calcula Z
  return eta + ( Y * np.exp(-eta) ) - 1 
  

# Calcula z0
z_i = calcula_Z(eta, Y)

# -------------------------------------------------------------------------------------------------------

# c) Utilizando W(0) e z(0) obtidos em (a) e (b), respectivamente,
# calcule a estimativa β(1) = (X'W(0)X)−1 X'W(0)z(0) referente a primeira iteração do IWLS.
# Mostre o script R e o resultado final de β(1).

# -------------------------------------------------------------------------------------------------------

# Define funções
def mat(A, B): 
  return np.matmul(A,B)

def t(A):
  return np.transpose(A)

def inv(A):
  return np.linalg.inv(A)


# Calcula beta_1
beta_1 = mat( mat( mat( inv( mat( mat( t(Xs), W) , Xs) ), t(Xs) ), W), z_i)

# Funcao estimativa de betas
def EMV(Y, Xs, beta0, it, tol):
  
  # Declara betas
  betas = pd.DataFrame({"it" : np.arange(it),
                        "beta0" : beta0[0],
                        "beta1" : beta0[1]})
  
  for r in np.arange(it) :
    
    # Declara beta da vez
    beta = np.array(betas.iloc[r,[1,2]])
    
    # Calcula eta0
    eta = mat(Xs, beta)

    # Calcula media inicial
    theta = np.exp(eta)
    
    # Calcula diagonal da matriz Wi
    W_i = calcula_W(eta, theta)
    
    # Declara matriz W
    W = np.eye(Y.shape[0])
    np.fill_diagonal(W, W_i)

    # Calcula z
    z_i = calcula_Z(eta, Y)
    
    # Calcula beta_0
    betas.iloc[r+1, 1 ] = mat( mat( mat( inv( mat( mat( t(Xs), W) , Xs) ), t(Xs) ), W), z_i)[0]
    
    # Calcula beta_1
    betas.iloc[r+1, 2 ] = mat( mat( mat( inv( mat( mat( t(Xs), W) , Xs) ), t(Xs) ), W), z_i)[1]
    
    # Calcula norma euclideana
    norma = np.sqrt(sum((abs(betas.iloc[r+1, [1,2]] - betas.iloc[r, [1,2]]))**2))
    
    # Condição de parada
    if(norma < tol):
      
      # Declara df resultado
      resultados = betas.iloc[:r+1,]
      
      # Iterrompe laço
      break
  
  # Retorno da função
  return(resultados)

# Declara beta0
beta0 = np.array([1,1])

# Executa EMV
sol = EMV(Y, Xs, beta0, 20, 10e-5) 

# calcula predicao
pred = np.exp(sol.iloc[-1,1] + sol.iloc[-1,2] * X)

# Df predicao
df_pred = pd.DataFrame({"Observado" : Y, "X" : X, "Estimado" : pred})

data = pd.melt(df_pred, id_vars='X', value_vars=['Observado', 'Estimado'], var_name = "Y", value_name = "Valor")

# Grafico de predicao
plt.close()
sns.scatterplot(data = data, x = "X", y = "Valor", hue = "Y")
plt.show()


# -------------------------------------------------------------------------------------------------------

# Simulação

# -------------------------------------------------------------------------------------------------------

# Declara Y
Y = np.sort(np.random.poisson(20, 25))[::-1]

# Declara X
X = np.sort(np.random.normal(0, 0.5, 25))

# Declara matriz de Xs
Xs = np.stack([np.transpose(np.repeat(1, X.shape)), np.transpose(X)], axis = 1)

# Declara beta 0
beta0 = np.array([1,1])

# Executa EMV
sol = EMV(Y, Xs, beta0, 40, 10e-5) 

# calcula predicao
pred = np.exp(sol.iloc[-1,1] + sol.iloc[-1,2] * X)

# Df predicao
df_pred = pd.DataFrame({"Observado" : Y, "X" : X, "Estimado" : pred})

data = pd.melt(df_pred, id_vars='X', value_vars=['Observado', 'Estimado'], var_name = "Y", value_name = "Valor")

# Grafico de predicao
plt.close()
sns.scatterplot(data = data, x = "X", y = "Valor", hue = "Y")
plt.show()

# -------------------------------------------------------------------------------------------------------

# d) Use os dados amostrais fornecidos no enunciado e aplique a função glm do R para estimar β0 e β1.
# Use o mesmo chute inicial sugerido em (a).
# Mostre a saída do comando summary e avalie a significancia dos coeficientes.
# Atenção! sua análise sobre a significância deve indicar claramente a estimativa do coeficiente.
# e o valor-p sob avaliacão.
# A questão não está pedindo para você realizar vários ajustes para encontrar o melhor modelo.
# O objetivo aqui e apenas ver se você sabe interpretar a saída computacional solicitada.


# -------------------------------------------------------------------------------------------------------

# Ajuste
ajuste = glm(Y ~ X1,
             family = poisson(link = "log"),
             data = df2,
             start = c(1,1) )

# Imprime resultados
summary(ajuste)

# Funcao de grafico
fun.gen1 <- function(X1) exp(ajuste$coef[1] + ajuste$coef[2] * X1)

# Funcao de grafico
fun.gen2 <- function(X1) exp(betas$beta0[15] + betas$beta1[15] * X1)

# Desenha gráfico
ggplot(data = df2, aes(x = X1, y = Y))+
  geom_point() +
  stat_function(fun = fun.gen1, linetype = "dotdash", size = 2, aes( color = "GLM")) +
  stat_function(fun = fun.gen2, linetype = "dotted", size = 2, aes( color = "IWLS")) +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              linetype = "dashed", aes( color = "ggplot"))+
  labs(title = "Comparação de Métodos MLG Poisson",
       color = "Método")
  
