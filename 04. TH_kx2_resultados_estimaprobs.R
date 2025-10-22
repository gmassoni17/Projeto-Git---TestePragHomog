library('tidyverse')
library('glue')
library('coda')
library('geometry')


path <- "~/' Gabriela Massoni/Doutorado/Tese/Hipóteses Pragmáticas"
setwd(path)
source(glue('{path}/Códigos/Teste Homogeneidade kx2/00. TH_kx2_funcoes_hpd.R'))
source(glue('{path}/Códigos/Teste Homogeneidade kx2/00. TH_kx2_funcoes_dist.R'))



#########################################################
############ Definir os parâmetros do modelo ############
#########################################################

confianca <- 0.8
k <- 3
a <- NULL
b <- NULL
n <- NULL
for(i in 1:k){a[i] <- 2}
for(i in 1:k){b[i] <- 2}
for(i in 1:k){n[i] <- 10}
N <- 100000 # número de vezes que vou gerar da posteriori


### Binomiais independentes (não H0)
# parâmetros para amostra (dif da região prag)

# chutes de theta1 e theta2
t1 <- c(0.1, 0.4, 0.3)
t2 <- c(0.9, 0.5, 0.2)
t3 <- c(0.5, 0.4, 0.2)
m <- NULL
for(i in 1:k){m[i] <- 500}
m_list <- m

caminho_dados_simu <- glue("./Teste Homogeneidade kx2/Dados/dados_simulacao_testeHomog_kx2_m_{paste(m_list, collapse = '_')}_k_{k}.rds")
## lendo dados gerados no script 3
dados_simu <- read_rds(caminho_dados_simu)


#############################################################################################
########### Passo 1: gerar pontos da posteriori
#############################################################################################

### Para um dado fixado (X_1, X_2), geramos N pontos da posteriore de theta
ii <- 5 # depois, fazer for (uma iteração para cada conjunto de dados observado)
dados <- dados_simu[ii,]

thetas_gerados <- purrr::map(1:length(m), ~{
  theta_post_gerado <- rbeta(N, a[.x] + as.numeric(dados[,.x]),
                             m[.x] + b[.x] - as.numeric(dados[,.x]))
  return(theta_post_gerado)
  })

# gráfico de densidade da posteriore gerada para theta_1 e theta_2
ggplot() +
  geom_density(aes(thetas_gerados[[1]]), col = 'aquamarine3', cex = 1) +
  geom_density(aes(thetas_gerados[[2]]), col = 'deeppink3', cex = 1) +
  geom_density(aes(thetas_gerados[[3]]), col = 'black', cex = 1) +
  theme_classic() +
  xlab('theta')


thetas <- tibble(theta1_gerado = thetas_gerados[[1]])
for(i in 2:length(m)){ # colunando os "k" thetas_gerados
  thetas <- thetas %>% 
    mutate(thetas_gerados_aux = thetas_gerados[[i]])
  names(thetas)[i] <- glue('theta{i}_gerado')
}
thetas # matriz com os valores gerados de (theta_1, theta_2, theta_3)


#############################################################################################
########### Passo 2: Calculando corte do HPD
#############################################################################################

counts <- dados
posterior_par_list <- lapply(1:k, function(j) {
  posterior(dados = as.integer(counts[, j]), 
            prior_par = prior_par_list[[j]],
            m = m_list[j]) %>% as.numeric()
  })
th <- hpd_post_cut(k = k,
                   posterior_par_list = posterior_par_list,
                   probability = confianca)
th #7.92
# Obs.: este corte refere-se ao log da probabilidade a posteriore



#############################################################################################
########### Passo 3: VERIFICANDO PERTENCIMENTO DO PONTO GERADO NA PRAGMÁTICA ESPECÍFICA 
#############################################################################################

# Obs: theta0_linha  =  (theta1_gerado + theta2_gerado + ... + thetak_gerado) / k

# ATENÇÃAAAAO: supondo n1=n1=...=nk - definido no início do problema
eps <- 0.1
dados_thetas_verificar <- thetas %>% # thetas gerados da posteriori
  mutate(sum_thetas = select(., starts_with("theta")) %>% rowSums(.,na.rm = TRUE)) %>% 
  mutate(th = th) %>% 
  mutate(theta0_linha = sum_thetas/ncol(thetas),
         delta_estrela = eps/(eps+n[1]) * (1/2 - theta0_linha),
         theta_verificar = theta0_linha + delta_estrela) # theta do centro da bolinha específica que eu tenho que verificar se o ponto gerado pertence
dados_thetas_verificar %>% glimpse()


###### (i) Calcular o "cubo" x1,x2,x3,y1,y2,y3 (caso k = 3): (retomar qual quadrado - DOC 7 - é o que o Helder desenhou)

dados_thetas_verificar <- dados_thetas_verificar %>%
  mutate(
    raio = sqrt(2*eps/n[1]*(theta0_linha + delta_estrela) *
                  (1 - (theta0_linha + delta_estrela))),
    limite_inf = theta_verificar - raio,
    limite_sup = theta_verificar + raio,
    indc_dentro = if_all(starts_with("theta"),
                    ~ .x >= limite_inf & .x <= limite_sup)
  )

# Calculo no caso k = 2
# dados_thetas_verificar <- dados_thetas_verificar %>% 
#   mutate(raio = sqrt(2*eps/n*(theta0_linha + delta_estrela)*(1- (theta0_linha + delta_estrela))),
#          x1 = theta_verificar - raio,
#          x2 = theta_verificar + raio,
#          y1 = theta_verificar - raio,
#          y2 = theta_verificar + raio) %>%
#   ## verificar se pertece
#   mutate(indc1_in_quadrado = theta1_gerado >= x1 & theta1_gerado <= x2,
#          indc2_in_quadrado = theta2_gerado >= y1 & theta2_gerado <= y2,
#          indc_in_quadrado = indc1_in_quadrado & indc2_in_quadrado)

###### (ii) Verificar se está dentro do raio
dados_thetas_verificar <- dados_thetas_verificar %>%
  mutate(across(ends_with('gerado'), ~ (.x-theta_verificar)^2, 
                .names = '{.col}_quadrado_menos_theta_vertificar')) %>% 
  mutate(sum_thetas_quadrado_menos_theta_vertificar = select(., ends_with('_quadrado_menos_theta_vertificar')) %>% 
           rowSums(.,na.rm = TRUE)) %>% 
  mutate(raiz_sum_thetas_quadrado_menos_theta_vertificar = sqrt(sum_thetas_quadrado_menos_theta_vertificar) <= raio) %>% 
  mutate(indc_in_prag = ifelse(indc_dentro,
                             ifelse(sqrt(sum_thetas_quadrado_menos_theta_vertificar) <= raio,1,0),0)) 
dados_thetas_verificar %>% glimpse()


#############################################################################################
########### Passo 4: ESTIMANDO P(pragmática) (a posteriori)
#############################################################################################

n_pontos_in_prag <- dados_thetas_verificar %>% filter(indc_in_prag==1) %>% nrow()
prob_hat_prag <- n_pontos_in_prag/nrow(dados_thetas_verificar)
prob_hat_prag
# Se dados_simu[1,] - 0% (X1 = 57, X2 = 435, X3 = 228, n1=n2=n3 = 10)
# Se dados_simu[2,] - 18,3% (X1 = 185, X2 = 245, X3 = 203) 
# Se dados_simu[5,] - 93,8% (X1 = 255, X2 = 269, X3 = 247) 
# Se dados_simu[6,] - 98,7% (X1 = 455, X2 = 450, X3 = 450, n = 10)


#############################################################################################
########### Passo 5: ESTIMAR CONJUNTO DO HPD 
#############################################################################################


# Seleção dos thetas do grid (assumindo nomes theta1, theta2, ..., thetak)
theta_cols <- paste0("theta", 1:k, "_gerado")
thetas_row <- thetas %>%
  select(all_of(theta_cols)) %>%
  as.data.frame()

conjunto_hpd_estimado <- thetas %>% # thetas gerados da posteriori 
  mutate(corte_hpd = th) %>% 
  mutate(log_poster = apply(thetas_row, 1, function(t_vec) {
                log_posterior(as.numeric(t_vec), posterior_par_list)})
    ) %>%
  mutate(hpd = log_poster > th) %>% 
  filter(hpd)

conjunto_hpd_estimado


## Carrega dados reticulado que formam a região pragmátcia 
dist <- "BP"
tipo <- "full"

if(dist == "BP"){
  if(tipo == 'simples'){
    dados_gerados <- readRDS(
      glue("{path}/Teste Homogeneidade kx2/Dados/dados_gerados_BP_simples_k_{k}_theta0_{theta0}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))
    tipo1 <- paste(tipo, glue('_theta0_{theta0}'), sep = "")
  }
  if(tipo == 'full'){
    dados_gerados <- readRDS(
      glue("{path}/Teste Homogeneidade kx2/Dados/dados_gerados_BP_full_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))
    tipo1 <- tipo
  }
}

# plot_prag(dados_regiao_prag) +  
#   geom_polygon(aes(x = theta1_gerado, y = theta2_gerado),
#                data = conjunto_hpd_estimado,
#                color = 'deeppink3', fill = 'deeppink3', alpha = 0.2)


#############################################################################################
########### Passo 6: ESTIMAR Prob da INTERSEÇÃO DE HPD e PRAG 
#############################################################################################


# quero ver a intersecção dos conjuntos (que formam o HPD e a região Pragmática)
# e calcular a prob dele


dados_regiao_prag <- dados_gerados %>%
  filter(prag) %>% 
  select(theta1, theta2, theta3)

# Conjunto 1: pontos que formam a região pragmática em 3D
# (colunas theta1, theta2, theta3)
regiao <- as.matrix(dados_regiao_prag[, c("theta1", "theta2", "theta3")])

# Conjunto 2: pontos a verificar
pontos <- as.matrix(conjunto_hpd_estimado[, c("theta1_gerado",
                                              "theta2_gerado",
                                              "theta3_gerado")])

# Construir o "polígono" 3D - conjunto 1 (na verdade, um poliedro convexo = convex hull)
hull <- convhulln(regiao, options = "FA")

# Verificar se cada ponto de 'pontos' está dentro do hull
dentro <- inhulln(hull, pontos)

# Adicionar ao tibble original
conjunto_hpd_estimado <- conjunto_hpd_estimado %>%
  mutate(dentro = dentro)

# Filtrar apenas os que estão dentro
pontos_dentro <- conjunto_hpd_estimado %>% filter(dentro)


# depois de calcular a coluna `dentro` (ponto está no polígono ou não):
p_intersec <- mean(conjunto_hpd_estimado$dentro)
p_intersec*100
# Se dados_simu[1,] - 0% (X1 = 57, X2 = 435, X3 = 228, n = 10)
# Se dados_simu[2,] - 39.9% (X1 = 185, X2 = 245, X3 = 203) 
# Se dados_simu[5,] - 100%  (X1 = 255, X2 = 269, X3 = 247) 
# Se dados_simu[5,] - 100% (X1 = 455, X2 = 450, X3 = 450, n = 10)
thetas

