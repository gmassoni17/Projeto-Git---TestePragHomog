library('tidyverse')
library('combinat')
library('glue')


path <- "~/' Gabriela Massoni/Doutorado/Tese/Hipóteses Pragmáticas"
setwd(path)

source(glue('{path}/Códigos/Teste Homogeneidade kx2/00. TH_kx2_funcoes_dist.R'))


###### Grid de valores

k <- 3 # número de subpops
n <- rep(10, k)
N <- 100 # tamanho do grid
eps <- 0.1

thetas <- seq(0.00000001, 0.99999999, length.out = N) # grid de thetas
## criar grid com todas as combinações de k thetas
aux <- expand_grid(thetas, thetas)
for(i in 3:(k)){
  aux <- expand_grid(aux, thetas)
  names(aux) <- glue('theta{1:i}')
}
grid_thetas <- aux
grid_thetas

# grid_thetas <- grid_thetas %>% 
#   mutate(theta3 = 0.5) %>% 
#   distinct()
# grid_thetas




############################ LER DADOS
theta0 <- 0.5

dados_gerados_BP_simplex <- readRDS(
  glue("{path}/Teste Homogeneidade kx2/Dados/dados_gerados_BP_simples_k_{k}_theta0_{theta0}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))

dados_gerados_BP_full <- readRDS(
  glue("{path}/Teste Homogeneidade kx2/Dados/dados_gerados_BP_full_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))



#######################################
###### PLOTAR GRÁFICOS NO R2 (projeção)
#######################################

## ESCOLHER QUAL DISTÂNCIA e TIPO
dist <- "BP"
tipo <- "full"

if(dist == "BP"){
  if(tipo == 'simples'){
    dados_gerados <- dados_gerados_BP_simplex
    tipo1 <- paste(tipo, glue('_theta0_{theta0}'), sep = "")
  }
  if(tipo == 'full'){
    dados_gerados <- dados_gerados_BP_full  
    tipo1 <- tipo
  }
}

# if(dist == "KL"){
#   if(tipo == 'simples'){
#     dados_gerados <- dados_gerados_KL_simplex
#     tipo1 <- paste(tipo, glue('_theta0_{theta0}'), sep = "")
#   }
#   if(tipo == 'full'){
#     dados_gerados <- dados_gerados_KL_full 
#     tipo1 <- tipo
#   }
# }

caminho_graf <- glue("./Teste Homogeneidade kx2/Graficos/{dist}/TesteHomog_k_{k}_{tipo1}_n_{paste(n, collapse = '_')}_eps_{eps}")
fig2 <- plot_prag(dados_gerados, col_3d = 'Viridis')
htmlwidgets::saveWidget(fig2, glue("{caminho_graf}.html"))
browseURL(glue("{caminho_graf}.html"))  


# ggsave(glue("{caminho_graf}"),
#        width = 5,
#        height = 5)


########### aqui (pra baixo, não adaptei)

#########################################################################################################
#########################################################################################################

## checando ponto teórico
if(tipo == 'simples'){
  k <- (n1+n2)*(theta0)*(1-theta0) ## constante calculada
  plot_prag(dados_regiao_prag) +
    ## plotar interseção da elipse(theta0) com r
    geom_point(aes(x = theta0 - sqrt((eps*k)/(n1^2+n2^2)),
                   y = theta0 + sqrt((eps*k)/(n1^2+n2^2)))) +
    geom_point(aes(x = theta0 + sqrt((eps*k)/(n1^2+n2^2)),
                   y = theta0 - sqrt((eps*k)/(n1^2+n2^2))))
  
  } 
if(tipo == 'full'){
  k <- (n1+n2)*thetas*(1-thetas)
  plot_prag(dados_regiao_prag) +
    ## plotar interseção da elipse(theta0) com r
    geom_line(aes(x = thetas - sqrt((eps*k)/(n1^2+n2^2)),
                   y = thetas + sqrt((eps*k)/(n1^2+n2^2)))) +
    geom_line(aes(x = thetas + sqrt((eps*k)/(n1^2+n2^2)),
                   y = thetas - sqrt((eps*k)/(n1^2+n2^2))))
  
}



## plotar g1 e g2

q <- eps/n
g1 <- (q + thetas - thetas*q + sqrt(q*(q - 4*thetas^2 + 4*thetas)))/(1 + q)
g2 <- (q + thetas - thetas*q - sqrt(q*(q - 4*thetas^2 + 4*thetas)))/(1 + q)

plot_prag(dados_regiao_prag) +
  geom_line(aes(x = thetas,
                 y = g1)) +
  geom_line(aes(x = thetas,
                y = g2))




################################

###### Comparando elipse centrada em theta0 e em theta' = theta0 + delta (hipótese simples)

### para theta0 = 0.2
delta <- 0.02
dados_gerados_BP_simplex_linha_deltapos <- gerar_infos(BP_diss,
                                                      thetas0 = theta0 + delta, ## muda aqui
                                                      grid_thetas,
                                                      eps = eps,
                                                      n1 = n1,
                                                      n2 = n2)
borda_sup <- dados_gerados_BP_simplex_linha_deltapos %>% 
 filter(prag) %>% 
 group_by(theta1) %>% 
 summarise(theta2 = max(theta2)) %>% 
 arrange(desc(theta1))
borda_inf <- dados_gerados_BP_simplex_linha_deltapos %>% 
 filter(prag) %>% 
 group_by(theta1) %>% 
 summarise(theta2 = min(theta2))
dados_regiao_prag_linha_deltapos <- borda_inf %>% bind_rows(borda_sup)

delta <- 0
dados_gerados_BP_simplex_linha_deltaneg <- gerar_infos(BP_diss,
                                                      thetas0 = theta0 + delta, ## muda aqui
                                                      grid_thetas,
                                                      eps = eps,
                                                      n1 = n1,
                                                      n2 = n2)
borda_sup <- dados_gerados_BP_simplex_linha_deltaneg %>% 
 filter(prag) %>% 
 group_by(theta1) %>% 
 summarise(theta2 = max(theta2)) %>% 
 arrange(desc(theta1))
borda_inf <- dados_gerados_BP_simplex_linha_deltaneg %>% 
 filter(prag) %>% 
 group_by(theta1) %>% 
 summarise(theta2 = min(theta2))
dados_regiao_prag_linha_deltaneg <- borda_inf %>% bind_rows(borda_sup)

## plot elipse de theta'
k <- (n1+n2)*(theta0)*(1-theta0) ## constante calculada
k1 <- n1^2/(eps*k)
k2 <- n2^2/(eps*k)
plot_prag(dados_regiao_prag) +
 ## Região Pragmática - delta Positivo
 geom_polygon(data = dados_regiao_prag_linha_deltapos,
              aes(x = theta1, y = theta2),
              col = 'deepskyblue3', fill = 'blue', alpha = 0.2) +
 ## Região Pragmática - delta Negativo
 # geom_polygon(data = dados_regiao_prag_linha_deltaneg,
  #            aes(x = theta1, y = theta2),
   #           col = 'darkolivegreen', fill = 'green', alpha = 0.2) +
 ## plotar centro 
 # geom_point(aes(x = theta0, y = theta0)) +
 # geom_point(aes(x = theta0 + delta, y = theta0 + delta))
 ## plotar reta perpendicular (r)
 geom_line(aes(x = thetas, y = 2*theta0 - thetas)) + 
  theme(axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)) +
  geom_text(aes(x = 0.99, y = 0.9, label = glue(" t1 = t2"))) +
  #   ## plotar reta perpendicular (r)
  geom_line(aes(x = thetas, y = 2*theta0 - thetas), col = 'deepskyblue3') +
  geom_text(aes(x = 0.9, y = 0, label = glue("r: 2t0 = t1 + t2")))

 ## plotar interseção da elipse(theta0) com r
 geom_point(aes(x = theta0 - sqrt((eps*k)/(n1^2+n2^2)),
                y = theta0 + sqrt((eps*k)/(n1^2+n2^2)))) +
 geom_point(aes(x = theta0 + sqrt((eps*k)/(n1^2+n2^2)),
                y = theta0 - sqrt((eps*k)/(n1^2+n2^2)))) +
 ## plotar interseção da elipse(theta1) com r
 geom_point(aes(x = 2*theta0 - ((k1*(theta0 - delta) + k2*(theta0 + delta))/(k1+k2) + sqrt(k1+k2 - 4*k1*k2*delta^2)/(k1+k2)),
                y = (k1*(theta0 - delta) + k2*(theta0 + delta))/(k1+k2) + sqrt(k1+k2 - 4*k1*k2*delta^2)/(k1+k2)
         )) +
 geom_point(aes(x = 2*theta0 - ((k1*(theta0 - delta) + k2*(theta0 + delta))/(k1+k2) - sqrt(k1+k2 - 4*k1*k2*delta^2)/(k1+k2)),
                y = (k1*(theta0 - delta) + k2*(theta0 + delta))/(k1+k2) - sqrt(k1+k2 - 4*k1*k2*delta^2)/(k1+k2)
 ))
 
 ggsave(glue("./Teste Homogeneidade 2x2/Graficos/{dist}/TesteHomog_2x2_{tipo1}_n1_{n1}_n2_{n2}_eps_{eps}_teste_retaperpend.png"),
        width = 5,
        height = 5)
 
# 
# ggsave(paste("./Gráficos/",
#              dist, 
#              "/TesteHomog_2x2_COMPARACAO_", 
#              "dist_", dist, 
#              "_n1_", n1, "_n2_", n2,
#              "_eps_", eps,
#              "_N_", N,
#              '_', tipo,
#              "_compararthetas",
#              ".png", sep=""),
#        width = 5,
#        height = 5)


############### checando se extremos da elipse coincidem com téorico (BP)
k <- (n1+n2)*(thetas)*(1-thetas) ## constante calculada
k1 <- n1^2/(eps^2*k)
k2 <- n2^2/(eps^2*k)


## plotar g1 para elipse centrada em theta' = theta0 + delta
delta <- -(k1-k2+4*k1*k2)/(8*k1*k2) +
  sqrt((k1-k2+4*k1*k2)^2 + 16*k1*k2*(k1+k2))/(8*k1*k2) ### delta que maximiza distância à reta r
prag_yl <- (k1*(thetas - delta) + k2*(thetas + delta))/(k1+k2) +
  sqrt(k1+k2 - 4*k1*k2*delta^2)/(k1+k2)
prag_xl <- 2*thetas - prag_yl

plot_prag(dados_regiao_prag) + 
  # geom_line(aes(prag_x, prag_y), col = 'blue') +
  geom_line(aes(prag_xl, prag_yl), col = 'red')


