library('tidyverse')
library('combinat')
library('glue')
library('geometry')
library('rgl')
# install.packages('sp')
# install.packages('sf')
library(sp)
library(sf)



path <- "~/' Gabriela Massoni/Doutorado/Tese/Hipóteses Pragmáticas"
setwd(path)
source(glue('{path}/Códigos/Teste Homogeneidade kx2/00. TH_kx2_funcoes_dist.R'))
source(glue('{path}/Códigos/Teste Homogeneidade kx2/00. TH_kx2_funcoes_hpd.R'))


confianca <- 0.8


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



# chutes de a,b,c,d (das betas a priori)
a1 <- b1 <- 1
a2 <- b2 <- 1
a3 <- b3 <- 1
prior_par_list <- list(c(a1, b1), c(a2, b2), c(a3, b3))

################ GERAR DADOS

### Binomiais independentes (não H0)
# parâmetros para amostra (dif da região prag)
m1 <- 500
m2 <- 500
m3 <- 500
m_list <- c(m1, m2, m3)

# chutes de theta1 e theta2
t1 <- c(0.1, 0.4, 0.3)
t2 <- c(0.9, 0.5, 0.2)
t3 <- c(0.5, 0.4, 0.2)

# sample_full <- tibble(x = rbinom(1, m1, t1[1]), 
#                       y = rbinom(1, m2, t2[1]),
#                       z = rbinom(1, m3, t3[1])) %>% 
#   bind_rows(tibble(x = rbinom(1, m1, t1[2]), 
#                    y = rbinom(1, m2, t2[2]),
#                    z = rbinom(1, m3, t3[2]))) %>% 
#   bind_rows(tibble(x = rbinom(1, m1, t1[3]), 
#                    y = rbinom(1, m2, t2[3]),
#                    z = rbinom(1, m3, t3[3])))
# 
# 
# ################ GERAR DADOS - Sob H0
# t1 <- t2 <- t3 <- c(0.1, 0.5, 0.9) 
# sample_h0 <- tibble(x = rbinom(1, m1, t1[1]), 
#                     y = rbinom(1, m2, t2[1]),
#                     z = rbinom(1, m3, t3[1])) %>% 
#   bind_rows(tibble(x = rbinom(1, m1, t1[2]), 
#                    y = rbinom(1, m2, t2[2]),
#                    z = rbinom(1, m3, t3[2]))) %>% 
#   bind_rows(tibble(x = rbinom(1, m1, t1[3]), 
#                    y = rbinom(1, m2, t2[3]),
#                    z = rbinom(1, m3, t3[3])))
# 
# dados <- bind_rows(sample_full, sample_h0)

## SALVAR DADOS 
caminho_dados_simu <- glue("./Teste Homogeneidade kx2/Dados/dados_simulacao_testeHomog_kx2_m_{paste(m_list, collapse = '_')}_k_{k}.rds")
# write_rds(dados, caminho_dados_simu)
dados <- read_rds(caminho_dados_simu)


############ Criar regiões HPDs

## criar um data frama de NA (a ser preenchido) com o número de colunas = k
hpd_labels <- tibble(idx = as.double(1:nrow(dados)), 
                    x = as.double(rep(NA, nrow(dados))), 
                    y = as.double(rep(NA, nrow(dados))), 
                    z = as.double(rep(NA, nrow(dados))))
hpd_grid = NULL


######## Grid de theta "inteligente" ###################
# gerar em um "retângulo" menor do que (0,1]^2 - onde se concentra a densidade post.

min <- 0
max <- 0.2
thetas <- seq(min, max, length.out = N) # grid de thetas
## criar grid com todas as combinações de k thetas
aux <- expand_grid(thetas, thetas)
for(i in 3:(k)){
  aux <- expand_grid(aux, thetas)
  names(aux) <- glue('theta{1:i}')
}
grid_thetas <- aux
grid_thetas

start <- Sys.time()
for(ii in 1:nrow(dados)){
  counts <- dados[ii,]
  # Geração dos parâmetros posteriores para cada beta
  posterior_par_list <- lapply(1:k, function(j) {
    posterior(dados = as.integer(counts[, j]), 
              prior_par = prior_par_list[[j]],
              m = m_list[j]) %>% as.numeric()
    })
  # Cálculo do corte HPD
  th <- hpd_post_cut(k = k,
                     posterior_par_list = posterior_par_list,
                     probability = confianca)
  
  # Seleção dos thetas do grid (assumindo nomes theta1, theta2, ..., thetak)
  theta_cols <- paste0("theta", 1:k)
  thetas_row <- grid_thetas %>%
    select(all_of(theta_cols)) %>%
    as.data.frame()
  
  # Filtragem do grid com base na densidade conjunta
  this_hpd_grid <- grid_thetas %>%
    mutate(color = apply(thetas_row, 1, function(t_vec) {
      log_posterior(as.numeric(t_vec), posterior_par_list) > th
    })) %>%
    filter(color)

  # #### Bordas
  # borda_sup <- this_hpd_grid %>% 
  #   filter(color) %>% 
  #   group_by(theta1) %>% 
  #   summarise(theta2 = max(theta2)) %>% 
  #   arrange(desc(theta1))
  # borda_inf <- this_hpd_grid %>% 
  #   filter(color) %>% 
  #   group_by(theta1) %>% 
  #   summarise(theta2 = min(theta2))
  
  reg_hpd <- this_hpd_grid %>% mutate(idx = ii)
  hpd_grid <- hpd_grid %>% bind_rows(reg_hpd)

  # extrai e normaliza os valores desejados (por padrão o primeiro de cada vetor)
  valores_normalizados <- sapply(posterior_par_list, function(p) {
    (p / sum(p))[1]
  })
  # salva no data.frame/matriz, incluindo o índice ii
  hpd_labels[ii, ] <- t(c(ii, valores_normalizados))
  print(ii)
}
end <- Sys.time()
end - start

hpd_grid <- hpd_grid %>% mutate(idx = as.factor(idx))



############################ LER DADOS (DA REGIÃOPRAGMÁTICA)

## ESCOLHER QUAL DISTÂNCIA e TIPO

dist <- "BP"
tipo <- "full"

theta0 <- 0.5

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


caminho_graf <- glue("./Teste Homogeneidade kx2/Graficos/{dist}/HPD_simulacao_testeHomog_{k}x2_{tipo1}_n_{paste(n, collapse = '_')}_eps_{eps}")
fig = add_hpd_grid(k
             ,plot_prag(dados_gerados, col_3d = 'Electric')
             ,hpd_grid
             ,hpd_labels
             ,col = "Greens")
htmlwidgets::saveWidget(fig, glue("{caminho_graf}.html"))
browseURL(glue("{caminho_graf}.html"))  


#################################################################
######### VERIFICANDO SE HPD ESTÁ CONTIDO NA PRAG ###############

dados_regiao_prag <- dados_gerados %>%
  filter(prag) %>% 
  select(theta1, theta2, theta3)
teste <- NULL
for(i in 1:nrow(dados)){
  # Seus dataframes
  df1 <- dados_regiao_prag %>% 
    bind_rows(dados_regiao_prag[1,]) # Polígono 1 (externo) - região pragmática
  df2 <- hpd_grid %>% filter(idx == i) %>% 
    bind_rows((hpd_grid %>% filter(idx == i))[1,])           # Polígono 2 (interno) - HPD
  df3 <- grid_thetas %>% bind_rows(grid_thetas[1,])          # ESPAÇO PARAMÉTRICO COMPLETO
  
  ##### Criar os dois polígonos (prag e hpd)
  prag <- st_polygon(list(as.matrix(df1[, c(1, 2, 3)])), dim = 'XYZ') %>% 
    st_sfc() %>% 
    st_sf()
  hpd <- st_polygon(list(as.matrix(df2[, c(1, 2, 3)])), dim = 'XYZ') %>% 
    st_sfc() %>% 
    st_sf()
  prag_convex <- st_convex_hull(prag)   # Calcular a envoltória convexa de ambos os polígonos
  hpd_convex <- st_convex_hull(hpd)   # Calcular a envoltória convexa de ambos os polígonos

  
  ##### CRIAR A "CASCA" DOS POLÍGONOS
  # 1. Calcula o convex hull dos pontos de prag
  hull_prag <- convhulln(df1[, c("theta1", "theta2", "theta3")], options = "FA")
  hull_faces_prag <- hull_prag$hull
  hull_hpd <- convhulln(df2[, c("theta1", "theta2", "theta3")], options = "FA")
  hull_faces_hpd <- hull_hpd$hull
  # 3. Agora, transformamos em sf POLYHEDRA (representação via pontos)
  coords_prag <- df1[, c("theta1", "theta2", "theta3")] |> as.matrix()
  coords_hpd  <- df2[, c("theta1", "theta2", "theta3")] |> as.matrix()
  
  #### Interseções entre a pragmática e o HPD
  intersecta <- st_intersects(prag_convex, hpd_convex, sparse = FALSE) # Verifica se há interseção entre os polígonos
  contidos <- geometry::inhulln(hull_prag, p = coords_hpd) # verifica se TODOS os pontos do HPD estão dentro da região convexa pragmática
  parcial <- intersecta & !all(contidos)  # Interseção mas não totalmente dentro
  
  # ### Plotando os polígonos
  
  coords <- as.matrix(df1[, c("theta1", "theta2", "theta3")]) %>% as.data.frame()
  hull_faces <- convhulln(coords)  # Retorna região convexa
  coords2 <- as.matrix(df2[, c("theta1", "theta2", "theta3")]) %>% as.data.frame()
  hull_faces2 <- convhulln(coords2)  # Retorna região convexa
  # Plot
  fig <- plot_ly() %>%
    # Adicionar a casca convexa - Prag
    add_trace(
      type = "mesh3d",
      x = coords[,1],
      y = coords[,2],
      z = coords[,3],
      i = hull_faces[,1] - 1,
      j = hull_faces[,2] - 1,
      k = hull_faces[,3] - 1,
      facecolor = rep("lightblue", nrow(hull_faces)),
      opacity = 0.5,
      name = "Casca Convexa"
    ) %>%
    # Adicionar a casca convexa - HPD
    add_trace(
      type = "mesh3d",
      x = coords2[,1],
      y = coords2[,2],
      z = coords2[,3],
      i = hull_faces2[,1] - 1,
      j = hull_faces2[,2] - 1,
      k = hull_faces2[,3] - 1,
      facecolor = rep("blue", nrow(hull_faces)),
      opacity = 0.5,
      name = "Casca Convexa"
    ) %>%
    layout(
      scene = list(
        xaxis = list(title = "theta1"),
        yaxis = list(title = "theta2"),
        zaxis = list(title = "theta3")
      )
    )
  caminho_graf <- glue("./Teste Homogeneidade kx2/Graficos/{dist}/HPD_dados_{i}")
  htmlwidgets::saveWidget(fig, glue("{caminho_graf}.html"))
  browseURL(glue("{caminho_graf}.html")) 
  
  # Exibe o resultado
  if (all(contidos)) {
    print(glue("O HPD {i} está contido na região pragmática (aceitação)."))
    aux <- 0
  } else if(parcial){
    print("Indefinido")
    aux <- 1/2
  } else {
    print(glue("O HPD {i} não está contido na região pragmática (rejeição)."))
    aux <- 1
  }
  print(aux)
  teste[i] <- aux
}


dados_result <- dados %>% 
  mutate(teste = teste) %>% 
  mutate(result = ifelse(teste == 1, 'Rejeita H0',
                         ifelse(teste == 1/2, 'Agnóstico', 
                                ifelse(teste == 0, 'Aceita H0', NA))))
dados_result
# A tibble: 6 × 5
# x     y     z teste result    
# <int> <int> <int> <dbl> <chr>     
# 1    57   435   228   1   Rejeita H0
# 2   185   245   203   0.5 Agnóstico 
# 3   162   105    90   0.5 Agnóstico 
# 4    55    41    37   0.5 Agnóstico 
# 5   255   269   247   0   Aceita H0 
# 6   455   450   450   0   Aceita H0 



