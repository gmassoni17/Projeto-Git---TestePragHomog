######################
#### DISTÂNCIAS
######################


#### Best Predicition Dissimilarity (BP)

BP_diss <- function(k,       # quantas populações temos (verificar se os vetores passados condizem com esse tamanho)
                    ts,      # valores dos parâmetros - deve ser um vetor com k valores (theta1, theta2,...,thetak)
                    theta0,  # valor fixado (simples) para a hipótese nula 
                    n        # tamanho das subpops - deve ser um vetor com k valores (n1,n2,...,nk)
){
  # if (length(ts) != k) stop(glue::glue("Erro: ts deve ter tamanho k = {k}"))
  if (length(n) != k) stop(glue::glue("Erro: n deve ter tamanho k = {k}"))
  
  # Cálculo da dissimilaridade
  thetas <- as.numeric(ts)  # Garante que ts seja numérico
  diff <- n * (theta0 - thetas)
  num <- sum(diff^2)
  den <- sum(n) * theta0 * (1 - theta0)
  
  return(as.numeric(num/den))
  # }  
}


############### REVEEEEER
# #### Kulback-Leiber Distance (KL)
# 
# KL_diss <- function(t1, t2,  # valores dos parâmetros
#                     theta0,  # valor fixado (simples) para a hipótese nula 
#                     n1 = 10, # número de observações da sub-pop 1
#                     n2 = 10  # número de observações da sub-pop 2
# ){
#   parte1 <- n1*t1*log((t1*(1-theta0))/(theta0*(1-t1))) + n2*t2*log((t2*(1-theta0))/(theta0*(1-t2)))
#   parte2 <- n1*log((1-t1)/(1-theta0)) + n2*log((1-t2)/(1-theta0))
#   KL <- parte1 + parte2
#   return(KL)
# }

######################
#### GERANDO GRID
######################

gerar_infos <- function(k,             # quantas populações temos (verificar se os vetores passados condizem com esse tamanho)
                        diss_function, # função que calcula dissimilaridade
                        thetas0,       # valores de thetas0 in (0,1) para testar se está em H0
                        grid_thetas,   # valores de thetas (no espaço amostral) para verificar dissimilaridade 
                        eps,           # epsilon aceitável (para "engordar" H0)
                        n              # VEOTR COM número de observações da sub-pops (n1,...,nk)
){
  
  theta0 = thetas0[1]
  
  # Convertendo diretamente grid_thetas para matriz
  grid_matrix <- as.matrix(grid_thetas)
  
  # Calcula dissimilaridade diretamente em toda a matriz
  diss_values <- apply(grid_matrix, 1, function(theta) diss_function(k, theta, theta0, n))
  
  # Criando o dataframe final com os resultados
  grid <- grid_thetas %>%
    mutate(diss = diss_values,
           prag = diss < eps)
  
  # JEITO ANTIGO
  # # cria uma base em que cada linha é uma lista com k valores vazios (a serem preenchidos)
  # grid <- grid_thetas %>%
  #   mutate(thetas = pmap(., ~ t(as.matrix(c(...))))) 
  # # Calcula dissimilaridade para cada theta em thetas0 de forma vetorizada
  # grid <- grid %>%
  #   mutate(diss = map_dbl(thetas, ~ diss_function(k, as.vector(.x), theta0, n)),
  #          prag = diss < eps)
  
  pb <- txtProgressBar(min = 0, 
                       max = length(thetas0), 
                       style = 3)
  
  # Calcular dissimilaridade para cada theta em thetas0
  grid_theta_prag <- grid 
  for(i in seq_along(thetas0)){ # para cada p a ser testado, calcula dissimilaridade
    t <- thetas0[i]
    
    # Calcula dissimilaridade para cada theta0 
    diss_values <- apply(grid_matrix, 1, function(theta) diss_function(k, theta, t, n))
    
    ## grid para usar no union
    aux <- grid_thetas %>%
      mutate(diss = diss_values,
             prag = diss < eps)
    grid_theta_prag$prag <- grid_theta_prag$prag | aux$prag # faz o UNION das regiões pragmáticas para cada p de H0
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(grid_theta_prag)
}



########################### REVEEEEEEEER!!!!

######################
#### GERANDO O PLOT
######################

plot_prag <- function(k, dados_gerados, thetas, col_3d){
  if(k <= 2){
    borda_sup <- dados_gerados %>% 
      filter(prag) %>% 
      group_by(theta1) %>% 
      summarise(theta2 = max(theta2)) %>% 
      arrange(desc(theta1))
    borda_inf <- dados_gerados %>% 
      filter(prag) %>% 
      group_by(theta1) %>% 
      summarise(theta2 = min(theta2))
    dados_regiao_prag <- borda_inf %>% bind_rows(borda_sup)
    
    fig <- ggplot() + 
      ## Região Pragmática
      geom_polygon(data = dados_regiao_prag,
                   aes(x = theta1, y = theta2),
                   col = 'deeppink3', fill = 'pink', alpha = 0.2) +
      theme_classic() +
      ## Região original de H0
      geom_line(aes(x = thetas, 
                    y = thetas) # indica H0: theta1 = theta2
      ) +
      # ## Ponto de H0 fixado
      # geom_point(aes(x = theta0, y = theta0)) +
      xlab(expression(theta[1])) +
      ylab(expression(theta[2]))
  }
  else if(k == 3){
    library(plotly)
    df = dados_gerados %>% 
      filter(prag) %>% 
      mutate(x = theta1, y = theta2, z = theta3) %>% 
      select(x, y, z)
    fig <- plot_ly(df, x = ~x, y = ~y, z = ~z, 
                   type = "scatter3d", 
                   mode = "markers",
                   marker = list(size = 4,
                                 color = ~z,
                                 colorscale = col_3d,  # Paletas: "Viridis", "Cividis", "Plasma", etc.
                                 opacity = 0.8,
                                 colorbar = list(title = "Z valor"))) %>% 
      layout(
        scene = list(
          xaxis = list(title = "X", range = c(0, 1)),
          yaxis = list(title = "Y", range = c(0, 1)),
          zaxis = list(title = "Z", range = c(0, 1))
        ) 
      ) %>% 
      layout(
        scene = list(
          xaxis = list(title = "Eixo X", showbackground = FALSE),
          yaxis = list(title = "Eixo Y", showbackground = FALSE),
          zaxis = list(title = "Eixo Z", showbackground = FALSE)
        ),
        title = "Gráfico da Região Pragmática",
        margin = list(l = 0, r = 0, b = 0, t = 30)
      ) 
    
  }  
  else{return(print('Não é possível plotar gráfico'))}
  
  return(fig)
  
}

