###############################################################################
################ FUNÇÃO PARA GERAR DADOS DA REGIÃO PRAGMÁTICA #################
###############################################################################

gerar_dados_prag <- function(k    # número de populações
                            ,n   # vetor de tamanho k com o tamanho das binomiais que gerar a pragmática
                            ,name_dist = 'BP' # função de distância preditiva
                            ,eps # parâmetro que define tamanho da distância preditiva "aceitável"
                            ,N   # tamanho do reticulado (em um eixo) - que será usado para desenhar a pragmática
                            ,path = '~/'# caminho raiz para salvar os dados gerados
                            ){
  if(name_dist == 'BP'){diss_function = BP_diss}
  
  ######################################
  ###### Gera thetas (reticulado) ######
  ######################################
  thetas <- seq(0.00000001, 0.99999999, length.out = N) # thetas no reticulado em um eixo
  ## criar grid com todas as combinações de k thetas
  aux <- expand_grid(thetas, thetas)
  if(k>2){
    for(i in 3:(k)){
      aux <- expand_grid(aux, thetas)
      names(aux) <- glue('theta{1:i}')
    }
  }else{
    names(aux) <- glue('theta{1:k}')
  }
  grid_thetas <- aux
  
  ######################################
  ############# Gera dados #############
  ######################################
  dados_gerados <- gerar_infos(k = k,
                               diss_function,
                               thetas0 = thetas,
                               grid_thetas,
                               eps = eps,
                               n = n)
  
  ######################################
  ###### Salvando dados e Gráficos #####
  ###################################### 
  # caminho_graf <- glue("{path}/DadosAplicacao/Test_Homog_{k}x2/{name_dist}/")
  # if (!dir.exists(caminho_graf)) {
  #   dir.create(caminho_graf, recursive = TRUE)
  # }
  caminho_graf <- path
  if(k == 2){
    fig <- plot_prag(k, dados_gerados, thetas) + 
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 14)) +
      geom_text(aes(x = 0.99, y = 0.9, label = glue(" t1 = t2"))) 
    ggsave(glue("{caminho_graf}/Graf_PRAG_TesteHomog_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.png"),
           width = 5,
           height = 5,
           plot = fig)
  }
  if(k > 2){
    fig <- plot_prag(k, dados_gerados, thetas, col_3d = 'Viridis')
    htmlwidgets::saveWidget(fig, glue("{caminho_graf}/Graf_PRAG_TesteHomog_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.html"))
  }
  saveRDS(dados_gerados,
          glue("{caminho_graf}/dados_gerados_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))
  
  return(list(dados_gerados, fig))  
}


###############################################################################
###### FUNÇÃO PARA APLICAR HPD - GRÁFICO E TESTE (Script 3 da Simulação) ######
###############################################################################
gerar_hpd <- function(k    # número de populações
                      ,n   # vetor de tamanho k com o tamanho das binomiais que gerar a pragmática
                      ,confianca = 0.8 # parâmetro para confiança do HPD
                      ,prior_par_list  # lista com parâmetros das prioris Beta (k pares (a_i, b_i))
                      ,m     # vetor de tamanho k com tamanho da binomial da função de verossimilhança
                      ,dados # base de dados que se deseja aplicar o HPD
                      ,name_dist = 'BP' # função de distância preditiva
                      ,eps # parâmetro que define tamanho da distância preditiva "aceitável"
                      ,N   # tamanho do reticulado (em um eixo) - que será usado para desenhar a pragmática
                      ,path = '~/'# caminho raiz para salvar os dados gerados
                      ){
  
  library('geometry')
  library('rgl')
  library('sp')
  library('sf')
  ######################################
  ###### Gera thetas (reticulado) ######
  ######################################
  thetas <- seq(0.00000001, 0.99999999, length.out = 5000) # thetas no reticulado em um eixo
  ## criar grid com todas as combinações de k thetas
  aux <- expand_grid(thetas, thetas)
  if(k>2){
    for(i in 3:(k)){
      aux <- expand_grid(aux, thetas)
      names(aux) <- glue('theta{1:i}')
    }
  }else{
    names(aux) <- glue('theta{1:k}')
  }
  grid_thetas <- aux
  
  #######################################
  ###### Gerar dados da pragmática ######
  #######################################
  retorno_prag = gerar_dados_prag(k = k, n = n,
                                  name_dist = name_dist,
                                  N = N, eps = eps, path = path)
  dados_gerados <- retorno_prag[[1]]
  fig_prag <- retorno_prag[[2]]
  
  
  ####################################################
  ###### Gera região HPD a partir do reticulado ######
  ####################################################
  hpd_labels <- tibble(id = as.double(1:nrow(dados)))
  for (j in 1:k) {
    var_name <- paste0("x_", j)
    hpd_labels[[var_name]] <- as.double(rep(NA, nrow(dados)))
  }
  hpd_grid = NULL
  for(ii in 1:nrow(dados)){
    counts <- dados[ii,]
    # Geração dos parâmetros posteriores para cada beta
    posterior_par_list <- lapply(1:k, function(j) {
      posterior(dados = as.integer(counts[, j]), 
                prior_par = prior_par_list[[j]],
                m = m[j]) %>% as.numeric()
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
      mutate(log_p = apply(thetas_row, 1, function(t_vec) {
        log_posterior(as.numeric(t_vec), posterior_par_list)
      }),
      color = apply(thetas_row, 1, function(t_vec) {
        log_posterior(as.numeric(t_vec), posterior_par_list) > th
      })) 
    
    reg_hpd <- this_hpd_grid %>% 
      mutate(idx = ii) %>%
      filter(color)
    hpd_grid <- hpd_grid %>% bind_rows(reg_hpd)
    
    # extrai e normaliza os valores desejados (por padrão o primeiro de cada vetor)
    valores_normalizados <- sapply(posterior_par_list, function(p) {
      (p / sum(p))[1]
    })
    # salva no data.frame/matriz, incluindo o índice ii
    hpd_labels[ii, ] <- t(c(ii, valores_normalizados))
    ## print(ii)
  }

  hpd_grid <- hpd_grid %>% mutate(idx = as.factor(idx))
  
  ajustar código

  ################################
  ###### Gráfico HPD + Prag ######
  ################################
  caminho_graf <- glue("{path}/DadosAplicacao/Test_Homog_{k}x2/{name_dist}/")
  if (!dir.exists(caminho_graf)) {
    dir.create(caminho_graf, recursive = TRUE)
  }
  
  fig = add_hpd_grid(k
                     ,fig_prag
                     ,hpd_grid
                     ,hpd_labels
                     ,my_color = "green4"
                     ,alph=0.2)
  if(k == 2){
    ggsave(glue("{caminho_graf}/Graf_HPD_TesteHomog_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.png"),
           width = 5,
           height = 5,
           plot = fig)
  }
  else if(k == 3){
    htmlwidgets::saveWidget(fig, glue("{caminho_graf}/Graf_HPD_TesteHomog_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.html"))
  }
  else{print(glue("Não é possível plotar gráfico - k > 3"))}

  ####################################################
  ###### Verifica se  HPD está contido na  Prag ######
  ####################################################
  dados_regiao_prag <- dados_gerados %>%
    filter(prag) %>% 
    select(starts_with("theta"))
  teste <- NULL
  for(i in 1:nrow(dados)){
    # Seus dataframes
    df1 <- dados_regiao_prag %>% 
      bind_rows(dados_regiao_prag[1,]) # Polígono 1 (externo) - região pragmática
    df2 <- hpd_grid %>% filter(idx == i) %>% 
      bind_rows((hpd_grid %>% filter(idx == i))[1,])           # Polígono 2 (interno) - HPD
    df3 <- grid_thetas %>% bind_rows(grid_thetas[1,])          # ESPAÇO PARAMÉTRICO COMPLETO
    
    
    # --- identifica dinamicamente as colunas theta ---
    theta_cols1 <- grep("^theta", names(df1), value = TRUE)
    theta_cols2 <- grep("^theta", names(df2), value = TRUE)
    
    ##### Criar os dois polígonos (prag e hpd)
    prag <- st_polygon(list(as.matrix(df1[, theta_cols1]))) %>%
      st_sfc(dim = paste0("XYZ"[seq_along(theta_cols1)])) %>%  # define dimensão (XYZ, XYM etc)
      st_sf()
    
    hpd <- st_polygon(list(as.matrix(df2[, theta_cols2]))) %>%
      st_sfc(dim = paste0("XYZ"[seq_along(theta_cols2)])) %>%
      st_sf()
    prag_convex <- st_convex_hull(prag)   # Calcular a envoltória convexa de ambos os polígonos
    hpd_convex <- st_convex_hull(hpd)   # Calcular a envoltória convexa de ambos os polígonos
    
    
    ##### CRIAR A "CASCA" DOS POLÍGONOS
    # 1. Calcula o convex hull dos pontos de prag
    coords_prag <- as.matrix(df1[, theta_cols1])
    coords_hpd  <- as.matrix(df2[, theta_cols2])
    
    hull_prag <- convhulln(coords_prag, options = "FA")
    hull_faces_prag <- hull_prag$hull
    
    hull_hpd <- convhulln(coords_hpd, options = "FA")
    hull_faces_hpd <- hull_hpd$hull
    
    #### Interseções entre a pragmática e o HPD
    intersecta <- st_intersects(prag_convex, hpd_convex, sparse = FALSE) # Verifica se há interseção entre os polígonos
    contidos <- geometry::inhulln(hull_prag, p = coords_hpd) # verifica se TODOS os pontos do HPD estão dentro da região convexa pragmática
    parcial <- intersecta & !all(contidos)  # Interseção mas não totalmente dentro
    
    # Exibe o resultado
    if (all(contidos)) {
      # print(glue("O HPD {i} está contido na região pragmática (aceitação)."))
      aux <- 0
    } else if(parcial){
      # print("Indefinido")
      aux <- 1/2
    } else {
      # print(glue("O HPD {i} não está contido na região pragmática (rejeição)."))
      aux <- 1
    }
    # print(aux)
    teste[i] <- aux
  }
  
  dados_result <- dados %>% 
    mutate(teste = teste) %>% 
    mutate(result = ifelse(teste == 1, 'Rejeita H0',
                           ifelse(teste == 1/2, 'Agnóstico', 
                                  ifelse(teste == 0, 'Aceita H0', NA))))
  
  return(list(resultado = dados_result, 
              hpd_grid = hpd_grid,
              grafico_hpd = fig,
              caminho_graf = glue("{caminho_graf}/Graf_HPD_TesteHomog_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.html"),
              parametros_entrada = list(k=k, 
                                        n = n,
                                        confianca = confianca,
                                        prior_par_list = prior_par_list,
                                        dados = dados,
                                        m = m,
                                        name_dist = name_dist,
                                        eps = eps,
                                        N = N,
                                        path = path)))  
  
}

###############################################################################
######### FUNÇÃO PARA CALCULAR PROBS PRAG+HPD (Script 4 da Simulação) #########
###############################################################################
estima_propb_prag_hpd <- function(k    # número de populações
                                  ,n   # vetor de tamanho k com o tamanho das binomiais que gerar a pragmática
                                  ,confianca = 0.8 # parâmetro para confiança do HPD
                                  ,prior_par_list  # lista com parâmetros das prioris Beta (k pares (a_i, b_i))

                                  ,dados_gerados_prag # dados gerados da região pragmática
                                  
                                  ,m     # vetor de tamanho k com tamanho da binomial da função de verossimilhança
                                  ,dados_full # base de dados que se deseja aplicar o HPD
                                  ,N     # número de vezes que vou gerar da posteriori
                                  ,eps   #
                                  ,path = '~/'# caminho raiz para salvar os dados gerados
){
  library('coda')
  library('geometry')
  
  p_intersec_full <- NULL
  prob_hat_prag_full <- NULL
  for(ii in 1:nrow(dados_full)){ # repete o procedimento para cada linha de dados

    #############################################################################################
    ########### Passo 1: gerar pontos da posteriori
    #############################################################################################
    a <- sapply(prior_par_list, `[`, 1)  # primeiros elementos
    b <- sapply(prior_par_list, `[`, 2)  # segundos elementos
    dados <- dados_full[ii,]
    thetas_gerados <- purrr::map(1:length(m), ~{
      theta_post_gerado <- rbeta(N, a[.x] + as.numeric(dados[,.x]),
                                 m[.x] + b[.x] - as.numeric(dados[,.x]))
      return(theta_post_gerado)
    })
    
    thetas <- tibble(theta1_gerado = thetas_gerados[[1]])
    for(i in 2:length(m)){ # colunando os "k" thetas_gerados
      thetas <- thetas %>% 
        mutate(thetas_gerados_aux = thetas_gerados[[i]])
      names(thetas)[i] <- glue('theta{i}_gerado')
    }

    
    #############################################################################################
    ########### Passo 2: Calculando corte do HPD
    #############################################################################################
    
    counts <- dados
    posterior_par_list <- lapply(1:k, function(j) {
      posterior(dados = as.integer(counts[, j]), 
                prior_par = prior_par_list[[j]],
                m = m[j]) %>% as.numeric()
    })
    th <- hpd_post_cut(k = k,
                       posterior_par_list = posterior_par_list,
                       probability = confianca)
    
    
    #############################################################################################
    ########### Passo 3: VERIFICANDO PERTENCIMENTO DO PONTO GERADO NA PRAGMÁTICA ESPECÍFICA 
    #############################################################################################
    
    # Obs: theta0_linha  =  (theta1_gerado + theta2_gerado + ... + thetak_gerado) / k
    
    # ATENÇÃAAAAO: supondo n1=n1=...=nk - definido no início do problema
    dados_thetas_verificar <- thetas %>% # thetas gerados da posteriori
      mutate(sum_thetas = select(., starts_with("theta")) %>% rowSums(.,na.rm = TRUE)) %>% 
      mutate(th = th) %>% 
      mutate(theta0_linha = sum_thetas/ncol(thetas),
             delta_estrela = eps/(eps+n[1]) * (1/2 - theta0_linha),
             theta_verificar = theta0_linha + delta_estrela) # theta do centro da bolinha específica que eu tenho que verificar se o ponto gerado pertence
    
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
    
    ###### (ii) Verificar se está dentro do raio
    dados_thetas_verificar <- dados_thetas_verificar %>%
      mutate(across(ends_with('gerado'), ~ (.x-theta_verificar)^2, 
                    .names = '{.col}_quadrado_menos_theta_vertificar')) %>% 
      mutate(sum_thetas_quadrado_menos_theta_vertificar = select(., ends_with('_quadrado_menos_theta_vertificar')) %>% 
               rowSums(.,na.rm = TRUE)) %>% 
      mutate(raiz_sum_thetas_quadrado_menos_theta_vertificar = sqrt(sum_thetas_quadrado_menos_theta_vertificar) <= raio) %>% 
      mutate(indc_in_prag = ifelse(indc_dentro,
                                   ifelse(sqrt(sum_thetas_quadrado_menos_theta_vertificar) <= raio,1,0),0)) 
    
    #############################################################################################
    ########### Passo 4: ESTIMANDO P(pragmática) (a posteriori)
    #############################################################################################
    
    n_pontos_in_prag <- dados_thetas_verificar %>% filter(indc_in_prag==1) %>% nrow()
    prob_hat_prag <- n_pontos_in_prag/nrow(dados_thetas_verificar)
    
    
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
    
    
    
    #############################################################################################
    ########### Passo 6: ESTIMAR Prob da INTERSEÇÃO DE HPD e PRAG 
    #############################################################################################
    

    # quero ver a intersecção dos conjuntos (que formam o HPD e a região Pragmática)
    # e calcular a prob dele
    dados_regiao_prag <- dados_gerados_prag %>%
      filter(prag) %>% 
      select(starts_with('theta'))
    # Conjunto 1: pontos que formam a região pragmática em 3D
    # (colunas theta1, ..., thetak)
    regiao <- as.matrix(dados_regiao_prag)
    # Conjunto 2: pontos a verificar
    pontos <- as.matrix(conjunto_hpd_estimado %>% select(ends_with('_gerado')))
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
    
    
    
    p_intersec_full[ii] <- p_intersec
    prob_hat_prag_full[ii] <- prob_hat_prag
    
  }
  dados_resultado_probs <- dados_full %>% 
    mutate(p_intersec_prag_hpd = p_intersec_full,
           prob_hat_prag = prob_hat_prag_full)
  
  return(dados_resultado_probs)
}



# name_dist = 'BP'
# k = 3
# start <- Sys.time()
# result_hpd <- gerar_hpd(k=k, 
#                         n = c(10,10,10),
#                         confianca = 0.8,
#                         prior_par_list = list(c(1, 1), c(1, 1), c(1, 1)),
#                         dados = as.tibble(t(matrix(c(45, 65, 100)))),
#                         m = c(50, 100, 200),
#                         name_dist = name_dist,
#                         eps = 0.1,
#                         N = 100,
#                         path = path)
# 
# end <- Sys.time()
# end-start
# ## Funcionou, mas demora 18 min pra rodar
# 
# 
# result_hpd$resultado
# browseURL(result_hpd$caminho_graf)  
# 
# 
# saveRDS(result_hpd, glue("{path}/DadosAplicacao/Test_Homog_{k}x2/{name_dist}/resultado_aplicacao_full.rds"))



# k=3
# n=c(10,10,10)
# confianca = 0.8 
# prior_par_list = list(c(1, 1), c(1, 1), c(1, 1))
# name_dist = "BP"
# eps = 0.1
# caminho_graf = glue("{path}/DadosAplicacao/Test_Homog_{k}x2/{name_dist}/")
# dados_gerados_prag = read_rds(glue("{caminho_graf}/dados_gerados_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))
# m = c(50, 100, 200)     # vetor de tamanho k com tamanho da binomial da função de verossimilhança
# dados_full = as.tibble(t(matrix(c(45, 95, 180)))) # base de dados que se deseja aplicar o HPD
# 
# start <- Sys.time()  
# estima_propb_prag_hpd(k
#                       ,n
#                       ,confianca
#                       ,prior_par_list
#                       ,dados_gerados_prag
#                       ,m
#                       ,dados_full
#                       ,N = 1000
#                       ,eps = 0.1)
#   
# end <- Sys.time()
# end-start  

  