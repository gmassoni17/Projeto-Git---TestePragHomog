######################
#### FUNÇÕES HPD
######################


# distribuições a priori:
# theta1 ~ Beta(a,b)
# theta2 ~ Beta(c,d)
# funcão calcula parâmetros das betas isoladas
posterior <- function(dados                 ## valores obervados de x ou y
                      ,prior_par = c(1, 2)  ## valores "iniciais" de a e b (ou c e d)
                      ,m               ## parametros "conhecidos" da Binomial ---- m = (m1,m2)
){
  aux <- tibble(a = prior_par[1] + dados, b = prior_par[2] + m - dados) # Posteriori é Beta(a + x, b + m - x)
  return(aux)
}

## A função a posteriori é uma multiplicação de betas
log_posterior = function(ts, # vetor com k valores de thetas
                         posterior_pars # lista com k pares de valores da priori (para calcular posteriori)
                         )
{
  if (length(ts) != length(posterior_pars)) {
    stop("O comprimento de ts deve ser igual ao comprimento de posterior_pars")
  }
  
  # Calcula a densidade beta para cada t e seus respectivos parâmetros
  beta_vals <- mapply(function(t, pars) {
    dbeta(t, pars[1], pars[2])
  }, ts, posterior_pars)
  
  return(log(prod(beta_vals)))
}


###############################################################################
################# Função de densidade posterior para \theta_i #################
###############################################################################

theta_posterior <- function(theta, X, n, a, b) {
  if (theta < 0 || theta > 1) {
    return(0)  # Retorna 0 se theta proposto estiver fora do intervalo [0, 1]
  }
  # Densidade posterior (não normalizada)
  return(dbeta(theta, a + X, b + n - X)) # beta da com parâmetros já levando em conta os dados
}



# The HPD is a region with parameters
# such that the posterior is greater than v.
# This function determines v such that
# the HPD has posterior probability of 'probability'.

hpd_post_cut = function(k,                  #
                        posterior_par_list, # lista com k vetores (a, b) para cada beta
                        probability = 0.95, # nível de confiança
                        B = 10^5            # tamanho da amostra gerada
                        )
{
  
  if(k != length(posterior_par_list)) {
    stop(glue("O comprimento de posterior_par_list deve ser igual a k = {k}"))
  }
  
  # Passo 1: gerar amostras para cada beta
  amostras <- lapply(posterior_par_list, function(par) {# lapply aplica a função a cada elemento de posterior_par_list
    rbeta(n = B, shape1 = par[1], shape2 = par[2])
  })
  
  # Passo 2: calcular densidade para cada conjunto de amostras
  densidades <- mapply(function(amostra, par) {
    dbeta(amostra, par[1], par[2])
    }, amostras, posterior_par_list)
  
  # Passo 3: multiplicar as densidades linha a linha e calcular log
  # densidades é uma matriz: cada coluna é uma variável, cada linha um B
  log_prod_densidades <- log(apply(densidades, 1, prod))
  
  # Passo 4: Calculo quantil "probability"% da multiplicação das densidades geradas
  quantil <- as.numeric(quantile(log_prod_densidades, 
                                 probs = 1 - probability))
  
  return(quantil)
}



## ADICIONAR GRÁFICO DO HPD
add_hpd_grid = function(k
                        ,plot
                        ,hpd_grid
                        ,hpd_labels
                        ,my_color = "green4" 
                        ,alph = 0.2)
{
  all_idx = unique(hpd_grid$idx)
  new_plot = plot
  if(k<=2){
    for(this_idx in all_idx){
      # print(this_idx)
      grid <- hpd_grid %>% filter(idx == this_idx)
     
      new_plot <- new_plot +
        ggplot2::geom_polygon(aes(x = theta1, y = theta2), data = grid,
                     color = my_color, fill = my_color, alpha = alph)
 
    }
    x <- hpd_grid %>% group_by(idx) %>% summarise(mean_t1 = mean(theta1))
    y <- hpd_grid %>% group_by(idx) %>% summarise(mean_t2 = mean(theta2))
    labels <- left_join(x, y, by = 'idx')
    return(
      new_plot +
        geom_text(aes(x = mean_t1, y = mean_t2, label = idx), data = labels)
    )
  }
  else if(k == 3){
    for(this_idx in all_idx){
        # print(this_idx)
        grid <- hpd_grid %>% filter(idx == this_idx)
        
        df <- grid %>% 
          mutate(x = theta1, y = theta2, z = theta3) %>% 
          select(x, y, z)
        hull_faces <- convhulln(as.matrix(df[, c("x", "y", "z")]))  # Retorna região convexa
        new_plot <- new_plot %>% 
          # add_trace(
          #   type = "mesh3d",
          #   x = as.matrix(df[, c("x")]),
          #   y = as.matrix(df[, c("y")]),
          #   z = as.matrix(df[, c("z")]),
          #   i = hull_faces[,1] - 1,
          #   j = hull_faces[,2] - 1,
          #   k = hull_faces[,3] - 1,
          #   facecolor = rep("green4", nrow(hull_faces)),
          #   opacity = 0.5,
          #   name = "Casca Convexa"
          # )
          add_trace(data = df,
                    x = ~x, y = ~y, z = ~z,
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 4,
                                  color = ~z,
                                  colorscale = my_color,
                                  opacity = 0.8,
                                  colorbar = list(title = "Z valor")))
    }
    x <- hpd_grid %>% group_by(idx) %>% summarise(mean_t1 = mean(theta1))
    y <- hpd_grid %>% group_by(idx) %>% summarise(mean_t2 = mean(theta2))
    z <- hpd_grid %>% group_by(idx) %>% summarise(mean_t3 = mean(theta3))
    labels <- left_join(x, y, z, by = 'idx')
    return(new_plot %>% 
            add_trace(
              data = labels,
              x = ~x,
              y = ~y,
              z = ~z,
              type = "scatter",
              mode = "text",
              text = ~idx,
              textposition = "top right",
              showlegend = FALSE
              )
    )
    
  }else{stop('Não é possível plotar gráfico')}
  
  
}
