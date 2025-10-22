library('tidyverse')
library('combinat')
library('glue')


###########################################################
######################### FUNÇÕES #########################
###########################################################
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
nrow(grid_thetas)

# teste
# grid_thetas <- grid_thetas %>%
#   mutate(theta3 = 0.5) %>%
#   distinct()
# grid_thetas

# names(grid_thetas) <- c('theta1','theta2' )

###########################################################
######################## APLICANDO ########################
###########################################################

#########################################
##################### theta0 FIXADO
#########################################
theta0 <- 0.5

start <- Sys.time()
dados_gerados_BP_simplex <- gerar_infos(k = k,
                                        BP_diss,
                                        thetas0 = theta0,
                                        grid_thetas,
                                        eps = eps,
                                        n = n)
end <- Sys.time()
end - start

saveRDS(dados_gerados_BP_simplex,
        glue("{path}/Teste Homogeneidade kx2/Dados/dados_gerados_BP_simples_k_{k}_theta0_{theta0}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))

dados_gerados_BP_simplex %>% 
  filter(prag) %>% 
  filter(round(theta1, 3) == 0.455,
         round(theta2, 3) == 0.444,
         round(theta3, 3) == 0.475)
# tem que resultar 0.0772




# dados_gerados_KL_simplex <- gerar_infos(KL_diss,
#                                         thetas0 = theta0,
#                                         grid_thetas,
#                                         eps = eps,
#                                         n1 = n1,
#                                         n2 = n2)
# saveRDS(dados_gerados_KL_simplex,
#         glue("{path}/Teste Homogeneidade 2x2/Dados/dados_gerados_KL_simples_theta0_{theta0}_n1_{n1}_n2_{n2}_eps_{eps}.rds"))
# 


#########################################
##################### theta0 in (0,1)
#########################################
thetas <- seq(0.00000001, 0.99999999, length.out = N) # grid de thetas

start <- Sys.time()
dados_gerados_BP_full <- gerar_infos(k = k,
                                     BP_diss,
                                     thetas0 = thetas,
                                     grid_thetas,
                                     eps = eps,
                                     n = n)
end <- Sys.time()
end - start

saveRDS(dados_gerados_BP_full,
        glue("{path}/Teste Homogeneidade kx2/Dados/dados_gerados_BP_full_k_{k}_n_{paste(n, collapse = '_')}_eps_{eps}.rds"))


# dados_gerados_KL_full <- gerar_infos(KL_diss,
#                                      thetas0 = thetas,
#                                      grid_thetas,
#                                      eps = eps,
#                                      n1 = n1,
#                                      n2 = n2)
# saveRDS(dados_gerados_KL_full,
#         glue("{path}/Teste Homogeneidade 2x2/Dados/dados_gerados_KL_full_n1_{n1}_n2_{n2}_eps_{eps}.rds"))
# 


