# //////////////////////////////////////////////////////////////////////////////#
#                                                                               #
#                                                                               #
#     PME3332 - Mecânica dos Fluidos : Noções, Laboratório e Aplicações (2022)  #
#     Primeiro Trabalho                                                         #
#                                                                               #
#     11808575 - Arthur Freitas Campos - arthurfc@usp.br                        #
#     11375225 - Ivan da Cruz Nunes de Moraes - ivandacruz1805@usp.br           #
#     11805593 - Pedro Ruiz Toniato - pedro.toniato@usp.br                      #
#                                                                               #
#                                                                               #
# //////////////////////////////////////////////////////////////////////////////#

#Functions
calcZs <- function(ids, n){
  result <- c()
  
  vectorAux <- c()
  for (i in ids){
    vectorAux <- append(vectorAux, as.double(substr(i, 3, 3)))
  }
  result <- append(result, mean(vectorAux))
  
  vectorAux <- c()
  for (i in ids){
    vectorAux <- append(vectorAux, as.double(substr(i, 4, 4)))
  }
  result <- append(result, sum(vectorAux) + 2)
  
  vectorAux <- c()
  for (i in ids){
    vectorAux <- append(vectorAux, as.double(substr(i, 5, 5)))
  }
  result <- append(result, sum(vectorAux) + 2)
  
  if(n==1){
    return(result[1])
  }
  else if(n==2){
    return(result[1]+result[2])
  }
  else if(n==3){
    return(result[1]+result[2]+result[3])
  }
}

calcHm <- function(ids, z3){
  result <- c()
  vectorAux <- c()
  for (i in ids){
    vectorAux <- append(vectorAux, as.double(substr(i, 6, 6)))
  }
  result <- append(result, sum(vectorAux))
  return(result[1]+z3+10)
}

clbrk <- function(relRoughness, re2){
  x <- (-2*log10(relRoughness/3.7 + 2.51/re2))**(-2)
  return (x)
}

hlnd <- function(relRoughness, re){
  x <- (-1.8*log10((relRoughness/3.7)**1.11 + 6.9/re))**(-2)
  return (x)
}

swJain <- function(relRoughness, re){
  x <- (-2*log10(relRoughness/3.7 + 5.74/(re**0.9)))**(-2)
  return (x)
}

converge <- function(solution, hj, z1, z2, z3, d, hm, la, lb, lc, leq, vsc, relRoughness, g){
  fcVc2 <- ((hj-z3)*2*g*d)/(lc+leq)
  fbVb2 <- ((hj-z2)*2*g*d)/lb
  re2_c <- (d/vsc)*sqrt(fcVc2)
  re2_b <- (d/vsc)*sqrt(fbVb2)
  fc <- clbrk(relRoughness, re2_c)
  fb <- clbrk(relRoughness, re2_b)
  Vc <- sqrt(fcVc2/fc)
  Vb <- sqrt(fbVb2/fb)
  Va <- Vc+Vb
  re_a <- Va*d/vsc
  fa <- hlnd(relRoughness, re_a)
  new_Hj <- z1+hm-fa*(la/d)*((Va**2)/(2*g))
  Qa <- Va*pi*d**2/4
  Qb <- Vb*pi*d**2/4
  Qc <- Vc*pi*d**2/4
  result <- c()
  result <- append(result, hj)
  result <- append(result, fa)
  result <- append(result, fb)
  result <- append(result, fc)
  result <- append(result, Qa)
  result <- append(result, Qb)
  result <- append(result, Qc)
  result <- append(result, new_Hj)
  solution[nrow(solution) + 1,] <- result
  
  return(solution)
}


  
#Data
d <- 5e-2
relRoughness <- 0.001
la <- 50
lb <- 100
lc <- 200
leq <- 50
rho <- 1000
vsc <- 1e-6
g <- 10
ids <- c("11808575", "11375225", "11805593")
z1 <- calcZs(ids, 1)
z2 <- calcZs(ids, 2)
z3 <- calcZs(ids, 3)
hm <- calcHm(ids, z3)

#solving
hj <- c(0)
fa <- c(0)
fb <- c(0)
fc <- c(0)
Qa <- c(0)
Qb <- c(0)
Qc <- c(0)
new_Hj <- c(0)
solution <-data.frame(hj, fa, fb, fc, Qa, Qb, Qc, new_Hj)
solution <- converge(solution, 41.8573275, z1, z2, z3, d, hm, la, lb, lc, leq, vsc, relRoughness, g)



