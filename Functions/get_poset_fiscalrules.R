get_poset_fiscalrules<-function(inputfile,index_variables){
  
  profiles_s <- inputfile %>%
    select(names(index_variables))%>%
    pop2prof()
  
  rmProfiles2 <- function(y, v, ...) { # Corrected function based on rmProfiles in parsec
    v <- which(v)
    y$profiles <- y$profiles[-v, ]
    y$freq <- y$freq[-v]
    return(y)
  }
  
  #must include leading 0
  #0 means no fiscal rule in place
  profiles_f <- var2prof(
    index_variables
  )
  
  profiles_f <- profiles_f %>% rmProfiles2(grepl("0", rownames(profiles_f$profiles)) & !grepl("0000", rownames(profiles_f$profiles))) # dropping all profiles with '0' in it, keeping '0000'
  
  
  poset_j <- merge(profiles_s, profiles_f, support = T)
  
  eval <- evaluation(poset_j, threshold = "0000", error = 10^(-3))
  return(eval)
}