
decay_instant <- function(t, hl, C_max){
     C_max * 2^(-t/hl)
 }

decay_instant_integral <- function(t, hl, C_max){
    -hl * C_max * 2^(-t/hl)/log(2)
 }

return_druginf <- function(){

    druginf <- list()

    druginf[['3TC']] <-
        c(IC_50 = 0.0298, slope = 1.15, c_max = 15.3, hl = 10, dosing = 2)
    druginf[['AZT']] <-
        c(IC_50 = 0.1823, slope = 0.85, c_max = 4.5, hl = 8.5, dosing = 2)
    druginf[['D4T']] <-
        c(IC_50 = 0.5524, slope = 1.13, c_max = 2.3, hl = 3.5, dosing = 2)
    druginf[['EFV']] <-
        c(IC_50 = 0.0035, slope = 1.69, c_max = 12.9, hl = 35.8, dosing = 1)
    druginf[['NFV']] <-
        c(IC_50 = 0.2360, slope = 1.88, c_max = 5.1, hl = 4.0, dosing = 3)
    druginf[['LPV']] <-
        c(IC_50 = 0.0380, slope = 2.1, c_max = 15.6, hl = 9.9, dosing = 2)
    druginf[['NVP']] <-
        c(IC_50 = 0.0490, slope = 1.49, c_max = 25.2, hl = 21.5, dosing = 1)
    
    return(druginf)
}


return_mutinf <- function(){
    mutinf <- list()


    #For 3TC, mutation options are: 
    #K65R, M184V
    mutinf[['3TC']] <-
        c(
            #s is taken from K65R
            s = 0.41, 
            #mu is equivalent 
            mu =  1.1 * (10^-5), 
            #rho and sigma are from M184V
            rho =  963, 
            sigma = -0.58)

    #For AZT, mutation options are: 
    #M184V, M41L, T215Y* (excluded, because multi-nuc)
    mutinf[['AZT']] <-
        c(
            #s is taken from M41L
            s = 0.17,
            #mu is taken from M184V
            mu = 1.1 * 10^(-5),
            #rho and sigma are taken from M41L
            rho = 2.2,
            sigma = 0.07)
    
    #For D4T, mutation options are: 
    #M41L, T215Y* (excluded, because multi-nuc)
    #all parameters taken from M41L
    mutinf[['D4T']] <-
        c(s = 0.17, 
          mu = 1.3 * 10^(-6),
#          rho = 1.0,
#          sigma = 0.07)
          rho = 1.08,
          sigma = -0.12)

    #For EFV, mutation options are:
    #G190S, K103N, Y181C
    mutinf[['EFV']] <-
        #s is taken from Y181C
        c(s = 0.26, 
          #mu is taken from G190S
          mu = 2.2 * 10^(-5), 
          #rho is taken from K103N
          rho = 85, 
          sigma = -0.17)

    #For NFV, mutation options are
    #D30N, L90M, 
    mutinf[['NFV']] <-
        #s is taken from D30N
        c(s = 0.27, 
          #mut is taken from D30N
          mu = 5.5 * 10^(-5), 
          #rho and sigma are taken from D30N
          rho = 2.3, 
          sigma = -0.29)

    #For LPV/r, mutation options are
    #I47A*, I47V, V32I, V82A, V82F,  
    mutinf[['LPV']] <-
        #s is taken from I47V
        c(s = 0.05, 
          #mut is taken from V32I
          mu = 4.1 * 10^(-5), 
          #rho and sigma are taken from I47V
          rho = 1.8, 
          sigma = -0.29)

    #For NVP, mutation options are
    #G190S, K103N, Y181C, #Y181I
    mutinf[['NVP']] <-
        #s is taken from K103N
        c(s = 0.3, 
          #mut is taken from G190S
          mu = 2.2 * 10^(-5), 
          #rho and sigma are taken from K103N
          rho = 94, 
          sigma = -0.15)

    return(mutinf)
    
}


