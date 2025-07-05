function pChoose = pChoose_betaKappa(pScene_1, pScene_2, beta, kappa)
       
    pChoose = kappa * 0.5 + (1-kappa) * ((1 + exp(-beta * log(pScene_1 / pScene_2)))^(-1));

end