#//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#//%%%%%%%%%%%  w0 and p0_z are RANDOM
#//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#SS <- S; TT <- T; NN_Tr <- N_Tr.tmp; nn_Tr <- n_Tr.tmp; NN_Tr_sum <- N_Tr_sum; BB_Tr <- B_Tr.tmp; hpara <- hyper; min_C <- Min_C; max_C <- Max_C; N_sam <- burn.in;
#SS <- S; TT <- T; NN_Tr <- N.tmp; nn_Tr <- n.tmp; NN_Tr_sum <- N_sum; BB_Tr <- rep(1, SS*TT); hpara <- hyper; min_C <- Min_C; max_C <- Max_C; N_sam <- burn.in;
fn.Tr.BayClone2 <- function(hpara, SS, TT, ind_B, MM_B, NN_Tr, nn_Tr, NN_Tr_sum, BB_Tr, min_C, max_C, N_sam)
{
    #HYPERPARAMETERS
    a <- hpara$a; b <- hpara$b;  #PHI
    a_w <- hpara$a_w; b_w <- hpara$b_w ## W^\star
    d0 <- hpara$d0; d <- hpara$d  #W
    
    #PRIOR FOR P_zO
    a_z0 <- hpara$a_z0; b_z0 <- hpara$b_z0 # FOR P_Z0
    
    #FOR L
    Q <- hpara$Q;
    alpha <- hpara$alpha
    beta <- hpara$beta
    gam <- hpara$gam
    
    ##FOR m
    kappa <- hpara$kappa
    
    #SAVE THE PROPOSALS FOR ALL C VALUES
    sam.all <- NULL
    
    sam.all$phi <- NULL
    sam.all$pi <- NULL
    
    sam.all$L <- NULL
    sam.all$Z <- NULL
    sam.all$A <- NULL
    
    sam.all$p0_z <- rep(NA, max_C-min_C+1)
    sam.all$w_star <- NULL
    sam.all$th <- NULL
    sam.all$w_tilde <- NULL
    sam.all$w <- NULL
    
    sam.all$m <- NULL
    sam.all$M <- NULL
    sam.all$p <- NULL

    #INITIALIZATION FOR EACH C
    for(i.c in min_C:max_C)
    {
        print(date())
        print(paste("i.c", i.c))
        set.seed(545754)
        
        ##L AND Z
        Z <- L <- array(0, dim=c(SS, i.c))
        
         #C--TO INITIALIZE
        p.ave <- apply(array(nn_Tr/NN_Tr, dim=c(SS, TT)), 1, mean)
        p.quant <- quantile(p.ave, probs=seq(0.10, 0.95, length.out=(i.c-1)))
        
        Z[,1] <- L[,1] <- 2 #CELL TYPE 0
        
        for(ii.c in 2:i.c) #EXCEPT CELL TYPE 0
        {
            Z[p.ave > p.quant[ii.c-1],ii.c] <- 2
        }
        
        ###CORRECTED MARCH-30TH
        for(ii.c in 2:i.c)
        for(i.s in 1:SS)
        {
            if(Z[i.s, ii.c]==Q)
            {
                L[i.s, ii.c] <- Q
            }else{
                L[i.s, ii.c] <- sample((Z[i.s, ii.c]:Q), 1, replace=FALSE)  ###If x has length 1, is numeric (in the sense of is.numeric) and x >= 1, sampling via sample takes place from 1:x. Note that this convenience feature may lead to undesired behaviour when x is of varying length in calls such as sample(x).
            }
        }
        
        ###  obey BB output if BB output is integer
        ## if (ind_B[i_s] == 1) {  // if BB output is not integer
        for(ii.c in 2:i.c)
        {
            L[ind_B==0, ii.c] <- (MM_B[1:SS])[ind_B==0]
        }
        
        for(ii.c in 2:i.c)
        {
            S.set <- (1:SS)[L[,ii.c] < Z[,ii.c]]
            
            for(ii.s in S.set)
            {
                Z[ii.s, ii.c] <- L[ii.s, ii.c]
            }
        }
    
        #COMPUTE A
        A <- array(NA, dim=c(i.c, (Q+1)))
        for(i.q in (0:Q))
        {
            A[,(i.q+1)] <- apply(L==i.q, 2, sum)
        }
        
        #PI
        ppi <- (A+5)/(SS + (Q+1)*5)
        
        #P0_Z
        p_z0 <- c(0.005, rep(1, i.c-1))

        #PHI
        phi <- rgamma(TT, a, b)
        
        # W^\STAR
        w_star <- rbeta(TT, a_w, b_w)
        
        #w_tilde AND W
        th <- w_tilde <- w <- array(NA, dim=c(TT, i.c))  #UN-NORMLAIZED WEIGHTS
        for(i_t in 1:TT) {
            th[i_t,] <- rgamma(i.c, c(d0, rep(d, i.c-1)), 1)
            w_tilde[i_t,] <- th[i_t,]/sum(th[i_t,])
            w[i_t,] <- (1-w_star[i_t])*w_tilde[i_t,]
        }
        
        
        #COMPUTE m, M AND P
        m <- array(NA, dim=c(SS, TT))
        M <- array(NA, dim=c(SS, TT))
        p <- array(NA, dim=c(SS, TT))
        for(i_s in 1:SS)
        {
            for(i_t in 1:TT)
            {
                m[i_s, i_t] <- sum(w_tilde[i_t,]*L[i_s,])
                M[i_s, i_t] <- sum(w[i_t,]*L[i_s,]) + w_star[i_t]*2 # last term is normal contamination
                p[i_s, i_t] <- sum(w[i_t,]*Z[i_s,]*p_z0)/M[i_s,i_t]
            }
        }
        
        #COPY CURRENT SAMPLE IN A DIFFERENT FORMAT
        A <- array(A, dim=c(1, i.c*(Q+1)))[1,]
        L <- array(L, dim=c(1, SS*i.c))[1,]
        Z <- array(Z, dim=c(1, SS*i.c))[1,]
        m <- array(M, dim=c(1, SS*TT))[1,]
        M <- array(M, dim=c(1, SS*TT))[1,]
        p <- array(p, dim=c(1, SS*TT))[1,]
        th <- array(th, dim=c(1, i.c*TT))[1,]
        w_tilde <- array(w_tilde, dim=c(1, i.c*TT))[1,]
        w <- array(w, dim=c(1, i.c*TT))[1,]
        ppi <- array(ppi, dim=c(1, i.c*(Q+1)))[1,]

        output <- .C("fn_CNV_MCMC_1", alpha=as.double(alpha/(i.c-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), a_w=as.double(a_w), b_w=as.double(b_w), d0=as.double(d0), d=as.double(d), aa0_z=as.double(a_z0), bb0_z=as.double(b_z0), kappa=as.double(kappa), CC=as.integer(i.c), ppi=as.double(ppi), phi=as.double(phi), L=as.double(L), Z=as.double(Z), A=as.double(A), p0_z=as.double(p_z0[1]), w_star=as.double(w_star), th=as.double(th), w_tilde=as.double(w_tilde), w=as.double(w), m=as.double(m), M=as.double(M), p=as.double(p), SS=as.integer(SS), TT=as.integer(TT), ind_B=as.integer(ind_B), M_B=as.double(MM_B), n=as.double(nn_Tr), N=as.double(NN_Tr), N_sum=as.double(NN_Tr_sum), BB=as.double(BB_Tr), NN_iter=as.integer(N_sam))
        
        sam.all$phi[[i.c-min_C+1]] <- output$phi
        sam.all$pi[[i.c-min_C+1]] <- output$ppi
        
        sam.all$L[[i.c-min_C+1]] <- output$L
        sam.all$Z[[i.c-min_C+1]] <- output$Z
        sam.all$A[[i.c-min_C+1]] <- output$A

        sam.all$p0_z[i.c-min_C+1] <- output$p0_z
        sam.all$w_star[[i.c-min_C+1]] <- output$w_star
        sam.all$th[[i.c-min_C+1]] <- output$th
        sam.all$w_tilde[[i.c-min_C+1]] <- output$w_tilde
        sam.all$w[[i.c-min_C+1]] <- output$w
        
        sam.all$m[[i.c-min_C+1]] <- output$m
        sam.all$M[[i.c-min_C+1]] <- output$M
        sam.all$p[[i.c-min_C+1]] <- output$p
    }

    return(sam.all)
}

#MCMC.sam <- fn.MCMC.CNV(sam.all.OR, hyper, S, T, N_Tr.tmp, n_Tr.tmp, N_Tr_sum, B_Tr.tmp, N_Te.tmp, n_Te.tmp, N_Te_sum, B_Te.tmp, N.tmp, n.tmp, N_sum, Min_C, Max_C, n.sam)
#NN_Te <- N_Te.tmp; nn_Te <- n_Te.tmp; NN_Te_sum <- N_Te_sum; BB_Te <- B_Te.tmp; NN <- N.tmp; nn <-  n.tmp; NN_sum <- N_sum; sam.all <- sam.all.OR; N_sam <- n_sam
fn.MCMC.BayClone2 <- function(sam.all, hpara, SS, TT, ind_B, MM_B, NN_Tr, nn_Tr, NN_Tr_sum, BB_Tr, NN_Te, nn_Te, NN_Te_sum, BB_Te, NN, nn, NN_sum, min_C, max_C, N_sam)
{
    #HYPERPARAMETERS
    r <- hpara$r
    
    a <- hpara$a; b <- hpara$b;  #PHI
    a_w <- hpara$a_w; b_w <- hpara$b_w ## W^\star
    d0 <- hpara$d0; d <- hpara$d  #W
    
    #PRIOR FOR P_zO
    a_z0 <- hpara$a_z0; b_z0 <- hpara$b_z0 # FOR P_Z0
    
    #FOR L
    Q <- hpara$Q;
    alpha <- hpara$alpha
    beta <- hpara$beta
    gam <- hpara$gam
    
    ##FOR m
    kappa <- hpara$kappa

    #SAVE MCMC SAMPLES
    sam <- NULL
    sam$C <- rep(0, N_sam)
    
    sam$L <- NULL
    sam$Z <- NULL
    
    sam$w_star <- array(NA, dim=c(N_sam, TT))
    sam$th <- NULL
    sam$w_tilde <- NULL
    sam$w <- NULL
    
    sam$phi <- array(NA, dim=c(N_sam, TT))
    sam$pi  <- NULL
    
    sam$p0_z <- rep(0, N_sam)
    
    sam$m <- array(NA, dim=c(N_sam, SS, TT))
    sam$M <- array(NA, dim=c(N_sam, SS, TT))
    sam$p <- array(NA, dim=c(N_sam, SS, TT))
    
    ##---randomly choose one value as a current value of C
    C_cur <- sample((min_C:max_C), 1)
    
    for(i_sam in 1:N_sam)
    {
        if((i_sam%%100)==0)
        {
            print(paste("i.sam=", i_sam))
            print(date())
        }
        
        ##UPDATE (C, X) THROUGH M-H###############################################################
        ###CURRENT C_CUR
        output_cur <- .C("fn_CNV_MCMC_2", alpha=as.double(alpha/(C_cur-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), a_w=as.double(a_w), b_w=as.double(b_w), d0=as.double(d0), d=as.double(d), aa0_z=as.double(a_z0), bb0_z=as.double(b_z0), kappa=as.double(kappa), CC=as.integer(C_cur), ppi=as.double(sam.all$pi[[C_cur-min_C+1]]), phi=as.double(sam.all$phi[[C_cur-min_C+1]]), L=as.double(sam.all$L[[C_cur-min_C+1]]), Z=as.double(sam.all$Z[[C_cur-min_C+1]]), A=as.double(sam.all$A[[C_cur-min_C+1]]), p0_z=as.double(sam.all$p0_z[[C_cur-min_C+1]]), w_star=as.double(sam.all$w_star[[C_cur-min_C+1]]), th=as.double(sam.all$th[[C_cur-min_C+1]]), w_tilde=as.double(sam.all$w_tilde[[C_cur-min_C+1]]), w=as.double(sam.all$w[[C_cur-min_C+1]]), m=as.double(sam.all$m[[C_cur-min_C+1]]), M=as.double(sam.all$M[[C_cur-min_C+1]]), p=as.double(sam.all$p[[C_cur-min_C+1]]), SS=as.integer(SS), TT=as.integer(TT), M_B=as.double(MM_B), n_tr=as.double(nn_Tr), N_tr=as.double(NN_Tr), N_tr_sum=as.double(NN_Tr_sum), ind_B=as.integer(ind_B), BB_tr=as.double(BB_Tr), n_te=as.double(nn_Te), N_te=as.double(NN_Te), N_te_sum=as.double(NN_Te_sum),  BB_te=as.double(BB_Te), loglike=as.double(0), NN_iter=as.integer(0))
        alpha_cur <- output_cur$loglike + (C_cur-2)*log(1-r)
        
        ###CURRENT C_PRO
        C_pro <- sample((min_C:max_C), 1)  ###uniform proposal
        output_pro <- .C("fn_CNV_MCMC_2", alpha=as.double(alpha/(C_pro-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), a_w=as.double(a_w), b_w=as.double(b_w), d0=as.double(d0), d=as.double(d), aa0_z=as.double(a_z0), bb0_z=as.double(b_z0), kappa=as.double(kappa), CC=as.integer(C_pro), ppi=as.double(sam.all$pi[[C_pro-min_C+1]]), phi=as.double(sam.all$phi[[C_pro-min_C+1]]), L=as.double(sam.all$L[[C_pro-min_C+1]]), Z=as.double(sam.all$Z[[C_pro-min_C+1]]), A=as.double(sam.all$A[[C_pro-min_C+1]]), p0_z=as.double(sam.all$p0_z[[C_pro-min_C+1]]), w_star=as.double(sam.all$w_star[[C_pro-min_C+1]]),th=as.double(sam.all$th[[C_pro-min_C+1]]), w_tilde=as.double(sam.all$w_tilde[[C_pro-min_C+1]]), w=as.double(sam.all$w[[C_pro-min_C+1]]), m=as.double(sam.all$m[[C_pro-min_C+1]]), M=as.double(sam.all$M[[C_pro-min_C+1]]), p=as.double(sam.all$p[[C_pro-min_C+1]]), SS=as.integer(SS), TT=as.integer(TT), ind_B=as.integer(ind_B), M_B=as.double(MM_B), n_tr=as.double(nn_Tr), N_tr=as.double(NN_Tr), N_tr_sum=as.double(NN_Tr_sum), BB_tr=as.double(BB_Tr), n_te=as.double(nn_Te), N_te=as.double(NN_Te), N_te_sum=as.double(NN_Te_sum), BB_te=as.double(BB_Te), loglike=as.double(0), NN_iter=as.integer(1))
        alpha_pro <- output_pro$loglike + (C_pro-2)*log(1-r)
        
        sam.all$phi[[C_pro-min_C+1]] <- output_pro$phi
        sam.all$pi[[C_pro-min_C+1]] <- output_pro$ppi
        
        sam.all$L[[C_pro-min_C+1]] <- output_pro$L
        sam.all$Z[[C_pro-min_C+1]] <- output_pro$Z
        sam.all$A[[C_pro-min_C+1]] <- output_pro$A
        
        sam.all$w_star[[C_pro-min_C+1]] <- output_pro$w_star
        sam.all$p0_z[C_pro-min_C+1] <- output_pro$p0_z
        sam.all$th[[C_pro-min_C+1]] <- output_pro$th
        sam.all$w_tilde[[C_pro-min_C+1]] <- output_pro$w_tilde
        sam.all$w[[C_pro-min_C+1]] <- output_pro$w
        
        sam.all$m[[C_pro-min_C+1]] <- output_pro$m
        sam.all$M[[C_pro-min_C+1]] <- output_pro$M
        sam.all$p[[C_pro-min_C+1]] <- output_pro$p
        
        #print(c(C_cur, C_pro, (alpha_pro - alpha_cur)))
        
        ##ACCEPT??
        if(log(runif(1)) < (alpha_pro - alpha_cur))
        {
            output_cur <- output_pro
            C_cur <- C_pro
        }
        
        ##UPDATE X (THE OTHER PARAMETERS BUT C) USING ALL DATA
        output_cur <- .C("fn_CNV_MCMC_1", alpha=as.double(alpha/(C_cur-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), a_w=as.double(a_w), b_w=as.double(b_w), d0=as.double(d0), d=as.double(d), aa0_z=as.double(a_z0), bb0_z=as.double(b_z0), kappa=as.double(kappa), CC=as.integer(C_cur), ppi=as.double(output_cur$ppi), phi=as.double(output_cur$phi), L=as.double(output_cur$L), Z=as.double(output_cur$Z), A=as.double(output_cur$A), p0_z=as.double(output_cur$p0_z), w_star=as.double(output_cur$w_star), th=as.double(output_cur$th), w_tilde=as.double(output_cur$w_tilde), w=as.double(output_cur$w), m=as.double(output_cur$m), M=as.double(output_cur$M), p=as.double(output_cur$p), SS=as.integer(SS), TT=as.integer(TT), ind_B=as.integer(ind_B), M_B=as.double(MM_B), n=as.double(nn), N=as.double(NN), N_sum=as.double(NN_sum), BB=as.double(rep(1, SS*TT)), NN_iter=as.integer(10))
        
        
        ####SAVE SAMPLES
        sam$C[i_sam] <- C_cur
        sam$L[[i_sam]] <- array(output_cur$L, dim=c(SS, C_cur))
        sam$Z[[i_sam]] <- array(output_cur$Z, dim=c(SS, C_cur))
        
        sam$w_star[i_sam,] <- output_cur$w_star
        sam$th[[i_sam]] <- array(output_cur$th, dim=c(TT, C_cur))
        sam$w_tilde[[i_sam]] <- array(output_cur$w_tilde, dim=c(TT, C_cur))
        sam$w[[i_sam]] <- array(output_cur$w, dim=c(TT, C_cur))
        
        sam$phi[i_sam,] <- output_cur$phi
        sam$pi[[i_sam]] <- array(output_cur$ppi, dim=c(C_cur, (Q+1)))
        
        sam$p0_z[i_sam] <- output_cur$p0_z
        
        sam$m[i_sam,,] <- array(output_cur$m, dim=c(SS, TT))
        sam$M[i_sam,,] <- array(output_cur$M, dim=c(SS, TT))
        sam$p[i_sam,,] <- array(output_cur$p, dim=c(SS, TT))
    }
    return(sam)
}



#//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#//%%%%%%%%%%%  w0 and p0_z are *FIXED**
#//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#SS <- S; TT <- T; NN_Tr <- N_Tr.tmp; nn_Tr <- n_Tr.tmp; NN_Tr_sum <- N_Tr_sum; BB_Tr <- B_Tr.tmp; hpara <- hyper; min_C <- Min_C; max_C <- Max_C; N_sam <- burn.in;
#SS <- S; TT <- T; NN_Tr <- N.tmp; nn_Tr <- n.tmp; NN_Tr_sum <- N_sum; BB_Tr <- rep(1, SS*TT); hpara <- hyper; min_C <- Min_C; max_C <- Max_C; N_sam <- burn.in;
fn.Tr.BayClone2.fixed.w0 <- function(hpara, SS, TT, ind_B, MM_B, NN_Tr, nn_Tr, NN_Tr_sum, BB_Tr, min_C, max_C, N_sam)
{
    #HYPERPARAMETERS
    a <- hpara$a; b <- hpara$b;  #PHI
    a_w <- hpara$a_w; b_w <- hpara$b_w ## W^\star
    d <- hpara$d  #W
    
    ## we fix p0 and w0 (for all t)
    fixed_p0_z <- hpara$p0_z
    fixed_w0 <- hpara$w0

    #PRIOR FOR P_zO
    a_z0 <- hpara$a_z0; b_z0 <- hpara$b_z0 # FOR P_Z0
    
    #FOR L
    Q <- hpara$Q;
    alpha <- hpara$alpha
    beta <- hpara$beta
    gam <- hpara$gam
    
    ##FOR m
    kappa <- hpara$kappa

    #SAVE THE PROPOSALS FOR ALL C VALUES
    sam.all <- NULL
    
    sam.all$phi <- NULL
    sam.all$pi <- NULL
    
    sam.all$L <- NULL
    sam.all$Z <- NULL
    sam.all$A <- NULL
    
    sam.all$p0_z <- rep(NA, max_C-min_C+1)
    sam.all$w_star <- NULL
    sam.all$th <- NULL
    sam.all$w_tilde <- NULL
    sam.all$w <- NULL
    
    sam.all$m <- NULL
    sam.all$M <- NULL
    sam.all$p <- NULL
    
    #INITIALIZATION FOR EACH C
    for(i.c in min_C:max_C)
    {
        print(date())
        print(paste("i.c", i.c))
        set.seed(545754)
        
        ##L AND Z
        Z <- L <- array(0, dim=c(SS, i.c))
        
        #C--TO INITIALIZE
        p.ave <- apply(array(nn_Tr/NN_Tr, dim=c(SS, TT)), 1, mean)
        p.quant <- quantile(p.ave, probs=seq(0.10, 0.95, length.out=(i.c-1)))
        
        Z[,1] <- L[,1] <- 2 #CELL TYPE 0
        
        for(ii.c in 2:i.c) #EXCEPT CELL TYPE 0
        {
            Z[p.ave > p.quant[ii.c-1],ii.c] <- 2
        }
        
        ###CORRECTED MARCH-30TH
        for(ii.c in 2:i.c)
        for(i.s in 1:SS)
        {
            if(Z[i.s, ii.c]==Q)
            {
                L[i.s, ii.c] <- Q
            }else{
                L[i.s, ii.c] <- sample((Z[i.s, ii.c]:Q), 1, replace=FALSE)  ###If x has length 1, is numeric (in the sense of is.numeric) and x >= 1, sampling via sample takes place from 1:x. Note that this convenience feature may lead to undesired behaviour when x is of varying length in calls such as sample(x).
            }
        }
        
        
        ###  obey BB output if BB output is integer
        ## if (ind_B[i_s] == 1) {  // if BB output is not integer
        for(ii.c in 2:i.c)
        {
            L[ind_B==0, ii.c] <- MM_B[ind_B==0]
        }
        
        #COMPUTE A
        A <- array(NA, dim=c(i.c, (Q+1)))
        for(i.q in (0:Q))
        {
            A[,(i.q+1)] <- apply(L==i.q, 2, sum)
        }
        
        #PI
        ppi <- (A+5)/(SS + (Q+1)*5)
        
        #P0_Z
        p_z0 <- c(fixed_p0_z, rep(1, i.c-1))
        
        #PHI
        phi <- rgamma(TT, a, b)
        
        # W^\STAR
        w_star <- rbeta(TT, a_w, b_w)
        
        #w_tilde AND W
        th <- w_tilde <- w <- array(-1, dim=c(TT, i.c))  #UN-NORMLAIZED WEIGHTS
        w_tilde[, 1] <- fixed_w0/(1-w_star)  ##backengineer w_tilde_t0
        w[,1] <- fixed_w0
        for(i_t in 1:TT) {
            th[i_t, -1] <- rgamma(i.c-1, d, 1)
            w_tilde[i_t, -1] <- th[i_t, -1]/sum(th[i_t, -1])*(1-w_tilde[i_t,1])
            w[i_t, -1] <- (1 - w_star[i_t] - fixed_w0)*w_tilde[i_t, -1]
        }
        
        #COMPUTE M AND P
        m <- array(NA, dim=c(SS, TT))
        M <- array(NA, dim=c(SS, TT))
        p <- array(NA, dim=c(SS, TT))
        for(i_s in 1:SS)
        {
            for(i_t in 1:TT)
            {
                m[i_s, i_t] <- sum(w_tilde[i_t,]*L[i_s,])
                M[i_s, i_t] <- sum(w[i_t,]*L[i_s,]) + w_star[i_t]*2 # last term is normal contamination
                p[i_s, i_t] <- sum(w[i_t,]*Z[i_s,]*p_z0)/M[i_s,i_t]
            }
        }
        
        #COPY CURRENT SAMPLE IN A DIFFERENT FORMAT
        A <- array(A, dim=c(1, i.c*(Q+1)))[1,]
        L <- array(L, dim=c(1, SS*i.c))[1,]
        Z <- array(Z, dim=c(1, SS*i.c))[1,]
        M <- array(M, dim=c(1, SS*TT))[1,]
        p <- array(p, dim=c(1, SS*TT))[1,]
        th <- array(th, dim=c(1, i.c*TT))[1,]
        w_tilde <- array(w_tilde, dim=c(1, i.c*TT))[1,]
        w <- array(w, dim=c(1, i.c*TT))[1,]
        ppi <- array(ppi, dim=c(1, i.c*(Q+1)))[1,]
        
        output <- .C("fn_CNV_MCMC_1_fixed_w0", alpha=as.double(alpha/(i.c-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), a_w=as.double(a_w), b_w=as.double(b_w), d=as.double(d), kappa=as.double(kappa), CC=as.integer(i.c), ppi=as.double(ppi), phi=as.double(phi), L=as.double(L), Z=as.double(Z), A=as.double(A), p0_z=as.double(fixed_p0_z), w_star=as.double(w_star), th=as.double(th), w_tilde=as.double(w_tilde), w=as.double(w), m=as.double(m), M=as.double(M), p=as.double(p), SS=as.integer(SS), TT=as.integer(TT), ind_B=as.integer(ind_B), M_B=as.double(MM_B), n=as.double(nn_Tr), N=as.double(NN_Tr), N_sum=as.double(NN_Tr_sum), BB=as.double(BB_Tr), NN_iter=as.integer(N_sam))
        
        sam.all$phi[[i.c-min_C+1]] <- output$phi
        sam.all$pi[[i.c-min_C+1]] <- output$ppi
        
        sam.all$L[[i.c-min_C+1]] <- output$L
        sam.all$Z[[i.c-min_C+1]] <- output$Z
        sam.all$A[[i.c-min_C+1]] <- output$A
        
        sam.all$w_star[[i.c-min_C+1]] <- output$w_star
        sam.all$th[[i.c-min_C+1]] <- output$th
        sam.all$w_tilde[[i.c-min_C+1]] <- output$w_tilde
        sam.all$w[[i.c-min_C+1]] <- output$w
        
        sam.all$m[[i.c-min_C+1]] <- output$m
        sam.all$M[[i.c-min_C+1]] <- output$M
        sam.all$p[[i.c-min_C+1]] <- output$p
        
    }
    
    return(sam.all)
}

#MCMC.sam <- fn.MCMC.CNV(sam.all.OR, hyper, S, T, N_Tr.tmp, n_Tr.tmp, N_Tr_sum, B_Tr.tmp, N_Te.tmp, n_Te.tmp, N_Te_sum, B_Te.tmp, N.tmp, n.tmp, N_sum, Min_C, Max_C, n.sam)
#NN_Te <- N_Te.tmp; nn_Te <- n_Te.tmp; NN_Te_sum <- N_Te_sum; BB_Te <- B_Te.tmp; NN <- N.tmp; nn <-  n.tmp; NN_sum <- N_sum; sam.all <- sam.all.OR
fn.MCMC.BayClone2.fixed.w0 <- function(sam.all, hpara, SS, TT, ind_B, MM_B, NN_Tr, nn_Tr, NN_Tr_sum, BB_Tr, NN_Te, nn_Te, NN_Te_sum, BB_Te, NN, nn, NN_sum, min_C, max_C, N_sam)
{
    #HYPERPARAMETERS
    r <- hpara$r
    
    a <- hpara$a; b <- hpara$b;  #PHI
    a_w <- hpara$a_w; b_w <- hpara$b_w ## W^\star
    d <- hpara$d  #W
    
    fixed_w0 <- hpara$w0
    fixed_p0_z <- hpara$p0_z
    
    #FOR L
    Q <- hpara$Q;
    alpha <- hpara$alpha
    beta <- hpara$beta
    gam <- hpara$gam
    
    ##FOR m
    kappa <- hpara$kappa

    #SAVE MCMC SAMPLES
    sam <- NULL
    sam$C <- rep(0, N_sam)
    
    sam$L <- NULL
    sam$Z <- NULL
    
    sam$w_star <- array(NA, dim=c(N_sam, TT))
    sam$th <- NULL
    sam$w_tilde <- NULL
    sam$w <- NULL
    
    sam$phi <- array(NA, dim=c(N_sam, TT))
    sam$pi  <- NULL
    
    sam$m <- array(NA, dim=c(N_sam, SS, TT))
    sam$M <- array(NA, dim=c(N_sam, SS, TT))
    sam$p <- array(NA, dim=c(N_sam, SS, TT))
    
    ##---randomly choose one value as a current value of C
    C_cur <- sample((min_C:max_C), 1)
    
    for(i_sam in 1:N_sam)
    {
        if((i_sam%%100)==0)
        {
            print(paste("i.sam=", i_sam))
            print(date())
        }
        
        ##UPDATE (C, X) THROUGH M-H###############################################################
        ###CURRENT C_CUR
        output_cur <- .C("fn_CNV_MCMC_2_fixed_w0", alpha=as.double(alpha/(C_cur-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), a_w=as.double(a_w), b_w=as.double(b_w), d=as.double(d), kappa=as.double(kappa), CC=as.integer(C_cur), ppi=as.double(sam.all$pi[[C_cur-min_C+1]]), phi=as.double(sam.all$phi[[C_cur-min_C+1]]), L=as.double(sam.all$L[[C_cur-min_C+1]]), Z=as.double(sam.all$Z[[C_cur-min_C+1]]), A=as.double(sam.all$A[[C_cur-min_C+1]]), p0_z=as.double(fixed_p0_z), w_star=as.double(sam.all$w_star[[C_cur-min_C+1]]), th=as.double(sam.all$th[[C_cur-min_C+1]]), w_tilde=as.double(sam.all$w_tilde[[C_cur-min_C+1]]), w=as.double(sam.all$w[[C_cur-min_C+1]]), m=as.double(sam.all$m[[C_cur-min_C+1]]), M=as.double(sam.all$M[[C_cur-min_C+1]]), p=as.double(sam.all$p[[C_cur-min_C+1]]), SS=as.integer(SS), TT=as.integer(TT), ind_B=as.integer(ind_B), M_B=as.double(MM_B), n_tr=as.double(nn_Tr), N_tr=as.double(NN_Tr), N_tr_sum=as.double(NN_Tr_sum), BB_tr=as.double(BB_Tr), n_te=as.double(nn_Te), N_te=as.double(NN_Te), N_te_sum=as.double(NN_Te_sum), BB_te=as.double(BB_Te), loglike=as.double(0), NN_iter=as.integer(0))
        alpha_cur <- output_cur$loglike + (C_cur-2)*log(1-r)
        
        ###CURRENT C_PRO
        C_pro <- sample((min_C:max_C), 1)  ###uniform proposal
        output_pro <- .C("fn_CNV_MCMC_2_fixed_w0", alpha=as.double(alpha/(C_pro-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), a_w=as.double(a_w), b_w=as.double(b_w), d=as.double(d), kappa=as.double(kappa), CC=as.integer(C_pro), ppi=as.double(sam.all$pi[[C_pro-min_C+1]]), phi=as.double(sam.all$phi[[C_pro-min_C+1]]), L=as.double(sam.all$L[[C_pro-min_C+1]]), Z=as.double(sam.all$Z[[C_pro-min_C+1]]), A=as.double(sam.all$A[[C_pro-min_C+1]]), p0_z=as.double(fixed_p0_z), w_star=as.double(sam.all$w_star[[C_pro-min_C+1]]),th=as.double(sam.all$th[[C_pro-min_C+1]]), w_tilde=as.double(sam.all$w_tilde[[C_pro-min_C+1]]), w=as.double(sam.all$w[[C_pro-min_C+1]]), m=as.double(sam.all$m[[C_pro-min_C+1]]), M=as.double(sam.all$M[[C_pro-min_C+1]]), p=as.double(sam.all$p[[C_pro-min_C+1]]), SS=as.integer(SS), TT=as.integer(TT), ind_B=as.integer(ind_B), M_B=as.double(MM_B), n_tr=as.double(nn_Tr), N_tr=as.double(NN_Tr), N_tr_sum=as.double(NN_Tr_sum), BB_tr=as.double(BB_Tr), n_te=as.double(nn_Te), N_te=as.double(NN_Te), N_te_sum=as.double(NN_Te_sum), BB_te=as.double(BB_Te), loglike=as.double(0), NN_iter=as.integer(1))
        alpha_pro <- output_pro$loglike + (C_pro-2)*log(1-r)
        
        sam.all$phi[[C_pro-min_C+1]] <- output_pro$phi
        sam.all$pi[[C_pro-min_C+1]] <- output_pro$ppi
        
        sam.all$L[[C_pro-min_C+1]] <- output_pro$L
        sam.all$Z[[C_pro-min_C+1]] <- output_pro$Z
        sam.all$A[[C_pro-min_C+1]] <- output_pro$A
        
        sam.all$w_star[[C_pro-min_C+1]] <- output_pro$w_star
        sam.all$th[[C_pro-min_C+1]] <- output_pro$th
        sam.all$w_tilde[[C_pro-min_C+1]] <- output_pro$w_tilde
        sam.all$w[[C_pro-min_C+1]] <- output_pro$w
        
        sam.all$m[[C_pro-min_C+1]] <- output_pro$m
        sam.all$M[[C_pro-min_C+1]] <- output_pro$M
        sam.all$p[[C_pro-min_C+1]] <- output_pro$p
        
        ##ACCEPT??
        if(log(runif(1)) < (alpha_pro - alpha_cur))
        {
            output_cur <- output_pro
            C_cur <- C_pro
        }
        
        ##UPDATE X (THE OTHER PARAMETERS BUT C) USING ALL DATA
        output_cur <- .C("fn_CNV_MCMC_1_fixed_w0", alpha=as.double(alpha/(C_cur-1)), beta=as.double(beta), gam=as.double(gam), QQ=as.integer(Q), a=as.double(a), b=as.double(b), a_w=as.double(a_w), b_w=as.double(b_w), d=as.double(d), kappa=as.double(kappa), CC=as.integer(C_cur), ppi=as.double(output_cur$ppi), phi=as.double(output_cur$phi), L=as.double(output_cur$L), Z=as.double(output_cur$Z), A=as.double(output_cur$A), p0_z=as.double(fixed_p0_z), w_star=as.double(output_cur$w_star), th=as.double(output_cur$th), w_tilde=as.double(output_cur$w_tilde), w=as.double(output_cur$w), m=as.double(output_cur$m), M=as.double(output_cur$M), p=as.double(output_cur$p), SS=as.integer(SS), TT=as.integer(TT), ind_B=as.integer(ind_B), M_B=as.double(MM_B), n=as.double(nn), N=as.double(NN), N_sum=as.double(NN_sum), BB=as.double(rep(1, SS*TT)), NN_iter=as.integer(10))
        
        
        ####SAVE SAMPLES
        sam$C[i_sam] <- C_cur
        sam$L[[i_sam]] <- array(output_cur$L, dim=c(SS, C_cur))
        sam$Z[[i_sam]] <- array(output_cur$Z, dim=c(SS, C_cur))
        
        sam$w_star[i_sam,] <- output_cur$w_star
        sam$th[[i_sam]] <- array(output_cur$th, dim=c(TT, C_cur))
        sam$w_tilde[[i_sam]] <- array(output_cur$w_tilde, dim=c(TT, C_cur))
        sam$w[[i_sam]] <- array(output_cur$w, dim=c(TT, C_cur))
        
        sam$phi[i_sam,] <- output_cur$phi
        sam$pi[[i_sam]] <- array(output_cur$ppi, dim=c(C_cur, (Q+1)))
        
        sam$m[i_sam,,] <- array(output_cur$m, dim=c(SS, TT))
        sam$M[i_sam,,] <- array(output_cur$M, dim=c(SS, TT))
        sam$p[i_sam,,] <- array(output_cur$p, dim=c(SS, TT))
    }
    return(sam)
}

#nn <- n; NN <- N; SS <- S; TT <- T; hpara <- hyper; b <- 0.9; min_C <- Min_C; max_C <- Max_C; burn_in <- burn.in; n_sam <- n.sam; MM_B <- M_B; ind_B <- Ind_B

BayClone2_B <- function(ind_B, MM_B, nn, NN, SS, TT, hpara, b, min_C, max_C, burn_in, n_sam, ind_w0_fixed)
{
    n_C <- max_C - min_C + 1
    
    #TRAINING
    B <- array(rbeta(SS*TT, b*1000, (1-b)*1000), dim=c(SS, TT))
    n_Tr <- nn*B
    n_Te <- nn - n_Tr
    
    N_Tr <- NN*B
    N_Te <- NN - N_Tr
    
    n_Tr.tmp <- array(n_Tr, dim=c(1, SS*TT))[1,]
    N_Tr.tmp <- array(N_Tr, dim=c(1, SS*TT))[1,]
    N_Tr_sum <- apply(N_Tr, 2, sum)
    
    n_Te.tmp <- array(n_Te, dim=c(1, SS*TT))[1,]
    N_Te.tmp <- array(N_Te, dim=c(1, SS*TT))[1,]
    N_Te_sum <- apply(N_Te, 2, sum)
    
    B_Tr.tmp <- array(B, dim=c(1, SS*TT))[1,]
    B_Te.tmp <- array(1-B, dim=c(1, SS*TT))[1,]
    
    MM_B.tmp <- array(MM_B, dim=c(1, SS*TT))[1,]
    
    n.tmp <- array(nn, dim=c(1, SS*TT))[1,]
    N.tmp <- array(NN, dim=c(1, SS*TT))[1,]
    N_sum <- apply(NN, 2, sum)

    
    if(ind_w0_fixed == 1)
    {
        #/// TRAINING
        sam.all.OR <- fn.Tr.BayClone2.fixed.w0(hpara, SS, TT, ind_B, MM_B.tmp, N_Tr.tmp, n_Tr.tmp, N_Tr_sum, B_Tr.tmp, min_C, max_C, burn_in)
        
        #/// MCMC
        MCMC.sam <- fn.MCMC.BayClone2.fixed.w0(sam.all.OR, hpara, SS, TT, ind_B, MM_B.tmp, N_Tr.tmp, n_Tr.tmp, N_Tr_sum, B_Tr.tmp, N_Te.tmp, n_Te.tmp, N_Te_sum, B_Te.tmp, N.tmp, n.tmp, N_sum, min_C, max_C, n_sam)
    }else{
        #/// TRAINING
        sam.all.OR <- fn.Tr.BayClone2(hpara, SS, TT, ind_B, MM_B.tmp, N_Tr.tmp, n_Tr.tmp, N_Tr_sum, B_Tr.tmp, min_C, max_C, burn_in)
        
        #/// MCMC
        MCMC.sam <- fn.MCMC.BayClone2(sam.all.OR, hpara, SS, TT, ind_B, MM_B.tmp, N_Tr.tmp, n_Tr.tmp, N_Tr_sum, B_Tr.tmp, N_Te.tmp, n_Te.tmp, N_Te_sum, B_Te.tmp, N.tmp, n.tmp, N_sum, min_C, max_C, n_sam)
    }
    
    #save(sam.all.OR, file = "training_sam.RData")
    return(MCMC.sam)
}


rep.row <- function(x,n){
    matrix(rep(x,each=n),nrow=n)
}


## MCMC_sam_post <- fn.BayClone2.post.fixed.w0(MCMC_sam, SS, TT, n_sam, w_cutoff, w0_fixed, p0_fixed)
## *FIXED* w0 and p0
## remove subclones with w < w_cutoff for all the samples
fn.BayClone2.post.fixed.w0 <- function(MCMC_sam, SS, TT, N_sam, w_cutoff, w0_fixed, p0_fixed)
{
    ## FIXED
    #> names(MCMC.sam)
    #[1] "C"       "w_star"  "phi"   "M"       "p"       "L"
    #[8] "Z"       "th"      "w_tilde" "w"       "pi"
    
    #SAVE MCMC SAMPLES after post-processing
    sam <- NULL
    sam$w_star <- MCMC_sam$w_star
    sam$phi <- MCMC_sam$phi
    
    sam$C <- rep(0, N_sam)
    sam$L <- NULL
    sam$Z <- NULL
    
    sam$th <- NULL
    sam$w_tilde <- NULL
    sam$w <- NULL
    
    sam$pi  <- NULL
    
    sam$m <- array(NA, dim=c(N_sam, SS, TT))
    sam$M <- array(NA, dim=c(N_sam, SS, TT))
    sam$p <- array(NA, dim=c(N_sam, SS, TT))
    
    for(i_sam in 1:N_sam)
    {
        C_cur <- MCMC_sam$C[i_sam]
        w_cur <- as.matrix(MCMC_sam$w[[i_sam]])
        ind_cut_w <- as.matrix(w_cur > w_cutoff, 2, sum)  ## 0 ==  all the samples have weight < cut-off
        
        if((C_cur > 2) & (sum(ind_cut_w[-1] == 0) > 0))  ## at least TWO subclones (except cell type 0) & NOT(all the samples > cutoff)
        {
            remove_c <- (2:C_cur)[ind_cut_w[-1]==0]
            C_new <- C_cur - length(remove_c)

            th_new <- as.matrix(MCMC_sam$th[[i_sam]])
            th_new <- th_new[,-remove_c]
            if(TT==1) th_new <- t(as.matrix(th_new))
            
            w_star_cur <- MCMC_sam$w_star[i_sam]
            w_tilde_new <- array(0, dim=c(TT, C_new))
            w_tilde_new[,1] <- w0_fixed/(1-w_star_cur)  ## backengineer
            for(i_t in 1:TT)
            {
                w_tilde_new[i_t, -1] <- th_new[i_t, -1]/sum(th_new[i_t, -1])*(1-w_tilde_new[i_t,1])
            }
            
            w_new <- w_tilde_new*(1- w_star_cur)

            L_new <- MCMC_sam$L[[i_sam]]
            L_new <- L_new[,-remove_c]
            
            Z_new <- MCMC_sam$Z[[i_sam]]
            Z_new <- Z_new[,-remove_c]
            
            pi_new <- MCMC_sam$pi[[i_sam]]
            pi_new <- pi_new[-remove_c,]
            
            ## have to recompute M and p
            m_new <- M_new <- p_new <- array(NA, dim=c(SS, TT))
            for(i_t in 1:TT)
            {
                w_tilde_tmp <- rep.row(w_tilde_new[i_t,], SS)
                m_new[,i_t] <- apply(w_tilde_tmp*L_new, 1, sum)
                
                w_tmp <- rep.row(w_new[i_t,], SS)
                M_new[,i_t] <- apply(w_tmp*L_new, 1, sum) + w_star_cur*2
                
                w_tmp[,1] <- w_tmp[,1]*p0_fixed
                p_new[,i_t] <- apply(w_tmp*Z_new, 1, sum)/M_new[,i_t]
            }
            
            ####SAVE SAMPLES
            sam$C[i_sam] <- C_new
            sam$L[[i_sam]] <- L_new
            sam$Z[[i_sam]] <- Z_new
            
            #sam$w_star[i_sam,] ## already saved
            sam$th[[i_sam]] <- th_new
            sam$w_tilde[[i_sam]] <- w_tilde_new
            sam$w[[i_sam]] <- w_new
            
            #sam$phi[i_sam,]  ## already saved
            sam$pi[[i_sam]] <- pi_new
            
            sam$m[i_sam,,] <- m_new
            sam$M[i_sam,,] <- M_new
            sam$p[i_sam,,] <- p_new
        }else{ ##if((C_cur > 2) & (sum(ind_cut_w[-1] == 0) > 0))  NOT
            ####SAVE SAMPLES
            sam$C[i_sam] <- MCMC_sam$C[[i_sam]]
            sam$L[[i_sam]] <- MCMC_sam$L[[i_sam]]
            sam$Z[[i_sam]] <- MCMC_sam$Z[[i_sam]]
            
            #sam$w_star[i_sam,] ## already saved
            sam$th[[i_sam]] <- MCMC_sam$th[[i_sam]]
            sam$w_tilde[[i_sam]] <- MCMC_sam$w_tilde[[i_sam]]
            sam$w[[i_sam]] <- MCMC_sam$w[[i_sam]]
            
            #sam$phi[i_sam,]  ## already saved
            sam$pi[[i_sam]] <- MCMC_sam$pi[[i_sam]]
            
            sam$m[i_sam,,] <- MCMC_sam$m[i_sam,,]
            sam$M[i_sam,,] <- MCMC_sam$M[i_sam,,]
            sam$p[i_sam,,] <- MCMC_sam$p[i_sam,,]
        } ##if((C_cur > 2) & (sum(ind_cut_w[-1] == 0) > 0))
        
    } ##for(i_sam in 1:N_sam)
    
    return(sam)
}

## MCMC_sam_post <- fn.BayClone2.post.fixed.w0(MCMC_sam, SS, TT, n_sam, w_cutoff, w0_fixed, p0_fixed)
## *RANDOM* w0 and p0
## remove subclones with w < w_cutoff for all the samples
fn.BayClone2.post.random.w0 <- function(MCMC_sam, SS, TT, N_sam, w_cutoff)
{
    
    ##RANDOM
    #> names(MCMC.sam)
    #[1] "C"       "w_star"  "phi"     "p0_z"    "M"       "p"       "L"
    #[8] "Z"       "th"      "w_tilde" "w"       "pi"
    
    #SAVE MCMC SAMPLES after post-processing
    sam <- NULL
    sam$w_star <- MCMC_sam$w_star
    sam$phi <- MCMC_sam$phi
    sam$p0_z <- MCMC_sam$p0_z
    
    sam$C <- rep(0, N_sam)
    sam$L <- NULL
    sam$Z <- NULL
    
    sam$th <- NULL
    sam$w_tilde <- NULL
    sam$w <- NULL
    
    sam$pi  <- NULL
    
    sam$m <- array(NA, dim=c(N_sam, SS, TT))
    sam$M <- array(NA, dim=c(N_sam, SS, TT))
    sam$p <- array(NA, dim=c(N_sam, SS, TT))
    
    ## w0 is fixed
    
    
    for(i_sam in 1:N_sam)
    {
        C_cur <- MCMC_sam$C[i_sam]
        w_cur <- as.matrix(MCMC_sam$w[[i_sam]])
        ind_cut_w <- as.matrix(w_cur > w_cutoff, 2, sum)  ## 0 ==  all the samples have weight < cut-off
        
        if(TT > 1)
        {
            ind_cut_w <- apply(ind_cut_w, 2, sum)  ## ind_cut_w==0 means subclone c has small w for all the samples
        }
        
        if((C_cur > 2) & (sum(ind_cut_w[-1] == 0) > 0))  ## at least TWO subclones (except cell type 0) & NOT(all the samples > cutoff)
        {
            remove_c <- (2:C_cur)[ind_cut_w[-1]==0]
            C_new <- C_cur - length(remove_c)
            
            th_new <- as.matrix(MCMC_sam$th[[i_sam]])
            th_new <- as.matrix(th_new[,-remove_c])  ### corrected 18th March Subhajit
            if(TT==1) th_new <- t(as.matrix(th_new))
            
            w_star_cur <- MCMC_sam$w_star[i_sam]
            w_tilde_new <- array(0, dim=c(TT, C_new))
            
            for(i_t in 1:TT)
            {
                w_tilde_new[i_t,] <- th_new[i_t,]/sum(th_new[i_t,])  #corrected -- J-- march/22/2016
            }
            
            w_new <- w_tilde_new*(1- w_star_cur)
            
            L_new <- MCMC_sam$L[[i_sam]]
            L_new <- L_new[,-remove_c]
            
            Z_new <- MCMC_sam$Z[[i_sam]]
            Z_new <- Z_new[,-remove_c]
            
            pi_new <- MCMC_sam$pi[[i_sam]]
            pi_new <- pi_new[-remove_c,]
            
            ## have to recompute m, M and p
            m_new <- M_new <- p_new <- array(NA, dim=c(SS, TT))
            for(i_t in 1:TT)
            {
                w_tilde_tmp <- rep.row(w_tilde_new[i_t,], SS)
                m_new[,i_t] <- apply(w_tilde_tmp*L_new, 1, sum)
                
                w_tmp <- rep.row(w_new[i_t,], SS)
                M_new[,i_t] <- apply(w_tmp*L_new, 1, sum) + w_star_cur*2
                
                w_tmp[,1] <- w_tmp[,1]*MCMC_sam$p0_z[i_sam]
                p_new[,i_t] <- apply(w_tmp*Z_new, 1, sum)/M_new[,i_t]
            }
            
            ####SAVE SAMPLES
            sam$C[i_sam] <- C_new
            sam$L[[i_sam]] <- L_new
            sam$Z[[i_sam]] <- Z_new
            
            #sam$w_star[i_sam,] ## already saved
			sam$th[[i_sam]] <- list(NULL)  ### corrected 18th March Subhajit
            sam$th[[i_sam]] <- th_new

			sam$w_tilde[[i_sam]] <- list(NULL)  ### corrected 18th March Subhajit
            sam$w_tilde[[i_sam]] <- w_tilde_new

			sam$w[[i_sam]] <- list(NULL)  ### corrected 18th March Subhajit
            sam$w[[i_sam]] <- w_new
            
            #sam$phi[i_sam,]  ## already saved
            sam$pi[[i_sam]] <- pi_new
            
            sam$m[i_sam,,] <- m_new
            sam$M[i_sam,,] <- M_new
            sam$p[i_sam,,] <- p_new
        }else{ ##if((C_cur > 2) & (sum(ind_cut_w[-1] == 0) > 0))  NOT
            ####SAVE SAMPLES
            sam$C[i_sam] <- MCMC_sam$C[[i_sam]]
            sam$L[[i_sam]] <- MCMC_sam$L[[i_sam]]
            sam$Z[[i_sam]] <- MCMC_sam$Z[[i_sam]]
            
            #sam$w_star[i_sam,] ## already saved
            sam$th[[i_sam]] <- MCMC_sam$th[[i_sam]]
            sam$w_tilde[[i_sam]] <- MCMC_sam$w_tilde[[i_sam]]
            sam$w[[i_sam]] <- MCMC_sam$w[[i_sam]]
            
            #sam$phi[i_sam,]  ## already saved
            sam$pi[[i_sam]] <- MCMC_sam$pi[[i_sam]]
            
            sam$m[i_sam,,] <- MCMC_sam$m[i_sam,,]
            sam$M[i_sam,,] <- MCMC_sam$M[i_sam,,]
            sam$p[i_sam,,] <- MCMC_sam$p[i_sam,,]
        } ##if((C_cur > 2) & (sum(ind_cut_w[-1] == 0) > 0))
        
    } ##for(i_sam in 1:N_sam)
    
    return(sam)
}


# MCMC_sam <- MCMC.sam; SS <- S; TT <- T; n_sam <- n.sam; w_cutoff <- 0.03; ind_fixed <- ind.fixed; w0_fixed <- hyper$w0; p0_fixed <- hyper$p0_z;

fn.BayClone2.post.process.B <- function(MCMC_sam, SS, TT, n_sam, w_cutoff, ind_fixed, w0_fixed, p0_fixed)
{
    
    ## FIXED
    #> names(MCMC.sam)
    #[1] "C"       "w_star"  "phi"   "M"       "p"       "L"
    #[8] "Z"       "th"      "w_tilde" "w"       "pi"
    
    ##RANDOM
    #> names(MCMC.sam)
    #[1] "C"       "w_star"  "phi"     "p0_z"    "M"       "p"       "L"
    #[8] "Z"       "th"      "w_tilde" "w"       "pi"
    
    if(ind_fixed == 1)
    {
        #print("fixed")
        ## fixed
        MCMC_sam_post <- fn.BayClone2.post.fixed.w0(MCMC_sam, SS, TT, n_sam, w_cutoff, w0_fixed, p0_fixed)
    }else{
        #print("random")
        ## random
        MCMC_sam_post <- fn.BayClone2.post.random.w0(MCMC_sam, SS, TT, n_sam, w_cutoff)
    }
    return(MCMC_sam_post)
}




rep.row <- function(x,n){
    matrix(rep(x,each=n),nrow=n)
}



## MCMC_sam <- MCMC.sam.post; SS <-  S; TT <- T; n_sam <- n.sam; ind_fixed <- ind.fixed; max_C <- Max_C; min_C <- Min_C
fn.BayClone2.point.B <- function(MCMC_sam, SS, TT, min_C, max_C, n_sam, ind_fixed)
{
    #POSTERIOR DIST OF C --- now consider the case no tumor subclone--march/22/16--J
    ##n_C <- max_C - min_C + 1
    post_dist <- rep(0, max_C)
    for(i.c in 1:max_C) post_dist[i.c] <- sum(MCMC_sam$C==i.c)/n_sam
    
    ###  finding a point estimate of L, Z, w
    CC <- which.max(post_dist)  ## posterior mode + includes cell type 0

    post_point <- NULL
    post_point$C_dist <- post_dist
    post_point$C_star <- CC
    
    sam.C <- (1:n_sam)[MCMC_sam$C==CC]
    
    if(CC > 1)  ## if we have more than 2 subclones (including subclone0)
    {
        ############  finding point estimate of L
        n.C <- sum(MCMC_sam$C==CC)
        D.C <- array(0, dim=c(n.C, n.C))
        
        ind_N <- factorial(CC-1)
        index_set <- array(NA, dim=c(ind_N, CC-1))
        tmp <- permn(CC-1)
        for(i.n in 1:ind_N)
        {
            index_set[i.n,] <- tmp[[i.n]]
        }
        
        index_set <- array(index_set-1, dim=c(1, ind_N*(CC-1)))[1,]
        
        print(paste("How many samples for point estimate?? ", n.C))
        
        for(i in 1:(n.C-1))
        {
            if((i%%100)==0)
            {
                print(paste("i.sam=", i))
                print(date())
            }
            
            for(j in (i+1):n.C)
            {
                i.sam <- sam.C[i]
                j.sam <- sam.C[j]
                
                L.i <- array(MCMC_sam$L[[i.sam]], dim=c(1, SS*CC))[1,]
                L.j <- array(MCMC_sam$L[[j.sam]], dim=c(1, SS*CC))[1,]
                
                output <- .C("fn_sum_Z_1", SS=as.integer(SS), CC=as.integer(CC-1), ind_NN=as.integer(ind_N), ind_set=as.integer(index_set), Z_1=as.integer(L.i), Z_2=as.integer(L.j), min_dist=as.double(0))
                
                D.C[i, j] <- D.C[j, i] <- output$min_dist
            }
        }
        
        sam.C <- (1:n_sam)[MCMC_sam$C==CC]
        sum.dist <- apply(D.C, 1, sum)
        ind.min <- which.min(sum.dist)
        ind.min <- sam.C[ind.min]
        
        post.L <- MCMC_sam$L[[ind.min]]
        post.Z <- MCMC_sam$Z[[ind.min]]
        post.w <- MCMC_sam$w[[ind.min]]
        post.w.tilde <- MCMC_sam$w_tilde[[ind.min]]
        
        
        post_point$L <- post.L
        post_point$Z <- post.Z
        post_point$w <- post.w
        post_point$w_tilde <- post.w.tilde
    }## if we have more than 2 subclones (including subclone0)
    
    ### we can still do the following even when no tumor subclone  -- 03/22/16---J
    post.w.star_dist <- as.matrix(MCMC_sam$w_star[sam.C,])
    
    if(ind_fixed == 0)
    {
        post.p0.z <- mean(MCMC_sam$p0_z[sam.C])
        post.p0.z_dist <- MCMC_sam$p0_z[sam.C]
        
        w0.C <- array(NA, dim=c(length(sam.C), TT))
        i <- 1
        for(i_iter in sam.C)
        {
            w0.C[i,] <- (MCMC_sam$w[[i_iter]])[,1]
            i <- i + 1
        }
    }
    
    if(TT > 1)
    {
        post.phi <- apply(MCMC_sam$phi[sam.C,], 2, mean)
    }else{
        post.phi <- mean(MCMC_sam$phi[sam.C,])
    }
    
    #####################################################
    N_st_post <- p_st_post <- M_st_post <- array(0, dim=c(SS, TT))
    m_st_post <- multi_st_post <- array(0, dim=c(SS, TT))
    ploidy_dist <-  array(0, dim=c(length(sam.C), TT))
    i_ind <- 1
    
    for(i_sam in sam.C)
    {
        phi_tmp <- MCMC_sam$phi[i_sam,]
        m_tmp <- as.matrix(MCMC_sam$m[i_sam,,])  ## tumor copy number == L weighted by w_tilde
        M_tmp <- as.matrix(MCMC_sam$M[i_sam,,])  ## L weighted by w
        p_tmp <- as.matrix(MCMC_sam$p[i_sam,,])
        
        m_st_post <- m_st_post + m_tmp
        M_st_post <- M_st_post + M_tmp
        p_st_post <- p_st_post + p_tmp
        
        Z_tmp <- as.matrix(MCMC_sam$Z[[i_sam]])
        p0_tmp <- MCMC_sam$p0_z[i_sam]
        w_tilde_tmp <- as.matrix(MCMC_sam$w_tilde[[i_sam]])
        
        for(i_t in 1:TT)
        {
            N_st_post[, i_t] <- N_st_post[, i_t] + M_tmp[,i_t]*phi_tmp[i_t]/2
            
            w_tilde_t <- w_tilde_tmp[i_t,]  ## w_tilde for sample t
            w_tilde_t[1] <- w_tilde_t[1]*p0_tmp  ## adjust for subclone 0
            multi_st_post[, i_t] <- multi_st_post[, i_t] + apply(Z_tmp*rep.row(w_tilde_t, SS), 1, sum)
        }
        
        ploidy_dist[i_ind, ] <- apply(M_tmp, 2, mean)
        i_ind <- i_ind + 1
    }
    
    multi_st_post <- multi_st_post/length(sam.C)
    m_st_post <- m_st_post/length(sam.C)
    M_st_post <- M_st_post/length(sam.C)
    p_st_post <- p_st_post/length(sam.C)
    N_st_post <- N_st_post/length(sam.C)
    
    post_point$w_star_dist <- post.w.star_dist
    
    post_point$phi <- post.phi
    post_point$m <- m_st_post
    post_point$M <- M_st_post
    post_point$ploidy_dist <- ploidy_dist
    post_point$p <- p_st_post
    post_point$N <- N_st_post
    post_point$m_multi <- multi_st_post
    
    if(ind_fixed == 0)
    {
        post_point$p0_z <- post.p0.z
        post_point$p0_z_dist <- post.p0.z_dist
        post_point$w0_dist <- w0.C
    }
    
    return(post_point)
}



## NN <- N; MM_B <- M_B; nn <- n; SS <- S
fn.preprocessing <- function(NN, MM_B, nn, SS)
{
    if(sd(MM_B) > 0)
    {
        tmp <- lm(NN~MM_B)
        coef <- summary(tmp)$coeff[,1]
        
        N_hat <- coef[1] + MM_B*coef[2]  ## estimated N from the linear regression
        u.b <- summary(tmp)$sigma*3 + N_hat
        l.b <- -summary(tmp)$sigma*3 + N_hat
        
    }else{
        u.b <- sd(NN)*3 + mean(NN)
        l.b <- -sd(NN)*3 + mean(NN)
    }
    
    ind_remove <- (((NN > u.b)|(NN < l.b)))
    snv_set <- (1:SS)[ind_remove]  ## set of snv's removed
    
    if(length(snv_set) > 0)
    {
        NN <- as.matrix(NN[-snv_set,])
        MM_B <- as.matrix(MM_B[-snv_set,])
        nn <- as.matrix(nn[-snv_set,])
    }
    
    print(paste("# of removed loci", length(snv_set)))
    
    return(list(N=NN, M_B=MM_B, n=nn, S=nrow(nn),removed_ind=((1:SS)%in%snv_set)))
}



##  NN <- NN[,i_t]; nn <- nn[,i_t]; #MM_B <- MM_B[,i_t]
fn.get.wstar <- function(NN, nn, MM_B)
{
    ##########################   1.  Prepare Data Set   #####################
    ###                                                                 ###
    ### Take the SNVs between the 40% percentile and 60% percentile of N's, and treat these SNVs as having copy number neutral. Then cluster n/N for these SNVs using DP.
    set.seed(919919)
    
    nn2 <- jitter(nn)
    NN2 <- jitter(NN)
    nN2 <- nn2/NN2

    keep.id <- (nN2 < 1 & NN2 > quantile(NN2, 0.01))

    nN2 <- nN2[keep.id]
    MM_B2 <- MM_B[keep.id]
    
    if(sum(MM_B2==2) > 0)
    {
        selnN <- nN2[MM_B2==2]
    }else{
        N40 <- quantile(NN2, 0.40)
        N60 <- quantile(NN2, 0.60)
        selnN <- nN2[(NN2 > N40)& (NN2 < N60)]
    }
    ########################################################################
    
    ##########################   2.  Set up Dirichlet process mixture #####
    ###                                                                 ###
    nburn <- 500
    nsave <- 1000
    nskip <- 5
    ndisplay <- 500
    mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
    
    priordp <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
    psiinv2=solve(diag(0.5,1)),
    nu1=4,nu2=4,tau1=1,tau2=100)  ### Used the prior4 setting in the DPdensity example, which gives noninformative priors
    
    fit <- DPdensity(y=selnN,prior=priordp,mcmc=mcmc,
    state=state,status=TRUE)  ### Fit an initial DP MCMC iteration
    
    #########    3. Trick DP to save clustering membership every 10 iterations for 100 samples  #######
    ### Note: DPpackage does not automatically save the MCMC samples for clustering memberships    ###
    
    state.save <- NULL
    input.state <- fit$state$ss
    
    #### Change start. March 5, 2016 YJ ########
    n.sim <- 300
    mcmc1 <- list(nburn=100, nsave=10, nskip=1, ndisplay=10)  ##
    
    for(i in 1:n.sim){
        fit=DPdensity(y=selnN, state=input.state,prior=priordp, mcmc=mcmc1, status=TRUE)  ### Change end. March 5, 2016 YJ ###
        input.state=fit$state$ss
        state.save <- cbind(state.save, input.state)
        cat(i, " ")
    } #for(i in 1:n.sim){
    
    
    #### 4. For each MCMC sample, use the right cluster and middle but <0.5 cluster to estimate w-star. The right cluster corresponds to homozygous somatic mutations, and middle but <0.5 cluster corresponds to heterozygous somatic mutations.
    
    right.clust <- NULL
    mid.clust <- NULL
    right.size <- NULL
    mid.size <- NULL
    
    for(i in 1:n.sim){
        curr.state <- state.save[,i]
        max.c <- max(curr.state)
        
        curr.mean <- NULL
        for(j in 1:max.c){
            curr.mean <- c(curr.mean, mean(selnN[curr.state==j]))
        }
        
        max.mean1 <- max(curr.mean)  ### location for homozygous
        
        n_cnt <- 2
        mid_loc <- NA
        
        while((n_cnt <= 8)&(is.na(mid_loc)))  ###  at least (0.25, 0.125)
        {
            if(max.mean1 > (1.0/n_cnt))
            {
                mid_loc <- 1.0/n_cnt
            }else{
                n_cnt <- n_cnt*2.0
            }
        }
        
        if((!is.na(max.mean1))&(!is.na(mid_loc)))
        {
            right.ix <- which(curr.mean == max.mean1)
            right.size <- c(right.size, length(selnN[curr.state %in% right.ix])) ## the size of the right cluster -- homozygous mutations
            right.clust <- c(right.clust, max.mean1)  ## the mean of the right cluster as an estimate for 1-wstar
            
            ind.tmp <- ((curr.mean - mid_loc)<=0)
            
            if(sum(ind.tmp) > 0)  ## nothing is below mid_loc
            {
                max.mean2 <- max(curr.mean[ind.tmp])   ### location for heterozygous
                
                mid.ix <- which(curr.mean== max.mean2)
                mid.size <- c(mid.size, length(selnN[curr.state %in% mid.ix])) ## the size of the middle cluster -- heterozygous mutations
                mid.clust <- c(mid.clust, max.mean2) ## the mean of the middle cluster as an estimate for wstar/2
            }else{  ### if(sum(ind.tmp) > 0)
                
                right.size[i] <- NA
                right.clust[i] <- NA
                mid.size <- c(mid.size, NA)
                mid.clust <- c(mid.clust, NA)
            }  ###   if(sum(ind.tmp) > 0)
        }else{
            
            right.size <- c(right.size, NA)
            right.clust <- c(right.clust, NA)
            mid.size <- c(mid.size, NA)
            mid.clust <- c(mid.clust, NA)
        } ##if((!is.na(max.mean1))&(!is.na(mid_loc)))
        
    } ##for(i in 1:n.sim){
    
    w <- mid.size/(mid.size+right.size)  ## w-star is estiamted as a weighted average of the means of the middle and right clusters, with weight proportional to the cluster size
    w.star <- w*mid.clust*2 + (1-w)*right.clust
    cat("Estimated Purity is", 1-median(w.star, na.rm = TRUE))
    
    return(1-median(w.star, na.rm = TRUE)) ## return the posterior median
}



fn.elicit.w.star.B <- function(MM_B, NN, nn, SS)
{
    loci_set <- (1:SS)[(MM_B > 1.8)&(MM_B < 2.2)]
    #length(loci_set) > 0.3*SS
    p_set <- (nn/NN)[loci_set,]
    
    p_set1 <- p_set[(p_set < 0.5)]
    p_set1_cutoff <- quantile(p_set1, 0.8)
    w_star_1 <- 1 - 2*mean(p_set1[p_set1_cutoff < p_set1])
    
    
    p_set2 <- p_set[(p_set > 0.5)]
    p_set2_cutoff <- quantile(p_set2, 0.8)
    w_star_2 <- 1 - mean(p_set2[p_set2_cutoff < p_set2])
    
    w_star <- (length(p_set1)*w_star_1 + length(p_set2)*w_star_2)/(length(p_set1) + length(p_set2))
    
    return(w_star)
    
}



fn.elicit.w.star <- function(NN, nn, SS, TT)
{
    w_star <- rep(NA, TT)
    
    for(i_t in 1:TT)
    {
        N_t <- NN[,i_t]
        n_t <- nn[,i_t]
        
        m_N <- median(N_t)
        loci_set <- (1:SS)[(N_t > m_N*0.9)&(N_t < m_N*1.1)]
        #length(loci_set) > 0.3*SS
        p_set <- (n_t/N_t)[loci_set]
        
        
        p_set1 <- p_set[(p_set < 0.5)]
        p_set1_cutoff <- quantile(p_set1, 0.8)
        w_star_1 <- 1 - 2*mean(p_set1[p_set1_cutoff < p_set1])
        
        p_set2 <- p_set[(p_set > 0.5)]
        p_set2_cutoff <- quantile(p_set2, 0.8)
        w_star_2 <- 1 - mean(p_set2[p_set2_cutoff < p_set2])
        
        w_star[i_t] <- (length(p_set1)*w_star_1 + length(p_set2)*w_star_2)/(length(p_set1) + length(p_set2))
    }
    return(w_star)
}


# QQ <- ceiling(max(M_B)); NN <- N; nn <- n; SS <- S; TT <- T; MM_B <- M_B; ind_LM <- ind.LM; ind_DP <- ind.DP
fn.elicit.hyper.B <- function(QQ, NN, nn, SS, TT, MM_B, ind_LM, ind_DP)
{
    
    ###################################
    #HYPER-PARAMETER
    hpara <- NULL
    
    #NUMBER OF CELL TYPES (GEOMETRIC DIST)
    hpara$r <- 0.2
    
    #PRIOR FOR cIBP
    if(QQ > 2)
    {
        hpara$Q <- QQ  #NUMBER OF COPIES -- q = 0, 1, 2, 3
    }else{
        hpara$Q <- 2  #NUMBER OF COPIES -- q = 0, 1, 2, 3
    }
    ##BETA-DIRICHLET
    hpara$alpha <- 1
    hpara$beta <- 10
    hpara$gam <- rep(0.5, hpara$Q)
    
    #PRIOR FOR PHI--TOTAL NUMBER OF READS IN A SAMPLE T
    ##  one common gamma
    if(sum(MM_B==2) > 0)
    {
        phi_tmp <- median(NN[MM_B==2])
    }else{
        phi_tmp <- median(NN)
    }
    hpara$b <- 50
    hpara$a <- phi_tmp*hpara$b
    
    #PRIOR FOR P_zO  --- if p0 is *random*
    hpara$a_z0 <- 0.3*100
    hpara$b_z0 <- 5*100
    
    #PRIOR FOR W_star
    hpara$a_w <- hpara$b_w <- rep(0, TT)
    
    for(i_t in 1:TT)
    {
        w_star_tmp <- rep(NA, 2)
        
        if(ind_LM == 1)
        {
            tmp <- summary(lm(NN[,i_t]~MM_B[,i_t]))$coeff[,1]
            w_star_tmp1 <- tmp[1]/phi_tmp
            w_star_tmp2 <- 1-2*tmp[2]/phi_tmp
            w_star_tmp[1] <- (w_star_tmp1 + w_star_tmp2)/2
            
            print(paste("w_star estimate from LM   ", w_star_tmp[1]))
        }
        
        if(ind_DP == 1)
        {
            ##  USE DP to estimate w_star
            w_star_tmp[2] <- fn.get.wstar(NN[,i_t], nn[,i_t], MM_B[,i_t])
            
            print(paste("w_star estimate from DP  ", w_star_tmp[2]))
        }
        
        w_star_tmp <- mean(w_star_tmp, na.rm=TRUE)
        print(c(i_t, w_star_tmp))
        
        hpara$a_w[i_t] <- SS*phi_tmp*w_star_tmp
        hpara$b_w[i_t] <- SS*phi_tmp*(1-w_star_tmp)
    }
    
    ###PRIOR FOR M-B
    hpara$kappa <- rep(64, TT)
    
    #PRIOR FOR W_tilde
    hpara$d0 <- 0.1  ## if w0 is *random*
    hpara$d <- 1
    
    ##  use these for the model with fixed w0 and p0
    hpara$w0 <- 0.02
    hpara$p0_z <- 0.05
    
    return(hpara)
}



# Point.est <- point.est; min_C <- Min_C; max_C <- Max_C; TT <- T; removed_ind <- removed.ind;
fn.gen.ICGC.output.B <- function(outTxt, sample_id, Point.est, min_C, max_C, TT, removed_ind, chrName, chrPos)
{
    
    #### March 08 Subhajit
    #### March 04 Juhee
    out = outTxt
    
    SS_or <- length(removed_ind)  ## SS before removing
    
    ### 0. Multiplicity
    tumor_CN <- multi <- rep(NA, SS_or)
    tumor_CN[removed_ind==0] <- Point.est$m
    multi[removed_ind==0] <- Point.est$m_multi
    
    ## cbind(chr, pos, tumor_CN, multi)  ### tumor_copy number & multiplicity  --- matrix of S*T (orignal S)   ----- !!!!!!!!!!!!!!!! Subhajit
    multiplicity_file_name = sprintf("%s/%s_multiplicity.txt",out,sample_id)
    fid0 <- file(multiplicity_file_name, "w")
    cat("chr\tpos\ttumor_copynumber\tmultiplicity\n",file=fid0)
    for(i in 1:SS_or)
    {
        cat(chrName[i],"\t",chrPos[i],"\t",tumor_CN[i],"\t",multi[i],"\n",file=fid0)
    }
    close(fid0)
    
    ### commented 2 lines below @subhajit
    ##w.st <- read.table(w_star_file_name)  ### what are these for???  --Juhee
    ##putiry_median = w.st[1] ### should we use this one ?? ### what are these for???  --Juhee
    
    ###  1. Purity and Ploidy
    purity_ploidy_file_name = sprintf("%s/%s_purity_ploidy.txt",out,sample_id)
    fid1 <- file(purity_ploidy_file_name, "w")
    cat(sample_id,"\t", file=fid1)
	purity_sample = 1.0 - mean(Point.est$w_star_dist)
    cat(purity_sample,"\t", file=fid1)
    cat(mean(Point.est$ploidy_dist),"\n", file=fid1)
    close(fid1)
    
#     #### 2-A-I The number of subclones
#     C_star <- Point.est$C_star - 1    #### @subhajit: should be -1 right ??
#     ## c(sample neame, C_star)   --- scalor   ----- !!!!!!!!!!!!!!!! Subhajit  @@@@ Juhee C_star-1 ??
#     subclone_cluster_file_name = sprintf("%s/%s_number_of_clusters.txt",out,sample_id)
#     fid2 <- file(subclone_cluster_file_name, "w")
#     cat(sample_id,"\t", file=fid2)
#     cat(C_star,"\n", file=fid2)
#     close(fid2)
#     
#     #### 2-A-II Distribuiton of the number of subclones
#     subclone_cluster_distr_file_name = sprintf("%s/%s_cluster_distributions.txt",out,sample_id)
#     fid3 <- file(subclone_cluster_distr_file_name, "w")
#     cat("number-of-subclones(C)","\t", file=fid3)
#     
#     c.post <- Point.est$C_dist  ## March-03 Juhee
#     #c.post.max <- max(c.post)  ### find the posterior mode of C
#     ## ind.c<-as.numeric(row.names(c.post))
#     for (i in min_C:(max_C-1))
#     {
#         j = i-1
#         cat(j,"\t",file=fid3)
#     }
#     j = max_C - 1
#     cat(j,"\n",file=fid3)
#     
#     cat(sample_id,"\t",file=fid3)
#     for (i in 1:(max_C-2))
#     {
#         cat(c.post[i],"\t",file=fid3)
#     }
#     cat(c.post[max_C-1],"\n",file=fid3)
#     close(fid3)
#     
    #cat("subclone #", ind.c[which(c.post==c.post.max)]-1, c.post.max, "\n",file = fid2)  ## find the posterior probability that C equals the posterior mode
    
    #     for (i in 1:(max_C-2))
    #          cat(Point.est$C_dist[i],"\t",file=fid2)
    #     cat(Point.est$C_dist[max_C-1],"\n",file=fid2)
    
    
    #### 2-B
    #### commented by Subhajit as yitan's output will be used instead
#     Z_star <- Point.est$Z
#     cnt_Z <- (Z_star > 0)
#     cnt_Z <- apply(cnt_Z, 2, sum)
#     cnt_Z <- cnt_Z[-1]
#     w <- Point.est$w
#     for(i_t in 1:TT) w[i_t, -1] <- w[i_t, -1]/sum(w[i_t, -1])
#     w[,1] <- 0
#     
#     for(i in 1:(C_star -1))
#     {
#         print(c(i, cnt_Z[i], w[1, i+1]))  ## sample 1 only    ----- !!!!!!!!!!!!!!!! Subhajit
#     }
#     
    
    #### 2-C  ----- !!!!!!!!!!!!!!!! Subhajit
    #### commented by Subhajit as yitan's output will be used instead
#     i_s <- 1
#     for(i in 1:SS_or)
#     {
#         if(removed_ind[i]==1)
#         {
#             print(c("chr", "pos", NA))
#         }else{
#           
#           for(i_c in 2:C_star)
#           {
#               if(Z_star[i_s, i_c]>0)
#               {
#                  print(c("chr", "pos", i_c-1))
#                  
#               }
#           }# for(i_c in 2:C_star)
#           i_s <- i_s + 1
#         }# if(removed_ind[i]==1)
#     }  #for(i in 1:SS_or)
#     
    
    
}

fn.make.figures.B <- function(outFig, sample_id, Point.est, min_C, max_C, ind_fixed, NN, nn, MM_B, TT, QQ)
{
    out = outFig
    
    fName = sprintf("%s_B-CC-plots-data.pdf",sample_id)
    pdf(paste(out,fName, sep="/"))
    par(mar=c(4.5, 4.5, 2.1, 2.1))
    hist(NN, xlab="N", main="Histogram of N")
    hist(nn/NN, xlab="n/N", main="Histogram of n/N")
    plot(nn/NN, NN, xlab="n/N", ylab="N", main="Scatterplot of n/N vs. N")
    abline(h=median(NN), lwd=2, col=2, lty=2)
    plot(MM_B, NN, ylab="N", xlab="M_B", main="Scatterplot of M_B vs. N")
    dev.off()
    
    CC_star <- Point.est$C_star
    
    fName = sprintf("%s_B-CC-plot-C.pdf",sample_id)
    pdf(paste(out,fName, sep="/"))
    par(mar=c(4.5, 4.5, 2.1, 2.1))
    plot((0:(max_C-1)), Point.est$C_dist, type="o", lwd=4, col=1, xlab="C", ylab="P(C|n, N)", main="", cex.axis=1.5, cex.lab=1.5)
    abline(v=(CC_star-1), lwd=2, lty=2, col=2)
    dev.off()
    
    fName = sprintf("%s_B-CC-hist-w-star.pdf",sample_id)
    pdf(paste(out, fName, sep="/"))
    par(mar=c(4.5, 4.5, 2.1, 2.1))
    hist(Point.est$w_star_dist, main="", xlab="w_star conditional on C^*")
    abline(v=mean(Point.est$w_star_dist), lwd=2, lty=2, col=2)
    dev.off()
    
    
    if(ind_fixed==0)
    {
        fName = sprintf("%s_B-CC-hist-w0-p0.pdf",sample_id)
        pdf(paste(out, fName,sep="/"))
        par(mar=c(4.5, 4.5, 2.1, 2.1))
        hist(Point.est$w0_dist, main="", xlab="random w0 conditional on C^*")
        abline(v=mean(Point.est$w0_dist), lwd=2, lty=2, col=2)
        
        hist(Point.est$p0_z_dist, main="", xlab="random p0 conditional on C^*")
        abline(v=mean(Point.est$p0_z_dist), lwd=2, lty=2, col=2)
        
        dev.off()
    }
    
    fName = sprintf("%s_B-CC-hist-ploidy.pdf",sample_id)
    pdf(paste(out, fName, sep="/"))
    par(mar=c(4.5, 4.5, 2.1, 2.1))
    for(i_t in 1:TT)
    {
        hist(Point.est$ploidy_dist[, i_t], xlab="ploidy", main=paste("sample", i_t))
        abline(v=mean(Point.est$ploidy_dist[,i_t]), lwd=2, lty=2, col=2)
    }
    dev.off()
    
    ###########HISTOGRAMS OF DIFFERENCES
    fName = sprintf("%s_B-CC-histogram-diff.pdf",sample_id)
    pdf(paste(out, fName,sep="/"))
    par(mar=c(4.5, 4.5, 2.1, 2.1))
    hist(Point.est$N - NN, main="", col="grey", cex.axis=1.5, cex.lab=1.5, lwd=4, xlab="N_st_hat - N_st")
    abline(v=0, lty=2, col=2, lwd=2)
    
    par(mar=c(4.5, 4.5, 2.1, 2.1))
    hist(Point.est$p - (nn/NN), main="", col="grey", cex.axis=1.5, cex.lab=1.5, lwd=4, xlab="p_st_hat - (n_st/N_st)", nclass=15)
    abline(v=0, lty=2, col=2, lwd=2)
    dev.off()
    
    
    if(CC_star > 2)
    {
        lmat = rbind(c(0,3),c(2,1),c(0,4))
        lwid = c(1.5,4)
        lhei = c(1.5,4,1.2)
        
        colnames(Point.est$L) <- (0:(CC_star-1))
        colnames(Point.est$Z) <- (0:(CC_star-1))
        
        
        #MATCHING THE COLOER KEYS --- MEANS (NOT ONE SAMPLE)
        fName = sprintf("%s_CC_L_Z.pdf",sample_id)
        pdf(paste(out, fName,sep="/"))
        tmp <- heatmap.2(cbind(Point.est$L[,-1], Point.est$Z[,-1]), trace="none", Colv=FALSE, Rowv=TRUE, dendrogram="none", scale="none", col=redgreen, colsep=(1:CC_star), sepcol=c("white", "white"), sepwidth=c(0.05, 0.1), main="",lmat = lmat, lwid = lwid, lhei = lhei)
        dev.off()
        
        #### fix me:  some problem in generating individual heatmap for L and Z
        
#         post.L.new <- Point.est$L[tmp$rowInd,]
#         post.Z.new <- Point.est$Z[tmp$rowInd,]
#         
#         rownames(post.L.new) <- tmp$rowInd
#         rownames(post.Z.new) <- tmp$rowInd
#         
#         
#         rg <- redgreen(100)      # the original color vector
#         min.c <- 0
#         max.c <- QQ
#         
#         rg.L <- rg[round((min(post.L.new)/(max.c - min.c)*100)):round(100*max(post.L.new)/(max.c - min.c))]
#         rg.Z <- rg[round((min(post.Z.new)/(max.c - min.c)*100)):round(100*max(post.Z.new)/(max.c - min.c))]
#         
#         fName = sprintf("%s_B-chk-heatmap-L.pdf",sample_id)
#         pdf(paste(out, fName,sep="/"))
#         heatmap.2(post.L.new[,-1], trace="none", Colv=FALSE, Rowv=FALSE, dendrogram="none", scale="none", col=rg.L, colsep=(1:(CC_star-1)), sepcol=c("white", "white"), sepwidth=c(0.05, 0.1), main="",lmat = lmat, lwid = lwid, lhei = lhei)
#         dev.off()
#         
#         fName = sprintf("%s_B-chk-heatmap-Z.pdf",sample_id)
#         pdf(paste(out, fName, sep="/"))
#         heatmap.2(post.Z.new[,-1], trace="none", Colv=FALSE, Rowv=FALSE, dendrogram="none", scale="none", col=rg.Z, colsep=(1:(CC_star-1)), sepcol=c("white", "white"), sepwidth=c(0.05, 0.1), main="",lmat = lmat, lwid = lwid, lhei = lhei)
#         dev.off()
        
        if(TT == 1)
        {
            if(ind_fixed==0)
            {
                max.w <- max(max(Point.est$w_tilde), max(Point.est$w))
                min.w <- min(min(Point.est$w_tilde), min(Point.est$w))
                
                fName = sprintf("%s_B-CC-plot-w.pdf",sample_id)
                pdf(paste(out, fName, sep="/"))
                par(mar=c(4.5, 4.5, 2.1, 2.1))
                plot((0:(CC_star-1)), Point.est$w_tilde, type="o", lwd=4, col=1, xlab="C", ylab="point w", main="", cex.axis=1.5, cex.lab=1.5, ylim=c(min.w, max.w))
                lines((0:(CC_star-1)), Point.est$w, type="o", lwd=4, col=2, lty=2)
                dev.off()
                
            }else{
                max.w <- max(max(Point.est$w_tilde[,-1]), max(Point.est$w[,-1]))
                min.w <- min(min(Point.est$w_tilde[,-1]), min(Point.est$w[,-1]))
                
                fName = sprintf("%s_B-CC-plot-w.pdf",sample_id)
                pdf(paste(out, fName, sep="/"))
                par(mar=c(4.5, 4.5, 2.1, 2.1))
                plot((1:(CC_star-1)), Point.est$w_tilde[, -1], type="o", lwd=4, col=1, xlab="C", ylab="point w", main="", cex.axis=1.5, cex.lab=1.5, ylim=c(min.w, max.w))
                lines((1:(CC_star-1)), Point.est$w[,-1], type="o", lwd=4, col=2, lty=2)
                dev.off()
            }
        }
        
    }else{
        print("# of subclones is 1. Can't make heatmaps!")
    }
    
}




