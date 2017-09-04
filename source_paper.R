
if (case == "setup") {
  cl_export = c("SL.gam3","screen.Main","screen10","screen6","SL.glmnet_1","SL.glmnet_2",
                "SL.glmnet_3","xgbFull","xgbMain","screen.Main",
                "xgb6","screen6","xgb10","screen10","rpartPrune","nnetFull",
                "nnetMain","screen.Main","earthFull","All","screen10","screen6",
                "earthMain","screen.Main","rangerFull","All","screen.Main",
                "ranger10","screen10","screen6","SL.glm","screen6","screen10",
                "SL.stepAIC","SL.hal","glm.mainint")
  
  SL.library1 = list(c("SL.gam3","screen.Main","screen6","screen10","All"),"SL.glmnet_1",
                     "SL.glmnet_2","SL.glmnet_3", c("SL.rpartPrune", "screen.Main"), 
                     c("nnetMain","screen.Main"), 
                     c("earthMain","screen.Main"),
                     c("SL.glm","screen.Main","screen6","screen10","All"),
                     "SL.stepAIC", c("SL.hal","screen.Main"),"SL.mean","glm.mainint")
  
  SL.library2 = list(c("SL.gam3","screen.Main","screen6","screen10","All"),
                     "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3",
                     c("SL.rpartPrune", "screen.Main"),"xgbFull",c("xgbMain","screen.Main"),
                     c("nnetMain","screen.Main"), c("earthMain","screen.Main"),
                     c("rangerFull","screen.Main"),
                     c("SL.glm","screen.Main","screen6","screen10","All"),
                     "SL.stepAIC", c("SL.hal","screen.Main"),"SL.mean","glm.mainint")
  
  SL.libraryG = list("nnetMain","SL.mean","SL.hal",
                     "SL.earth","SL.glm","SL.step.interaction",
                     "SL.glm.interaction")
} else {
  if (case == "LRcase2a") {
    g0 = g0_linear
    Q0 = Q0_trig1
    testdata=gendata(1e6, g0 = g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    Qform = paste(colnames(gendata(1,g0 = g0, Q0 = Q0))[2:5], collapse = "+")
    Qform = paste0("Y ~ A*(", Qform, ")")
    
    # SL.library = c("glm.mainint")
    # SL.libraryG = c("SL.glm")
    
    if (!resultsGotten) {
      cl = makeCluster(no.cores, type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,cl_export)
      
      # run this on a 24 core node
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                  .errorhandling = "remove")%dopar%
                  {sim_hal(n, g0, Q0, gform=formula("A~."), 
                           Qform = formula(Qform), V=10)}
      
      
    }
    
    results = data.matrix(data.frame(do.call(rbind, ALL)))
    results = results_LRcase2a
    B = nrow(results)
    
    type= c(rep("tmle with lr initial",B),rep("lr initial",B))
    types = c("tmle with lr initial","lr initial")
    inds = c(1,22)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    colors = c("red","blue")
    
    plotdf = data.frame(ests=ests, type=type)
    
    ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle("blip variance sampling distributions", subtitle=
                "tmle with logistic regression initial estimates")
    ggover = ggover+geom_vline(xintercept = var0,color="black")
    assign(paste0("gg_",case), ggover)
    assign(paste0("results_",case), results)
  }
  
  if (case == "LRcase2b"){
    g0 = g0_linear
    Q0 = Q0_trig
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    Qform = paste(colnames(gendata(1,g0 = g0, Q0 = Q0))[2:5], collapse = "+")
    Qform = paste0("Y ~ A*(", Qform, ")")
    
    # SL.library = c("glm.mainint")
    # SL.libraryG = c("SL.glm")
    
    if (!resultsGotten) {
      cl = makeCluster(no.cores, type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,cl_export)
      
      # run this on a 24 core node
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                  .errorhandling = "remove")%dopar%
                  {sim_hal(n, g0, Q0, gform=formula("A~."), 
                           Qform = formula(Qform), V=10)}
      
    }
    
    results = data.matrix(data.frame(do.call(rbind, ALL)))
    B = nrow(results)
    type= c(rep("tmle with lr initial",B),rep("lr initial",B))
    types = c("tmle with lr initial","lr initial")
    inds = c(1,22)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    colors = c("red","blue")
    
    plotdf = data.frame(ests=ests, type=type)
    
    ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle("blip variance sampling distributions", subtitle=
                "tmle with logistic regression initial estimates")
    ggover = ggover+geom_vline(xintercept = var0,color="black")
    assign(paste0("gg_",case), ggover)
    assign(paste0("results_",case), results)
  }
  
  ####
  ####
  ####
  # hal
  if (case == "HALcase2a") {
    
    g0 = g0_linear
    Q0 = Q0_trig1
    testdata=gendata(1e6, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    SL.library = list(c("SL.hal","screen.Main"))
    SL.libraryG = c("SL.glm")
    
    if (!resultsGotten) {
      cl = makeCluster(no.cores, type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,cl_export)
      # run this on a 24 core node
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                  .errorhandling = "remove")%dopar%
                  {sim_hal(n, g0 = g0, Q0 = Q0, gform = formula("A~."))}
      
      
    }
    
    results = data.matrix(data.frame(do.call(rbind, ALL)))
    B = nrow(results)
    
    varind = c("1step tmle HAL" = 1,"init est HAL" = 22)
    
    type= c(rep("TMLE HAL",B), rep("initial est HAL",B))
    types = c("TMLE HAL", "initial est HAL")
    inds = varind
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    
    inds = inds[order(types)]
    colors = c("red","blue")
    
    plotdf = data.frame(ests=ests, type=type)
    
    ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle("blip variance sampling distributions", subtitle=
                "tmle with hal initial estimates")
    ggover = ggover+geom_vline(xintercept = var0,color="black")+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])
    
    capt = "Truth is at black vline."
    ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", 
                          fontface = "plain",colour = "black", size = 10, angle = 0, 
                          lineheight = 0.9))
    assign(paste0("gg_",case), ggover)
    
    coverage = vapply(varind[1], 
                      FUN = function(x) cov.check(results, var0, x),  
                      FUN.VALUE = 1, USE.NAMES = TRUE)
    
    SE_true = sd(results[,varind[1]])*sqrt(B-1)/sqrt(B)
    ci_true = data.frame(results[,varind[1]], 
                         results[,varind[1]] - 1.96*SE_true,
                         results[,varind[1]] + 1.96*SE_true)
    cov_true = cov.check(ci_true, var0, 1)
    coverage = c(coverage, cov_true)
    names(coverage)[2] = paste0(names(varind)[1], " TRUE VAR")
    
    performance.sig = t(apply(results[,varind], 2, perf,var0))
    rownames(performance.sig) = names(varind)
    
    assign(paste0("results_",case), results)
    assign(paste0("performance.sig_",case), performance.sig)
    assign(paste0("coverage_",case), coverage)
  }
  
  if (case == "HALcase2b") {
    g0 = g0_linear
    Q0 = Q0_trig
    testdata=gendata(1e6, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    SL.library = list(c("SL.hal","screen.Main"))
    SL.libraryG = c("SL.glm")
    
    if (!resultsGotten) {
      cl = makeCluster(no.cores, type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,cl_export)
      # run this on a 24 core node
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                  .errorhandling = "remove")%dopar%
                  {sim_hal(n, g0 = g0, Q0 = Q0, gform = formula("A~."))}
      
      
    }
    
    results = data.matrix(data.frame(do.call(rbind, ALL)))
    B = nrow(results)
    
    varind = c("1step tmle HAL" = 1,"init est HAL" = 22)
    
    type= c(rep("TMLE HAL",B), rep("initial est HAL",B))
    types = c("TMLE HAL", "initial est HAL")
    inds = varind
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    
    inds = inds[order(types)]
    colors = c("red","blue")
    
    plotdf = data.frame(ests=ests, type=type)
    
    ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle("blip variance sampling distributions", subtitle=
                "tmle with hal initial estimates")
    ggover = ggover+geom_vline(xintercept = var0,color="black")+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])
    
    capt = "Truth is at black vline."
    ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", 
                          fontface = "plain",colour = "black", size = 10, angle = 0, 
                          lineheight = 0.9))
    assign(paste0("gg_",case), ggover)
    
    coverage = vapply(varind[1], 
                      FUN = function(x) cov.check(results, var0, x),  
                      FUN.VALUE = 1, USE.NAMES = TRUE)
    
    SE_true = sd(results[,varind[1]])*sqrt(B-1)/sqrt(B)
    ci_true = data.frame(results[,varind[1]], 
                         results[,varind[1]] - 1.96*SE_true,
                         results[,varind[1]] + 1.96*SE_true)
    cov_true = cov.check(ci_true, var0, 1)
    coverage = c(coverage, cov_true)
    names(coverage)[2] = paste0(names(varind)[1], " TRUE VAR")
    
    performance.sig = t(apply(results[,varind], 2, perf,var0))
    rownames(performance.sig) = names(varind)
    
    assign(paste0("results_",case), results)
    assign(paste0("performance.sig_",case), performance.sig)
    assign(paste0("coverage_",case), coverage)  
    }
  
  if (case == "case2a") {
    
    g0 = g0_linear
    Q0 = Q0_trig1
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    SL.library = SL.library1
    SL.libraryG = list("SL.glm")
    
    if (!resultsGotten) {
      cl = makeCluster(no.cores, type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,cl_export)
      
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                  .errorhandling = "remove")%dopar%
                  {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                          SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE
                  )}
      }
    
    results = data.matrix(data.frame(do.call(rbind, ALL)))
    B = nrow(results)
    
    varind = c("1step tmle" = 1,"simultaneous tmle" = 7, "init est" = 37)
    ateind = c("1step tmle" = 28,"simultaneous tmle" = 25, "init est" = 39)
    
    performance.sig = t(apply(results[,varind], 2, perf,var0))
    performance.ate = t(apply(results[,ateind], 2, perf,ATE0))
    
    rownames(performance.sig) = names(varind)
    rownames(performance.ate) = names(ateind)
    
    coverage = c(cov.check(results, var0, 1),
                 cov.simul(results, c(var0, ATE0), c(7,25)),
                 cov.check(results, ATE0, 28), cov.check(results, ATE0, 34))
    
    # getting coveage using the true variance of the estimator
    dd = data.frame(psi = results[,1],l = results[,1]-1.96*sqrt(performance.sig[1,1]),
                    r = results[,1]+1.96*sqrt(performance.sig[1,1]), psi = results[,28],
                    l = results[,28]-1.96*sqrt(performance.ate[1,1]),
                    r = results[,28]+1.96*sqrt(performance.ate[1,1]))
    
    cov.sig.1step = cov.check(dd, var0, 1)
    cov.ate = cov.check(dd, ATE0, 4)
    cov = c(coverage, cov.sig.1step, cov.ate)
    names(cov) = c("TMLE Blip Variance SL1", "Simultaneous TMLE SL1", "TMLE ATE SL1",
                   "TMLE ATE LR", "TMLE Blip Var using true Var", "TMLE ATE using true Var")
    coverage = data.frame(coverage = cov)
    rownames(coverage) = names(cov)
    
    # getting superlearner results
    LL = 0
    for (i in 1:length(SL.library)) {
      if (length(SL.library[[i]]) > 1) {
        LL = LL + length(SL.library[[i]])-1} else
        {LL = LL + 1}
    }
    SL_results = data.frame(colMeans(results[,65:(64+LL)]))
    rownames(SL_results) = colnames(results)[65:(64+LL)]
    
    type = c(rep("TMLE SL1",B), rep("init. est SL1",B),
             rep("TMLE LR",B), rep("init est. LR",B))
    types = c("simul TMLE SL1","iter TMLE SL1","init. est SL1","init est. LR")
    inds = c(25, 39, 34, 40)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("red", "blue","green","orange")
    
    ateests = data.frame(ests=ests,type=type)
    ggover = ggplot(ateests,aes(ests, fill=type,color=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("ATE sampling distributions, ", case))
    ggover = ggover+geom_vline(xintercept = ATE0,color="black")+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    capt = paste0("Truth is at black vline.\n",
                  "\ninitial est using SuperLearner library 1 is excellent", 
                  "\nTMLE for ATE using IC approx for variance covers at ", 
                  100*round(coverage[3,1],3),"%",
                  "\ninit est. LR which is Logistic Reg. with main terms and interactions ",
                  "\nplug-in clearly biased but TMLE corrects it and it covers at ",
                  100*round(coverage[4,],3),"%\n")
    ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", 
                          fontface = "plain",colour = "black", size = 10, angle = 0, 
                          lineheight = 0.9))
    assign(paste0("gg_ATE", case), ggover)
    
    
    type = c(rep("TMLE SL1",B), rep("init est SL1",B), 
             rep("init est LR",B), rep("TMLE LR",B))
    types = c("TMLE SL1","init est SL1","init est LR", "TMLE LR")
    inds = c(1, 37, 38, 16)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("blue","green","orange","red")
    varests = data.frame(ests=ests,type=type)
    
    ggover2 = ggplot(varests,aes(ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("Blip Variance sampling distributions, ", case))
    ggover2 = ggover2+geom_vline(xintercept = var0,color="black")+
      theme(plot.title = element_text(size=12), 
            plot.subtitle = element_text(size=10))+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    
    cap = paste0("truth is black line\n",
                 "tmle SL1, which used Superlearner Library 1 for initial ests\n", 
                 "attains near nominal coverage at ", 100*round(coverage[1,1],3),"%\n", 
                 "and debiases initial SuperLearner lib 1 estimate.\n", 
                 "tmle LR used logistic regression with main terms and interactions\n",
                 "has very biased initial estimates and leads to bad targeting for computing\n",
                 "the accompanying TMLE.")
    ggover2=ggdraw(add_sub(ggover2,cap, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 10, angle = 0, 
                           lineheight = 0.9))
    assign(paste0("gg_BV", case), ggover2)
    
    assign(paste0("results_",case), results)
    assign(paste0("performance.sig_",case), performance.sig)
    assign(paste0("performance.ate_",case), performance.ate)
    assign(paste0("coverage_",case), coverage)
    assign(paste0("SL_results_",case), SL_results)
  }
  
  if (case == "combo_LRandSL1case2a") {
    
    g0 = g0_linear
    Q0 = Q0_trig1
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    SL.library = SL.library1
    
    SL.libraryG = list("SL.glm")
    
    varind = c("1step tmle LR" = 1,"1step tmle SL1" = 5,
               "init est LR" = 22,"init est tmle SL1" = 41)
    
    performance.sig = lapply(results[varind],perf,var0)
    performance.sig = t(as.data.frame(performance.sig))
    rownames(performance.sig) = names(varind)
    
    coverage_tmle = cov.check(results_case2a, var0, 1)
    coverage_tmle
    coverage_tmlesimul = cov.simul(results_case2a, c(ATE0, var0), c(25,7))
    coverage_tmlesimul
    coverage_LR = cov.check(results_LRcase2a, var0, 1)
    coverage_LR
    coverage_LRsimul = cov.simul(results_LRcase2a, c(ATE0, var0), c(16,7))
    coverage_LRsimul
    coverage = c(coverage_halglm, coverage_hal, .434, NA, NA)
    
    cov.check(results_LRcase2a, ATE0, 16)
    MSE_cov = cbind(performance.sig, coverage)
    MSE_cov
    # getting superlearner results
    LL = 0
    for (i in 1:length(SL.library)) {
      if (length(SL.library[[i]]) > 1) {
        LL = LL + length(SL.library[[i]])-1} else
        {LL = LL + 1}
    }
    SL_results = data.frame(colMeans(results_halglm[,65:(65+LL-1)]))
    rownames(SL_results) = colnames(results)[65:(65+LL-1)]
    
    LG=0
    for (i in 1:length(SL.libraryG)) {
      if (length(SL.libraryG[[i]]) > 1) {
        LG = LG + length(SL.libraryG[[i]])-1} else
        {LG = LG + 1}
    }
    SL_resultsG = data.frame(colMeans(results_halglm[,(65+LL):(65+LL+LG-1)]))
    rownames(SL_results) = colnames(results)[(65+LL):(65+LL+LG-1)]
    
    B = nrow(results_halglm)
    B1 = nrow(results_hal)
    type = c(rep(names(varind)[1],B), rep(names(varind)[2],B1), rep(names(varind)[3],B),
             rep(names(varind)[4],B), rep(names(varind)[5],B1), rep(names(varind)[6],B))
    # types = c("1step TMLE LR","1step TMLE HAL","1step TMLE HAL+glm")
    types = names(varind)
    inds = varind
    ests = unlist(lapply(inds, FUN = function(x) results[[x]]))
    inds = inds[order(types)]
    
    colors = c("blue","green","orange","red", "purple", "yellow")
    varests = data.frame(ests=ests,type=type)
    
    ggover2 = ggplot(varests,aes(ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("Blip Variance sampling distributions, ", case))
    ggover2 = ggover2+geom_vline(xintercept = var0,color="black")+
      theme(plot.title = element_text(size=12), 
            plot.subtitle = element_text(size=10))+
      geom_vline(xintercept=mean(results[[inds[1]]]),color = colors[1])+
      geom_vline(xintercept=mean(results[[inds[2]]]),color = colors[2])+
      geom_vline(xintercept=mean(results[[inds[3]]]),color = colors[3])+
      geom_vline(xintercept=mean(results[[inds[4]]]),color = colors[4])+
      geom_vline(xintercept=mean(results[[inds[5]]]),color = colors[5])+
      geom_vline(xintercept=mean(results[[inds[6]]]),color = colors[6])
    cap = paste0("truth at black line.\n",
                 "tmle LR uses glm with interactions for outcome model and glm for\n",
                 "treatment mechanism initial estimates. tmle LR CI's cover at", 
                 round(MSE_cov[3,4] ,1),"\n",
                 "hal tmle uses highly adaptive lasso for initial estimates of both\n",
                 "outcome and treatment mech initial estimates and these cover\n",
                 "at", round(MSE_cov[3,4] ,2),
                 "\ntmle hal+glm SL uses a SuperLearner with hal and glm for\n",
                 "outcome and treatment mechanism initial estimates and cover at ",
                 round(MSE_cov[3,4] ,3),"\n",
                 "All coverage here is using infuence curve approx for inference.")
    ggover2=ggdraw(add_sub(ggover2,cap, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 10, angle = 0, 
                           lineheight = 0.9))
    assign(paste0("results_",case), results)
    assign(paste0("gg_BV",case), ggover2)
    assign(paste0("MSE_cov_",case), MSE_cov)
    assign(paste0("SL_results_",case), SL_results) 
    assign(paste0("SL_resultsG_",case), SL_resultsG) 
  }
  
  if (case == "case2bSL1") {
    g0 = g0_linear
    Q0 = Q0_trig
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    SL.library = SL.library1
    SL.libraryG = list("SL.glm")
    
    if (!resultsGotten) {
      cl = makeCluster(no.cores, type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,cl_export)
      
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                  .errorhandling = "remove")%dopar%
                  {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                          SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE
                  )}
      results = data.matrix(data.frame(do.call(rbind, ALL)))
    }
    
    B = nrow(results)
    
    varind = c("1step tmle" = 1,"simultaneous tmle" = 7, "init est" = 37)
    ateind = c("1step tmle" = 28,"simultaneous tmle" = 25, "init est" = 39)
    
    performance.sig = t(apply(results[,varind], 2, perf,var0))
    performance.ate = t(apply(results[,ateind], 2, perf,ATE0))
    
    rownames(performance.sig) = names(varind)
    rownames(performance.ate) = names(ateind)
    
    coverage = c(cov.check(results, var0, 1),
                 cov.simul(results, c(var0, ATE0), c(7,25)),
                 cov.check(results, ATE0, 28), cov.check(results, ATE0, 34))
    
    # getting coveage using the true variance of the estimator
    dd = data.frame(psi = results[,1],l = results[,1]-1.96*sqrt(performance.sig[1,1]),
                    r = results[,1]+1.96*sqrt(performance.sig[1,1]), psi = results[,28],
                    l = results[,28]-1.96*sqrt(performance.ate[1,1]),
                    r = results[,28]+1.96*sqrt(performance.ate[1,1]))
    
    cov.sig.1step = cov.check(dd, var0, 1)
    cov.ate = cov.check(dd, ATE0, 4)
    cov = c(coverage, cov.sig.1step, cov.ate)
    names(cov) = c("TMLE Blip Variance SL1", "Simultaneous TMLE SL1", "TMLE ATE SL1",
                   "TMLE ATE LR", "TMLE Blip Var using true Var", "TMLE ATE using true Var")
    coverage = data.frame(coverage = cov)
    rownames(coverage) = names(cov)
    
    # getting superlearner results
    LL = 0
    for (i in 1:length(SL.library)) {
      if (length(SL.library[[i]]) > 1) {
        LL = LL + length(SL.library[[i]])-1} else
        {LL = LL + 1}
    }
    SL_results = data.frame(colMeans(results[,65:(64+LL)]))
    rownames(SL_results) = colnames(results)[65:(64+LL)]
    
    type = c(rep("TMLE SL1",B), rep("init. est SL1",B),
             rep("TMLE LR",B), rep("init est. LR",B))
    types = c("simul TMLE SL1","iter TMLE SL1","init. est SL1","init est. LR")
    inds = c(25, 39, 34, 40)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("red", "blue","green","orange")
    
    ateests = data.frame(ests=ests,type=type)
    ggover = ggplot(ateests,aes(ests, fill=type,color=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("ATE sampling distributions, ", case))
    ggover = ggover+geom_vline(xintercept = ATE0,color="black")+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    capt = paste0("Truth is at black vline.\n",
                  "\ninitial est using SuperLearner library 1 is excellent", 
                  "\nTMLE for ATE using IC approx for variance covers at ", 
                  100*round(coverage[3,1],3),"%",
                  "\ninit est. LR which is Logistic Reg. with main terms and interactions ",
                  "\nplug-in clearly biased but TMLE corrects it and it covers at ",
                  100*round(coverage[4,],3),"%\n")
    ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", 
                          fontface = "plain",colour = "black", size = 10, angle = 0, 
                          lineheight = 0.9))
    assign(paste0("gg_ATE",case), ggover)
    
    
    type = c(rep("TMLE SL1",B), rep("init est SL1",B), rep("init est LR",B))
    types = c("TMLE SL1","init est SL1","init est LR")
    inds = c(1, 37, 38)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("blue","green","orange","red")
    varests = data.frame(ests=ests,type=type)
    
    ggover2 = ggplot(varests,aes(ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("Blip Variance sampling distributions, ", case))
    ggover2 = ggover2+geom_vline(xintercept = var0,color="black")+
      theme(plot.title = element_text(size=12), 
            plot.subtitle = element_text(size=10))+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])
    cap = paste0("truth is black line\n",
                 "tmle SL1, which used Superlearner Library 1 for initial ests\n", 
                 "attains near nominal coverage at ", 100*round(coverage[1,1],3),"%\n", 
                 "and debiases initial SuperLearner lib 1 estimate.\n", 
                 "init est LR used logistic regression with main terms and interactions\n",
                 "plug-in estimator and is a disaster, which TMLE cannot help.\n")
    ggover2=ggdraw(add_sub(ggover2,cap, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 10, angle = 0, 
                           lineheight = 0.9))
    
    assign(paste0("gg_BV",case), ggover2)
    assign(paste0("results_",case), results)
    assign(paste0("performance.sig_",case), performance.sig)
    assign(paste0("performance.ate_",case), performance.ate)
    assign(paste0("coverage_",case), coverage)
    assign(paste0("SL_results_",case), SL_results)
  }
  
  if (case == "case2bCVSL2") {
    B = nrow(results)
    varind = c("1step tmle" = 1,"simultaneous tmle" = 7, "init est" = 37)
    ateind = c("1step tmle" = 28,"simultaneous tmle" = 25, "init est" = 39)
    
    performance.sig = t(apply(results[,varind], 2, perf,var0))
    performance.ate = t(apply(results[,ateind], 2, perf,ATE0))
    
    rownames(performance.sig) = names(varind)
    rownames(performance.ate) = names(ateind)
    
    coverage = c(cov.check(results, var0, 1),
                 cov.simul(results, c(var0, ATE0), c(7,25)),
                 cov.check(results, ATE0, 28), cov.check(results, ATE0, 34))
    
    # getting coveage using the true variance of the estimator
    dd = data.frame(psi = results[,1],l = results[,1]-1.96*sqrt(performance.sig[1,1]),
                    r = results[,1]+1.96*sqrt(performance.sig[1,1]), psi = results[,28],
                    l = results[,28]-1.96*sqrt(performance.ate[1,1]),
                    r = results[,28]+1.96*sqrt(performance.ate[1,1]))
    
    cov.sig.1step = cov.check(dd, var0, 1)
    cov.ate = cov.check(dd, ATE0, 4)
    cov = c(coverage, cov.sig.1step, cov.ate)
    names(cov) = c("TMLE Blip Variance SL2", "Simultaneous TMLE SL2", "TMLE ATE SL2",
                   "TMLE ATE LR", "TMLE Blip Var using true Var", "TMLE ATE using true Var")
    coverage = data.frame(coverage = cov)
    rownames(coverage) = names(cov)
    
    # getting superlearner results
    LL = 0
    for (i in 1:length(SL.library)) {
      if (length(SL.library[[i]]) > 1) {
        LL = LL + length(SL.library[[i]])-1} else
        {LL = LL + 1}
    }
    SL_results = data.frame(colMeans(results[,65:(64+LL)]))
    rownames(SL_results) = colnames(results)[65:(64+LL)]
    
    type = c(rep("TMLE SL2",B), rep("init. est SL2",B),
             rep("TMLE LR",B), rep("init est. LR",B))
    types = c("simul TMLE SL2","iter TMLE SL2","init. est SL2","init est. LR")
    inds = c(25, 39, 34, 40)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("red", "blue","green","orange")
    
    ateests = data.frame(ests=ests,type=type)
    ggover = ggplot(ateests,aes(ests, fill=type,color=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("ATE sampling distributions, ", case))
    ggover = ggover+geom_vline(xintercept = ATE0,color="black")+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    capt = paste0("Truth is at black vline.\n",
                  "\ninitial est using SuperLearner library 2 is excellent", 
                  "\nTMLE for ATE using IC approx for variance covers at ", 
                  100*round(coverage[3,1],3),"%",
                  "\ninit est. LR which is Logistic Reg. with main terms and interactions ",
                  "\nplug-in clearly biased but TMLE corrects it and it covers at ",
                  100*round(coverage[4,],3),"%\n")
    ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", 
                          fontface = "plain",colour = "black", size = 10, angle = 0, 
                          lineheight = 0.9))
    assign(paste0("gg_ATE",case), ggover)
    
    
    type = c(rep("TMLE SL2",B), rep("init est SL2",B), rep("init est LR",B))
    types = c("TMLE SL2","init est SL2","init est LR")
    inds = c(1, 37, 38)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("blue","green","orange","red")
    varests = data.frame(ests=ests,type=type)
    
    ggover2 = ggplot(varests,aes(ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("Blip Variance sampling distributions, ", case))
    ggover2 = ggover2+geom_vline(xintercept = var0,color="black")+
      theme(plot.title = element_text(size=12), 
            plot.subtitle = element_text(size=10))+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])
    cap = paste0("truth is black line\n",
                 "tmle SL2, which used Superlearner Library 2 for initial ests\n", 
                 "is normally dist, hasn't many outliers, covers at ", 100*round(coverage[1,1],3),"%\n", 
                 "and slightly debiases initial estimate using Superlearner Lib 2 \n",
                 "which defies the donsker condition but such is not needed for\n",
                 "CV-TMLE. The TMLE coverage using the True variance was ",
                 100*round(coverage[5,1], 3), "%, \nso the IC approx was a bit anticonservative.\n",
                 "init est LR used logistic regression with main terms and interactions\n",
                 "plug-in estimator and is a disaster, which TMLE cannot help.\n")
    ggover2=ggdraw(add_sub(ggover2,cap, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 10, angle = 0, 
                           lineheight = 0.9))
    
    assign(paste0("gg_BV",case), ggover2)
    assign(paste0("results_",case), results)
    assign(paste0("performance.sig_",case), performance.sig)
    assign(paste0("performance.ate_",case), performance.ate)
    assign(paste0("coverage_",case), coverage)
    assign(paste0("SL_results_",case), SL_results) 
  }
  
  if (case == "case2bSL2") { 
    
    g0 = g0_linear
    Q0 = Q0_trig
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    SL.library = SL.library2
    SL.libraryG = list("SL.glm")
    
    if (!resultsGotten) {
      cl = makeCluster(no.cores, type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,cl_export)
      
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                  .errorhandling = "remove")%dopar%
                  {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                          SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE
                  )}
      results = data.matrix(data.frame(do.call(rbind, ALL)))
    }
    
    B = nrow(results)
    
    varind = c("1step tmle" = 1,"simultaneous tmle" = 7, "init est" = 37)
    ateind = c("1step tmle" = 28,"simultaneous tmle" = 25, "init est" = 39)
    
    performance.sig = t(apply(results[,varind], 2, perf,var0))
    performance.ate = t(apply(results[,ateind], 2, perf,ATE0))
    
    rownames(performance.sig) = names(varind)
    rownames(performance.ate) = names(ateind)
    
    coverage = c(cov.check(results, var0, 1),
                 cov.simul(results, c(var0, ATE0), c(7,25)),
                 cov.check(results, ATE0, 28), cov.check(results, ATE0, 34))
    
    # getting coveage using the true variance of the estimator
    dd = data.frame(psi = results[,1],l = results[,1]-1.96*sqrt(performance.sig[1,1]),
                    r = results[,1]+1.96*sqrt(performance.sig[1,1]), psi = results[,28],
                    l = results[,28]-1.96*sqrt(performance.ate[1,1]),
                    r = results[,28]+1.96*sqrt(performance.ate[1,1]))
    
    cov.sig.1step = cov.check(dd, var0, 1)
    cov.ate = cov.check(dd, ATE0, 4)
    cov = c(coverage, cov.sig.1step, cov.ate)
    names(cov) = c("TMLE Blip Variance SL2", "Simultaneous TMLE SL2", "TMLE ATE SL2",
                   "TMLE ATE LR", "TMLE Blip Var using true Var", "TMLE ATE using true Var")
    coverage = data.frame(coverage = cov)
    rownames(coverage) = names(cov)
    
    # getting superlearner results
    LL = 0
    for (i in 1:length(SL.library)) {
      if (length(SL.library[[i]]) > 1) {
        LL = LL + length(SL.library[[i]])-1} else
        {LL = LL + 1}
    }
    SL_results = data.frame(colMeans(results[,65:(64+LL)]))
    rownames(SL_results) = colnames(results)[65:(64+LL)]
    
    type = c(rep("TMLE SL2",B), rep("init. est SL2",B),
             rep("TMLE LR",B), rep("init est. LR",B))
    types = c("simul TMLE SL2","iter TMLE SL2","init. est SL2","init est. LR")
    inds = c(25, 39, 34, 40)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("red", "blue","green","orange")
    
    ateests = data.frame(ests=ests,type=type)
    ggover = ggplot(ateests,aes(ests, fill=type,color=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("ATE sampling distributions, ", case))
    ggover = ggover+geom_vline(xintercept = ATE0,color="black")+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    capt = paste0("Truth is at black vline.\n",
                  "\ninitial est using SuperLearner library 2 is excellent", 
                  "\nTMLE for ATE using IC approx for variance covers at ", 
                  100*round(coverage[3,1],3),"%",
                  "\ninit est. LR which is Logistic Reg. with main terms and interactions ",
                  "\nplug-in clearly biased but TMLE corrects it and it covers at ",
                  100*round(coverage[4,],3),"%\n")
    ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", 
                          fontface = "plain",colour = "black", size = 10, angle = 0, 
                          lineheight = 0.9))
    gg_ATEcase2bSL2 = ggover
    
    
    type = c(rep("TMLE SL2",B), rep("init est SL2",B), rep("init est LR",B))
    types = c("TMLE SL2","init est SL2","init est LR")
    inds = c(1, 37, 38)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("blue","green","orange","red")
    varests = data.frame(ests=ests,type=type)
    
    ggover2 = ggplot(varests,aes(ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("Blip Variance sampling distributions, ", case))
    ggover2 = ggover2+geom_vline(xintercept = var0,color="black")+
      theme(plot.title = element_text(size=12), 
            plot.subtitle = element_text(size=10))+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])
    cap = paste0("truth is black line\n",
                 "tmle SL2, which used Superlearner Library 2 for initial ests\n", 
                 "is skewed, has many outliers and covers at ", 100*round(coverage[1,1],3),"%\n", 
                 "and does not debias initial estimate due to overfit outcome regressions\n",
                 "which defy the donsker condition--need CV-TMLE!\n",
                 "init est LR used logistic regression with main terms and interactions\n",
                 "plug-in estimator and is a disaster, which TMLE cannot help.\n")
    ggover2=ggdraw(add_sub(ggover2,cap, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 10, angle = 0, 
                           lineheight = 0.9))
    assign(paste0("gg_BV",case), ggover2)
    assign(paste0("results_",case), results)
    assign(paste0("performance.sig_",case), performance.sig)
    assign(paste0("performance.ate_",case), performance.ate)
    assign(paste0("coverage_",case), coverage)
    assign(paste0("SL_results_",case), SL_results) 
  }
  if (case == "combo_case2b") {
    g0 = g0_linear
    Q0 = Q0_trig
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    coverage_case2b = c(cov.check(results_case2bSL1, var0, c(1)),
                        cov.check(results_case2bSL2, var0, c(1)),
                        cov.check(results_case2bCVSL2, var0, c(1)))
    
    coverage_simulcase2b = c(cov.simul(results_case2bSL1, c(var0,ATE0), c(7,25)),
                             cov.simul(results_case2bSL2, c(var0,ATE0), c(7,25)),
                             cov.simul(results_case2bCVSL2, c(var0,ATE0), c(7,25)))
    
    coverage_case2b = c(coverage_case2b, coverage_simulcase2b, 
                        rep(NA,3))[c(1,4,7,2,5,7,3,6,7)]
    
    MSE_case2b = rbind(t(apply(results_case2bSL1[,c(1,7,37)],2,perf,var0)),
                       MSE_of = t(apply(results_case2bSL2[,c(1,7,37)],2,perf,var0)),
                       MSE_cv = t(apply(results_case2bCVSL2[,c(1,7,37)],2,perf,var0)))
    
    MSE_cov = cbind(MSE_case2b, coverage_case2b)
    rownames(MSE_cov) = c("TMLE SL1","Simul. TMLE SL1","init est SL1",
                          "TMLE SL2","Simul. TMLE SL2","init est SL2",
                          "CV-TMLE SL2","Simul. CV-TMLE SL2","init est CV-SL1")
    
    # df_trueSL1 = data.frame(results_case2bSL1[,1], 
    #                      results_case2bSL1[,1]-1.96*sqrt(MSE_cov[1,1]),
    #                      results_case2bSL1[,1]+1.96*sqrt(MSE_cov[1,1]))
    # df_trueSL2 = data.frame(results_case2bSL1[,1], 
    #                         results_case2bSL1[,1]-1.96*sqrt(MSE_cov[4,1]),
    #                         results_case2bSL1[,1]+1.96*sqrt(MSE_cov[4,1]))
    # df_trueCVSL2 = data.frame(results_case2bSL1[,1], 
    #                         results_case2bSL1[,1]-1.96*sqrt(MSE_cov[7,1]),
    #                         results_case2bSL1[,1]+1.96*sqrt(MSE_cov[7,1]))
    # 
    # cov_true = unlist(lapply(list(df_trueSL1, df_trueSL2, df_trueCVSL2), 
    #        FUN = function(x) cov.check(x, var0, c(1))))
    # names(cov_true) = c("TMLE SL1", "TMLE SL2", "CV-TMLE SL2")
    
    B1 = nrow(results_case2bSL1)
    B2 = nrow(results_case2bSL2)
    B3 = nrow(results_case2bCVSL2)
    
    type = c(rep("tmle SL1",B1), rep("tmle SL2",B2),
             rep("cv-tmle SL2",B3))
    types = c("tmle SL1","tmle SL2",
              "cv-tmle SL2")
    colors = c("red","green","blue")
    ests = c(results_case2bSL1[,1], results_case2bSL2[,1], 
             results_case2bCVSL2[,1])
    plotdf = data.frame(ests=ests, type=type)
    
    ggover = ggplot(plotdf, aes(x=ests, fill = type, color = type))+geom_density(alpha=.5)
    ggover = ggover + scale_fill_manual(values=colors)
    ggover = ggover + scale_color_manual(values=colors)
    ggover = ggover + geom_vline(xintercept = var0, color= "black")+
      ggtitle("The virtue of cv-tmle, case 2b",
              subtitle = "cv-tmle maintains normality, eliminates skewing, bad outliers")+
      theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))+
      geom_vline(xintercept = mean(as.numeric(results_case2bCVSL2[,1])), color = colors[1])+
      geom_vline(xintercept = mean(as.numeric(results_case2bSL1[,1])), color = colors[2])+
      geom_vline(xintercept = mean(as.numeric(results_case2bSL2[,1])), color = colors[3])
    
    cap = paste0("truth is black line\n",
                 "tmle SL2, which used Superlearner Library 2 for initial ests\n", 
                 "is skewed, has many outliers and covers at ", 100*round(MSE_cov[4,4],3),"%\n", 
                 "has over double the MSE of the other two and is less efficient.\n",
                 "tmle SL1 uses a non-overfitting SuperLearner and covers near nominally at ", 
                 100*round(MSE_cov[1,4],3), "%\n",
                 "cv-tmle SL2 uses the overfitting SuperLearner library which defies the donsker\n",
                 "condition for initial ests, but cv-tmle does not require the donsker condition.\n",
                 "cv-tmle has lowest MSE, no skewing and covers respectably at ",
                 100*round(MSE_cov[8,4],3),"%.")
    ggover=ggdraw(add_sub(ggover,cap, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 10, angle = 0, 
                           lineheight = 0.9))
    
    gg_cvadvert = ggover
  }
  if (case == "case4_hal") {
    g0 = g0_1
    Q0 = Q0_1
    cl = makeCluster(no.cores, type = "SOCK")
    registerDoSNOW(cl)
    clusterExport(cl,cl_export)
    
    ALL_hal=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                    .errorhandling = "remove")%dopar%
                    {sim_hal(n, g0 = g0, Q0 = Q0)}
    
    results_hal = data.matrix(data.frame(do.call(rbind, ALL_hal)))
  }
  if (case == "case4_halglm"){
    g0 = g0_1
    Q0 = Q0_1
    
    SL.library = list(c("glm.mainint", "screen.Main"), 
                      c("SL.hal", "screen.Main"))
    SL.libraryG = list("SL.glm", "SL.hal")
    
    cl = makeCluster(detectCores(), type = "SOCK")
    registerDoSNOW(cl)
    clusterExport(cl,cl_export)
    
    ALL_halglm=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                       .errorhandling = "remove")%dopar%
                       {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                               SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE
                       )}
    results_halglm = data.matrix(data.frame(do.call(rbind, ALL_halglm)))
  }
  
  if (case == "case4"){
    g0 = g0_1
    Q0 = Q0_1
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    SL.library = list(c("glm.mainint", "screen.Main"), 
                      c("SL.hal", "screen.Main"))
    SL.libraryG = list("SL.glm", "SL.hal")
    
    varind = c("1step tmle LR" = 20,"1step tmle HAL" = 1, "1step tmle HAL+glm" = 5,
               "init est LR" = 42,"init est HAL" = 4, "init est HAL+glm" = 41)
    
    performance.sig = lapply(results[varind],perf,var0)
    performance.sig = t(as.data.frame(performance.sig))
    rownames(performance.sig) = names(varind)
    
    coverage_halglm = cov.check(results_halglm, var0, c(16,1))
    coverage_halglm 
    coverage_hal = cov.check(results_hal, var0, 1)
    coverage_hal
    coverage = c(coverage_halglm, coverage_hal, .434, NA, NA)
    
    MSE_cov = cbind(performance.sig, coverage)
    MSE_cov
    # getting superlearner results
    LL = 0
    for (i in 1:length(SL.library)) {
      if (length(SL.library[[i]]) > 1) {
        LL = LL + length(SL.library[[i]])-1} else
        {LL = LL + 1}
    }
    SL_results = data.frame(colMeans(results_halglm[,65:(65+LL-1)]))
    rownames(SL_results) = colnames(results)[65:(65+LL-1)]
    
    LG=0
    for (i in 1:length(SL.libraryG)) {
      if (length(SL.libraryG[[i]]) > 1) {
        LG = LG + length(SL.libraryG[[i]])-1} else
        {LG = LG + 1}
    }
    SL_resultsG = data.frame(colMeans(results_halglm[,(65+LL):(65+LL+LG-1)]))
    rownames(SL_results) = colnames(results)[(65+LL):(65+LL+LG-1)]
    
    B = nrow(results_halglm)
    B1 = nrow(results_hal)
    type = c(rep(names(varind)[1],B), rep(names(varind)[2],B1), rep(names(varind)[3],B),
             rep(names(varind)[4],B), rep(names(varind)[5],B1), rep(names(varind)[6],B))
    # types = c("1step TMLE LR","1step TMLE HAL","1step TMLE HAL+glm")
    types = names(varind)
    inds = varind
    ests = unlist(lapply(inds, FUN = function(x) results[[x]]))
    inds = inds[order(types)]
    
    colors = c("blue","green","orange","red", "purple", "yellow")
    varests = data.frame(ests=ests,type=type)
    
    ggover2 = ggplot(varests,aes(ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("Blip Variance sampling distributions, ", case))
    ggover2 = ggover2+geom_vline(xintercept = var0,color="black")+
      theme(plot.title = element_text(size=12), 
            plot.subtitle = element_text(size=10))+
      geom_vline(xintercept=mean(results[[inds[1]]]),color = colors[1])+
      geom_vline(xintercept=mean(results[[inds[2]]]),color = colors[2])+
      geom_vline(xintercept=mean(results[[inds[3]]]),color = colors[3])+
      geom_vline(xintercept=mean(results[[inds[4]]]),color = colors[4])+
      geom_vline(xintercept=mean(results[[inds[5]]]),color = colors[5])+
      geom_vline(xintercept=mean(results[[inds[6]]]),color = colors[6])
    cap = paste0("truth at black line.\n",
                 "tmle LR uses glm with interactions for outcome model and glm for\n",
                 "treatment mechanism initial estimates. tmle LR CI's cover at", 
                 round(MSE_cov[3,4] ,1),"\n",
                 "hal tmle uses highly adaptive lasso for initial estimates of both\n",
                 "outcome and treatment mech initial estimates and these cover\n",
                 "at", round(MSE_cov[3,4] ,2),
                 "\ntmle hal+glm SL uses a SuperLearner with hal and glm for\n",
                 "outcome and treatment mechanism initial estimates and cover at ",
                 round(MSE_cov[3,4] ,3),"\n",
                 "All coverage here is using infuence curve approx for inference.")
    ggover2=ggdraw(add_sub(ggover2,cap, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 10, angle = 0, 
                           lineheight = 0.9))
    assign(paste0("results_",case), results)
    assign(paste0("gg_BV",case), ggover2)
    assign(paste0("MSE_cov_",case), MSE_cov)
    assign(paste0("SL_results_",case), SL_results) 
    assign(paste0("SL_resultsG_",case), SL_resultsG) 
  }
  
  if (case == "case3") {
    
    g0 = g0_1
    Q0 = Q0_2
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    SL.library = SL.library1
    SL.libraryG = SL.libraryG
    
    if (!resultsGotten) {
      cl = makeCluster(no.cores, type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,cl_export)
      
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                  .errorhandling = "remove")%dopar%
                  {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                          SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE
                  )}
      results = data.matrix(data.frame(do.call(rbind, ALL)))
    }
    
    B = nrow(results)
    
    varind = c("1step tmle" = 1,"simultaneous tmle" = 7, "init est" = 37)
    ateind = c("1step tmle" = 28,"simultaneous tmle" = 25, "init est" = 39)
    
    performance.sig = t(apply(results[,varind], 2, perf,var0))
    performance.ate = t(apply(results[,ateind], 2, perf,ATE0))
    
    rownames(performance.sig) = names(varind)
    rownames(performance.ate) = names(ateind)
    
    coverage = c(cov.check(results, var0, 1),
                 cov.simul(results, c(var0, ATE0), c(7,25)),
                 cov.check(results, ATE0, 28), cov.check(results, ATE0, 34),
                 cov.check(results, var0, 16))
    
    # getting coveage using the true variance of the estimator
    dd = data.frame(psi = results[,1],l = results[,1]-1.96*sqrt(performance.sig[1,1]),
                    r = results[,1]+1.96*sqrt(performance.sig[1,1]), psi = results[,28],
                    l = results[,28]-1.96*sqrt(performance.ate[1,1]),
                    r = results[,28]+1.96*sqrt(performance.ate[1,1]))
    
    cov.sig.1step = cov.check(dd, var0, 1)
    cov.ate = cov.check(dd, ATE0, 4)
    cov = c(coverage, cov.sig.1step, cov.ate)
    names(cov) = c("TMLE Blip Variance SL1", "Simultaneous TMLE SL1", "TMLE ATE SL1",
                   "TMLE ATE LR", "TMLE Blip Variance LR", 
                   "TMLE Blip Var using true Var", 
                   "TMLE ATE using true Var")
    coverage = data.frame(coverage = cov)
    rownames(coverage) = names(cov)
    
    # getting superlearner results
    LL = 0
    for (i in 1:length(SL.library)) {
      if (length(SL.library[[i]]) > 1) {
        LL = LL + length(SL.library[[i]])-1} else
        {LL = LL + 1}
    }
    SL_results = data.frame(colMeans(results[,65:(64+LL)]))
    rownames(SL_results) = colnames(results)[65:(64+LL)]
    
    type = c(rep("TMLE SL1",B), rep("init. est SL1",B),
             rep("TMLE LR",B), rep("init est. LR",B))
    types = c("simul TMLE SL1","iter TMLE SL1","init. est SL1","init est. LR")
    inds = c(25, 39, 34, 40)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("red", "blue","green","orange")
    
    ateests = data.frame(ests=ests,type=type)
    ggover = ggplot(ateests,aes(ests, fill=type,color=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("ATE sampling distributions, ", case))
    ggover = ggover+geom_vline(xintercept = ATE0,color="black")+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    capt = paste0("Truth is at black vline.\n",
                  "\ninitial est using SuperLearner library 1 is excellent", 
                  "\nTMLE for ATE using IC approx for variance covers at ", 
                  100*round(coverage[3,1],3),"%", 
                  "\ninit est. LR which is Logistic Reg. with main terms and interactions ",
                  "\nplug-in unbiased and TMLE keeps it as is, covering at ",
                  100*round(coverage[4,],3),"%\n")
    ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", 
                          fontface = "plain",colour = "black", size = 10, angle = 0, 
                          lineheight = 0.9))
    gg_ATEcase3 = ggover
    
    
    type = c(rep("TMLE SL1",B), rep("init est SL1",B), rep("init est LR",B),
             rep("TMLE LR",B))
    types = c("TMLE SL1","init est SL1","init est LR","TMLE LR")
    inds = c(1, 37, 38, 16)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("blue","green","orange","red")
    varests = data.frame(ests=ests,type=type)
    
    ggover2 = ggplot(varests,aes(ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("Blip Variance sampling distributions, ", case),
              subtitle = "Recovering both treatment mech and outcome")
    ggover2 = ggover2+geom_vline(xintercept = var0,color="black")+
      theme(plot.title = element_text(size=12), 
            plot.subtitle = element_text(size=10))+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    cap = paste0("truth is black line\n",
                 "tmle SL1, which used Superlearner Library 1 for initial ests\n", 
                 "attains very good coverage at ", 100*round(coverage[1,1],3),"%\n", 
                 "and debiases initial SuperLearner lib 1 estimate.\n", 
                 "Simultaneous TMLE covers both ATE and blip var at ", 
                 100*round(coverage[2,1],3),
                 "%\n",
                 "init est LR, which used logistic regression with main terms and\n",
                 "interactions, and corresponding TMLE and more biased.\n",
                 "This underscores the importance of machine learning in making\n",
                 "the initial estimates for both propensity score and outcome here.")
    ggover2=ggdraw(add_sub(ggover2,cap, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 10, angle = 0, 
                           lineheight = 0.9))
    
    assign(paste0("gg_BV",case), ggover2)
    assign(paste0("results_",case), results)
    assign(paste0("performance.sig_",case), performance.sig)
    assign(paste0("performance.ate_",case), performance.ate)
    assign(paste0("coverage_",case), coverage)
    assign(paste0("SL_results_",case), SL_results) 
  }
  
  if (case == "noise"|case =="noise_neg"){
    g0 = g0_linear
    Q0 = Q0_noise
    testdata=gendata_noise(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    rate = 1/3
    if (case == "noise") {
      biasQ = function(A,W1,W2,W3,W4,n,rate)  {
        n^-rate*1.5*(-.2+1.5*A+.2*W1+1*W2-A*W3+ 1*W4)
      }
      # hist(with(truth,biasQ(A,W1,W2,W3,W4,n=n,rate=rate)))
      sdQ = function(A,W1,W2,W3,W4,n,rate) {
        (n^-rate)*.8*(abs(3.5+.5*W1+.15*W2+.33*W3*W4-W4))
      }
      N=20000} else {
        biasQ = function(A,W1,W2,W3,W4,n,rate) { 
          -n^-rate*.8*(+.2+1.5*A+.2*W1+1*W2+5*A*W3^2+ 1*W4)
        }
        
        sdQ = function(A,W1,W2,W3,W4,n,rate) {
          (n^-rate)*.8*(abs(3.5+.5*W1+.15*W2+.33*W3*W4-W4))
        }
        N=50000}
    cl = makeCluster(no.cores, type = "SOCK")
    registerDoSNOW(cl)
    L = list()
    i=1
    sizes = seq(250,N,250)
    for (n in sizes){
      rate = 1/3
      print(i)
      time = proc.time()
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm", "Simulations"))%dopar%
      {simBlipvar(n = n, rate = rate, g0 = g0, Q0 = Q0, biasQ, sdQ)}
      
      L[[i]] = ALL
      print(proc.time()-time)
      i=i+1
    }
    res_noise = vapply(1:length(sizes), FUN = function(x) {
      res = getRes(L[[x]],B, ATE0=ATE0, var0=var0)
      bias_init = res[[2]][2,2]
      bias_tmle = res[[2]][1,2]
      mse_init = res[[2]][2,3]
      mse_tmle = res[[2]][1,3]
      coverage = res[[3]][2]
      return(c(bias_init, bias_tmle, mse_init, mse_tmle, coverage))
    }, FUN.VALUE = c(1,1,1,1,1))
    res_noise = t(res_noise)
    plotdf_biasmse = data.frame(n = rep(sizes,2), 
                                bias = c(res_noise[,1], res_noise[,2]), 
                                mse = c(res_noise[,3], res_noise[,4]),
                                type = c(rep("init", length(sizes)), 
                                         rep("tmle", length(sizes))))
    
    capt=paste0("sample size n=250, blip bias has L2 norm = O_p(n^{-1/3})",
                "\no_p(n^{-.25}) as required")
    p250 = ggdraw(add_sub(getRes(L[[1]],1000, ATE0=ATE0, var0=var0)$ggover2,capt, 
                          x= .05, y = 0.5, hjust = 0, vjust = 0.5, 
                          vpadding = grid::unit(1, "lines"), 
                          fontfamily = "", fontface = "plain",
                          colour = "black", size = 9, angle = 0, lineheight = 0.9))
    capt=paste0("sample size n=1000, blip bias has L2 norm = O_p(n^{-1/3})",
                "\no_p(n^{-.25}) as required")
    p1000 = ggdraw(add_sub(getRes(L[[4]],1000, ATE0=ATE0, var0=var0)$ggover2,capt, 
                           x= .05, y = 0.5, hjust = 0, vjust = 0.5, 
                           vpadding = grid::unit(1, "lines"), 
                           fontfamily = "", fontface = "plain",
                           colour = "black", size = 9, angle = 0, lineheight = 0.9))
    capt=paste0("sample size n=5000, blip bias has L2 norm = O_p(n^{-1/3})",
                "\no_p(n^{-.25}) as required")
    p5000 = ggdraw(add_sub(getRes(L[[20]],1000, ATE0=ATE0, var0=var0)$ggover2,capt, 
                           x= .05, y = 0.5, hjust = 0, vjust = 0.5, 
                           vpadding = grid::unit(1, "lines"), 
                           fontfamily = "", fontface = "plain",
                           colour = "black", size = 9, angle = 0, lineheight = 0.9))
    capt=paste0("sample size n=10000, blip bias has L2 norm = O_p(n^{-1/3})",
                "\no_p(n^{-.25}) as required")
    p10000 = ggdraw(add_sub(getRes(L[[40]],1000, ATE0=ATE0, var0=var0)$ggover2,capt, 
                            x= .05, y = 0.5, hjust = 0, vjust = 0.5, 
                            vpadding = grid::unit(1, "lines"), 
                            fontfamily = "", fontface = "plain",
                            colour = "black", size = 9, angle = 0, lineheight = 0.9))
    
    ml=marrangeGrob(list(p250,p1000,p5000,p10000),ncol=2,nrow=2, 
                    widths = c(3.5,3.5),
                    heights = c(1,1))
    
    if (case=="noise"){
      truth = gendata_noise(1e6, g0, Q0)
      
      grobs = lapply(c(1,4,10,20), FUN = function(x) {
        coverage = res_noise[x,5]
        noise_analysis(sizes[x],1/3,truth,coverage=coverage,
                       Q0, biasQ, sdQ)[[5]]
      })
      
      ml1=marrangeGrob(grobs,ncol=2,nrow=2, widths = c(3.5,3.5),
                       heights=c(1,1))
    } else {
      truth = gendata_noise(1e6, g0, Q0)
      grobs = lapply(c(1,4,10,160), FUN = function(x) {
        coverage = res_noise[x,5]
        noise_analysis(sizes[x],1/3,truth,coverage=coverage,
                       Q0, biasQ, sdQ)[[5]]
      })
      
      ml1=marrangeGrob(grobs,ncol=2,nrow=2, widths = c(3.5,3.5),
                       heights=c(1,1))
    }
    
    
  }
  if (case == "wells") {
    
    a = seq(0,15,.5)
    b = seq(0,15,.5)
    truevars = c()
    trueates = c()
    gg_wells = list()
    oc_list = list()
    for (i in 1:31){
      Q0 = function (A, W1, W2, W3, W4) 
      {
        plogis(.14*(2* A + W1 + a[i]*A * W1 - b[i]*A*W2 + W2 - W3+ W4))
      }
      testdata=gendata(1000000, g0=g0_1, Q0 = Q0)
      blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
      truevars = c(truevars, var(blip_true))
      trueates = c(trueates, mean(blip_true))
      oc_list = append(oc_list,Q0)
    }
    
    for (nn in c(250,500,1000)){
      cl = makeCluster(detectCores(), type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,c("a","b","oc_list"))
      for (i in seq(1,31,2)){
        print(i)
        B = 1000
        n = nn
        Q0 = oc_list[[i]]
        SL.library =
          ALL = foreach(rep=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations"),
                        .export = c("a","b","i"))%dopar%
                        {sim_lr(n, g0 = g0_linear, Q0 = Q0, 
                                formQ = formula("Y~A*(W1 + W2) +W3 +W4"),
                                formG = formula("A~."))}
        
        results = do.call(rbind,ALL)
        cnames = colnames(results)
        results = apply(results,2,as.numeric)
        colnames(results) = cnames
        
        var0 = truevars[i]
        ATE0 = trueates[i]
        cov = vapply(c(1,4),FUN = function(x){
          covs = results[,x+1]<=var0&results[,x+2]>=var0
          mean(covs)
        }, FUN.VALUE = 1)
        cov_simulsig = results[,8] <= var0 & results[,9] >= var0
        cov_simulATE = results[,11] <= ATE0 & results[,12] >= ATE0
        cov_simul = mean(cov_simulsig * cov_simulATE)
        cov_ate = mean(results[,14] <= ATE0 & results[,15] >= ATE0)
        cov = c(cov, cov_simul, cov_ate)
        coverage = rbind(coverage,cov)
        
        MSE = apply(results[,c(1,4,7,16)], 2, perf, var0)
        MSE = t(MSE)
        
        bias1 = MSE[,2]
        bias = rbind(bias, bias1)
        
        type = c(rep("1step tmle",1000), 
                 rep("initial",1000))
        types = c("1step tmle", "initial")
        
        types = types[c(1,16)]
        inds = c(1,16)[order(types)]
        types = types[order(types)]
        colors = c("red", "blue","green","orange")
        rbind(colors, types)
        
        ests = c(results[,1],results[,16])
        plotdf = data.frame(ests=ests, type=type)
        
        gg_hal = ggplot(plotdf, aes(x=ests, fill = type, color = type))+geom_density(alpha=.5)
        gg_hal = gg_hal + scale_fill_manual(values=colors)
        gg_hal = gg_hal + scale_color_manual(values=colors)
        gg_hal = gg_hal + geom_vline(xintercept = var0, color= "black")+
          theme(plot.title = element_text(size=12))+
          ggtitle(paste0("tmle sampling dists \nwell-spec models, n=",nn))+
          geom_vline(xintercept = mean(as.numeric(results[,inds[1]])), color = colors[1])+
          geom_vline(xintercept = mean(as.numeric(results[,inds[2]])), color = colors[2])
        caption = paste0("truth at black line, true blip variance = ",round(var0,6),
                         "\ntmle bias is ",round(MSE[2,2],6))
        gg_hal=ggdraw(add_sub(gg_hal,caption,x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                              vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                              colour = "black", size = 10, angle = 0, lineheight = 0.9))
        gg_hal
        gg_wells[[count]] = gg_hal
        count = count+1
        
      }
    }
    data = gendata(1e6, g0_linear, Q0_trig1)
    pscores = with(data, g0_linear(W1,W2,W3,W4))
    df = data.frame(pscores=pscores, type = rep("p", 1e6))
    gg_pscoresWell = ggplot(df, aes(x = pscores, fill = type)) + geom_density()+
      scale_fill_manual(values=c("blue"))+
      scale_color_manual(values=c("blue"))+
      theme(legend.position="none")
    gg_pscoresWell = ggdraw(add_sub(gg_pscoresWell,"no positivity issues here", 
                                    x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                                    vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                                    colour = "black", size = 10, angle = 0, lineheight = 0.9))
  }
  if (case%in%list("CI_LRcase2a", "CI_LRcase2b", "CI_LRcase3", "CI_LRcase4")) {
    if(case=="CI_LRcase2a"){
      g0 = g0_linear
      Q0 = Q0_trig1
    }
    if(case=="CI_LRcase2b"){
      g0 = g0_linear
      Q0 = Q0_trig
    }
    if(case=="CI_LRcase3"){
      g0 = g0_1
      Q0 = Q0_2
    }
    if(case=="CI_LRcase4"){
      g0 = g0_1
      Q0 = Q0_1
    }
    
    # get the truth
    testdata=gendata(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    
    cl = makeCluster(detectCores(), type = "SOCK")
    registerDoSNOW(cl)
    clusterExport(cl,cl_export)
    
    ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                .errorhandling = "remove")%dopar%
                {LR.BVinference(n, g0 = g0, Q0 = Q0)}
    
    resultsLR = data.matrix(data.frame(do.call(rbind, ALL)))
    coverageLR = c(cov.check(resultsLR, var0, 1),
                   cov.check(resultsLR, ATE0, 7),
                   cov.simul(resultsLR, c(ATE0,var0), c(10,4)))
    
    names(coverageLR) = c("blip variance", "ATE", "simultaneous")
    
    assign(paste0("coverageLR_",case),coverageLR)
  }
  
  if (case == "case2b_2G") {
    B = nrow(results)
    
    varind = c("1step tmle" = 1,"simultaneous tmle" = 7, "init est" = 37)
    ateind = c("1step tmle" = 28,"simultaneous tmle" = 25, "init est" = 39)
    
    performance.sig = t(apply(results[,varind], 2, perf,var0))
    performance.ate = t(apply(results[,ateind], 2, perf,ATE0))
    
    rownames(performance.sig) = names(varind)
    rownames(performance.ate) = names(ateind)
    
    coverage = c(cov.check(results, var0, 1),
                 cov.simul(results, c(var0, ATE0), c(7,25)),
                 cov.check(results, ATE0, 28))
    
    # getting coveage using the true variance of the estimator
    dd = data.frame(psi = results[,1],l = results[,1]-1.96*sqrt(performance.sig[1,1]),
                    r = results[,1]+1.96*sqrt(performance.sig[1,1]), psi = results[,28],
                    l = results[,28]-1.96*sqrt(performance.ate[1,1]),
                    r = results[,28]+1.96*sqrt(performance.ate[1,1]))
    
    cov.sig.1step = cov.check(dd, var0, 1)
    cov.ate = cov.check(dd, ATE0, 4)
    cov = c(coverage, cov.sig.1step, cov.ate)
    names(cov) = c("TMLE Blip Variance", "Simultaneous TMLE", "TMLE ATE",
                   "TMLE Blip Var using true Var", "TMLE ATE using true Var")
    coverage = data.frame(coverage = cov)
    rownames(coverage) = names(cov)
    
    # getting superlearner results
    LL = 0
    for (i in 1:length(SL.library)) {
      if (length(SL.library[[i]]) > 1) {
        LL = LL + length(SL.library[[i]])-1} else
        {LL = LL + 1}
    }
    SL_results = data.frame(colMeans(results[,65:(64+LL)]))
    rownames(SL_results) = colnames(results)[65:(64+LL)]
    
    type = c(rep("iter TMLE SL1",B), rep("init. est SL1",B),
             rep("TMLE LR",B), rep("init est. LR",B))
    types = c("simul TMLE SL1","iter TMLE SL1","init. est SL1","init est. LR")
    inds = c(25, 39, 34, 40)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("red", "blue","green","orange")
    
    ateests = data.frame(ests=ests,type=type)
    ggover = ggplot(ateests,aes(ests, fill=type,color=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("ATE sampling distributions, ", case))
    ggover = ggover+geom_vline(xintercept = ATE0,color="black")+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    ggover_ATEcase2b2G = ggover
    
    
    type = c(rep("TMLE SL1",B), rep("init est SL1",B), rep("init est LR",B))
    types = c("TMLE SL1","init est SL1","init est LR")
    inds = c(1, 37, 38)
    ests = unlist(lapply(inds, FUN = function(x) results[,x]))
    inds = inds[order(types)]
    
    colors = c("blue","green","orange","red")
    
    varests = data.frame(ests=ests,type=type)
    
    ggover2 = ggplot(varests,aes(ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("Blip Variance sampling distributions, ", case))
    ggover2 = ggover2+geom_vline(xintercept = var0,color="black")+
      theme(plot.title = element_text(size=12), 
            plot.subtitle = element_text(size=10))+
      geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])
    
    ggover_BVcase2b2G = ggover2
  }
  
  if (case == "example"){
    # get rid of a few with missing data
    nas = apply(wcgs, 1, FUN = function(x) any(is.na(x)))
    bads  = which(nas)
    goods = which(!nas)
    
    #select covariates as per the paper
    data = wcgs[goods, c(2:7,9:11)]
    # assign typical names to treatment and outcome
    colnames(data)[8:9] = c("A","Y")
    
    #####
    #####
    X = data
    X$Y = NULL
    W = X
    W$A = NULL
    A = X$A
    
    mainform = paste0(paste(colnames(data)[1:6],"+",collapse=""),colnames(data)[7])
    
    squares = paste0(paste0("I(",colnames(data)[1:7]),"^2)")
    
    squares = paste0(paste(squares[1:6],"+",collapse=""),squares[7])
    
    mainsq = paste0(mainform,"+",squares)
    
    mainsq.int = paste0("Y~A*(",mainsq,")")
    mainsq.int = formula(mainsq.int) 
    
    X = model.matrix(mainsq.int,data)
    X = as.data.frame(X[,-1])
    colnames(X)[2:ncol(X)] = paste0("X",2:ncol(X))
    head(X)
    X1 = X0 = X
    X0$A = 0
    X1$A = 1
    Y = data$Y
    newdata = rbind(X,X1,X0)
    
    SL.library = list(c("nnetMain","screen.Main"),c("nnetMain1","screen.Main"),
                      c("earthFull", "screen6", "screen12","All"),
                      c("SL.earth", "screen.Main"),c("xgboostFull","screen12","All"),
                      c("xgboostMain", "screen.Main"),
                      c("gamFull", "screen6", "screen12","All"), 
                      c("SL.gam", "screen.Main"),"SL.rpartPrune",
                      c("rangerMain", "screen.Main"), c("rpartMain", "screen.Main"),
                      "SL.stepAIC", c("SL.glm","screen6", "screen12","screen.Main","All"),
                      c("SL.hal", "screen.Main"), "SL.mean")
    
    SL.libraryG = list("nnetMainG","nnetMainG1","SL.earth","SL.rpartPrune",
                       "SL.gam", "rpartMain", "SL.step.interaction", "SL.glm", 
                       "SL.hal", "SL.mean","xgboostG")
    
    # SL.library = SL.libraryG = c("SL.mean","SL.glm")
  }
}


