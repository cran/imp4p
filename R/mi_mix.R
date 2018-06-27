
####################################
#
#Function to perform multiple imputation
#
####################################

mi.mix=function (tab, tab.imp, prob.MCAR, conditions, repbio = NULL,
                 reptech = NULL, nb.iter = 3, nknn = 15, weight = 1, selec = "all",
                 siz = 500, ind.comp = 1, methodi = "mle", q = 0.95, progress.bar = TRUE, details= FALSE)
{
  if (is.null(repbio)) {
    repbio = as.factor(1:length(conditions))
  }
  if (is.null(reptech)) {
    reptech = as.factor(1:length(conditions))
  }
  if (selec>=nrow(tab)){
    selec=nrow(tab)-1;
  }
  if (selec == "all") {
    if (siz > nrow(tab)) {
      siz = nrow(tab) - 1;
    }
  }
  else {
    if (siz > selec) {
      siz = selec - 1;
    }
  }
  data_imp = array(NA, dim = c(nrow(tab), ncol(tab), nb.iter))
  l.NA = matrix(0, nrow(tab), ncol(tab))
  if (progress.bar == TRUE) {
    cat(paste("\n Iterations: \n"))
  }
  iter = 1
  while ((iter <= nb.iter)) {
    if (progress.bar == TRUE) {
      cat(paste("\n", iter, "/", nb.iter, " - "))
    }
    tab.mvs = as.matrix(tab)
    l.MCAR = matrix(mapply(prob.MCAR, FUN = rbinom, size = 1,
                           n = 1), nrow(prob.MCAR), ncol(prob.MCAR))
    l.NA[which(is.na(tab.mvs))] = 1
    l.MNAR = l.NA - l.MCAR
    l.MCAR = l.NA * l.MCAR
    l.MNAR = l.NA * l.MNAR
    tab.mvs[which(l.MCAR == 1)] = tab.imp[which(l.MCAR ==
                                                  1)]
    data_imp[, , iter] = tab.mvs
    tab.mvs.imp = impute.igcda(tab = tab.mvs, tab.imp = tab.imp,
                               conditions = conditions, q = q)
    if (progress.bar == TRUE) {
      cat(paste("Imputation MNAR OK - \n"))
    }
    rna = apply(l.NA, 1, sum)
    rna = which(rna > 0)
    if (progress.bar == TRUE) {
      cat(paste("Imputation MCAR in progress - \n"))
    }
    for (i in 1:length(rna)) {
      tab.mod = tab.mvs.imp
      tab.mod.imp = tab.mod
      tab.mod[rna[i], ] = tab.mvs[rna[i], ]
      if (methodi == "mle") {
        nb_cond = length(levels(conditions))
        lab = (1:nrow(tab.mod))[-rna[i]]
        sel = selec
        if (selec == "all") {
          sel = nrow(tab.mod) - 1
        }
        list.select = sample(lab, size =  max(sel, min(siz,nrow(tab.mod) - 1)), replace = FALSE)
        list.select = c(list.select, rna[i])
        tab.mod.imp2 = NULL
        nb_rep = rep(0, nb_cond)
        k = 1
        for (it in 1:nb_cond) {
          nb_rep[it] = sum((conditions == levels(conditions)[it]))
          tab.mod.imp2 = cbind(tab.mod.imp2, impute.wrapper.MLE(tab.mod[list.select,
                                                                        (k:(k + nb_rep[it] - 1))]))
          k = k + nb_rep[it]
        }
        tab.mod.imp[list.select, ] = tab.mod.imp2
      }
      else {
        lab = (1:nrow(tab.mod))[-rna[i]]
        sel = selec
        if (selec == "all") {
          sel = nrow(tab.mod) - 1
        }
        #        print(lab)
        #       print(max(sel, min(siz,nrow(tab.mod) - 1)))
        list.select = sample(lab, size = max(sel, min(siz,nrow(tab.mod) - 1)), replace = FALSE)
        list.select = c(list.select, rna[i])
        tab.mod.imp[list.select, ] = impute.slsa(tab.mod[list.select,
                                                         ], conditions = conditions, repbio = repbio,
                                                 reptech = reptech, nknn = nknn, weight = weight,
                                                 selec = selec, progress.bar = FALSE, ind.comp = ind.comp)
      }
      data_imp[rna[i], , iter] = tab.mod.imp[list.select,
                                             ][length(list.select), ]
      if (progress.bar == TRUE) {
        if (i%/%floor(0.01 * length(rna)) == i/floor(0.01 *
                                                     length(rna))) {
          cat(c((i%/%(0.01 * length(rna)))), "% - ")
        }
      }
    }
    iter = iter + 1
  }
  if (details == TRUE){
      data_fin = list(
      imputed.matrix=apply(data_imp, 1:2, mean),
      sd.imputed.matrix=apply(data_imp, 1:2, sd),
      all.imputed.matrices=data_imp);
  }else{
      data_fin = apply(data_imp, 1:2, mean)
  }
  return(data_fin)
}
