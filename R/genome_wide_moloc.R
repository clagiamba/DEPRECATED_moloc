#' Genome-wide co-localization with three traits
#'
#' @title coloc.eqtl.biom
#' @param listData, a list of data frames in this order: gwas, eqtl, mqtl
#' The data frames need to have columns: "SNP" or "CHR", "POS"; "BETA", "SE"; 
#'                                        If eqtl/mqtl: "ProbeID"
#'                                        Optionally, "A1", "A2" if want to match alleles
#' @param bed, a data frame with "CHR", "START", "STOP", and 
#'             "ProbeID" matching either eqtl or mqtl
#' @param outfolder, path of folder where the results will be written to
#' @param prefix, character, name to give to the output file
#' @param save.SNP.info, logical, should info of each SNP be saved? 
#'        If this option is true, beware that it will takes up a lot of time and space.
#' @param cores, number. See foreach package.
#' @param have_alleles, logical. If TRUE, matches alleles to be the same as gwas data.
#' @param bychrpos, logical. If TRUE, uses "CHR", "POS" columns 
#'        as SNP id to match across data frames;
#'        If TRUE, make sure not to have another column called "SNP".
#' @param prior_var, numeric vector specifying any number of values;
#'        These will be used to average the ABF computation across different prior variances.
#' @param priors, numeric vectors with three numbers; 
#'        First number is the prior of associaion of a SNP with any of the traits;
#'        Second number is the prior of a SNP being associated with two traits;
#'        Third number is the prior of a SNP being associated with three traits.
#' @param min_nsnps, number. The minimum number of matching SNPs across 
#'        the three data frames to consider. 
#' @return results in outfolder: text file with colocalization results for each locus; 
#'         text file with removed snps; 
#'         if save.SNP.info = TRUE, a folder called "coloc.output.perSNP"
#'         will save a file per locus with SNP information such as ABF (one row per matching SNP)
#' @examples
#' listData = list(gwas=biom.df, eqtl=eqtl.df, mqtl=methyl.df)
#' 
#' @export
#' @author Claudia Giambartolomei
coloc.eqtl.biom.mqtl <- function(listData, bed, outfolder = "test", prefix = "pref", save.SNP.info=FALSE, cores=20, have_alleles=TRUE, bychrpos=TRUE, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06, 1e-07), min_nsnps = 50){

  ########
  merge_results  <- function(a, b) {
    if(is.null(a) & is.null(b)) {
        return(NULL)
    } else if (is.null(a)){
        return(b)
    } else if (is.null(b)){
        return(a)
    } else {
        return(rbind(a,b))
    }
  }
  ########
  # You need the suggested package for this function    
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Pkg foreach needed for this function to work. Please install it.",
      call. = FALSE)
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Pkg doParallel needed for this function to work. Please install it.",
      call. = FALSE)
  }

  ########## 
  outfname = paste(outfolder, prefix, '_summary.tab', sep='')
  if (!file.exists(outfolder)) dir.create(outfolder)
  
  if (bychrpos) {
    chrpos = lapply(listData, function(x) x$chrpos = paste(x$CHR, x$POS, sep=":"))
    listData = Map(cbind, listData, SNP = chrpos)
    listData <- lapply(listData, function(df) {df[c("SNP")] <- lapply(df[c("SNP")], as.character); df})
  }
  
  for (i in 1:length(listData)) {
    if (length(listData[[i]]$BETA[listData[[i]]$BETA<0])>0) log=TRUE  else log=FALSE # if there are negative value, then it is a logOR?
    if (!log) message("Dataset ",  i, " should be a case-control: CHECK!!")
    # before taking the log must remove the SNPs with beta = 0
    if (!log) (listData[[i]] = subset(listData[[i]], BETA != 0))
    if (!log) (listData[[i]]$BETA = log(listData[[i]]$BETA))
  }
                         
  if (!all(c("CHR", "START", "STOP") %in% names(bed))) stop("Bed file is missing info")
  bed$CHR <- gsub("chr", "", bed$CHR)

  qtl_index <- which(unlist(lapply(listData, function(x) any(bed$ProbeID %in% x$ProbeID))))
  qtl_defined_regions <- ("ProbeID" %in% names(bed)) && (length(qtl_index)>0)
  if (qtl_defined_regions) {
    message("Regions defined around the QTL dataset number ", qtl_index[1])
    if (length(qtl_index)>1) stop("Please specify unique ProbeIDs for the bed and matching expression dataset - different from other datasets")
    eqtl.df <- listData[[qtl_index]]
    bed = bed[bed$ProbeID %in% unique(eqtl.df$ProbeID),]
    eqtl.dfByProbe = split(seq(nrow(eqtl.df)), eqtl.df$ProbeID)
  }
  # Apart from main QTL, find others so can find the loops to go through
  qtl_no_main = listData
  qtl_no_main[[qtl_index]] <- NULL
  qtl_datasets = qtl_no_main[which(unlist(lapply(qtl_no_main, function(x) length(x$ProbeID)>0)))]

  ##################################################### now start the loop
  # Now go over all regions that overlap between eQTL table and input.data
  message("Running in parallel")
  registerDoParallel(cores=cores)

  message("Looping through ", nrow(bed), " regions")

  res.all <- data.frame()
  removed_snp_list = data.frame()

  res.all  <-  foreach(i=1:nrow(bed), .combine=merge_results) %dopar% {
      data <- listData
      if (qtl_defined_regions) {
          ProbeID = as.character(bed$ProbeID[i]) ##the character bit is important for probe names that are numbers
          message("Looping through ProbeID: ", ProbeID)
          region.eqtl = eqtl.df[eqtl.dfByProbe[[as.character(ProbeID)]],]
          data[[qtl_index]] <- region.eqtl
       } 
       
       pos.start = bed$START[i]
       pos.end = bed$STOP[i]
       chrom = bed$CHR[i]
       
       listRegion = list()
       for (i in 1:length(data)) {
         region = data[[i]]
         matches <- which(region$CHR==chrom & region$POS > pos.start & region$POS < pos.end )
         x <- region[matches, ]
         listRegion[[length(listRegion)+1]] <- x
       }

       matches_in_dfs = all(unlist(lapply(listRegion, function(x) nrow(x)>0)))
       if (!matches_in_dfs) next("Empty matches")
       if (matches_in_dfs) {

       region.methyl.all = listRegion[[3]]
       if ("ProbeID" %in% names(region.methyl.all)) {
            region.methyl.all$ProbeID = as.character(region.methyl.all$ProbeID)
            methyl.dfByProbe = split(seq(nrow(region.methyl.all)), f=region.methyl.all$ProbeID)
          } else {
            methyl.dfByProbe = list(region.methyl.all)
          }

       res.out = data.frame()
       for (j in 1:length(methyl.dfByProbe)){
           data2 <- listRegion
           if ("ProbeID" %in% names(region.methyl.all)) {
             ProbeIDmet = as.character(names(methyl.dfByProbe)[j])
             message("Looping through ProbeIDmet: ", ProbeIDmet)
             region.methyl = region.methyl.all[methyl.dfByProbe[[as.character(ProbeIDmet)]],]
             data2[[3]] <- region.methyl
           }
           
           try = list()
           for (i in 1:length(data2)) {
              region = data2[[i]]
              matches <- which(region$CHR==chrom & region$POS > pos.start & region$POS < pos.end )
              x <- region[matches, ]
              try[[length(try)+1]] <- x
            }
           listRegion2 = try

           # Out of each data frame:
           # 1. Remove duplicated IDs
           # 2. keep only common SNPs in all data
           # if (all(c("A1","A2") %in% names(data))) {
           # duplSNPs <- lapply(listRegion2, function(x) remove_dupl(data=x, snpcol="SNP"))
           # try = lapply(duplSNPs, function(x) {
           try = list()
           for (d in 1:length(listRegion2)) {
                  x = remove_dupl(data=listRegion2[[d]], snpcol="SNP")
                  if (is.na(x[[2]]$Marker_removed[1])) {
                     removed_snp_list = rbind(removed_snp_list, x[[2]])
                  }
                  try[[length(try)+1]] <- x[[1]]
            }
           listRegion2 = try

           listRegion2 <- lapply(listRegion2, function(x) x[x$SNP %in% Reduce(intersect, Map("[[", listRegion2, "SNP")), ])

           if (nrow(listRegion2[[1]])>0) {

           listRegion2 <- lapply(listRegion2, function(df){
                                 df[order(df$SNP),]
                                 })
           # Check that the alleles match
           # consider A/G T/C, match by strand (complement) 
           if (have_alleles){
               # listRegion2 <- lapply(listRegion2, function(x) sapply(x[,"A1", "A2"], as.character))
               listRegion2 <- lapply(listRegion2, function(df) {df[c("A1", "A2")] <- lapply(df[c("A1", "A2")], as.character); df})
               listRegion2 <- lapply(listRegion2, function(df) {df[c("A1", "A2")] <- lapply(df[c("A1", "A2")], toupper); df})
               listRegion2 <- lapply(listRegion2, change_indels)
               match_correct2 = (listRegion2[[1]]["A1"] ==listRegion2[[2]]["A1"]) & (listRegion2[[1]]["A2"]== listRegion2[[2]]["A2"])
               match_flip2 = (listRegion2[[1]]["A1"] == listRegion2[[2]]["A2"]) & (listRegion2[[1]]["A2"] == listRegion2[[2]]["A1"])
           if (any(which(match_flip2)>0)) {
               listRegion2[[2]][which(match_flip2), "A1"]=listRegion2[[1]][which(match_flip2), "A1"]
               listRegion2[[2]][which(match_flip2), "A2"]=listRegion2[[1]][which(match_flip2), "A2"]
               listRegion2[[2]][which(match_flip2), "BETA"]=-listRegion2[[2]][which(match_flip2), "BETA"]
           }
           # eQTL
           match_comp_one2 = (listRegion2[[1]]["A1"] == complement_snp(listRegion2[[2]]["A1"])) & (listRegion2[[1]]["A2"]== complement_snp(listRegion2[[2]]["A2"]))
           match_comp_two2 = (listRegion2[[1]]["A1"] == complement_snp(listRegion2[[2]]["A2"])) & (listRegion2[[1]]["A2"] == complement_snp(listRegion2[[2]]["A2"]))
           snp_allele_match2 = match_flip2 | match_correct2 | match_comp_one2 | match_comp_two2
           print(listRegion2[[2]][!snp_allele_match2,])
           # methyl
           match_correct3 = (listRegion2[[1]]["A1"] ==listRegion2[[3]]["A1"]) & (listRegion2[[1]]["A2"]== listRegion2[[3]]["A2"])
           match_flip3 = (listRegion2[[1]]["A1"] == listRegion2[[3]]["A2"]) & (listRegion2[[1]]["A2"] == listRegion2[[3]]["A1"])
           if (any(which(match_flip3)>0)) {
               listRegion2[[3]][which(match_flip3), "A1"]=listRegion2[[1]][which(match_flip3), "A1"]
               listRegion2[[3]][which(match_flip3), "A2"]=listRegion2[[1]][which(match_flip3), "A2"]
               listRegion2[[3]][which(match_flip3), "BETA"]=-listRegion2[[3]][which(match_flip3), "BETA"]
           }
           # mQTL
           match_comp_one3 = (listRegion2[[1]]["A1"] == complement_snp(listRegion2[[3]]["A1"])) & (listRegion2[[1]]["A2"]== complement_snp(listRegion2[[3]]["A2"]))
           match_comp_two3 = (listRegion2[[1]]["A1"] == complement_snp(listRegion2[[3]]["A2"])) & (listRegion2[[1]]["A2"] == complement_snp(listRegion2[[3]]["A2"]))
           snp_allele_match3 = match_flip3 | match_correct3 | match_comp_one3 | match_comp_two3
           print(listRegion2[[3]][!snp_allele_match3,])

          if (sum(!snp_allele_match2)>0 | sum(!snp_allele_match3)>0) {
             if (sum(!snp_allele_match2)) {
             removed_snp_list = rbind(removed_snp_list, data.frame(Marker_removed=listRegion2[[2]]["SNP"][!snp_allele_match2], reason="Alleles do not match"))
             listRegion2[[2]] = listRegion2[[2]][which(snp_allele_match2),]
             }
             if (sum(!snp_allele_match3)) {
             removed_snp_list = rbind(removed_snp_list, data.frame(Marker_removed=listRegion2[[3]]["SNP"][!snp_allele_match3], reason="Alleles do not match"))
             listRegion2[[3]] = listRegion2[[3]][which(snp_allele_match3),]
             }
             listRegion2 <- lapply(listRegion2, function(x) x[x$SNP %in% Reduce(intersect, Map("[[", listRegion2, "SNP")), ])
             listRegion2 <- lapply(listRegion2, function(df){
                                 df[order(df$SNP),]
                                 })
            }
          } # end of have_alleles


           nsnps = nrow(listRegion2[[1]])
           message("Number of SNPs matching ", nsnps)
           if (nsnps>min_nsnps) {
              moloc = moloc_test(listRegion2, overlap=FALSE, prior_var, priors)       
              ppa = as.data.frame(signif(t(moloc$priors_lkl_ppa[4]), digits=2)) # configuration PPA using set priors
              bf = as.data.frame(signif(t(moloc$priors_lkl_ppa[2]), digits=2)) # configuration BF before priors
              names(bf) = paste("bf.",  names(bf), sep="")
              best.colnames =  rownames(moloc$best_snp)
              best.snp.PPA= as.data.frame(signif(t(moloc$best_snp$coloc_ppas), digits=2))
              names(best.snp.PPA) = paste("best.snp.PPA.",  best.colnames, sep="")
              best.snp= as.data.frame(t(as.character(moloc$best_snp$best.snp.coloc)))
              names(best.snp) = paste("best.snp.",  best.colnames, sep="")
              if (any(is.na(ppa))) stop("Moloc gives missing values for ", ProbeID, " ", ProbeIDmet)

              best.snp.betas.df =c()
              for (d in 1:length(listRegion2)) {
                 best.snp.betas =c()
                 for (s in (as.character(moloc$best_snp$best.snp.coloc))) {
                 b <- signif(listRegion2[[d]]$BETA[listRegion2[[d]]$SNP==s], digits=3)
                best.snp.betas = c(best.snp.betas, b)
                }
                best.snp.betas = paste(best.snp.betas, collapse=",")
                best.snp.betas.df = c(best.snp.betas.df,  best.snp.betas)
              }
              names(best.snp.betas.df) = paste("best.snp.betas.", c("GWAS", "eQTL", "mQTL"), sep="")

              best.snp.se.df =c()
              for (d in 1:length(listRegion2)) {
                 best.snp.se =c()
                 for (s in (as.character(moloc$best_snp$best.snp.coloc))) {
                 b = signif(listRegion2[[d]]$SE[listRegion2[[d]]$SNP==s], digits=3)
                 best.snp.se = c(best.snp.se, b)
                 }
                best.snp.se = paste(best.snp.se, collapse=",")
                best.snp.se.df = c(best.snp.se.df,  best.snp.se)
              }
              names(best.snp.se.df) = paste("best.snp.se.", c("GWAS", "eQTL", "mQTL"), sep="")

              nsnps = moloc$nsnps
              # Find min pval and best snps from the list of all the traits
              minpi = paste(lapply(listRegion2, function(x) signif(min(x$PVAL), digits=2)), collapse=",")
              minpi.snp = paste(lapply(listRegion2, function(x) x$SNP[which.min(x$PVAL)]), collapse=",")
         

              d <- letters[1:length(listRegion2)]
              configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
              configs <- do.call(paste0, configs_cases)
              coloc_configs <- configs[nchar(configs)>1]

              if (any(as.numeric(ppa[coloc_configs])>0.5)) {
              message("Found a colocalized signal!!!")
              if (save.SNP.info) {
                coloc.out = paste(outfolder, "/coloc.output.perSNP/", sep="")
                if (!file.exists(coloc.out)) dir.create(coloc.out)
                # dt <- do.call("rbind", lapply(listRegion2, "[", c("SNP", "CHR", "POS", "BETA", "SE", "PVAL")))
                ABF <- data.frame(ABF=adjust_bfs(listRegion2, overlap=FALSE, prior_var))
                ABF$SNP=rownames(ABF)
                dt <- rbindlist(listRegion2, idcol= "id", fill=TRUE)
                dt$id = gsub("1", "gwas", dt$id)
                dt$id = gsub("2", "eqtl", dt$id)
                dt$id = gsub("3", "mqtl", dt$id)
                dt = data.frame(dt)
                dt = merge(dt, ABF, by="SNP")
           write.table(x=dt, file=paste(coloc.out, ProbeID,'_', ProbeIDmet, '_results.tab', sep=''),row.names = FALSE, quote = FALSE, sep = '\t')

              }
              }

              names(ppa) = paste("PPA.",  names(ppa), sep="")

              res = cbind(ProbeID = ProbeID, ProbeIDmet = ProbeIDmet, CHR = chrom, START=pos.start, STOP=pos.end, nsnps, minpi=minpi, minpi.snp, bf, ppa, best.snp, best.snp.betas, best.snp.se, best.snp, best.snp.PPA)

              res.out = rbind(res.out,res)

           } # if nsnps>50 # after finding common alelles if have_alleles = TRUE
           } # if nrow(listRegion2[[1]])>0  # first matching of SNPs across data frames
       } # for each j (ProbeIDmet) #*# TRY PUTTING THIS AFTER?

      if(nrow(res.out)==0){
        return(NULL)
       }
      return(res.out)

       } #  matches_in_dfs

    write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
   }

  write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
  if (nrow(removed_snp_list)>0) write.table(x = removed_snp_list, file= paste(outfolder, prefix, '_removed.snps', sep=''), row.names = FALSE, quote = FALSE, sep = '\t')
  return(list(res.all, bed))

}
