CpGisl.level.info = function(set = ILS.hypo.drift5, isls = Islands.hypo.drift5, 
                             cpgs = CpGs.hypo.drift5.isl, dat, full = T, prom = T) {
  # input set: island names of interest
  # input isls and cpgs: vectors of the same length with CpG names and associated island names
  # input data: methylation data 
  # input str1 and str2: optional, to constrain selection of CpGs for mean calculation
  
  len = length(set) # number of Islands
  Mvals = sd.Mvals = matrix(0,ncol=ncol(dat),nrow=len)
  genes = list()
  islands = vector()
  Prom = rep(T,len)
  
  ii = 1
  for (i in 1:len) {
    
    if(full==T) {
      cpgsi = manifestDATA$Name[manifestDATA$UCSC_CpG_Islands_Name == set[i]]  # all cpgs on the island
    } else {
      cpgsi = cpgs[isls == set[i]]  # cpgs provided on the island
    }
    
    # remove shelves !!! already checked
    dum = manifestDATA[cpgsi,"Relation_to_UCSC_CpG_Island"]
    cpgsi = cpgsi[!grepl("Shelf",dum)]  # & dum!="N_Shore" & dum!="S_Shore"] 
    
    # check exclusions (only informative islands etc. )
    dum = manifestDATA[cpgsi,"UCSC_RefGene_Group"]
    tmp = cpgsi[grepl("TSS",dum)] # & !grepl("Body",dum) & !grepl("3'UTR",dum)]
    if(prom == T) { # only promoters islands are included, else skip 
      if(length(tmp)==0) {next} else {cpgsi = tmp}
    }
    
    if(length(tmp)==0) Prom[ii]=FALSE
    
    # see which cpgs are in dat, if any .... 
    idum = na.omit(match(cpgsi,rownames(dat)))
    len.idum = length(idum)
    if(len.idum > 0) {
      if(len.idum ==1) {
        Mvals[ii,]=dat[idum,]; sd.Mvals[ii,]=0
      } else {
        Mvals[ii,] = apply(dat[idum,],2,mean); sd.Mvals[ii,] = apply(dat[idum,],2,sd)
      }
    
      # get gene name for island
      genes[[ii]] = dum = unique(unlist(strsplit(manifestDATA[cpgsi,"UCSC_RefGene_Name"],split=';')))
      islands[ii] = set[i]
    
      ii = ii+1
    }
  }
  
  rownames(Mvals) = rownames(sd.Mvals) = set
  colnames(Mvals) = colnames(sd.Mvals) = colnames(dat)
  Mvals = Mvals[1:(ii-1),]; sd.Mvals = sd.Mvals[1:(ii-1),]

  tmp = lapply(genes,strsplit,';')
  tmp = lapply(tmp,unlist)
  genes = lapply(tmp,unique)
  Prom = Prom[1:(ii-1)]
  
  return(list(Mvals=Mvals, sd.Mvals=sd.Mvals, genes=genes, islands=islands, prom=Prom, ncheck = ii-1))
}


