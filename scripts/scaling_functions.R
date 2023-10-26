## The purpose of this script is to outline the allometric scaling
## functions for macromolecules as determined by empirical scaling studies
## Finkel et al, 2016 in the main text references
## (Protein, Lipids, Carbs)
## And biophysical scaling calculations from Kempes et al 2016
## (Genome, Ribosomes)
## Macromolecular - inversion of convert_to_genome
genome_to_volume<-function(bp){
  vol<-exp(log(bp/(6.3*1e26))/1.16) # function of bp
  return(vol)
}

convert_to_protein<-function(volume){
  out<-10^(1.09+0.83*log10(volume/100)) # Finkel 2016 scaling um^3 -> pg
  return(out)
}

convert_to_genome<-function(volume){
  out<-(6.3*1e26)*(volume^1.16) #scaling m^3 to bp
  out<-out*(617.96/(6.02*1e23))*1e12 #scaling bp to pg
  return(out)
}

convert_to_ribo<-function(volume){
  out<-1.54e-7*(volume^(0.73)) #cell vol to m^3 ribo
  out<-out*(4.32/3.04)*1e18 # m^3 ribo to pg ribo
  return(out)
}

convert_to_lipid<-function(volume){
  #out<-6.31e-12*(volume^(0.76))
  out<-10^(0.8+0.8*log10(volume/100)) # Finkel 2016 scaling um^3 -> pg
  return(out)
}

convert_to_carb<-function(volume){
  out<-10^(0.72+0.93*log10(volume/100)) # Finkel 2016 scaling um^3 -> pg
  return(out)
}

## All of these functions take total pg macromolecule as the input
## and return pg C,N,P as output as defined by the ratios 
## from La Roche Redfield Revisited cited in the main text
protein_to_cnp<-function(pg){
  c<-0.53*pg
  n<-0.16*pg
  p<-0*pg
  return(c(c,n,p))
}

carbs_to_cnp<-function(pg){
  c<-0.4*pg
  n<-0*pg
  p<-0*pg
  return(c(c,n,p))
}

genome_to_cnp<-function(pg){
  c<-0.36*pg
  n<-0.16*pg
  p<-0.095*pg
  return(c(c,n,p))
}

lipid_to_cnp<-function(pg){
  phospho<-mean(c(5/65,50/65)) # average proportion phospholipids
  ### 42.3% 
  c<-0.76*pg*(1-phospho)+(0.64*pg*phospho)
  n<-0.008*pg*phospho
  p<-0.043*phospho*pg
  return(c(c,n,p))
}

ribo_to_cnp<-function(pg){
  prot_ribo<-1.66/4.32 #proportion of ribosomal mass that is protein vs rna
  c<-0.53*pg*prot_ribo+(0.34*pg*(1-prot_ribo))
  n<-0.16*pg*prot_ribo+(0.155*pg*(1-prot_ribo))
  p<-0.091*pg*(1-prot_ribo)
  return(c(c,n,p))
}


