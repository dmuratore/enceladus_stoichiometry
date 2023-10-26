### The purpose of this script is to run the analyses for the
## Macromolecular scaling from the methanogen genomic database
## cited in the main text, described in Methods
## Loading libraries
library(tidyverse)
library(MASS)
source('scripts/scaling_functions.R')
## Read in data
methanogen_table<-readxl::read_xlsx('data/dsac048_suppl_supplementary_tables.xlsx',
                                    skip=1)
psychrophiles<-methanogen_table %>%
  filter(`Assembly Level`=='Complete Genome',
         (`temperature group`=='psychrotolerant'|`T min`<=0)) %>%
  drop_na(`T min`)

## Modeling the rank-genome size distribution
rank_model<-lm(log10(`Genome Size [bp]`)~rank(1/`Genome Size [bp]`),
               data=methanogen_table %>% filter(`Assembly Level`=='Complete Genome',
                                                !is.na(`temperature group`)))

## Creating Figure 2a
genome_size_dist<-ggplot(methanogen_table %>% filter(`Assembly Level`=='Complete Genome',
                                   !is.na(`temperature group`)))+
  geom_point(aes(x=rank(1/`Genome Size [bp]`),
                 y=`Genome Size [bp]`,
                 col=`temperature group`),
             size=2)+
  theme_bw()+
  stat_function(fun=function(x) 10^(6.65-0.0116*x),
                linetype='dashed')+
  scale_color_viridis_d(name='Temperature Tolerance',
                        labels=c('Mesophilic',
                                 'Psychrotolerant',
                                 'Thermotolerant'))+
  xlab('Rank Genome Size')+
  ylab('Genome Size [bp]')

## Fitting a gamma distribution to the genome length distribution to impute
## A smooth size distribution
full_methanogens<-methanogen_table %>% filter(`Assembly Level`=='Complete Genome')
out<-fitdistr(log10(full_methanogens$`Genome Size [bp]`),
         'gamma')
## Writing to output
methanogen_table %>% filter(`Assembly Level`=='Complete Genome') %>%
  write.csv(file='data/methanogen_genomes.csv')

## Estimating cell volumes from genome length
with_vol<-methanogen_table %>% filter(`Assembly Level`=='Complete Genome') %>%
  mutate(est_volume=genome_to_volume(`Genome Size [bp]`)*10^18) %>%
  mutate(pg_prot=convert_to_protein(est_volume),
         pg_genome=convert_to_genome(est_volume/10^18),
         pg_ribo=convert_to_ribo(est_volume/10^18),
         pg_lipid=convert_to_lipid(est_volume),
         pg_carb=convert_to_carb(est_volume))

## Fitting a gamma distribution to cell volumes and plotting Figure 2B
vol_dist<-MASS::fitdistr(log10(with_vol$est_volume),'gamma')

better_vol<-ggplot(with_vol)+
  geom_freqpoly(aes(x=log10(est_volume)),stat='density',
                linewidth=2)+
  stat_function(fun=dgamma,args=list(shape=vol_dist$estimate[1],
                                     rate=vol_dist$estimate[2]),
                aes(x=log10(est_volume)),
                linetype='dashed')+
  theme_bw()+
  ylab('Relative Abundance')+
  xlab('Log10 Cell Volume [um^3]')

## Using the gamma distribution to impute a smooth cell size distribution
theoretical_psd<-data.frame(vol=seq(0.1,0.8,by=0.01)) %>%
  mutate(dn=dgamma(vol,shape=vol_dist$estimate[1],
                   rate=vol_dist$estimate[2])) %>%
  mutate(vol=(10^vol)) %>%
  mutate(pg_prot=convert_to_protein(vol),
         pg_genome=convert_to_genome(vol/1e18),
         pg_ribo=convert_to_ribo(vol/1e18),
         pg_lipid=convert_to_lipid(vol),
         pg_carb=convert_to_carb(vol))

## Plotting macromolecular components as a function of cell volume
ggplot(theoretical_psd)+
  geom_line(aes(x=vol,y=pg_prot*dn,
                 col='Protein'))+
  geom_line(aes(x=vol,y=pg_genome*dn,
                col='Genome'))+
  geom_line(aes(x=vol,y=pg_ribo*dn,
                col='Ribosome'))+
  geom_line(aes(x=vol,y=pg_lipid*dn,
                col='Lipid'))+
  geom_line(aes(x=vol,y=pg_carb*dn,
                col='Carbohydrate'))+
  theme_bw()+
  xlab('Volume [um^3]')+
  scale_y_log10(name='Weighted Molecule Abundance')+
  scale_color_manual(name='Macromolecule',
                     labels=c('Carbohydrate',
                              'Genome',
                              'Lipid',
                              'Protein',
                              'Ribosome'),
                     values=c('brown4',
                              'magenta2',
                              'darkgreen',
                              'midnightblue',
                              'purple2'))

## Summing over the cell size distribution
community_total<-c('Carbs'=sum(theoretical_psd$pg_carb*theoretical_psd$dn),
                   'Genome'=sum(theoretical_psd$pg_genome*theoretical_psd$dn),
                   'Lipid'=sum(theoretical_psd$pg_lipid*theoretical_psd$dn),
                   'Protein'=sum(theoretical_psd$pg_prot*theoretical_psd$dn),
                   'Ribo'=sum(theoretical_psd$pg_ribo*theoretical_psd$dn))


community_frame<-data.frame(molecule=names(community_total),
                            value=community_total)

## Writing output to data
theoretical_psd %>%
  rename('vol_um3'=vol,
         'relative_frequency'=dn) %>%
  write.csv('data/psd_with_macromolecules.csv')

## Doing stoichiometric conversions of macromolecules to C,N,P
theoretical_stoich<-rbind(carbs_to_cnp(community_total[1]),
                          genome_to_cnp(community_total[2]),
                          lipid_to_cnp(community_total[3]),
                          protein_to_cnp(community_total[4]),
                          ribo_to_cnp(community_total[5]))
colnames(theoretical_stoich)<-c('pg_C',
                                'pg_N',
                                'pg_P')
rownames(theoretical_stoich)<-c('carbs','genome','lipid',
                                'protein','ribosome')
theoretical_stoich<-data.frame(theoretical_stoich) %>%
  rownames_to_column('macromolecule')

## Preparing for plotting
protein_by_size<-matrix(protein_to_cnp(theoretical_psd$pg_prot),ncol=3)
genome_by_size<-matrix(genome_to_cnp(theoretical_psd$pg_genome),ncol=3)
ribo_by_size<-matrix(ribo_to_cnp(theoretical_psd$pg_ribo),ncol=3)
lipid_by_size<-matrix(lipid_to_cnp(theoretical_psd$pg_lipid),ncol=3)
carb_by_size<-matrix(carbs_to_cnp(theoretical_psd$pg_carb),ncol=3)

element_by_size<-cbind(protein_by_size,
                       genome_by_size,
                       ribo_by_size,
                       lipid_by_size,
                       carb_by_size)

colnames(element_by_size)<-paste0(rep(c('prot','gen','ribo','lip','carb'),
                                      each=3),'_',
                                  rep(c('C','N','P'),5))

psd_with_element<-cbind(theoretical_psd,element_by_size) %>%
  mutate(total_C=prot_C+gen_C+ribo_C+lip_C+carb_C,
         total_N=prot_N+gen_N+ribo_N+lip_N+carb_N,
         total_P=prot_P+gen_P+ribo_P+lip_P+carb_P)

##Plotting Figure 2C
element_dist<-ggplot(psd_with_element)+
  geom_line(aes(x=vol^(1/3)*0.75/pi,y=total_C*dn,
                col='C'))+
  geom_line(aes(x=vol^(1/3)*0.75/pi,y=total_N*dn,
                col='N'))+
  geom_line(aes(x=vol^(1/3)*0.75/pi,y=total_P*dn,
                col='P'))+
  scale_y_log10()+
  xlab('Cell Radius [um]')+
  ylab('pg DW Element')+
  theme_bw()+
  scale_color_brewer(name='Element',palette='Set1')

psd_props<-psd_with_element %>%
  pivot_longer(c(contains('pg'),contains('_')),
               names_to='pool',values_to='value')

psd_mac<-psd_props %>% filter(str_detect(pool,'pg'))

psd_element<-psd_props %>% filter(str_detect(pool,'_.$')) %>%
  mutate(element=gsub('^.*_','',pool),
         mac=gsub('_.*$','',pool)) %>%
  filter(!str_detect(pool,'total'))

macro_scale<-ggplot(psd_element)+
  geom_col(aes(x=vol^(1/3)*0.75/pi,
               y=value,
               fill=factor(mac,
                           levels=c('prot','lip','carb','gen','ribo'),
                           labels=c('Protein','Lipid','Carbohydrate','Genome',
                                    'Ribosome'))))+
  facet_wrap(~element,ncol=1,scale='free_y')+
  scale_fill_manual(name='Macromolecule',
                    labels=c('Carbohydrate',
                             'Genome',
                             'Lipid',
                             'Protein',
                             'Ribosome')[c(4,3,1,5,2)],
                    values=c('brown4',
                             'magenta2',
                             'darkgreen',
                             'midnightblue',
                             'purple2')[c(4,3,1,5,2)])+
  theme_bw()+
  ylab('Mass Element [pg]')+
  xlab('Cell radius [um]')+
  scale_x_log10()+
  theme(strip.background=element_blank(),
        strip.text=element_text(size=14,face='bold'))



## Plotting Figure 2D
np_allo<-ggplot(psd_with_element)+
  geom_line(aes(x=vol^(1/3)*0.75/pi,y=total_N/total_P))+
  theme_bw()+
  xlab('Cell Radius [um]')+
  ylab('Macromolecular N:P Ratio')

np_dist<-ggplot(psd_with_element)+
  geom_line(aes(y=dn,x=(total_N/total_P)))+
  theme_bw()+
  xlab('Macromolecular N:P Ratio')+
  ylab('Expected Abundance')
  
## Assembling completed figure
library(patchwork)

full_f2<-(genome_size_dist+theme(legend.position='bottom')|better_vol)/(macro_scale+theme(legend.position='bottom')|np_dist)+
  plot_annotation(tag_levels = 'A')
ggsave('figures/f2.pdf',full_f2,scale=1.75)


