## Conducting dissolved:particulate scaling analysis
## and generating results for Figure 4 and extended data
#loading libraries
library(tidyverse)
library(lmodel2)
library(ggExtra)
library(BGLR)
source('scripts/scaling_functions.R')
source('scripts/figure3_4b_script.R')

## Reading in data
stoich_data<-data.table::fread('data/Stoichiometries-Moore.csv')

# Converting umol/kg to umol:umol C to compare with cell quota data
stoich_data<-stoich_data %>%
  mutate(dissolved_rat=`Mean ocean concentration(umol kg)`/2247) %>%
  drop_na(`Mean ocean concentration(umol kg)`) %>%
  arrange(desc(`Phytoplankton quota(mol:mol C)`))

## Adding additional data about the elements (position on the periodic table)
# and # redox states
stoich_data$block<-c('P','P','P','S','P','P',
                     'S','S','D','S','D','D',
                     'D','D','D','D','D')

stoich_data$n_redox<-c(9,4,3,2,5,4,2,2,3,2,4,2,2,2,2,3,3)

## Simulating out all combinations of elemental ratios
cell_quotas<-c()
dissolved_stoich<-c()
element_combos<-c()
k<-1

for(i in 1:nrow(stoich_data)){
  for(j in 1:nrow(stoich_data)){
    if(j<=i){
      next
    }
    element_combos[k]<-paste0(stoich_data$Element[i],'_',stoich_data$Element[j])
    dissolved_stoich[k]<-stoich_data$`Mean ocean concentration(umol kg)`[i]/stoich_data$`Mean ocean concentration(umol kg)`[j]
    cell_quotas[k]<-stoich_data$`Phytoplankton quota(mol:mol C)`[i]/stoich_data$`Phytoplankton quota(mol:mol C)`[j]
    k<-k+1
  }
}
## Checking calculations worked
plot(cell_quotas,dissolved_stoich,log='xy')
hist(log10(cell_quotas),breaks=20)
hist(log10(dissolved_stoich),breaks=20)

## Generating a table for plotting of all modeled ratios
ratio_analysis<-data.frame(combo=element_combos,
                           main_element=gsub('_.*$','',element_combos),
                           secondary_element=gsub('^.*_','',element_combos),
                           bio=cell_quotas,
                           abio=dissolved_stoich)
write.csv(ratio_analysis,'data/dissolved_particulate_ratios.csv')

## Trying different power law estimation methods 
## to make sure estimates converge to similar results
ols_model<-lm(log10(cell_quotas)~log10(dissolved_stoich))
typeii_out<-lmodel2(log10(bio)~log10(abio),range.x='interval',
                    range.y='interval',data=ratio_analysis,nperm=1e4)
bayes_linear<-BLR(log10(cell_quotas),XR=log10(dissolved_stoich))
bayes_noeffect<-BLR(log10(cell_quotas))

## Creating figure 4A
joint_dist<-ggplot(ratio_analysis)+
  geom_density_2d_filled(aes(y=log10(bio),x=log10(abio)))+
  geom_abline(slope=1,intercept=0,col='white',size=2,linetype='dashed')+
  geom_abline(slope=bayes_linear$bR,
              intercept=bayes_linear$mu,
              col='white',
              size=2)+
  geom_point(aes(x=log10(abio),y=log10(bio)),
             shape=1,size=1.5,col='white')+
  ggrepel::geom_label_repel(data=filter(ratio_analysis,combo %in% c('C_N','C_P',
                                                     'N_P','P_Fe',
                                                     'N_Fe','C_Fe',
                                                     'N_Mg','C_Mg','P_Mg')),
             aes(x=log10(abio),
                 y=log10(bio),
                 label=paste0(main_element,'/',secondary_element)))+
  scale_y_continuous(expand=c(0,0),
                     name='Log10 Phytoplankton Stoichiometry',
                     breaks=seq(0,7,by=1))+
  scale_x_continuous(expand=c(0,0),
                     name='Log10 Seawater Stoichiometry',
                     breaks=seq(-6,10,by=2))+
  scale_fill_viridis_d(option='B',name='Probability Density')+
  theme_bw()+
  theme(legend.position='bottom')


## Updating for marginal distributions
out2<-ggExtra::ggMarginal(joint_dist,aes(x=log10(abio),y=log10(bio)),type='histogram')

## Estimating ranges for Enceladus based on N:P ratios of 10-10^6
sample_intercepts<-rnorm(1000,bayes_linear$mu,bayes_linear$SD.mu)
sample_slopes<-rnorm(1000,bayes_linear$bR,sqrt(bayes_linear$varBR))
errors<-rnorm(1000,0,sqrt(bayes_linear$varE))
stoichs<-seq(1,6,by=0.05)
outcome_mat<-matrix(NA,nrow=1000,ncol=length(stoichs))
for(i in 1:length(stoichs)){
  outcome_mat[,i]<-errors+sample_intercepts+stoichs[i]*sample_slopes
}

plot(stoichs,colMeans(outcome_mat))

search_density_est<-density(log10(full_par_frame$result))
endpoint<-search_density_est$y[which(round(search_density_est$x,2) %in% c(1.78,2.80))]   
min_x<-which(round(search_density_est$x,2)==1.78)
max_x<-which(round(search_density_est$x,2)==2.8)

## Probability of a simulated ecosystem scaling like Enceladus stoichiometries
p_mass<-sum(search_density_est$y[min_x:max_x])/sum(search_density_est$y)
print(p_mass)

## Creating figure 4B
density_with_annotation<-ggplot()+
  geom_ribbon(aes(search_density_est$x[min_x:max_x],
                ymax=search_density_est$y[min_x:max_x],
                ymin=0.001),
                col='cyan3',
            fill='cyan3',
            alpha=0.25)+
  geom_density(aes(x=log10(result)),
               col='black',
               size=1.5,
               data=full_par_frame)+
  geom_vline(aes(xintercept=log10(16)),
             col='red4',
             linetype='dashed',size=1.5)+
  geom_ribbon(aes(x=log10(5:50),
                  ymin=0.001,
                  ymax=0.4),
              col='red4',fill='red4',
              alpha=0.25)+
  geom_point(aes(x=c(1.78,2.8),y=endpoint),shape=4,
             size=6,col='cyan3',
             stroke=2)+
  geom_text(aes(x=-1,y=0.375,label='Earth Bulk N:P'),
            col='red4',fontface='bold')+
  geom_text(aes(x=-0.85,y=0.35,label='Enceladus-Like\nConditions'),
            col='cyan3',fontface='bold')+
  scale_y_continuous(name='Probability Density Over Simulations',limits=c(0,0.4),
                     expand=c(0,0))+
  theme_bw()+
  xlab('Log10 Ecosystem N:P')

library(patchwork)
## Composing final figure 4
full_fig4<-(patchwork::wrap_elements(out2)|density_with_annotation)+
  plot_annotation(tag_levels='A')
ggsave('figures/f4.pdf',full_fig4,width=unit(10,'in'),height=unit(8,'in'))
