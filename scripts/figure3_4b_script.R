## Doing a MCMC type sampling of w/exhaustive grid search of 2 things:
## N:P power law exponents e[-1,1] (based on the fact that empirical is close to 0)
## PSDs prior conditioned on terrestrial PSD shapes using Gamma distributions
## Gamma will represent power laws when shape is extreme,
## Increasing PSD when shape is extreme the other way
## 'Humps' for intermediate sizes in other regimes
library(metR)
library(tidyverse)
## Start by picking a characteristic average 'N:P' to constrain models
## Using Redfield as a starting point, spread around this for a
## 'Redfield ratio' spanning 1 to 100. 

redfield<-seq(1,100,by=0.5)

## Now we want to sample from possible exponents for the scaling
## given that from our analysis of methanogen genomes from around
## Earth, the average N:P scaling exponent was approximately 0.03,
## We conservatively say the exponents exist in the range of -1 to 1

alpha<-seq(-1,1,by=0.05)

## Now based on our calculation of the average volume for the
## methanogen dataset from genome scaling, we end up with an average
## volume on the order of 10^-18 m^3, or ~1 um^3. 
## We will use this as our 'average volume' to associate with our 'redfield'
## characteristic number
char_vol<-1e-18

## This lets us calculate the intercept by doing the following operation
gamma0<-c()
k<-1
for(i in 1:length(redfield)){
  for(j in 1:length(alpha)){
    gamma0[k]<-exp(log(redfield[i])-alpha[j]*log(char_vol))
    k<-k+1
    print(k)
  }
}

## Now we've made our exhaustive search of power laws for potential
## N:P scaling relationships
## We can turn this into a data frame for easier multiplication 
## down the line

scaling_distribution<-data.frame(g=gamma0,
                                 a=rep(alpha,length(redfield)),
                                 mid=rep(redfield,each=length(alpha)))

## Now we're going to simulate potential size distributions 
## Assume everything is a power law size distribution and vary
## the exponent from 0 to -3 and simulate from a volume range from
## 10^-20 m^3 to 10^-17 m^3
psd_exp<-seq(-3,0,by=0.05)
psd_range<-10^seq(-20,-17,by=0.05)
psd_matrix<-matrix(data=NA,nrow=length(psd_exp),
                   ncol=length(psd_range))
for(i in 1:length(psd_exp)){
  psd_matrix[i,]<-psd_range^psd_exp[i]/sum(psd_range^psd_exp[i])
}

## Now we make a function to integrate our characteristic scaling
## with the power-law PSD

integrate_stoich<-function(psd,vol,a,g){
  scales<-g*vol^a
  total<-sum(psd*scales)
  return(total)
}

## Now we sample from the parameter space
out<-c()
k<-1
for(i in 1:nrow(scaling_distribution)){
  for(j in 1:nrow(psd_matrix)){
    out[k]<-integrate_stoich(psd=psd_matrix[j,],
                             vol=psd_range,
                             a=scaling_distribution$a[i],
                             g=scaling_distribution$g[i])
    k<-k+1
    print(k/(nrow(scaling_distribution)*nrow(psd_matrix)))
  }
}

## Generating figure 3
full_par_frame<-data.frame(result=out,
                           alpha=rep(scaling_distribution$a,
                                 each=length(psd_exp)),
                           gamma=rep(scaling_distribution$g,
                                 each=length(psd_exp)),
                           mid=rep(scaling_distribution$mid,
                                   each=length(psd_exp)),
                           exp=rep(psd_exp,
                                   nrow(scaling_distribution)))

hold_alpha<-ggplot(full_par_frame %>% group_by(exp,alpha) %>%
         summarize(mean_result=10^mean(log10(result))))+
  geom_raster(aes(x=exp,y=alpha,fill=mean_result))+
  scale_fill_viridis_c(name='N:P',trans='log10')+
  theme_bw()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab('Cell Size-Abundance Scaling Exponent')+
  ylab('N:P Scaling Exponent')+
  stat_contour2(aes(x=exp,y=alpha,z=mean_result,
                    label=after_stat(level)),
                breaks=c(1,5,10,50,100,500,1000))+
  geom_point(aes(x=-1,y=0.03),size=4,col='white')+
  ggrepel::geom_label_repel(data=data.frame(),
                            aes(x=-1,y=0.03,label='Earth Conditions'),
                            size=4)

hold_mid<-ggplot(full_par_frame %>% group_by(exp,mid) %>%
         summarize(mean_result=10^mean(log10(result))))+
  geom_raster(aes(x=exp,y=mid,fill=mean_result))+
  scale_fill_viridis_c(name='N:P',trans='log10')+
  theme_bw()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab('Cell Size-Abundance Scaling Exponent')+
  ylab('N:P Scaling Intercept\n(N:P For 1 um Diameter Cell)')+
  stat_contour2(aes(x=exp,y=mid,z=mean_result,
                    label=after_stat(level)),
                breaks=c(1,5,10,50,100,500,1000))+
  geom_point(aes(x=-1,y=16),
             size=4,col='white')+
  ggrepel::geom_label_repel(data=data.frame(),
                            aes(x=-1,y=16,label='Earth Conditions'),
                            size=4)

hold_alpha/hold_mid+plot_annotation(tag_levels='A')

ggsave('figures/f3.pdf',scale=1.25)

## Writing extended data files
write.csv(full_par_frame,'data/parameter_search_simulations.csv')
