### Results of chemostat modeling
## Loading libraries
library(tidyverse)
library(metR)
library(cmocean)

## Reading in results from mathematica chemostat simulations
chemostat_biomass_data<-data.table::fread('data/Energy-P-Biomass-Heatmap.csv')
chemostat_biomass_data<-chemostat_biomass_data %>%
  mutate(p_molar=`P Concentration (mol/m^3)`*1e-3)

## Creating F1
biomass_heatmap<-ggplot(chemostat_biomass_data)+
  geom_raster(aes(x=p_molar,y=`Relative Energy Availability (fraction of Earth)`,
                  fill=`Steady State Biomass (cells/ml)`))+
  geom_contour2(aes(x=p_molar,
                    y=`Relative Energy Availability (fraction of Earth)`,
                    z=`Steady State Biomass (cells/ml)`,
                    label=..level..),
                breaks=10^seq(3,14,by=1),
                fontface='bold',
                label_colour='white',
                color='white')+
  scale_y_log10(expand=c(0,0),name='Relative Growth Efficiency')+
  scale_x_log10(expand=c(0,0),name='Dissolved PO4 Concentration [M]')+
  cmocean::scale_fill_cmocean(name='deep',trans='log10',direction=-1)+
  geom_hline(linetype='dashed',col='white',yintercept=0.68)+
  geom_label(aes(x=1e-3,y=0.68,label='Fast-Growing E. coli'))+
  geom_hline(linetype='dashed',col='white',yintercept=5.56*1e-7)+
  geom_label(aes(x=1e-3,y=5.56*1e-7,label='Marine Sediment Bacteria'))+
  guides(fill=guide_colorbar(title='Steady State Cells/mL',
                             barwidth=unit(0.5,'npc')))+
  theme_bw()+
  theme(text=element_text(size=14),
        legend.position='bottom',
        axis.text.x=element_text(size=8))

## Creating extended data figure 1
mortality_plot<-ggplot(chemostat_biomass_data)+
  geom_point(aes(y=1/`Flow (Mortality) Rate (same units as growth rate)`,
                    x=`Steady State Biomass (cells/ml)`,
                 col=p_molar))+
  scale_y_log10(name='Cell Turnover Time')+
  scale_x_log10(name='Steady State Cells/mL')+
  cmocean::scale_color_cmocean(name='thermal',trans='log10',direction=1)+
  guides(color=guide_colorbar(title='Dissolved PO4\nConcentration [M]',
                             barwidth=unit(0.35,'npc')))+
  theme_bw()+
  theme(text=element_text(size=14),
        legend.position='bottom',
        axis.text.x=element_text(size=14),
        legend.text=element_text(size=8))

ggsave('figures/f1.pdf',scale=1,biomass_heatmap)
ggsave('figures/ed_f1.pdf',scale=1,mortality_plot)

