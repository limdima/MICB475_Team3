NOTE: Currently, the plots for disease course and treatment status have the same colours, but I will change the colours for the treatment status! Disease course colours will stay

***FOR DISEASE COURSE***
- notes to be added




***FOR TREATMENT STATUS***

***Alpha diversity: (alpha_diversity_plots.png)***

- Shannon : richness and evenness. comparing this metric between treatment status's indicate all groups have simular balance of species, mostly comparable in number of species and how evenly distributed the species are. 

- Simpson : also measure richness and evenness, but gives more weight to the most abundant species (values closer to zero indicate higher diversity). Groupings are simular for dominant species' presence and evenness. 

- Chao1 : richness (focuses on unseen/rare species in the samples). Potential unobserved or rare species (species present in low frequencies) are simular across groups. However, high Chao1 suggests each group has high richness (perhaps complex communities), but since this is shared between groups it shows that each treatment status results in a shared community complexity.

- Observed : simular number of directly counted species. shows straightforward count of total unique species observed in the sample.


***Beta diversity: (beta_plots.png)***

- LEFT: Weighted unifrac - measures differences in community composiiton (takes into account phylogenetic relationships and abundance). Groupings are simular in composition and the evolutionary relationships.
- *weighted unifrac: R2 is 0.01868, p-value is 0.002*

- RIGHT: Bray Curtis - focus on species abundance; there are not clear differences on the PCoA plot suggesting the groups have simular composition in terms of presence and abundance. Groupings are populated by simular species.
    - low % variance - variability is spread across many dimensions
    - *bray: R2 is 0.01198, p-value is 0.001*
 
