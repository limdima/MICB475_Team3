Running Experimental Aim 1
- rarify_phyloseq.R code is used to rarify phyloseq object and remove the control group from disease course (since only MS sub types are compared for alpha/beta analysis)
- Alpha Beta Diversity.R runs alpha and beta metrics and creates plots


***PLOT NOTES***
***Alpha diversity:***

- Shannon : richness and evenness. comparing this metric between disease courses's indicate all groups have simular balance of species, mostly comparable in number of species and how evenly distributed the species are. 

- Simpson : also measure richness and evenness, but gives more weight to the most abundant species (values closer to zero indicate higher diversity). Groupings are simular for dominant species' presence and evenness. 

- Chao1 : richness (focuses on unseen/rare species in the samples). Potential unobserved or rare species (species present in low frequencies) are simular across groups. However, high Chao1 suggests each group has high richness (perhaps complex communities), but since this is shared between groups it shows that each disease course results in a shared community complexity.

- Observed : simular number of directly counted species. shows straightforward count of total unique species observed in the sample.


***Beta diversity:***

- LEFT: Weighted unifrac - measures differences in community composiiton (takes into account phylogenetic relationships and abundance). Groupings are simular in composition and the evolutionary relationships.
- *weighted unifrac: R2 is 0.01868, p-value is 0.002*

- RIGHT: Bray Curtis - focus on species abundance; there are not clear differences on the PCoA plot suggesting the groups have simular composition in terms of presence and abundance. Groupings are populated by simular species.
    - low % variance - variability is spread across many dimensions
    - *bray: R2 is 0.01198, p-value is 0.001*
 
