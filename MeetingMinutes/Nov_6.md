# Meeting Agenda
1. Review progress from last week (new filtered phyloseq object, new analyses)
2. Should we/do we know how to use git in R for version and directory control?
3. Questions, any troubleshooting, objectives for next week

# Meeting Notes
**November 13 will be a zoom meeting at 2 pm**

* For next week, bring legend and short interpretation of the results + try aim 3?
* To stay on track, by November 20th acquire all of the data
* Create a readme for each aim directory and write any results or some notes
* Move all .qza files onto the local computer as a backup in case of issues with the server
  
* Stick with the same colour scheme/theme for all figures (choose colours that make sense):
  * ColorBrewer2
  * colorblind safe
  * categorical/qualitative colors
* Add a figure legends

### Aim 1
* PCoA plot:
  *  The "U" pattern may represent something important (Hans recalls hearing something in the past but is unsure)
  * Confirm significance and if the "golden standard" was used for beta diversity
  * check the richness of the samples and explore which microbial populations were lost
* Alpha Diversity:
  * add more details about the significance ("ns")
  * consider testing other metrics of diversity
  * consider makeing a violin plot
* Consider a bee swarm plot for the shannon_div.png box plot

### Aim 2
* Core Microbiome
  * Try to find out abundance of unique microbial species identified through core microbiome analysis (so they present lower/higher between MS types?)
* Heatmap could be used as a supplemental figure to explain why the specific abundance threshold was set. (Could alternatively do bar charts). 

* Indicator species analysis
  * do not worry too much about the statistics
  * low abundance does not correlate with relevance: look into species even if they are not abundant
* core microbiome could instead be used for indicator species analysis
  
* DESeq2
  * look into top 3-5 most and least abundant genera
