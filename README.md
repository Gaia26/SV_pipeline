# SV_pipeline

## Collapse Analysis 
*Collapse_Analysis.py* is a python script needed to:
*  count how many variants are present in each gene;
*  count in how many samples each gene is seen;
*  annotate the genes.

It takes in input three file:
- *report*, a *.csv* or *.tsv* file format containing information exctracted from the annotated VCF outputted from the pipeline;
- *hp0002205*, a *.csv* or *.tsv* file format downloaded from HPO; collection of genes involved with recurrent respiratory infections;
- *covid_panel*, a *.csv* or *.tsv* file format downloaded from PanelApp; collection of genes involved in the COVID-19 research.
