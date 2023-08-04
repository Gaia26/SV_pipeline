# SV_pipeline

## Collapse Analysis 
*Collapse_Analysis.py* is a python script needed to:
*  count how many variants are present in each gene;
*  count in how many samples each gene is seen;
*  annotate the genes with list of genes involved in COVID19 or chronic respiratory infections from NCBI, PanelApp and Human Phenotype Ontology.
Before performing the collapse analysis, the script is intended to manipulate the file coming from the annotation with VEP and performing the following operations:
*  filtering based on MAF;
*  split based on SVTYPE;
*  transposition of the dataframe on the specific columns, the ones corresponding to samples and genotypes.



The inputs needed are:
*  _report_tab_: file report coming from VEP (as _661_covid_patients_finalreport.csv_);
*  _file_patients_ = excel file containing all info on patients (as _661_samples_info_complete.xlsx_);
*   _hpo_: list of genes involved in chronic respiratory infections coming from human phenotype ontology (as _Recurrent_respiratory_infections.csv_);
*  _covid_tab_: list of genes involved in COVID19 from previous researches coming from PanelApp (as _COVID-19 research.tsv_);
*  _covid_ncbi_tab_:  list of genes involved in COVID19 from previous researches coming from NCBI (as _gene_result.txt_).
The last three files are necessary to perform the annotation at the end of the collapse.
  
To run the script use the following code:
_python Collapse_Analysis.py 661_covid_patients_finalreport.csv 661_samples_info_complete.xlsx Recurrent_respiratory_infections.csv COVID-19 research.tsv gene_result.txt_
