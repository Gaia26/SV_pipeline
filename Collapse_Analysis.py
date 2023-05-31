import pandas as pd
import sys

report = sys.argv[1] # your report file. Informations exctracted from the annotated VCF
hp0002205 = sys.argv[2] # recurrent respiratory infections table for annotations from HPO
covid_panel = sys.argv[3] # COVID-19 panel from PanelApp

tab = pd.read_csv(report, delimiter = "\t")
respiratory_infections = pd.read_csv(hp0002205, delimiter = "\t")
covid = pd.read_csv(covid_panel,  delimiter = "\t")

# create a new column Max Allele Frequency that takes the max value of the SV_overlap_AF
for i in range(len(tab)):
    allele_freq = tab.loc[i, 'SV_overlap_AF'].split(",")[0]
    max_allele_freq = max(allele_freq.split("&"))
    tab.loc[i, 'Max Allele Frequency'] = max_allele_freq

# move the new column to fifth
tab_columns = tab.columns.tolist()
last_column = tab_columns[-1]
df_without_last = tab.drop(columns=[last_column])
fifth_column = tab_columns[4]
df_new = pd.concat([df_without_last.iloc[:, :4], tab[last_column], df_without_last.iloc[:, 4:]], axis=1)


# filter the report table for SVLEN and Max Allele Frequency
filtered_tab = df_new[(df_new['SVLEN'] <= 50000) & ((df_new['Max Allele Frequency'] == '.') | (df_new['Max Allele Frequency'].apply(lambda x: x != '.' and float(x) <= 0.0001)))].reset_index(drop=True)

# 1. creating table for each SVTYPE
inv_df = filtered_tab[filtered_tab['SVTYPE'] == 'INV']
INV_reshape = pd.melt(inv_df, id_vars=inv_df.columns[:36], value_vars=inv_df.columns[36:], var_name='sample_name', value_name='genotype')
INV = INV_reshape.drop(INV_reshape[INV_reshape['genotype'] == ' ./.'].index)
column_names_inv = INV.columns.tolist()
new_column_order_inv = column_names_inv[-2:] + column_names_inv[:-2]
INV = INV[new_column_order_inv].reset_index(drop=True)


del_df = filtered_tab[filtered_tab['SVTYPE'] == 'DEL']
DEL_reshape = pd.melt(del_df, id_vars=del_df.columns[:36], value_vars=del_df.columns[36:], var_name='sample_name', value_name='genotype')
DEL = DEL_reshape.drop(DEL_reshape[DEL_reshape['genotype'] == ' ./.'].index)
column_names_del = DEL.columns.tolist()
new_column_order_del = column_names_del[-2:] + column_names_del[:-2]
DEL = DEL[new_column_order_del].reset_index(drop=True)


ins_df = filtered_tab[filtered_tab['SVTYPE'] == 'INS']
INS_reshape = pd.melt(ins_df, id_vars=ins_df.columns[:36], value_vars=ins_df.columns[36:], var_name='sample_name', value_name='genotype')
INS = INS_reshape.drop(INS_reshape[INS_reshape['genotype'] == ' ./.'].index)
column_names_ins = INS.columns.tolist()
new_column_order_ins = column_names_ins[-2:] + column_names_ins[:-2]
INS = INS[new_column_order_ins].reset_index(drop=True)


dup_df = filtered_tab[filtered_tab['SVTYPE'] == 'DUP']
DUP_reshape = pd.melt(dup_df, id_vars=dup_df.columns[:36], value_vars=dup_df.columns[36:], var_name='sample_name', value_name='genotype')
DUP = DUP_reshape.drop(DUP_reshape[DUP_reshape['genotype'] == ' ./.'].index)
column_names_dup = DUP.columns.tolist()
new_column_order_dup = column_names_dup[-2:] + column_names_dup[:-2]
DUP = DUP[new_column_order_dup].reset_index(drop=True)

# 2. Define two functions for counting how many variants are present in each gene and in which samples
def counter_gene(tab):
    dict_gene = {}
    seen_elements = []
    for i in range(len(tab)):
        gene_list = tab.loc[i, 'SYMBOL'].split(",")
        chrom = tab.loc[i,'CHROM']
        pos = tab.loc[i,'POS']
        end = tab.loc[i,'END']
        if (chrom, pos, end, gene_list) not in seen_elements: # if the variant in the gene has already been seen, do not consider it again
            for gene in gene_list:
                if gene != '.':
                    if gene not in dict_gene:
                        dict_gene[gene] = 0
                    dict_gene[gene] += 1
        seen_elements.append((chrom, pos, end, gene_list))
    return dict_gene
        
    
def tracker_sample(tab):
    dict_sample = {}
    for i in range(len(tab)):
        sample = (tab.loc[i, 'sample_name'])
        gene_list = tab.loc[i, 'SYMBOL'].split(",")
        for gene in gene_list:
            if gene != '.':
                if gene not in dict_sample:
                    dict_sample[gene] = sample
                else:
                    if sample not in dict_sample[gene]:
                        val = dict_sample[gene]
                        val = val + ', ' + sample
                        dict_sample[gene] = val
                        
        
    return dict_sample

# 4. Define a function for counting in how many samples the variants in the genes have been seen
def sample_count(tab):
    tab['Samples_Count'] = ''
    for i in range(len(tab)):
        sample_list = tab.loc[i,'Samples'].split(",")
        count = len(sample_list)
        tab.loc[i, 'Samples_Count'] =  count
    return tab


# 5. applying the functions to each table and creating a new dataframe with the following structure:
# SYMBOL COUNT SVTYPE SAMPLE
inv_gene = counter_gene_alt(INV)
inv_sample = tracker_sample(INV)
inv_table = pd.DataFrame(list(inv_gene.items()), columns = ['SYMBOL', 'Count'])
inv_table['SVTYPE'] = 'INV'
inv_table['Samples'] = list(inv_sample.values())
inv_complete_tab = sample_count(inv_table)
    

del_gene = counter_gene(DEL)
del_sample = tracker_sample(DEL)
del_table = pd.DataFrame(list(del_gene.items()), columns = ['SYMBOL', 'Count'])
del_table['SVTYPE'] = 'DEL'
del_table['Samples'] = list(del_sample.values())
del_complete_tab = sample_count(del_table)

dup_gene = counter_gene(DUP)
dup_sample = tracker_sample(DUP)
dup_table = pd.DataFrame(list(dup_gene.items()), columns = ['SYMBOL', 'Count'])
dup_table['SVTYPE'] = 'DUP'
dup_table['Samples'] = list(dup_sample.values())
dup_complete_tab = sample_count(dup_table)

ins_gene = counter_gene(INS)
ins_sample = tracker_sample(INS)
ins_table = pd.DataFrame(list(ins_gene.items()), columns = ['SYMBOL', 'Count'])
ins_table['SVTYPE'] = 'INS'
ins_table['Samples'] = list(ins_sample.values())
ins_complete_tab = sample_count(ins_table)




# 6. annotate the genes that are involved with respiratory infections from hp:0002205
def annotation_respiratory_infections(tab, tab_disease):
    tab['Disease'] = ''
    for i in range(len(tab)):
        gene = tab.loc[i, 'SYMBOL']
        for j in range(len(tab_disease)):
            gene_association = tab_disease.loc[j, 'GENE_SYMBOL']
            if gene == gene_association:
                disease = tab_disease.loc[j, 'DISEASE_NAME']
                if tab.loc[i, 'Disease'] == '':
                    tab.loc[i, 'Disease'] = disease
                else:
                    if disease != tab.loc[i, 'Disease']:
                        if disease not in tab.loc[i, 'Disease'].split(","):
                            tab.loc[i, 'Disease'] = tab.loc[i, 'Disease'] + ',' + disease
    return tab


inv_annotated_res_inf = annotation_respiratory_infections(inv_complete_tab,respiratory_infections)
del_annotated_res_inf = annotation_respiratory_infections(del_complete_tab, respiratory_infections)
dup_annotated_res_inf = annotation_respiratory_infections(dup_complete_tab, respiratory_infections)                
ins_annotated_res_inf = annotation_respiratory_infections(ins_complete_tab,respiratory_infections)


# 7. annotate the genes that are involved with respiratory infections from covid list
def covid_annot(tab, tab_disease):
    for i in range(len(tab)):
        gene = tab.loc[i, 'SYMBOL']
        for j in range(len(tab_disease)):
            gene_association = tab_disease.loc[j, 'Gene Symbol']
            if gene == gene_association:
                if tab.loc[i, 'Disease'] == '':
                    tab.loc[i, 'Disease'] = 'COVID-19 research'
                else:
                    tab.loc[i, 'Disease'] = tab.loc[i, 'Disease'] + ',' + 'COVID-19 research'
    return tab

inv_annot_res_covid = covid_annot(inv_annotated_res_inf, covid)
dup_annot_res_covid = covid_annot(dup_annotated_res_inf, covid)
del_annot_res_covid = covid_annot(del_annotated_res_inf, covid)
ins_annot_res_covid = covid_annot(ins_annotated_res_inf, covid)

# concatenate together the now annotated tables of each SVTYPE 
concatenated_df = pd.concat([inv_annot_res_covid, dup_annot_res_covid, del_annot_res_covid, ins_annot_res_covid])
concatenated_df = concatenated_df.sort_values('SYMBOL')

#export the dataframes
inv_annot_res_covid.to_csv('INV_661_samples_report.csv', index=False, sep='\t')
inv_annot_res_covid.to_excel('INV_661_samples_report.xlsx', index=False)
ins_annot_res_covid.to_csv('INS_661_samples_report.csv', index=False, sep='\t')
ins_annot_res_covid.to_excel('INS_661_samples_report.xlsx', index=False)
dup_annot_res_covid.to_csv('DUP_661_samples_report.csv', index=False, sep='\t')
dup_annot_res_covid.to_excel('DUP_661_samples_report.xlsx', index=False)
del_annot_res_covid.to_csv('DEL_661_samples_report.csv', index=False, sep='\t')
del_annot_res_covid.to_excel('DEL_661_samples_report.xlsx', index=False)

concatenate_df.to_csv('SVs_661_samples_report.csv', index=False, sep='\t')
concatenate_df.to_excel('SVs_661_samples_report.xlsx', index=False)



