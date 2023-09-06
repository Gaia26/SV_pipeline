import pandas as pd
import numpy as np
import sys


report_tab = sys.argv[1] # report file coming from VEP
tab = pd.read_csv(report_tab, delimiter = "\t")
file_patients = sys.argv[2] # excel file containing all info on patients 
info_file = pd.read_excel(file_patients)

hpo = sys.argv[3] # Recurrent_respiratory_infections.csv
covid_tab = sys.argv[4] # panelApp
covid_ncbi_tab = sys.argv[5] # ncbi
respiratory_infections = pd.read_csv(hpo, delimiter = "\t")
covid = pd.read_csv(covid_tab,  delimiter = "\t")
covid_ncbi = pd.read_csv(covid_ncbi_tab, delimiter= "\t")

# select the MAF
for i in range(len(tab)):
    allele_freq = tab.loc[i, 'SV_overlap_AF'].split(",")[0]
    values = allele_freq.split("&")
    float_list = []
    
    for value in values:
        if value == ".":
            float_list.append(float('nan'))
        elif 'e' in value.lower():  # Check if value contains scientific notation
            float_list.append(float(value))
        else:
            float_list.append(float(value))
    
    max_allele_freq = max(float_list)
    tab.loc[i, 'MAF'] = max_allele_freq

tab_columns = tab.columns.tolist()
last_column = tab_columns[-1]
df_without_last = tab.drop(columns=[last_column])
fifth_column = tab_columns[5]
# Create a new DataFrame with the last column inserted at the fifth position
df_new = pd.concat([df_without_last.iloc[:, :5], tab[last_column], df_without_last.iloc[:, 5:]], axis=1)

# filter for MAF either as a point or smaller than the threshold
filtered_tab = df_new[((np.isnan(df_new['MAF'])) |
                      ((~np.isnan(df_new['MAF'])) & (df_new['MAF'] <= 0.001999)))].reset_index(drop=True)

# creating table for each svtype
sample_1 = '001DWB'
index_df = inv_df.columns.get_loc(sample_1) # get  the index of the column 001DWB from which the following reshape start

inv_df = filtered_tab[filtered_tab['SVTYPE'] == 'INV']
# from column index_df on there are the samples with their genotypes and that part of the dataframe must be transposed.
INV_reshape = pd.melt(inv_df, id_vars=inv_df.columns[:index_df], value_vars=inv_df.columns[index_df:], var_name='sample_name', value_name='genotype')
INV = INV_reshape.drop(INV_reshape[INV_reshape['genotype'] == ' ./.'].index)
column_names_inv = INV.columns.tolist()
new_column_order_inv = column_names_inv[-2:] + column_names_inv[:-2]
INV = INV[new_column_order_inv].reset_index(drop=True)

del_df = filtered_tab[filtered_tab['SVTYPE'] == 'DEL']
DEL_reshape = pd.melt(del_df, id_vars=del_df.columns[:index_df], value_vars=del_df.columns[index_df:], var_name='sample_name', value_name='genotype')
DEL = DEL_reshape.drop(DEL_reshape[DEL_reshape['genotype'] == ' ./.'].index)
column_names_del = DEL.columns.tolist()
new_column_order_del = column_names_del[-2:] + column_names_del[:-2]
DEL = DEL[new_column_order_del].reset_index(drop=True)


ins_df = filtered_tab[filtered_tab['SVTYPE'] == 'INS']
INS_reshape = pd.melt(ins_df, id_vars=ins_df.columns[:index_df], value_vars=ins_df.columns[index_df:], var_name='sample_name', value_name='genotype')
INS = INS_reshape.drop(INS_reshape[INS_reshape['genotype'] == ' ./.'].index)
column_names_ins = INS.columns.tolist()
new_column_order_ins = column_names_ins[-2:] + column_names_ins[:-2]
INS = INS[new_column_order_ins].reset_index(drop=True)


dup_df = filtered_tab[filtered_tab['SVTYPE'] == 'DUP']
DUP_reshape = pd.melt(dup_df, id_vars=dup_df.columns[:index_df], value_vars=dup_df.columns[index_df:], var_name='sample_name', value_name='genotype')
DUP = DUP_reshape.drop(DUP_reshape[DUP_reshape['genotype'] == ' ./.'].index)
column_names_dup = DUP.columns.tolist()
new_column_order_dup = column_names_dup[-2:] + column_names_dup[:-2]
DUP = DUP[new_column_order_dup].reset_index(drop=True)


# Remove Bellaria
inv_wgs_wps_df = INV[~INV['sample_name'].str.startswith(('WGS', 'WP2'))].reset_index(drop=True)
ins_wgs_wps_df = INS[~INS['sample_name'].str.startswith(('WGS', 'WP2'))].reset_index(drop=True)
del_wgs_wps_df = DEL[~DEL['sample_name'].str.startswith(('WGS', 'WP2'))].reset_index(drop=True)
dup_wgs_wps_df = DUP[~DUP['sample_name'].str.startswith(('WGS', 'WP2'))].reset_index(drop=True)

# remove variants that have svlen between 50 and 100
filtered_inv = inv_wgs_wps_df[~inv_wgs_wps_df['SVLEN'].between(50, 100)].reset_index(drop=True)
filtered_ins = ins_wgs_wps_df[~ins_wgs_wps_df['SVLEN'].between(50, 100)].reset_index(drop=True)
filtered_dup = dup_wgs_wps_df[~dup_wgs_wps_df['SVLEN'].between(50, 100)].reset_index(drop=True)
filtered_del = del_wgs_wps_df[~del_wgs_wps_df['SVLEN'].between(-100, -50)].reset_index(drop=True)


# add the WHOCPS from the excel file containing all info from the patients
filtered_inv['WHOCPS'] = filtered_inv['sample_name'].map(info_file.set_index('Sample Name')['WHOCPS'])
filtered_del['WHOCPS'] = filtered_del['sample_name'].map(info_file.set_index('Sample Name')['WHOCPS'])
filtered_dup['WHOCPS'] = filtered_dup['sample_name'].map(info_file.set_index('Sample Name')['WHOCPS'])
filtered_ins['WHOCPS'] = filtered_ins['sample_name'].map(info_file.set_index('Sample Name')['WHOCPS'])

# exporting the dataframes to use them for further analysis on R
filtered_inv.to_csv('INV_vep.csv', index=False, sep='\t')
filtered_ins.to_csv('INS_vep.csv', index=False, sep='\t')
filtered_del.to_csv('DEL_vep.csv', index=False, sep='\t')
filtered_dup.to_csv('DUP_vep.csv', index=False, sep='\t')



#function counter_gene to count how many different SV are found in each gene
def counter_gene(tab):
    dict_gene = {}
    seen_elements = []
    for i in range(len(tab)):
        gene_list = tab.loc[i, 'SYMBOL'].split(",")
        chrom = tab.loc[i,'CHROM']
        pos = tab.loc[i,'POS']
        end = tab.loc[i,'END']
        if (chrom, pos, end, gene_list) not in seen_elements:
            for gene in gene_list:
                if gene not in dict_gene:
                    dict_gene[gene] = 0
                dict_gene[gene] += 1
        seen_elements.append((chrom, pos, end, gene_list))
    return dict_gene






def tracker_sample(tab):
    dict_sample = {}
    dict_whocps = {}
    for i in range(len(tab)):
        sample = (tab.loc[i, 'sample_name'])
        whocps = tab.loc[i, 'WHOCPS']
        gene_list = tab.loc[i, 'SYMBOL'].split(",")
        for gene in gene_list:
            if gene != '.':
                if gene not in dict_sample:
                    dict_sample[gene] = sample
                    dict_whocps[gene] = whocps
                else:
                    if sample not in dict_sample[gene]:
                        val = dict_sample[gene]
                        val = val + ', ' + sample
                        dict_sample[gene] = val
                        val2 = dict_whocps[gene]
                        val2 = val2 + ', ' + whocps
                        dict_whocps[gene] = val2
    return dict_sample, dict_whocps










# def tracker_sample(tab):
#     dict_sample = {}
#     for i in range(len(tab)):
#         sample = (tab.loc[i, 'sample_name'])
#         gene_list = tab.loc[i, 'SYMBOL'].split(",")
#         for gene in gene_list:
#             if gene not in dict_sample: 
#                 dict_sample[gene] = sample
#             else:
#                 if sample not in dict_sample[gene]:
#                     val = dict_sample[gene]
#                     val = val + ', ' + sample
#                     dict_sample[gene] = val
#     return dict_sample        

           


# count how many samples are present in the newly formed Sample column
def sample_count(tab):
    tab['Samples_Count'] = ''
    for i in range(len(tab)):
        sample_list = tab.loc[i,'Samples'].split(",")
        count = len(sample_list)
        tab.loc[i, 'Samples_Count'] =  count
    return tab


inv_gene = counter_gene(filtered_inv)
inv_sample, inv_whocps = tracker_sample(filtered_inv)
inv_table = pd.DataFrame(list(inv_gene.items()), columns = ['SYMBOL', 'Count', 'WHOCPS'])
inv_table['SVTYPE'] = 'INV'
inv_table['Samples'] = list(inv_sample.values())
inv_table['WHOCPS'] = list(inv_whocps.values())
inv_complete_tab = sample_count(inv_table)


del_gene = counter_gene(filtered_del)
del_sample = tracker_sample(filtered_del)
del_table = pd.DataFrame(list(del_gene.items()), columns = ['SYMBOL', 'Count'])
del_table['SVTYPE'] = 'DEL'
del_table['Samples'] = list(del_sample.values())
del_complete_tab = sample_count(del_table)


dup_gene = counter_gene(filtered_dup)
dup_sample = tracker_sample(filtered_dup)
dup_table = pd.DataFrame(list(dup_gene.items()), columns = ['SYMBOL', 'Count'])
dup_table['SVTYPE'] = 'DUP'
dup_table['Samples'] = list(dup_sample.values())
dup_complete_tab = sample_count(dup_table)


ins_gene = counter_gene(filtered_ins)
ins_sample = tracker_sample(filtered_ins)
ins_table = pd.DataFrame(list(ins_gene.items()), columns = ['SYMBOL', 'Count'])
ins_table['SVTYPE'] = 'INS'
ins_table['Samples'] = list(ins_sample.values())
ins_complete_tab = sample_count(ins_table)

# check if the genes present are found in the list of genes taken from hpo for respiratory infections
word_to_add = 'Chronic Respiratory Infections'
inv_complete_tab['Disease'] = ''
ins_complete_tab['Disease'] = ''
dup_complete_tab['Disease'] = ''
del_complete_tab['Disease'] = ''
inv_complete_tab.loc[((inv_complete_tab['SYMBOL'].isin(respiratory_infections['GENE_SYMBOL']))) & (~inv_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
ins_complete_tab.loc[((ins_complete_tab['SYMBOL'].isin(respiratory_infections['GENE_SYMBOL']))) & (~ins_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
dup_complete_tab.loc[((dup_complete_tab['SYMBOL'].isin(respiratory_infections['GENE_SYMBOL']))) & (~dup_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
del_complete_tab.loc[((del_complete_tab['SYMBOL'].isin(respiratory_infections['GENE_SYMBOL']))) & (~del_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add

# add "COVID-19 research if the genes are involved with it from info taken from panelApp and NCBI"

word_to_add = 'COVID-19 research'
inv_complete_tab.loc[((inv_complete_tab['SYMBOL'].isin(covid_ncbi['Symbol'])) | (inv_complete_tab['SYMBOL'].isin(covid_ncbi['Aliases']))) & (~inv_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
ins_complete_tab.loc[((ins_complete_tab['SYMBOL'].isin(covid_ncbi['Symbol'])) | (ins_complete_tab['SYMBOL'].isin(covid_ncbi['Aliases']))) & (~ins_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
dup_complete_tab.loc[((dup_complete_tab['SYMBOL'].isin(covid_ncbi['Symbol'])) | (dup_complete_tab['SYMBOL'].isin(covid_ncbi['Aliases']))) & (~dup_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
del_complete_tab.loc[((del_complete_tab['SYMBOL'].isin(covid_ncbi['Symbol'])) | (del_complete_tab['SYMBOL'].isin(covid_ncbi['Aliases']))) & (~del_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add

inv_complete_tab.loc[((inv_complete_tab['SYMBOL'].isin(covid['Gene Symbol']))) & (~inv_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
ins_complete_tab.loc[((ins_complete_tab['SYMBOL'].isin(covid['Gene Symbol']))) & (~ins_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
dup_complete_tab.loc[((dup_complete_tab['SYMBOL'].isin(covid['Gene Symbol']))) & (~dup_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add
del_complete_tab.loc[((del_complete_tab['SYMBOL'].isin(covid['Gene Symbol']))) & (~del_complete_tab['Disease'].str.contains(word_to_add, na=False)), 'Disease'] += ' ' + word_to_add

# export the tables 
inv_complete_tab.to_csv('INV_COVID_samples_report.csv', index=False, sep='\t')
inv_complete_tab.to_excel('INV_COVID_samples_report.xlsx', index=False)
ins_complete_tab.to_csv('INS_COVID_samples_report.csv', index=False, sep='\t')
ins_complete_tab.to_excel('INS_COVID_samples_report.xlsx', index=False)
dup_complete_tab.to_csv('DUP_COVID_samples_report.csv', index=False, sep='\t')
dup_complete_tab.to_excel('DUP_COVID_samples_report.xlsx', index=False)
del_complete_tab.to_csv('DEL_COVID_samples_report.csv', index=False, sep='\t')
del_complete_tab.to_excel('DEL_COVID_samples_report.xlsx', index=False)



