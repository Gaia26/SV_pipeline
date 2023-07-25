






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

# function tracker_sample is to check in which patients the variants are found
def tracker_sample(tab):
    dict_sample = {}
    for i in range(len(tab)):
        sample = (tab.loc[i, 'sample_name'])
        gene_list = tab.loc[i, 'SYMBOL'].split(",")
        for gene in gene_list:
            if gene not in dict_sample: 
                dict_sample[gene] = sample
            else:
                if sample not in dict_sample[gene]:
                    val = dict_sample[gene]
                    val = val + ', ' + sample
                    dict_sample[gene] = val
    return dict_sample

# count how many samples are present in the newly formed Sample column
def sample_count(tab):
    tab['Samples_Count'] = ''
    for i in range(len(tab)):
        sample_list = tab.loc[i,'Samples'].split(",")
        count = len(sample_list)
        tab.loc[i, 'Samples_Count'] =  count
    return tab


inv_gene = counter_gene(filtered_inv)
inv_sample = tracker_sample(filtered_inv)
inv_table = pd.DataFrame(list(inv_gene.items()), columns = ['SYMBOL', 'Count'])
inv_table['SVTYPE'] = 'INV'
inv_table['Samples'] = list(inv_sample.values())
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



