from bokeh.plotting import figure, output_file, save, ColumnDataSource
from bokeh.models import BoxAnnotation, Range1d, CustomJS, LabelSet, ColumnDataSource, NumeralTickFormatter, RadioButtonGroup, CustomJSFilter, CDSView, Button
from bokeh.layouts import gridplot, row, column, widgetbox
from bokeh.palettes import YlOrRd
from bokeh.io import output_notebook
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn, Slider, Spinner, TextInput, CheckboxGroup
from bokeh.models.tools import HoverTool
from bokeh.embed import json_item

from more_itertools import sort_together

import numpy as np
import pandas as pd
# import scipy.stats as st
import os
import json
import math
import datetime
import time

# Turns a set of values into a confidence interval
# def points_to_ci(values_list:list, alpha=0.95):
#     lower_ci,upper_ci = st.t.interval(alpha=alpha, df=len(values_list)-1, loc=np.mean(values_list), scale=st.sem(values_list))
#     return lower_ci, upper_ci

def slider_spinner(start: float=0.0, end: float=1.0,value: float=0.5, step: float=0.01, title: str="Slider", name:str=None, width=300):
    spinner = Spinner(title=" ", low = start,high = end, step = step, value = value, width = 60, name = name + '_spinner', margin=(2,0,0,0),height=32)
    slider = Slider(start=start, end=end, value=value, step=step, title=title, name=name, width=max(120,width-60), margin=(5,10,5,10))
    spinner.js_link('value', slider, 'value')
    slider.js_link('value',spinner, 'value')
    return slider, spinner

# Filter excessive CNV points #

def filter_cnv_points(cnvloh_dna_dfs,num_stdev=0.5,hard_filtering:bool=False):
    cnv_x_array = []
    cnv_y_array = []
    cnv_filter_array = []
    
    bin_sd_dict = {}
    bin_mean_dict = {}
    for i in range(0,len(cnvloh_dna_dfs['cnv_bins_df']['BIN'])):
        bin_id = f"{cnvloh_dna_dfs['cnv_bins_df']['CHROM'][i]}:{int(cnvloh_dna_dfs['cnv_bins_df']['BIN'][i])}"
        bin_sd_dict[bin_id] = cnvloh_dna_dfs['cnv_bins_df']['STDEV'][i]
        bin_mean_dict[bin_id] = cnvloh_dna_dfs['cnv_bins_df']['MEAN'][i]
    #print(bin_sd_dict)
    for i in range(0,len(cnvloh_dna_dfs['cnv_df']['LOG2'])):
        point_bin = f"{cnvloh_dna_dfs['cnv_df']['CHROM'][i]}:{int(cnvloh_dna_dfs['cnv_df']['BIN'][i])}"
        point_x = cnvloh_dna_dfs['cnv_df']['MIDPOINT'][i]
        point_y = cnvloh_dna_dfs['cnv_df']['LOG2'][i]
        if point_bin in bin_mean_dict:
            cutoff_max = (bin_mean_dict[point_bin] + abs(bin_sd_dict[point_bin])*num_stdev)-0.05
            cutoff_min = (bin_mean_dict[point_bin] - abs(bin_sd_dict[point_bin])*num_stdev)+0.05
            if point_y >= cutoff_max or point_y <= cutoff_min:
                filter_status = 0
            else:
                filter_status = 1
            if filter_status == 0 or not hard_filtering:
                cnv_x_array.append(point_x)
                cnv_y_array.append(point_y)
                cnv_filter_array.append(filter_status)
        else:
            cnv_x_array.append(point_x)
            cnv_y_array.append(point_y)
            cnv_filter_array.append(0)
    return cnv_x_array, cnv_y_array, cnv_filter_array

# CNV gain/loss color for helper plot

cnv_col_list = ["blue","dodgerblue","lightblue","lightgray","palevioletred","firebrick","red"]
loh_col_list = YlOrRd[7]

def calc_color_cnv_list(cnv_values:list, color_list:list=cnv_col_list, cnv_cutoff:float=0.25, gene_list:list = None, gene_id:list = None):
    if gene_list is None:
        color_choice_list = ["transparent"]*len(cnv_values)
    else:
        if gene_list is str:
            gene_list = set(gene_list.split(","))
        else:
            gene_list = set(gene_list)
        
        color_choice_list = []
        
        for i in range(0,len(cnv_values)):
            if gene_id is not None and gene_id[i] in gene_list:
                color_choice_list.append(cnv_values[i], color_list, cnv_cutoff)
            else:
                color_choice_list.append("transparent")
    return color_choice_list

def calc_col_cnv(cnv_value:float,color_list:list=cnv_col_list,cnv_cutoff:float=0.25):
    # CNV color for helper plot
    if not cnv_cutoff:
        cnv_cutoff = 0
    try:
        if abs(cnv_value) >= cnv_cutoff:
            if cnv_value > 1.58:
                color_choice = 0
            elif cnv_value > 0.85:
                color_choice = 1
            elif cnv_value > cnv_cutoff:
                color_choice = 2
            elif cnv_value < -2:
                color_choice = 6
            elif cnv_value < -0.85:
                color_choice = 5
            elif cnv_value < -1*cnv_cutoff:
                color_choice = 4
            else:
                color_choice = 3
            color_choice = color_list[color_choice]
        else:
            color_choice = "transparent"
    except ValueError:
        color_choice = "transparent"
    return color_choice

# LOH color for LOH plot

def calc_col_loh(loh_value:float,color_list:list=loh_col_list):
    try:
        if loh_value < 0.48:
            loh_calc = (0.5-loh_value)*2
            index_value = min(max(0,round((1-loh_calc)*(len(color_list)))),len(color_list)-1)
            color_choice = color_list[index_value]
        else:
            color_choice="transparent"
    except ValueError:
        color_choice = "transparent"
    return color_choice

# Takes a genes table file and subsets it down to just the genes in the provided genes sublist. If there is not a genes sublist, just returns the same table. ##
def subset_genes_table(genes_table, gene_sublist):
    try:
        if gene_sublist is not None: # Blank subset list = show no genes
            if type(gene_sublist) is list and len(gene_sublist) == 0:
                genes_table_subset = None
            else: # Non-blank subset list = show listed genes
                genes_table_subset = genes_table[[str(a) in gene_sublist for a in genes_table['ENTREZ_ID'].tolist()]]
        else: # No subset list = show all genes
            genes_table_subset = genes_table
    except Exception as e: # Fail and show no genes
        print(str(e))
        genes_table_subset = None
    return genes_table_subset

### Create data source objects ###

def create_sourceseg_datasource(cnvloh_dna_dfs, source_name, num_stdev:float=0.5, abs_adj_list:list=None, fill_alpha:float=1.0, genes_table_dict:dict=None,genome_info:dict=None,cnv_cutoff=0.25, loh_cutoff = 0.42, baseline_adj_list:list=None):
    label_pos_adjustment_top = 0
    label_pos_adjustment_bottom = 0

    seg_index_local = []
    loh_status_list = []
    cnvloh_status_list = []
    
    for i in range(0,len(list(cnvloh_dna_dfs['seg_df']['CHROM']))):
        local_index = 1
        for j in list(cnvloh_dna_dfs['seg_df'][cnvloh_dna_dfs['seg_df']['CHROM'] == cnvloh_dna_dfs['seg_df']['CHROM'][i]]['START']):
            if j < list(cnvloh_dna_dfs['seg_df']['START'])[i]:
                local_index += 1
        seg_index_local.append(local_index)
        cnvloh_status = cnvloh_to_string(cnvloh_dna_dfs['seg_df']['LOG2'][i], cnvloh_dna_dfs['seg_df']['LOH'][i], cnv_cutoff, loh_cutoff)
        if cnvloh_dna_dfs['seg_df']['LOH'][i] <= loh_cutoff:
            loh_status = 'YES'
        else:
            loh_status = 'NO'
        loh_status_list.append(loh_status)
        cnvloh_status_list.append(cnvloh_status)
        
    if baseline_adj_list is None:
        baseline_adj_list = [0]*len(list(cnvloh_dna_dfs['seg_df']['CHROM']))
    
    source_seg_data = {
        'seg_index':list(cnvloh_dna_dfs['seg_df']['INDEX']),
        'seg_index_local':seg_index_local,
        'seg_per_chrom':np.array([list(cnvloh_dna_dfs['seg_df']['CHROM']).count(i) for i in list(cnvloh_dna_dfs['seg_df']['CHROM'])]),
        'seg_chrom':list(cnvloh_dna_dfs['seg_df']['CHROM']),
        'seg_x0':np.array(cnvloh_dna_dfs['seg_df']['START']),
        'seg_x1':np.array(cnvloh_dna_dfs['seg_df']['END']),
        'seg_length':np.array(cnvloh_dna_dfs['seg_df']['END']) - np.array(cnvloh_dna_dfs['seg_df']['START']),
        'cnv_y_stored':np.array(cnvloh_dna_dfs['seg_df']['LOG2']),
        'cnv_y':np.array(cnvloh_dna_dfs['seg_df']['LOG2']),
        'cnv_points':np.array(cnvloh_dna_dfs['seg_df']['LOG2_MARKS']),
        'cnv_cov_depth': np.array(cnvloh_dna_dfs['seg_df']['LOG2_READS']) * 150 / (np.array(cnvloh_dna_dfs['seg_df']['END']) - np.array(cnvloh_dna_dfs['seg_df']['START'])),
        'cnv_reads':np.array(cnvloh_dna_dfs['seg_df']['LOG2_READS']),
        'cnv_y_height':np.array(cnvloh_dna_dfs['seg_df']['LOG2_STDEV']).astype(float)*num_stdev*2,
        'cnv_y_height_stored':np.array(cnvloh_dna_dfs['seg_df']['LOG2_STDEV']).astype(float)*num_stdev*2,
        'seg_loh_x':(np.array(cnvloh_dna_dfs['seg_df']['START'])+np.array(cnvloh_dna_dfs['seg_df']['END']))/2,
        'seg_loh_width':(np.array(cnvloh_dna_dfs['seg_df']['END'])-np.array(cnvloh_dna_dfs['seg_df']['START'])),
        'seg_loh_y':np.array([0.5]*len(cnvloh_dna_dfs['seg_df']['START'])),
        'seg_loh_baf':np.array(cnvloh_dna_dfs['seg_df']['LOH']),
        'seg_loh_height':np.array(1-cnvloh_dna_dfs['seg_df']['LOH']*2),
        'seg_loh_status':np.array(loh_status_list),
        'seg_loh_height_stored':np.array(1-cnvloh_dna_dfs['seg_df']['LOH']*2),
        'seg_loh_color':np.array([calc_col_loh(i) for i in cnvloh_dna_dfs['seg_df']['LOH']]),
        'seg_loh_line':np.array(["grey"]*len(cnvloh_dna_dfs['seg_df']['LOH'])),
        'seg_alpha':np.array([fill_alpha]*len(cnvloh_dna_dfs['seg_df']['LOH'])),
        'cnv_stdev': np.array(cnvloh_dna_dfs['seg_df']['LOG2_STDEV']),
        'cnv_num_copy':( 2 ** (np.array(cnvloh_dna_dfs['seg_df']['LOG2'])+1)) / (1+np.array(baseline_adj_list)),
        'cnvloh_status': (np.array(cnvloh_status_list)),
        'cnv_y_display':np.array(cnvloh_dna_dfs['seg_df']['LOG2']),
        'cnv_label_alpha':(np.array([0 for x in cnvloh_status_list])),
        'cnv_label_pos':(np.array([-1000.0 for x in cnvloh_status_list])),
        'cnv_label_text':np.array(cnvloh_dna_dfs['seg_df']['LOG2']),
        'seg_baseline_adj':np.array(baseline_adj_list)
        }
    for x in range(0,len(source_seg_data['cnv_y_display'])):
        if source_seg_data['cnv_y_display'][x] > 1.97 or source_seg_data['cnv_y_display'][x] < -1.97:
            source_seg_data['cnv_y_display'][x] = max(-1.97,min(1.97,source_seg_data['cnv_y_display'][x]))
            if source_seg_data['cnv_y'][x] > 0:
                source_seg_data['cnv_label_pos'][x] = 1.7 - label_pos_adjustment_top
                if label_pos_adjustment_top == 0.0:
                    label_pos_adjustment_top = 0.2
                else:
                    label_pos_adjustment_top = 0.0
            else:
                source_seg_data['cnv_label_pos'][x] = -2.0 + label_pos_adjustment_bottom
                if label_pos_adjustment_bottom == 0.0:
                    label_pos_adjustment_bottom = 0.2
                else:
                    label_pos_adjustment_bottom = 0.0
            source_seg_data['cnv_label_text'][x] = round(source_seg_data['cnv_label_text'][x],1)
            source_seg_data['cnv_label_alpha'][x] = 1
            source_seg_data['cnv_num_copy'][x] = source_seg_data['cnv_num_copy'][x] / (1+source_seg_data['seg_baseline_adj'][x])
    print(source_seg_data)
            
    if abs_adj_list:
        source_seg_data['seg_x0_absolute']=np.array(cnvloh_dna_dfs['seg_df']['START'])+np.array(abs_adj_list)
        source_seg_data['seg_x1_absolute']=np.array(cnvloh_dna_dfs['seg_df']['END'])+np.array(abs_adj_list)
        source_seg_data['seg_loh_x_absolute']=(np.array(cnvloh_dna_dfs['seg_df']['START'])+np.array(cnvloh_dna_dfs['seg_df']['END']))/2+np.array(abs_adj_list)
    if genes_table_dict:
        #print('genes_table_dict found in create_sourceseg_datasource')
        genes_dict_list = []
        seg_ids = list(cnvloh_dna_dfs['seg_df']['ID'])
        #print('seg_ids:')
        #print(seg_ids)
        for i in seg_ids:
            if i in genes_table_dict:
                genes_dict_list.append(genes_table_dict[i])
            else:
                genes_dict_list.append({})
        source_seg_data['genes_dict']=genes_dict_list
    if genome_info:
        all_cytobands = []
        for i in range(0,len(source_seg_data['seg_x0'])):
            seg_chrom = cnvloh_dna_dfs['seg_df']['CHROM'][i]
            seg_start = cnvloh_dna_dfs['seg_df']['START'][i]
            seg_end = cnvloh_dna_dfs['seg_df']['END'][i]
            if seg_chrom in genome_info:
                chrom_info = genome_info[seg_chrom]
            else:
                chrom_info = genome_info
            all_cytobands.append(get_cytobands(seg_chrom,seg_start,seg_end,chrom_info))
        source_seg_data['cytobands']=np.array(all_cytobands, dtype=object)
    for i in source_seg_data:
        if type(source_seg_data[i]) is np.ndarray:
            source_seg_data[i] = list(source_seg_data[i])
    source_seg = ColumnDataSource(data=source_seg_data, name=source_name)
    return source_seg

### Create gene table objects ###

def filter_nan_values(d):
    return {k: v for (k, v) in d.items() if v is not np.NaN}

def table_to_dict(table, key_1, key_2, exclude_columns:list=None, remove_blanks:bool=True, string_keys_only:bool=False):
    start_time = time.time()
    table_dict = {}
    if exclude_columns:
        table = table.drop(columns=exclude_columns)
    if string_keys_only:
        table[key_1] = table[key_1].astype(str)
        table[key_2] = table[key_2].astype(str)
    for i in list(set(table[key_1])):
        sub_table = table[table[key_1] == i].copy()
        sub_table.set_index(key_2, drop=True, inplace=True)
        sub_table = sub_table[~sub_table.index.astype(str).duplicated(keep='first')]
        table_dict[i] = sub_table.to_dict('index')
        if remove_blanks:
            for j in list(table_dict[i].keys()):
                for k in list(table_dict[i][j].keys()):
                    if table_dict[i][j][k] is None or table_dict[i][j][k] is np.nan or table_dict[i][j][k] == "" or str(table_dict[i][j][k]) == 'nan':
                        table_dict[i][j].pop(k)
    elapsed_time = time.time() - start_time
    print(f'\tTable to Dict Runtime: {elapsed_time:.2f} seconds.')
    return table_dict

def add_phenotype_gene_status(genes_info_dict:dict, phenotype_genes_list:list):
    phenotype_genes_list = [str(x) for x in phenotype_genes_list]
    # Updates dict in place
    for i in genes_info_dict:
        for j in genes_info_dict[i]:
            if str(j) in phenotype_genes_list:
                genes_info_dict[i][j]['phenotype_gene']=True
            else:
                genes_info_dict[i][j]['phenotype_gene']=False

# Takes number of copies and LOH status and converts to a string
def cnvloh_to_string(log2:float=0, baf:float=0.5, cnv_cutoff:float=0.25, loh_cutoff:float=0.42):
    if abs(log2) >= cnv_cutoff:
        num_copies = 2 ** (1+log2)
    else:
        num_copies = 2
    
    if baf < loh_cutoff:
        is_loh = True
    else:
        is_loh = False
    
    if num_copies >= 5.0:
        copy_status = 'Amplification'
        if not is_loh:
            copy_status = copy_status + ' (Balanced)'
    elif num_copies > 2.6:
        copy_status = 'Gain'
        if num_copies > 3.5:
            if is_loh:
                copy_status = copy_status + ' (Unbalanced)'
            else:
                copy_status = copy_status + ' (Balanced)'
        else:
            if not is_loh:
                copy_status = copy_status + ' (!Anomalous LOH)'
    elif num_copies > 2.2:
        copy_status = 'Weak Gain'
        if is_loh:
            copy_status = copy_status + ' (LOH Supported)'
    elif num_copies < 0.5:
        copy_status = 'Biallelic Loss'
    elif num_copies < 1.4:
        copy_status = 'Loss'
        if not is_loh:
            copy_status = copy_status + ' (!Anomalous LOH)'
    elif num_copies < 1.8:
        copy_status = 'Weak Loss'
        if is_loh:
            copy_status = copy_status + ' (LOH Supported)'
    else:
        copy_status = 'Copy-Neutral'
        if is_loh:
            copy_status = copy_status + ' LOH'
    return copy_status

def get_cytobands(chromosome:str, seg_start:int, seg_end:int, chrom_info_dict:dict):
    bands_list = []
    for k in chrom_info_dict['cytobands']:
            band_name = k['name']
            cyto_start = k['start']
            cyto_end = k['end']
            if not ((cyto_start < seg_start and cyto_end < seg_start) or (cyto_start > seg_end and cyto_end > seg_end)):
                bands_list.append(band_name)
    return bands_list


def make_gene_table_data(genes_table, cnv_cutoff:float=0, loh_cutoff:float=0.5, genes_list:list=None, focals_dict:dict=None, sd_cutoff:float=2.0, focal_effects_datasource=None):

    print("######### Genes Table ##########")
    print(genes_table)
    print("######### Genes Table End ##########")

    print("######### Genes List ##########")
    print(genes_list)
    print("######### Genes List End ##########")
    
    count_top = 0
    count_middle = 0
    count_bottom = 0
    
    # if (genes_table) is not None and not genes_table.empty and not (type(genes_list) is list and len(genes_list) == 0):
    if (genes_table) is not None and not genes_table.empty:
        genes_table.sort_values(by=['START'])

        # Default UI state is to show labels for phenotype genes only.  Here we check that
        # either we have not been passed a genes_list or that the Entrez Id is in genes_list
        # before deciding whether the alpha is 1 (shown) or 0 (hidden)
        #if genes_list is None:
        label_alpha=[1 for i in list(genes_table['loh_allele_fraction'])]
        #else:
        #    # Get the list of entrez ids in our table
        #    entrez_ids = list(genes_table['ENTREZ_ID'])
        #    # Instantiate a list of zeroes to the same length
        #    label_alpha=[0] * len(entrez_ids)
        #    # Check each entrez id, and if it's in genes_list, flip the value to 1
        #    for i, entrez_id in enumerate(entrez_ids):
        #        if str(entrez_id) in genes_list:
        #            label_alpha[i] = 1

        gene_table_data = dict(
            SYMBOL=list(genes_table['SYMBOL']),
            ENTREZ_ID=list(genes_table['ENTREZ_ID']),
            CHR=list(genes_table['CHR']),
            START=list(genes_table['START']),
            END=list(genes_table['END']),
            seg_mean=list(genes_table['seg_mean']),
            # gene_log2 begins as the segment log2, but later we will check for focal events
            gene_log2=list(genes_table['seg_mean']),
            gene_x=list((np.array(genes_table['START'])+np.array(genes_table['END']))/2),
            gene_width=list((np.array(genes_table['END'])-np.array(genes_table['START']))),
            gene_y=list(np.array([0]*len(genes_table['START']))),
            gene_y_label=list(np.array([-1]*len(genes_table['START']))),
            gene_height=list(np.array([1]*len(genes_table['START']))),
            gene_cnv_color=list(np.array(calc_color_cnv_list(list(genes_table['seg_mean']),gene_list=genes_list, gene_id=list(genes_table['ENTREZ_ID'])))),
            seg_baf=list(genes_table['loh_allele_fraction']),
            seg_loh=[0 for i in list(genes_table['loh_allele_fraction'])],
            seg_loh_status=['' for i in list(genes_table['loh_allele_fraction'])],
            gene_loh_color=['transparent' for i in list(genes_table['loh_allele_fraction'])],
            label_alpha=label_alpha
        )

        if focals_dict:
            if not sd_cutoff:
                sd_cutoff = 1
            entrez_id_dict={}
            
            for i in range(0,len(focals_dict)):
                log2 = focals_dict[i]['weighted_average_log2_copy_ratio']
                delta_sd = focals_dict[i]['event_stdevs_from_segment']
                if focals_dict[i]['focal_event_type'] == 'exon':
                    focal_type = 'exon'
                    entrez_id = focals_dict[i]['events'][0]['event_info']['entrez_id']
                else:
                    focal_type = 'gene'
                    entrez_id = focals_dict[i]['events'][0]['event_id']
            
                if entrez_id in entrez_id_dict:
                    if focal_type in entrez_id_dict[entrez_id]:
                        if abs(delta_sd) > abs(entrez_id_dict[entrez_id][focal_type]['delta_sd']):
                            entrez_id_dict[entrez_id][focal_type]={'log2':log2, 'delta_sd':delta_sd}
                elif genes_list is not None and entrez_id in genes_list:
                    entrez_id_dict[entrez_id]={focal_type:{'log2':log2, 'delta_sd':delta_sd}}
                    
            
            gene_table_data['focal_log2_gene'] = []
            gene_table_data['focal_sd_gene']= []
            gene_table_data['focal_sd_exon'] = []
            gene_table_data['focal_log2_exon'] = []
            gene_table_data['focal_log2'] = []
            gene_table_data['focal_status'] = []

            for j in range(0,len(gene_table_data['ENTREZ_ID'])):
                entrez_id = str(gene_table_data['ENTREZ_ID'][j])
                
                focal_status_list = []
                focal_log2 = 0
                
                if entrez_id in entrez_id_dict and 'gene' in entrez_id_dict[entrez_id]:
                    focal_log2_gene = entrez_id_dict[entrez_id]['gene']['log2']
                    focal_sd_gene = entrez_id_dict[entrez_id]['gene']['delta_sd']
                    if abs(focal_sd_gene) >= sd_cutoff:
                        focal_status_list.append('gene')
                        focal_log2 = focal_log2_gene
                else:
                    focal_log2_gene = 0
                    focal_sd_gene = 0
                    
                if entrez_id in entrez_id_dict and 'exon' in entrez_id_dict[entrez_id]:
                    focal_log2_exon = entrez_id_dict[entrez_id]['exon']['log2']
                    focal_sd_exon = entrez_id_dict[entrez_id]['exon']['delta_sd']
                    if abs(focal_sd_exon) >= sd_cutoff:
                        focal_status_list.append('exon')
                        if abs(focal_sd_exon) > (focal_sd_gene):
                            focal_log2 = focal_log2_exon
                else:
                    focal_log2_exon = 0
                    focal_sd_exon = 0
                
                if len(focal_status_list) > 0:
                    focal_status = 'Focal: ' + "+".join(focal_status_list)
                    # if abs(focal_log2) > abs(gene_table_data['gene_log2'][j]):
                    gene_table_data['gene_log2'][j] = focal_log2
                else:
                    focal_status = ''
                
                gene_table_data['focal_log2_gene'].append(focal_log2_gene)
                gene_table_data['focal_sd_gene'].append(focal_sd_gene)
                gene_table_data['focal_sd_exon'].append(focal_log2_exon)
                gene_table_data['focal_log2_exon'].append(focal_log2_gene)
                gene_table_data['focal_log2'].append(focal_log2)
                gene_table_data['focal_status'].append(focal_status)
                
                if genes_list is None:
                    gene_table_data['gene_cnv_color'][j] = "transparent"
                else:
                    gene_table_data['gene_cnv_color'][j] = calc_col_cnv(gene_table_data['gene_log2'][j])
                
                # Focal reversals get a special color
                if focal_effects_datasource is not None:
                    #print('focal_effects_datasource.data[entrez_id]: ' + entrez_id)
                    
                    if entrez_id in focal_effects_datasource.data['dict'][0]:
                        print(focal_effects_datasource.data['dict'][0][entrez_id])

                    if entrez_id in focal_effects_datasource.data['dict'][0] and 'Revers' in focal_effects_datasource.data['dict'][0][entrez_id]['status'] and gene_table_data['gene_cnv_color'][j] != "transparent":
                        gene_table_data['gene_cnv_color'][j] = 'mediumorchid'
                
                #print(f'{gene_table_data["ENTREZ_ID"][j]} {gene_table_data["SYMBOL"][j]}: {gene_table_data["focal_log2"][j]} {gene_table_data["gene_cnv_color"][j]}')
        else:
            print('No focals dict')

            
        if loh_cutoff:
            for i in range(0,len(gene_table_data['seg_baf'])):
                gene_table_data['seg_loh'][i] = 2*(0.5 - gene_table_data['seg_baf'][i])
                if gene_table_data['seg_baf'][i] <= loh_cutoff:
                    gene_table_data['seg_loh_status'][i] = 'YES'
                    if (genes_list is not None and (genes_list is list or genes_list is str) and ("all" in genes_list or (gene_table_data['ENTREZ_ID'][i] is not None and gene_table_data['ENTREZ_ID'][i] in genes_list))):
                        gene_table_data['gene_loh_color'][i] = 'gold'
                    else:
                        gene_table_data['gene_loh_color'][i]="transparent"
                else:
                    gene_table_data['seg_loh_status'][i] = 'NO'
                    gene_table_data['gene_loh_color'][i] = 'transparent'
        
        for i in range(0,len(gene_table_data['seg_baf'])):
            if gene_table_data['gene_cnv_color'][i] == "transparent" and gene_table_data['gene_loh_color'][i] == "transparent":
                gene_table_data['gene_y'][i] = -1
                gene_table_data['gene_height'][i] = 0
                gene_table_data['gene_y_label'][i] = -1
            elif gene_table_data['gene_cnv_color'][i] != "transparent" and gene_table_data['gene_loh_color'][i] == "transparent":
                gene_table_data['gene_y'][i] = 0.3
                gene_table_data['gene_height'][i] = 0.6
                gene_table_data['gene_y_label'][i] = 0.50 - ((count_top % 5.5) * 0.060)
                count_top += 1
            elif gene_table_data['gene_cnv_color'][i] == "transparent" and gene_table_data['gene_loh_color'][i] != "transparent":
                gene_table_data['gene_y'][i] = -0.3
                gene_table_data['gene_height'][i] = 0.6
                gene_table_data['gene_y_label'][i] = -0.33 - ((count_bottom % 5.5) * 0.060)
                count_bottom += 1
            else:
                gene_table_data['gene_y'][i] = 0
                gene_table_data['gene_height'][i] = 1.2
                gene_table_data['gene_y_label'][i] = 0.28 - ((count_middle % 10.5) * 0.072)
                count_middle += 1
            
        if cnv_cutoff or loh_cutoff:
            if not cnv_cutoff:
                cnv_cutoff = 0
            if not loh_cutoff:
                loh_cutoff = 0.5
            for i in range(0,len(gene_table_data['SYMBOL'])):
                if gene_table_data['gene_log2'][i] < cnv_cutoff and gene_table_data['gene_log2'][i] > cnv_cutoff*-1 and gene_table_data['seg_baf'][i] > loh_cutoff:
                    for f in list(gene_table_data.keys()):
                        gene_table_data[f][i]=""
        y_slots = np.array(range(0,11))/10 - 0.5
    else:
        gene_table_data = dict(
            SYMBOL=[],
            ENTREZ_ID=[],
            CHR=[],
            START=[],
            END=[],
            seg_mean=[],
            gene_log2=[],
            gene_x=[],
            gene_width=[],
            gene_y=[],
            gene_y_label=[],
            gene_height=[],
            gene_cnv_color=[],
            seg_baf=[],
            seg_loh=[],
            seg_loh_status=[],
            gene_loh_color=[],
            label_alpha = []
        )

    return gene_table_data

def focals_json_to_datasources(focals_dict:dict, genome_info:dict, filter_to_chrom:str=None,genes_list:list=None,filter_by_list:bool=True,abs_position_dict:dict=None,stdev_cutoff:float=None,log2_cutoff:float=None,gene_color:str='#FFD700',exon_color:str='#CC6677', male_adjustment:bool=False):
    print('focals_json_to_datasources genes_list:')
    if (genes_list is None):
        print('genes_list is None')
    else:
        print(','.join(genes_list))

    focals_dict_genes = {'event_id':[], 'chrom':[],'start':[],'end':[],'log2':[],'parent_log2':[],'delta_log2':[],'focal_copies':[],'delta_copies':[],'text_status':[],'delta_sd':[],'event_midpoint':[], 'symbol':[], 'entrez_id':[],'color':[],'cytoband':[],'display_y':[],'gene':[],'baseline_adj':[]}
    focals_dict_exons = {'event_id':[], 'chrom':[],'start':[],'end':[],'log2':[],'parent_log2':[],'delta_log2':[],'focal_copies':[],'delta_copies':[],'text_status':[],'delta_sd':[],'event_midpoint':[], 'symbol':[], 'entrez_id':[],'mane_status':[], 'transcript_id':[], 'transcripts_affected':[],'exon':[], 'transcript_text':[], 'color':[], 'cytoband':[], 'display_y':[],'gene':[],'baseline_adj':[]}
    focals_dict_stored = {'event_type':[], 'event_id':[], 'chrom':[],'start':[],'end':[],'log2':[],'parent_log2':[],'delta_log2':[],'focal_copies':[],'delta_copies':[],'text_status':[],'delta_sd':[],'event_midpoint':[], 'symbol':[], 'entrez_id':[],'mane_status':[], 'transcript_id':[], 'transcripts_affected':[], 'exon':[],'transcript_text':[], 'color':[], 'cytoband':[], 'display_y':[],'gene':[], 'baseline_adj':[]}
    
    focal_effects_dict={}
    
    if not (type(genes_list) is list and len(genes_list) == 0):
        if genes_list is not None:
            genes_list = [str(x) for x in genes_list]

        if not filter_by_list:
            focals_dict_genes['phenotype_gene']=[]
            focals_dict_exons['phenotype_gene']=[]
            focals_dict_stored['phenotype_gene']=[]
        
        for i in focals_dict:
            focal_event_type = i['focal_event_type']
            event_start = i['event_start']
            event_end = i['event_end']
            delta_sd = i['event_stdevs_from_segment']
            chromosome = i['chromosome']
            gene = i['GENE']
        
            if chromosome in genome_info:
                
                if male_adjustment and ('x' in chromosome.lower() or 'y' in chromosome.lower()):
                    baseline_adj = 1
                else:
                    baseline_adj = 0
        
                cytobands = []
                for j in genome_info[chromosome]['cytobands']:
                    cyto_start = j['start']
                    cyto_end = j['end']
                    if (event_start >= cyto_start and event_start <= cyto_end) or (event_end <= cyto_end and event_end >= cyto_start):
                        cytobands.append(j['name'])
        
                log2 = i['weighted_average_log2_copy_ratio']
                parent_log2 = i['parent_log2_copy_ratio']
                delta_log2 = log2 - parent_log2
                if abs(log2) <= log2_cutoff:
                    focal_copy_number = 2
                else:
                    focal_copy_number = 2 ** (1+log2)
                focal_copy_number = focal_copy_number / (1+baseline_adj)
                delta_copy_number = focal_copy_number - ((2 ** (1+parent_log2))/(baseline_adj+1))
                event_id = i['events'][0]['event_id']
                entrez_id = gene['ENTREZ_ID']
                if 'symbol' in i['events'][0]['event_info']:
                    symbol = i['events'][0]['event_info']['symbol']
                else:
                    symbol = 'Obsolete Gene #' + entrez_id
                print('Focal event entrez id: ' + entrez_id)

                # Flag the focal event as phenotype-related if it is in the provided list of phenotype genes
                if genes_list is not None and entrez_id in genes_list:
                    gene['PHENOTYPE_GENE'] = True

                # if 'entrez_id' in i['events'][0]['event_info']:
                #     entrez_id = str( i['events'][0]['event_info']['entrez_id'])
                # else:
                #     entrez_id = event_id
                event_midpoint = (event_start + event_end)/2
                if abs_position_dict:
                    event_midpoint += abs_position_dict[chromosome]
        
                if log2_cutoff and (parent_log2 > log2_cutoff or parent_log2 < (-1*log2_cutoff)):
                    #Additional focal on top of gene CNV
                    if not log2 > log2_cutoff and not log2 < (-1*log2_cutoff):
                        # Parent is CNV, focal on its own is not
                        text_status = 'Focal Re-Normalization'
                    elif delta_log2 > log2_cutoff:
                        if delta_copy_number > 3:
                            text_status = 'Focal Amplification'
                        else:
                            text_status = 'Focal Gain'
                        if parent_log2 > log2_cutoff:
                            text_status = f'Additional {text_status}'
                        elif parent_log2 < (-1*log2_cutoff):
                            text_status = f'{text_status} (Reversal)'
                    elif delta_log2 < (-1 * log2_cutoff):
                        if focal_copy_number < (0.2 / (1 + baseline_adj)):
                            text_status = 'Focal Complete Loss'
                        else:
                            text_status = 'Focal Loss'
                        if parent_log2 > log2_cutoff:
                            text_status = f'{text_status} (Reversal)'
                        elif parent_log2 < (-1*log2_cutoff):
                            text_status = f'Additional {text_status}'
                    else:
                        text_status = 'Indeterminate Focal Event'
                else:
                    if log2 > log2_cutoff:
                        if focal_copy_number > 5:
                            text_status = 'Focal Amplification'
                        else:
                            text_status = 'Focal Gain'
                    elif log2 < (-1*log2_cutoff):
                        if focal_copy_number < (0.2 / (1+baseline_adj)):
                            text_status = 'Focal Biallelic Loss'
                        else:
                            text_status = 'Focal Loss'
                    else:
                        text_status = 'Indeterminate Focal Event'
                focal_effects_dict[entrez_id] = {'log2':log2,'status':text_status}
            
                dict_link = focals_dict_stored
                if not filter_to_chrom or str(chromosome) == str(filter_to_chrom):
                    if focal_event_type == 'exon':
                        mane_list = []
                        transcript_list = []
                        exon_list = []
                        transcript_text_list = []
                        sorting_scores = []
                        slice = 0
                        for j in i['events']:
                            event_info = j['event_info']
                            if 'mane_status' in event_info:
                                mane_status = event_info['mane_status']
                            else:
                                mane_status = ""
                            mane_list.append(mane_status)
                
                            if 'transcript_id' in event_info:
                                transcript_id = event_info['transcript_id']
                            else:
                                transcript_id = ""
                            transcript_list.append(transcript_id)
                    
                            if 'exon_number' in event_info:
                                exon_number = event_info['exon_number']
                            else:
                                exon_number = ""
                            exon_list.append(exon_number)
                    
                            transcript_string = f'{transcript_id} Exon #{exon_number}'
                            if mane_status == 'MANE Plus Clinical':
                                transcript_string = f'<i>{transcript_string}</i>'
                            transcript_text_list.append(transcript_string)
                    
                            if mane_status == 'MANE Select':
                                sorting_scores.append(0)
                                slice += 1
                            elif mane_status == 'MANE Plus Clinical':
                                sorting_scores.append(1)
                                slice += 1
                            elif transcript_id.startswith('NM'):
                                sorting_scores.append(2)
                            elif transcript_id.startswith('X'):
                                sorting_scores.append(4)
                            else:
                                sorting_scores.append(3)
                        
                        sorting_scores, mane_list, transcript_list, exon_list, transcript_text_list = (sort_together([sorting_scores, mane_list, transcript_list, exon_list, transcript_text_list]))
                
                        transcript_text = "<br>\n".join(transcript_text_list[0:max(slice,1)])
                        transcripts_affected = len(transcript_list)
            
                    if len(entrez_id) > 0 and (genes_list is None or genes_list is [] or entrez_id in genes_list):
                        dict_link['event_type'].append(focal_event_type)
                        dict_link['chrom'].append(chromosome)
                        dict_link['start'].append(event_start)
                        dict_link['event_midpoint'].append(event_midpoint)
                        dict_link['end'].append(event_end)
                        dict_link['log2'].append(log2)
                        dict_link['parent_log2'].append(parent_log2)
                        dict_link['delta_log2'].append(delta_log2)
                        dict_link['delta_sd'].append(delta_sd)
                        dict_link['symbol'].append(symbol)
                        dict_link['event_id'].append(event_id)
                        dict_link['entrez_id'].append(entrez_id)
                        dict_link['focal_copies'].append(focal_copy_number)
                        dict_link['delta_copies'].append(delta_copy_number)
                        dict_link['text_status'].append(text_status)
                        dict_link['cytoband'].append(cytobands)
                        dict_link['baseline_adj'].append(baseline_adj)
                        print(f'Adding gene:')
                        print(str(gene))
                        dict_link['gene'].append(gene)
                        if log2 > 1.97:
                            display_y = 1.97
                        elif log2 < -1.97:
                            display_y = -1.97
                        else:
                            display_y = log2
                        dict_link['display_y'].append(display_y)
                        if focal_event_type == 'exon':
                            dict_link['mane_status'].append(mane_list)
                            dict_link['transcript_id'].append(transcript_list)
                            dict_link['transcripts_affected'].append(transcripts_affected)
                            dict_link['exon'].append(exon_list)
                            dict_link['transcript_text'].append(transcript_text)
                            dict_link['color'].append(exon_color)
                        else:
                            dict_link['mane_status'].append([])
                            dict_link['transcript_id'].append([])
                            dict_link['transcripts_affected'].append('')
                            dict_link['exon'].append([])
                            dict_link['transcript_text'].append("")
                            dict_link['color'].append(gene_color)
    
        for i in range(0,len(focals_dict_stored['event_type'])):
            if (genes_list is not None and focals_dict_stored['entrez_id'][i] in genes_list):
                if focals_dict_stored['event_type'][i] == 'exon':
                    dict_link = focals_dict_exons
                else:
                    dict_link = focals_dict_genes

                if (not stdev_cutoff or abs(focals_dict_stored['delta_sd'][i]) >= stdev_cutoff):
                    for j in list(dict_link.keys()):
                        dict_link[j].append(focals_dict_stored[j][i])

    print(focals_dict_stored)

    stored_focal_datasource = ColumnDataSource(data=focals_dict_stored, name='stored_focals')
    gene_focal_datasource = ColumnDataSource(data=focals_dict_genes, name='focal_genes')
    exon_focal_datasource = ColumnDataSource(data=focals_dict_exons, name='focal_exons')
    focal_effects_datasource = ColumnDataSource(data={'dict':[focal_effects_dict]}, name='focal_effects')
    return gene_focal_datasource, exon_focal_datasource, stored_focal_datasource, focal_effects_datasource


### Functions to generate custom JS to update plots. ###
def generate_callbackJS(source_seg, source_cnv, source_cnv_bins, source_loh, source_colors, cnv_slider, loh_slider, cell_slider, offset_slider, alpha_slider, genome:bool=False):
    sliders_dict={'cnv_slider':(cnv_slider,'cnv_cutoff'),'loh_slider':(loh_slider,'loh_cutoff'),'cell_slider':(cell_slider,'cell'),'offset_slider':(offset_slider,'offset'),'alpha_slider':(alpha_slider,'alpha')}
    sliders = []
    for i in sliders_dict.keys():
        if sliders_dict[i][0]:
            sliders.append(f'const {sliders_dict[i][1]} = {i}.value;')
        else:
            sliders.append(f'const {sliders_dict[i][1]} = 1;')
    sliders = "\n".join(sliders)
    functions = """
function cnvloh_to_status(log2, baf, cnv_cutoff, loh_cutoff) {
    var num_copies;
    var is_loh;
    var copy_status;
    if (Math.abs(log2) >= cnv_cutoff){num_copies = 2 ** (1+log2);} else {num_copies = 2;}
    if (baf <= loh_cutoff || 1-baf <= loh_cutoff){is_loh = true;} else {is_loh = false;}
        
    if (num_copies >= 5.0){
        copy_status = 'Amplification';
        if (is_loh == false){copy_status = copy_status + ' (Balanced)';}
        }
    else if (num_copies > 2.6){
        copy_status = 'Gain'
        if (num_copies > 3.5){
            if (is_loh == true){copy_status = copy_status + ' (Unbalanced)';} else {copy_status = copy_status + ' (Balanced)';}
            }
        else {if (is_loh == false){copy_status = copy_status + ' (! Anomalous LOH)';}}
        }
    else if (num_copies > 2.2){
        copy_status = 'Weak Gain';
        if (is_loh == true){copy_status = copy_status + ' (LOH Supported)';}
        }
    else if (num_copies < 0.5){copy_status = 'Biallelic Loss';}
    else if (num_copies < 1.4){
        copy_status = 'Loss';
        if (is_loh == false){copy_status = copy_status + ' (!Anomalous LOH)';}
        }
    else if (num_copies < 1.8){
        copy_status = 'Weak Loss';
        if (is_loh == true){copy_status = copy_status + ' (LOH Supported)';}
        }
    else {
        copy_status = 'Copy-Neutral';
        if (is_loh == true){copy_status = copy_status + ' LOH';}
        }
    return copy_status
    }
    """
    
    constants_pt_1="""
const data = source.data;
const bins_data = source_bins.data;
const colors_data = source_colors.data;
const loh_colors = colors_data['loh_colors'];

const seg_y = data['cnv_y'];
const seg_y_stored = data['cnv_y_stored'];
const seg_num_copy = data['cnv_num_copy'];
const seg_y_height = data['cnv_y_height'];
const seg_y_height_stored = data['cnv_y_height_stored'];
const seg_baseline_adj = data['seg_baseline_adj'];
const seg_loh = data['seg_loh_height'];
const seg_loh_stored = data['seg_loh_height_stored'];
const seg_loh_baf = data['seg_loh_baf'];
const seg_loh_status = data['seg_loh_status'];
const seg_loh_color = data['seg_loh_color'];
const seg_alpha = data['seg_alpha'];
const seg_cnvloh_status = data['cnvloh_status'];
const cnv_y_display = data['cnv_y_display'];
const cnv_label_alpha = data['cnv_label_alpha'];
const cnv_label_pos = data['cnv_label_pos'];
const cnv_label_text = data['cnv_label_text'];

const bin_y = bins_data['cnv_y'];
const bin_height = bins_data['cnv_bin_height'];
const bin_sd = bins_data['cnv_bin_sd'];
const bin_y_stored = bins_data['cnv_y_stored'];
const bin_height_stored = bins_data['cnv_bin_height_stored'];
const bin_sd_stored = bins_data['cnv_bin_sd_stored'];
"""
    code_pt_1="""
    console.log('Begin segments code_pt_1');
    const startTime = new Date();

let label_pos_adjustment_top = 0;
let label_pos_adjustment_bottom = 0
    
for (var i = 0; i < seg_y.length; i++) {
    seg_y[i] = Math.log2(Math.max(0,( Math.pow(2,(seg_y_stored[i]+offset+1))-(2*(1-cell)) )/cell )) - 1;
    seg_y_height[i] = Math.log2(Math.max(0,( Math.pow(2,(seg_y_height_stored[i]+1))-(2*(1-cell)) )/cell )) - 1;
    seg_num_copy[i] = Math.pow(2,seg_y[i]+1) / (1 + seg_baseline_adj[i]);
    var loh_calc = Math.min(1,Math.max(0,(0.5-(0.5-seg_loh_stored[i]/cell))));
    var baf_calc = 0.5 - loh_calc/2;
    seg_loh[i] = loh_calc;
    seg_loh_baf[i] = baf_calc;
    seg_loh_color[i] = loh_colors[Math.round((1-loh_calc)*loh_colors.length)];
    seg_alpha[i] = alpha;
    if (seg_loh_baf[i] <= loh_cutoff){seg_loh_status[i] = "YES";} else {seg_loh_status[i] = "NO";}
    seg_cnvloh_status[i] = cnvloh_to_status(seg_y[i],seg_loh_baf[i],cnv_cutoff,loh_cutoff);
    cnv_y_display[i] = Math.min(1.97,Math.max(-1.97,seg_y[i]));
    cnv_label_text[i] = String(seg_y[i].toFixed(1));
    if (seg_y[i] > 0){cnv_label_pos[i] = Bokeh.documents[0].get_model_by_name('cnv_plot').y_range.end-0.3-label_pos_adjustment_top; if(label_pos_adjustment_top == 0){label_pos_adjustment_top = 0.2;}else{label_pos_adjustment_top = 0;}} else {cnv_label_pos[i] = Bokeh.documents[0].get_model_by_name('cnv_plot').y_range.start + label_pos_adjustment_bottom; if(label_pos_adjustment_bottom == 0){label_pos_adjustment_bottom = 0.2;}else{label_pos_adjustment_bottom = 0;}}
    if (cnv_y_display[i] != seg_y[i]){cnv_label_alpha[i] = 1;} else {cnv_label_alpha[i] = 0;}
    }
    
for (var i = 0; i < bin_y.length; i++){
    bin_y[i] = Math.log2(Math.max(0,( Math.pow(2,(bin_y_stored[i]+offset+1))-(2*(1-cell)) )/cell )) - 1;
    bin_height[i] = Math.log2(Math.max(0,( Math.pow(2,(bin_height_stored[i]+1))-(2*(1-cell)) )/cell )) - 1;
    bin_sd[i] = bin_sd_stored[i] / cell;
    }
    
    console.log('End segments code_pt_1');
    const endTime = new Date();
    console.log('Elapsed time: ' + (endTime - startTime) + ' ms')
    
"""

    change_pt_1="""
source.change.emit();
source_bins.change.emit();
"""
    if not genome:
        constants_pt_2="""
const cnv_data = source_cnv.data;
const loh_data = source_loh.data;
const cnv_x = cnv_data['cnv_x'];
const cnv_y = cnv_data['cnv_y'];
const cnv_y_stored = cnv_data['cnv_y_stored'];
const cnv_alpha = cnv_data['cnv_alpha'];
const loh_y = loh_data['loh_y'];
const loh_y_stored = loh_data['loh_y_stored'];
const loh_alpha = loh_data['loh_alpha'];

"""
        code_pt_2="""
for (var i = 0; i < cnv_y.length; i++) {
    cnv_y[i] = Math.log2(Math.max(0,( Math.pow(2,(cnv_y_stored[i]+offset+1))-(2*(1-cell)) )/cell )) - 1;
    cnv_alpha[i] = alpha;
    }
for (var i = 0; i < loh_y.length; i++) {
    const y_calc = ((loh_y_stored[i] - 0.5)/cell) + 0.5;
    loh_y[i] = Math.min(1,Math.max(0,y_calc));
    loh_alpha[i] = alpha;
    }
"""
        change_pt_2="""
source_cnv.change.emit();
source_loh.change.emit();
"""
        code = "".join([sliders,functions,constants_pt_1,constants_pt_2,code_pt_1,code_pt_2,change_pt_1,change_pt_2])
        callback = CustomJS(args=dict(
        source=source_seg, source_cnv=source_cnv, source_bins=source_cnv_bins, source_loh=source_loh, source_colors=source_colors,
        cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider,
        offset_slider=offset_slider, alpha_slider=alpha_slider),code=code)
    else:
        code = "".join([sliders,functions,constants_pt_1,code_pt_1,change_pt_1])
        callback = CustomJS(args=dict(
        source=source_seg, source_bins=source_cnv_bins, source_colors=source_colors,
        cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider,
        offset_slider=offset_slider, alpha_slider=alpha_slider),code=code)
    return callback

def generate_callbackJS_genelist(source_genetable_stored, genelist_field):
    callback_genelist = CustomJS(args=dict(
        source_stored=source_genetable_stored,
        genelist_field = genelist_field),
                        code="""
        console.log("Genelist was updated:");
        console.log(genelist_field.value_input);
                        
        const genelist_string = genelist_field.value_input;
        var genes_set = undefined;
        
        if (!(!genelist_string || genelist_string == '' || genelist_string == '[]'))
        {
            genes_set = new Set(genelist_string.replace('[','').replace(']','').split(','));
        }
        console.log(genes_set);
        
        const data_stored = source_stored.data;
        
        const seg_mean_stored = data_stored['seg_mean'];
        const chr_stored = data_stored['CHR'];
        const symbol_stored = data_stored['SYMBOL'];
        const entrez_stored = data_stored['ENTREZ_ID'];
        const gene_start_stored = data_stored['START'];
        const gene_end_stored = data_stored['END'];
        const gene_x_stored = data_stored['gene_x'];
        const gene_width_stored = data_stored['gene_width'];
        const gene_y_stored = data_stored['gene_y'];
        const gene_height_stored = data_stored['gene_height'];
        const gene_cnv_color_stored = data_stored['gene_cnv_color'];
        const seg_baf_stored = data_stored['seg_baf'];
        const seg_loh_stored = data_stored['seg_loh'];
        var stored_label_alpha = data_stored['label_alpha'];
        
        for (var i = 0; i < seg_mean_stored.length; i++) {
            if (genes_set == null || genes_set.has("all") || (genes_set.has(entrez_stored[i].toString()) && !genes_set.has("none")))
            {stored_label_alpha[i] = 1;} else {stored_label_alpha[i] = 0;}
        }
    """)
    return callback_genelist

def generate_callbackJS_gene_table(source_genetable, source_genetable_stored, source_focal_effects, source_colors, cnv_slider, loh_slider, cell_slider, offset_slider, genelist_field):
    callback_table = CustomJS(args=dict(
        source=source_genetable, source_stored=source_genetable_stored, source_colors=source_colors, 
        source_focal_effects=source_focal_effects,
        cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider,
        offset_slider=offset_slider,genelist_field=genelist_field),
                        code="""
function get_cnv_color(log2,cnv_cutoff) {
    var color_index;
    if (log2 > 1.56){color_index = 0;}
    else if (log2 > 0.85){color_index = 1;}
    else if (log2 > cnv_cutoff){color_index = 2;}
    else if (log2 < -2){color_index = 6;}
    else if (log2 < -0.85){color_index = 5;}
    else if (log2 < -1*cnv_cutoff){color_index = 4;}
    else {color_index = 3;}
    return color_index;
}
        
        console.log('Attempting gene table update');
        const startTime = new Date();

        const data = source.data;
        const data_stored = source_stored.data;
        const data_focal_effects = source_focal_effects.data;
    
        const colors_data = source_colors.data;
        const cnv_cutoff = cnv_slider.value;
        const loh_cutoff = loh_slider.value;
        const cell = cell_slider.value;
        const offset = offset_slider.value;
    
        const cnv_colors = colors_data['cnv_colors'];
        const seg_mean = data['seg_mean'];
        const gene_log2 = data['gene_log2'];
        const chr = data['CHR'];
        const symbol = data['SYMBOL'];
        const entrez_id = data['ENTREZ_ID'];
        const gene_start = data['START'];
        const gene_end = data['END'];
        const gene_x = data['gene_x'];
        const gene_width = data['gene_width'];
        const gene_y = data['gene_y'];
        const gene_y_label = data['gene_y_label'];
        const gene_height = data['gene_height'];
        const gene_cnv_color = data['gene_cnv_color'];
        const gene_loh_color = data['gene_loh_color'];
        const seg_baf = data['seg_baf'];
        const seg_loh = data['seg_loh'];
        const seg_loh_status = data['seg_loh_status'];
        const gene_label_alpha = data['label_alpha'];
    
        const seg_mean_stored = data_stored['seg_mean'];
        const chr_stored = data_stored['CHR'];
        const symbol_stored = data_stored['SYMBOL'];
        const entrez_stored = data_stored['ENTREZ_ID'];
        const gene_start_stored = data_stored['START'];
        const gene_end_stored = data_stored['END'];
        const gene_x_stored = data_stored['gene_x'];
        const gene_width_stored = data_stored['gene_width'];
        const gene_y_stored = data_stored['gene_y'];
        const gene_height_stored = data_stored['gene_height'];
        const gene_cnv_color_stored = data_stored['gene_cnv_color'];
        const seg_baf_stored = data_stored['seg_baf'];
        const seg_loh_stored = data_stored['seg_loh'];
        const stored_label_alpha = data_stored['label_alpha'];
        
        var count = 0;
        var count_top = 0;
        var count_middle = 0;
        var count_bottom = 0;
        
        for (var i = 0; i < seg_mean_stored.length; i++) {
            var cnv_calc = Math.log2(Math.max(0,( Math.pow(2,(seg_mean_stored[i] + offset + 1))-(2*(1-cell)) )/cell )) - 1;
            var calc_gene_log2 = cnv_calc;
            var loh_calc = Math.min(1,Math.max(0,(0.5-(0.5-seg_loh_stored[i]/cell))));
            var baf_calc = 0.5 - loh_calc/2;
            var loh_status = "NO";
            var loh_color = "transparent";
            var cnv_color = "transparent";
            var gene_symbol = symbol_stored[i];
            var gene_id = entrez_stored[i];
            var cnv_calc_color = 3;
            var seg_listed = stored_label_alpha[i];
            if (baf_calc <= loh_cutoff){loh_status = "YES"; loh_color = "gold";};
            if ((seg_listed == 1) && (cnv_calc <= cnv_cutoff*-1 || cnv_calc >= cnv_cutoff || loh_status == 'YES' || gene_id in data_focal_effects)) {
                chr[count] = chr_stored[i];
                symbol[count] = symbol_stored[i];
                entrez_id[count] = entrez_stored[i];
                gene_start[count] = gene_start_stored[i];
                gene_end[count] = gene_end_stored[i];
                gene_x[count] = gene_x_stored[i];
                gene_width[count] = gene_width_stored[i];
                gene_loh_color[count] = loh_color;
                if (gene_id in data_focal_effects){
                    // console.log(`${gene_symbol} (${gene_id}) is focal.`);
                    
                    // Find the most extreme focal log2 for this gene, which becomes the gene's reported log2 value.
                    // This code sorts the log2 values by their absolute value in descending order and selects the first entry.
                    calc_gene_log2 = data_focal_effects[gene_id]['log2'].sort((a, b) => Math.abs(b) - Math.abs(a))[0];
                    
                    if (data_focal_effects[gene_id]['status'][0].includes('Revers')){
                        cnv_color = "mediumorchid";
                    } else {
                        cnv_calc_color = get_cnv_color(calc_gene_log2,cnv_cutoff);
                        cnv_color = cnv_colors[cnv_calc_color];
                    }
                    if (loh_status == "YES"){gene_y[count] = 0; gene_height[count] = 1.2; gene_y_label[count] = 0.28 - ((count_middle % 10.5)*0.072); count_middle ++;}
                    else {gene_y[count] = 0.25; gene_height[count] = 0.5; gene_y_label[count] = 0.5 - ((count_top % 5.5)*0.06); count_top ++;}
                    //console.log(`${cnv_calc_color}: ${cnv_color}`);
                }
                else if (cnv_calc <= cnv_cutoff*-1 || cnv_calc >= cnv_cutoff){
                    cnv_calc_color = get_cnv_color(cnv_calc,cnv_cutoff);
                    cnv_color = cnv_colors[cnv_calc_color];
                    if (loh_status == "YES"){gene_y[count] = 0; gene_height[count] = 1.2; gene_y_label[count] = 0.28 - ((count_middle % 10.5)*0.072); count_middle ++;}
                    else {gene_y[count] = 0.3; gene_height[count] = 0.6; gene_y_label[count] = 0.50 - ((count_top % 5.5)*0.06); count_top ++;}
                } else {gene_y[count] = -0.3; gene_height[count] = 0.6; gene_y_label[count] = -0.33 - ((count_bottom % 5.5)*0.06); count_bottom ++;}
                seg_mean[count] = cnv_calc;
                gene_log2[count] = calc_gene_log2;
                gene_cnv_color[count] = cnv_color;
                seg_baf[count] = baf_calc;
                seg_loh[count] = loh_calc;
                seg_loh_status[count] = loh_status;
                gene_label_alpha[count] = stored_label_alpha[i];
                count ++;
            }
        }
        for (var j = count; j < seg_mean_stored.length; j++) {
                seg_mean[j] = '';
                chr[j] = '';
                symbol[j] = '';
                entrez_id[j] = '';
                gene_start[j] = '';
                gene_end[j] = '';
                gene_x[j] = '';
                gene_width[j] = '';
                gene_y[j] = '';
                gene_height[j] = '';
                gene_cnv_color[j] = '';
                seg_baf[j] = '';
                seg_loh[j] = '';
                gene_label_alpha[j] = 0;
                gene_log2[j] = '';
        }
    source.change.emit();
    console.log('Gene table update succeeded');
    const endTime = new Date();
    console.log('Elapsed time: ' + (endTime - startTime) + ' ms')
    """)
    return callback_table


def generate_callbackJS_focals(source_focals_stored, source_focals_genes, source_focals_exons, source_focal_effects, cnv_slider, cell_slider, offset_slider, sd_slider, genelist_field):
    callback_table = CustomJS(args=dict(
        source_genes=source_focals_genes, source_exons=source_focals_exons, source_stored=source_focals_stored, source_focal_effects=source_focal_effects,
        cnv_slider=cnv_slider, cell_slider=cell_slider,
        offset_slider=offset_slider,sd_slider=sd_slider,
        genelist_field=genelist_field),
                        code="""
        function focal_to_status(log2, parent_log2, cnv_cutoff) {
            var num_copies;
            var parent_num_copies;
            var delta_num_copies;
            var delta_log2;
            var text_status;
            if (Math.abs(log2) >= cnv_cutoff){num_copies = 2 ** (1+log2);} else {num_copies = 2;}
            if (Math.abs(parent_log2) >= cnv_cutoff){num_copies = 2 ** (1+parent_log2);} else {parent_num_copies = 2;}
            delta_log2 = log2 - parent_log2;
            delta_num_copies = num_copies - parent_num_copies;
            
            if (Math.abs(parent_log2) > cnv_cutoff){
                if (Math.abs(log2) < cnv_cutoff){text_status = 'Focal Re-Normalization';}
                else if (delta_log2 > cnv_cutoff){
                    if (delta_num_copies > 3){text_status = 'Focal Amplification';}
                    else {text_status = 'Focal Gain';}
                    if (parent_log2 > cnv_cutoff){text_status = 'Additional ' + text_status;}
                    else if (parent_log2 < (-1*cnv_cutoff)){text_status = text_status + ' (Reversal)';}
                }
                else if (delta_log2 < (-1 * cnv_cutoff)){
                    if (num_copies < 0.2) {text_status = 'Focal Complete Loss';}
                    else {text_status = 'Focal Loss';}
                    if (parent_log2 < (-1*cnv_cutoff)){text_status = 'Additional ' + text_status;}
                    else if (parent_log2 > cnv_cutoff){text_status = text_status + ' (Reversal)';}
                }
                else {text_status = 'Indeterminate Focal Event';}
            } 
            else {
                if (log2 > cnv_cutoff){
                    if (num_copies > 5){text_status = 'Focal Amplification';}
                    else {text_status = 'Focal Gain';}
                }
                else if (log2 < -1*cnv_cutoff){
                    if (num_copies < 0.2){text_status = 'Focal Biallelic Loss';}
                    else {text_status = 'Focal Loss';}
                }
                else {text_status = 'Indeterminate Focal Event';}
            }
        return text_status;
        }
            
            
        console.log('Attempting focals data update');
        const startTime = new Date();

        const data_genes = source_genes.data;
        const data_exons = source_exons.data;
        const data_stored = source_stored.data;
        const data_focal_effects = source_focal_effects.data;
        console.log(data_focal_effects);
        
        const cell = cell_slider.value;
        const offset = offset_slider.value;
        const cnv_cutoff = cnv_slider.value;
        const sd_cutoff = sd_slider.value;
        
        const genelist_string = genelist_field.value_input;
        var genes_set = undefined;
        
        if (!(!genelist_string || genelist_string == '' || genelist_string == '[]'))
        {
            genes_set = new Set(genelist_string.replace('[','').replace(']','').split(','));
        }
        
        const type_stored = data_stored['event_type'];
        const id_stored = data_stored['event_id'];
        const chr_stored = data_stored['chrom'];
        const symbol_stored = data_stored['symbol'];
        const start_stored = data_stored['start'];
        const end_stored = data_stored['end'];
        const focal_x_stored = data_stored['event_midpoint'];
        const focal_y_stored = data_stored['log2'];
        const parent_y_stored = data_stored['parent_log2'];
        const delta_sd_stored = data_stored['delta_sd'];
        const entrez_stored = data_stored['entrez_id'];
        const mane_stored = data_stored['mane_status'];
        const transcript_stored = data_stored['transcript_id'];
        const transcript_text_stored = data_stored['transcript_text'];
        const transcripts_affected_stored = data_stored['transcripts_affected'];
        const exon_stored = data_stored['exon'];
        const color_stored = data_stored['color'];
        const cytobands_stored = data_stored['cytoband'];
        const gene_stored = data_stored['gene'];
        const baseline_adj_stored = data_stored['baseline_adj'];
        
        const id_gene = data_genes['event_id'];
        const chr_gene = data_genes['chrom'];
        const symbol_gene = data_genes['symbol'];
        const start_gene = data_genes['start'];
        const end_gene = data_genes['end'];
        const focal_x_gene = data_genes['event_midpoint'];
        const parent_y_gene = data_genes['parent_log2'];
        const focal_y_gene = data_genes['log2'];
        const delta_y_gene = data_genes['delta_log2'];
        const delta_sd_gene = data_genes['delta_sd'];
        const entrez_gene = data_genes['entrez_id'];
        const focal_copy_gene = data_genes['focal_copies'];
        const delta_copy_gene = data_genes['delta_copies'];
        const color_gene = data_genes['color'];
        const text_status_gene = data_genes['text_status'];
        const cytobands_gene = data_genes['cytoband'];
        const gene_gene = data_genes['gene'];
        const display_y_gene = data_genes['display_y'];
        const baseline_adj_gene = data_genes['baseline_adj'];
        
        const id_exon = data_exons['event_id'];
        const chr_exon = data_exons['chrom'];
        const symbol_exon = data_exons['symbol'];
        const start_exon = data_exons['start'];
        const end_exon = data_exons['end'];
        const focal_x_exon = data_exons['event_midpoint'];
        const focal_y_exon = data_exons['log2'];
        const parent_y_exon = data_exons['parent_log2'];
        const delta_y_exon = data_exons['delta_log2'];
        const delta_sd_exon = data_exons['delta_sd'];
        const entrez_exon = data_exons['entrez_id'];
        const focal_copy_exon = data_exons['focal_copies'];
        const delta_copy_exon = data_exons['delta_copies'];
        const mane_exon = data_exons['mane_status'];
        const transcript = data_exons['transcript_id'];
        const transcripts_affected = data_exons['transcripts_affected'];
        const transcript_text = data_exons['transcript_text'];
        const exon = data_exons['exon'];
        const color_exon = data_exons['color'];
        const text_status_exon = data_exons['text_status'];
        const cytobands_exon = data_exons['cytoband'];
        const gene_exon = data_exons['gene'];
        const display_y_exon = data_exons['display_y'];
        const baseline_adj_exon = data_exons['baseline_adj'];
        
        var g = 0;
        var e = 0;
        
        var cnv_calc;
        var parent_cnv_calc
        var parent_cnv_calc;
        var text_status;
        var display_y;
        
        Object.keys(data_focal_effects).forEach(key => delete data_focal_effects[key]);
        for (var i = 0; i < type_stored.length; i++) {
            if ( Math.abs(delta_sd_stored[i]) >= sd_cutoff && (!genes_set || genes_set.has('all') || genes_set.has(entrez_stored[i]))){
                cnv_calc = Math.log2(Math.max(0,( Math.pow(2,(focal_y_stored[i] + offset + 1))-(2*(1-cell)) )/cell )) - 1;
                if (cnv_calc > Bokeh.documents[0].get_model_by_name('cnv_plot').y_range.end-0.05){
                    display_y = Bokeh.documents[0].get_model_by_name('cnv_plot').y_range.end-0.05;
                    } else if (cnv_calc < Bokeh.documents[0].get_model_by_name('cnv_plot').y_range.start+0.05){
                    display_y = Bokeh.documents[0].get_model_by_name('cnv_plot').y_range.start+0.05;}
                    else {display_y = cnv_calc;}
                parent_cnv_calc = Math.log2(Math.max(0,( Math.pow(2,(parent_y_stored[i] + offset + 1))-(2*(1-cell)) )/cell )) - 1;
                text_status = focal_to_status(cnv_calc,parent_cnv_calc,cnv_cutoff);
                if (!(entrez_stored[i] in data_focal_effects)){data_focal_effects[entrez_stored[i]]={'log2':[],'status':[]};}
                data_focal_effects[entrez_stored[i]]['log2'].push(cnv_calc);
                data_focal_effects[entrez_stored[i]]['status'].push(text_status);
                if ( type_stored[i] == 'exon' ){
                    id_exon[e] = id_stored[i];
                    chr_exon[e] = chr_stored[i];
                    symbol_exon[e] = symbol_stored[i];
                    start_exon[e] = start_stored[i];
                    end_exon[e] = end_stored[i]
                    baseline_adj_exon[e] = baseline_adj_stored[i];
                    focal_x_exon[e] = focal_x_stored[i];
                    focal_y_exon[e] = cnv_calc;
                    parent_y_exon[e] = parent_cnv_calc;
                    delta_y_exon[e] = cnv_calc - parent_cnv_calc;
                    if (Math.abs(cnv_calc) <= cnv_cutoff) {
                        focal_copy_exon[e] = 2;
                    } else {
                        focal_copy_exon[e] = (2 ** (1+cnv_calc));
                    }
                    focal_copy_exon[e] = focal_copy_exon[e];
                    delta_copy_exon[e] = focal_copy_exon[e] - (2 ** (1+parent_cnv_calc));
                    focal_copy_exon[e] = focal_copy_exon[e] / (1 + baseline_adj_exon[e]);
                    delta_copy_exon[e] = delta_copy_exon[e] / (1 + baseline_adj_exon[e]);
                    delta_sd_exon[e] = delta_sd_stored[i];
                    entrez_exon[e] = entrez_stored[i];
                    mane_exon[e] = mane_stored[i];
                    transcript[e] = transcript_stored[i];
                    transcripts_affected[e] = transcripts_affected_stored[i];
                    transcript_text[e] = transcript_text_stored[i];
                    exon[e] = exon_stored[i];
                    color_exon[e] = color_stored[i];
                    text_status_exon[e] = text_status;
                    cytobands_exon[e] = cytobands_stored[i];
                    gene_exon[e] = gene_stored[i];
                    display_y_exon[e] = display_y;
                    e ++;
                } else {
                    id_gene[g] = id_stored[i];
                    baseline_adj_gene[g] = baseline_adj_stored[i];
                    chr_gene[g] = chr_stored[i];
                    symbol_gene[g] = symbol_stored[i];
                    start_gene[g] = start_stored[i];
                    end_gene[g] = end_stored[i];
                    focal_x_gene[g] = focal_x_stored[i];
                    focal_y_gene[g] = cnv_calc;
                    parent_y_gene[g] = parent_cnv_calc;
                    delta_y_gene[g] = cnv_calc - parent_cnv_calc;
                    if (Math.abs(cnv_calc) <= cnv_cutoff) {
                        focal_copy_gene[g] = 2;
                    } else {
                        focal_copy_gene[g] = (2 ** (1+cnv_calc));
                    }
                    focal_copy_gene[g] = focal_copy_gene[g];
                    delta_copy_gene[g] = focal_copy_gene[g] - (2 ** (1+parent_cnv_calc));
                    focal_copy_gene[g] = focal_copy_gene[g] / (1+baseline_adj_gene[g]);
                    delta_copy_gene[g] = delta_copy_gene[g] / (1+baseline_adj_gene[g]);
                    delta_sd_gene[g] = delta_sd_stored[i];
                    entrez_gene[g] = entrez_stored[i];
                    color_gene[g] = color_stored[i];
                    text_status_gene[g] = text_status;
                    cytobands_gene[g] = cytobands_stored[i];
                    gene_gene[g] = gene_stored[i];
                    display_y_gene[g] = display_y;
                    g++;
                }
            }
        }
        for (var h = g; h < id_gene.length; h++ ){
            id_gene[h] = null;
            chr_gene[h] = null;
            symbol_gene[h] = null;
            start_gene[h] = null;
            end_gene[h] = null;
            focal_x_gene[h] = null;
            focal_y_gene[h] = null;
            delta_sd_gene[h] = null;
            entrez_gene[h] = null;
            color_gene[h] = 'transparent';
        }
        for (var f = e; f < id_exon.length; f++ ){
            id_exon[f] = null;
            chr_exon[f] = null;
            symbol_exon[f] = null;
            start_exon[f] = null;
            end_exon[f] = null;
            focal_x_exon[f] = null;
            focal_y_exon[f] = null;
            delta_sd_exon[f] = null;
            entrez_exon[f] = null;
            mane_exon[f] = null;
            transcript[f] = null;
            exon[f] = null;
            color_exon[f] = 'transparent';
        }

        //console.log('Focal genes after update:');
        //console.log(data_genes);
        //console.log(source_genes);

        //console.log('Focal exons after update:');
        //console.log(data_exons);
        //console.log(source_exons);

        console.log('Firing focal genes change');
        source_genes.change.emit();
        console.log('Firing focal exons change');
        source_exons.change.emit();
        //console.log('Focal source update succeeded');
        const endTime = new Date();
        console.log('Elapsed time: ' + (endTime - startTime) + ' ms')
    """)
    return callback_table

def generate_cellcalcJS(source_seg, cnv_slider, loh_slider, cell_slider, offset_slider, calcmode_radio, cell_calc_button):
    sliders_dict={'cnv_slider':(cnv_slider,'cnv_cutoff'),'loh_slider':(loh_slider,'loh_cutoff'),'cell_slider':(cell_slider,'cell'),'offset_slider':(offset_slider,'offset')}
    sliders = []
    for i in sliders_dict.keys():
        if sliders_dict[i][0]:
            sliders.append(f'const {sliders_dict[i][1]} = {i}.value;')
        else:
            sliders.append(f'const {sliders_dict[i][1]} = 1;')
    sliders = "\n".join(sliders)
    functions = """
function loh_to_baf(loh_input){
    var baf_result = 0.5 - loh_input;
    return baf_result;
    }
    
function log2_to_copy(log2_input){
    var copy_number = Math.max(0, Math.pow(2,log2_input+1) );
    return copy_number;
    }

function copy_to_log2(copy_num_input){
    var log2_calc = Math.log2(Math.max(copy_num_input),0.01) - 1;
    return log2_calc
    }
    
function cell_adj_copy(copy_num, cell_adj){
    var adj_copy = Math.max(0,((copy_num - 2) / cell_adj)+2);
    return adj_copy
    }

function cell_adj_baf(baf_input, cell_adj){
    var deviation = Math.abs(0.5 - baf_input);
    var adj_deviation = deviation / cell_adj;
    var adj_baf = Math.min(1, Math.max(0, 0.5 - adj_deviation));
    return adj_baf;
    }

function closest_copy_num(calc_copy, has_loh){
    var nearest_copy = Math.round(calc_copy);
    if (!has_loh && nearest_copy % 2 != 0)
    { 
        if (nearest_copy < 2)
        {nearest_copy = 0;}
        else {nearest_copy = nearest_copy + 1;}
    }
    return nearest_copy;
}

function get_loh_diff(copy_num, seg_baf)
{
    let percent_list = [];
    let diff_list = [];
    let genotype_list = [];
    
    let int_copy = Math.round(copy_num);
    
    for (let i = 0; i <= Math.ceil(int_copy/2); i++)
    {
        if (int_copy == 0)
        {
            percent_list.push(0.50);
            diff_list.push(Math.abs( 0.50 - seg_baf));
            genotype_list.push("-");
        }
        else
        {
            percent_list.push(i/(int_copy+0.0001));
            diff_list.push(Math.abs(  (i/(int_copy+0.0001)) - seg_baf));
            genotype_list.push( "A".repeat( int_copy-i ) + "B".repeat(i));
            
            percent_list.push(1 - (i/ (int_copy+0.0001)));
            diff_list.push(Math.abs(  (1 - (i/ (int_copy+0.0001))) - seg_baf));
            genotype_list.push( "A".repeat( i ) + "B".repeat(int_copy - i));
        }
    }

    let min_index = 0;
    for (let i = 1; i < diff_list.length; i++)
    {
        if (diff_list[i] < diff_list[min_index])
        {
            min_index = i;
        }
    }
    
    return {
        "diffScore":diff_list[min_index], 
        "closestBAF": percent_list[min_index], 
        "closestGenotype": genotype_list[min_index]
    };
}

function segment_badness(seg_log2, seg_baf, cnv_cutoff, loh_cutoff, cell_adj, offset) {
    var copy_num = log2_to_copy(seg_log2 + offset);
    var adj_copy_num = cell_adj_copy(copy_num, cell_adj);
    var adj_baf = cell_adj_baf(seg_baf, cell_adj);
    var has_loh = (adj_baf <= loh_cutoff);
    var has_cnv = (Math.abs(copy_to_log2(adj_copy_num)) >= cnv_cutoff);
    var copy_distance = adj_copy_num - closest_copy_num(adj_copy_num, has_loh);
    var loh_diff = get_loh_diff(adj_copy_num, adj_baf);
    var diff_score = loh_diff.diffScore;

    var badness = (Math.pow(copy_distance, 2) + Math.pow(10*diff_score,2)) / Math.pow(cell_adj,2);
    if (!has_cnv && !has_loh)
    {
        badness = badness / 2;
    }
    // console.log(badness);
    return badness;
}

function badness_of_fit(log2_list, loh_list, seg_len, cnv_cutoff, loh_cutoff, cell, offset){
    var total_badness = 0
    var seg_badness = 0;
    for (let i = 0; i < log2_list.length; i++) 
    {
        seg_badness = segment_badness(log2_list[i], loh_to_baf(loh_list[i]), cnv_cutoff, loh_cutoff, cell, offset);
        total_badness += seg_badness * Math.min(50000000,seg_len[i])/50000000;
    }
    return total_badness;
}

function calc_best_fit_cell(log2_list, baf_list, seg_len, cnv_cutoff, loh_cutoff, offset, show_logging = false){
    var cellularities = [];
    var fit_scores = [];
    var cell_iter = 1;
    var score = 99999;
    var best_fit_score = 99999;
    
    for (let i = 0; i < 800; i++)
    {
        cell_iter = 1 - i/1000;
        cellularities.push(cell_iter);
        score = badness_of_fit(log2_list, baf_list, seg_len, cnv_cutoff, loh_cutoff, cell_iter, offset);
        fit_scores.push(score);
        if (show_logging){console.log(`${cellularities[cellularities.length - 1]*100}% : ${fit_scores[fit_scores.length - 1]}`);}
    }
    
    best_fit_score = Math.min.apply(null,fit_scores); 
    
    for (let i = 0; i < fit_scores.length; i++)
    {
        if (best_fit_score == fit_scores[i])
        {
            return [cellularities[i], fit_scores[i]];
        }
    }
    return [undefined, undefined];
}

function calc_best_fit_baseline(log2_list, baf_list, seg_len, cnv_cutoff, loh_cutoff, cell, show_logging = false){
    var baselines = [];
    var fit_scores = [];
    
    var baseline_shift = 0;
    var score = 99999;
    
    for (let i = 0; i <= 100; i++)
    {
        baseline_shift = i/100;
        
        baselines.push(-1 * baseline_shift);
        score = badness_of_fit(log2_list, baf_list, seg_len, cnv_cutoff, loh_cutoff, cell, -1 * baseline_shift) / Math.pow(1.01 - baseline_shift,0.33);
        fit_scores.push(score);
        
        baselines.push(baseline_shift);
        score = badness_of_fit(log2_list, baf_list, seg_len, cnv_cutoff, loh_cutoff, cell, baseline_shift) / Math.pow(1.01 - baseline_shift,0.33);
        fit_scores.push(score);
        if (show_logging){console.log(`${baselines[baselines.length - 1]} : ${fit_scores[fit_scores.length - 1]}`);}
    }
    
    var best_fit_score = Math.min.apply(null,fit_scores); 
    
    for (let i = 0; i < fit_scores.length; i++)
    {
        if (best_fit_score == fit_scores[i])
        {
            return [baselines[i], fit_scores[i]];
        }
    }
    return [undefined, undefined];
}

function input_lists_exclude_haploid(log2_list, baf_list, seg_len, is_haploid){
    var log2_list_cleaned = [];
    var baf_list_cleaned = [];
    var seg_len_cleaned = [];

    for (let i = 0; i < is_haploid.length; i++)
    {
    console.log(is_haploid[i]);
        if (is_haploid[i] != 1)
        {
            log2_list_cleaned.push(log2_list[i]);
            baf_list_cleaned.push(baf_list[i]);
            seg_len_cleaned.push(seg_len[i]);
        }
    }
    
    return [log2_list_cleaned, baf_list_cleaned, seg_len_cleaned];
}

function calc_best_fit_combined(log2_list, baf_list, seg_len, cnv_cutoff, loh_cutoff, show_logging = false){

    var baselines = [];
    var cellularities = [];
    var fit_scores = [];
    
    for (let i = 0; i <= 100; i++)
    {
        var baseline_shift = i/100;
        
        var fit_results_down = calc_best_fit_cell(log2_list, baf_list, seg_len, cnv_cutoff, loh_cutoff, -1 * baseline_shift);
        baselines.push(baseline_shift * -1);
        cellularities.push(fit_results_down[0]);
        fit_scores.push(fit_results_down[1] / Math.pow(1.01-baseline_shift,0.333));
        if (show_logging){console.log(`(${baselines[baselines.length - 1]},${cellularities[cellularities.length - 1]*100}%) : ${fit_scores[fit_scores.length - 1]}`);}
        
        var fit_results_up = calc_best_fit_cell(log2_list, baf_list, seg_len, cnv_cutoff, loh_cutoff, baseline_shift);
        baselines.push(baseline_shift);
        cellularities.push(fit_results_up[0]);
        fit_scores.push(fit_results_up[1] / Math.pow(1.01-baseline_shift,0.333));
        if (show_logging){console.log(`(${baselines[baselines.length - 1]},${cellularities[cellularities.length - 1]*100}%) : ${fit_scores[fit_scores.length - 1]}`);}
    }
    
    var best_fit_score = Math.min.apply(null,fit_scores); 
    
    for (let i = 0; i < fit_scores.length; i++)
    {
        if (best_fit_score == fit_scores[i])
        {
            return [cellularities[i], baselines[i], fit_scores[i]];
        }
    }
    console.log("No best fit returned.");
    return [undefined, undefined, undefined];
}
"""
    
    constants="""
const data = source.data;

console.log(source.data);

const calc_mode = calcmode_radio.active;

console.log(data);

console.log(data['seg_baseline_adj']);

const cleaned_list_results = input_lists_exclude_haploid(data['cnv_y_stored'], data['seg_loh_height_stored'], data['seg_length'], data['seg_baseline_adj']);

console.log(cleaned_list_results);

const seg_y = cleaned_list_results[0];
const seg_loh = cleaned_list_results[1];
const seg_len = cleaned_list_results[2];



var fit_results = undefined;
"""
    code="""
    cell_calc_button.label = 'Adjusting...'
    cell_calc_button.disabled = true
    calcmode_radio.disabled = true

    setTimeout(() => {
        if (calc_mode == 0)
        {
            console.log("Recalculating in mode: Cellularity only.");
            fit_results = calc_best_fit_cell(seg_y, seg_loh, seg_len, cnv_cutoff, loh_cutoff, offset, true);
            console.log(`Best fit returned: ${(fit_results[0] * 100).toFixed(1)}% (Score: ${fit_results[1].toFixed(2)}).`);
            cell_slider.value = fit_results[0];
        } else if (calc_mode == 1){
            console.log("Recalculating in mode: Baseline only.");
            fit_results = calc_best_fit_baseline(seg_y, seg_loh, seg_len, cnv_cutoff, loh_cutoff, cell, true);
            console.log(`Best fit returned: ${(fit_results[0]).toFixed(2)} (Score: ${fit_results[1].toFixed(2)}).`);
            offset_slider.value = fit_results[0];
        } else if (calc_mode == 2){
            fit_results = calc_best_fit_combined(seg_y, seg_loh, seg_len, cnv_cutoff, loh_cutoff, true);
            console.log(`Best fit returned: Cellularity ${(fit_results[0] * 100).toFixed(1)}%, Baseline ${fit_results[1].toFixed(2)} (Score: ${fit_results[2].toFixed(2)}).`);
            cell_slider.value = fit_results[0];
            offset_slider.value = fit_results[1];
        } else {
            console.log("Unsupported mode selected (somehow). Cannot perform auto-adjustment.")
        }

        cell_calc_button.label = 'Auto-Adjustment'
        calcmode_radio.disabled = false
        cell_calc_button.disabled = false
    }, 1) 
"""
    code = "".join([sliders,functions,constants,code])
    callback = CustomJS(args=dict(
    source=source_seg,cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider,
    offset_slider=offset_slider,calcmode_radio=calcmode_radio,cell_calc_button=cell_calc_button),code=code)

    return callback

# Make ideograms

def plot_chrom_ideogram(chrom_info_dict:dict,x_range=None):
    chrom_length = int(chrom_info_dict['length'])

    stain_to_color={'gneg':'white','gpos25':'lightgrey','gpos50':'grey','gpos75':'darkgrey','gpos100':'black','acen':'lightcoral','gvar':'powderblue','stalk':'lemonchiffon'}

    cytobands=[]
    cyto_start=[]
    cyto_end=[]
    cyto_x=[]
    cyto_width=[]
    cyto_y=[]
    cyto_height=[]
    cyto_stain=[]
    cyto_color=[]

    for i in chrom_info_dict['cytobands']:
        cytobands.append(i['name'])
        cyto_start.append(i['start'])
        cyto_end.append(i['end'])
        cyto_x.append((i['start'] + i['end'])/2)
        cyto_width.append(i['end'] - i['start'])
        cyto_y.append(0.5)
        cyto_height.append(1)
        cyto_stain.append(i['stain'])
        cyto_color.append(stain_to_color[i['stain']])
    
    source_ideogram = ColumnDataSource(data=dict(
        cytoband = cytobands,
        start = cyto_start,
        end = cyto_end,
        pos_x = cyto_x,
        pos_y = cyto_y,
        width = cyto_width,
        height = cyto_height,
        stain = cyto_stain,
        color = cyto_color,
        ))

    if not x_range:
        x_range=(0,chrom_length)
        
    ideogram_fig = figure(x_range=x_range, y_range=(0,1), y_axis_label="Ideo",toolbar_location=None,tools=TOOLS,frame_height=25)
    ideogram_fig.xaxis.visible = False
    ideogram_fig.yaxis.major_label_text_color="#FFFFFE"
    ideogram_fig.yaxis.formatter=NumeralTickFormatter(format="0.00")
    ideogram_fig.rect(source=source_ideogram,x="pos_x",y="pos_y",width="width",height="height",line_color="black", fill_color="color", name="cytobands")
    
    hover_ideogram = HoverTool(tooltips=f"""
        <b>@cytoband</b><br>
        Region: @start - @end<br>
        Stain: @stain
    """,
    names=['cytobands'])
    ideogram_fig.add_tools(hover_ideogram)
    return ideogram_fig

def plot_chrom_headers(abs_pos_dict:dict, x_range=None):
    genome_length = abs_pos_dict['end']
    
    chroms = list(abs_pos_dict.keys())
    label_pos = []
    for i in range(0,len(chroms)):
        if chroms[i] == 'end':
            label_pos.append(abs_pos_dict['end']*1.1)
        else:
            label_pos.append((abs_pos_dict[chroms[i]]+abs_pos_dict[chroms[i+1]])/2)
            chroms[i] = chroms[i].replace('chr','')
    
    source_chrom_header = ColumnDataSource(data=dict(
        label=chroms,
        label_pos=label_pos,
        label_y=[0.1]*len(abs_pos_dict)
        ))
    
    if not x_range:
        x_range=(0,genome_length)
        
    chrom_header = figure(x_range=x_range, y_range=(0,1), y_axis_label="Chr",toolbar_location=None,frame_height=20, tools=TOOLS)
    chrom_header.xaxis.visible = False
    chrom_header.yaxis.major_label_text_color="#FFFFFE"
    chrom_header.yaxis.formatter=NumeralTickFormatter(format="0.00")
    chrom_header.rect(x=(genome_length/2), y=0.5, width=genome_length*1.05, height=1.1, line_color=None, fill_color="white")
    
    chrom_labels = LabelSet(x='label_pos', y='label_y', text='label',
                      source=source_chrom_header, render_mode='canvas', text_align="center")
    chrom_header.add_layout(chrom_labels)
    
    return chrom_header

# Tools to display for Bokeh plots
TOOLS="pan,xzoom_in,xzoom_out,xwheel_zoom,xbox_zoom,reset,tap,save"

### TESTING FOCALS ###

def get_focals(focals_csv,cnv_cutoff:float=0.25,sd_cutoff:float=3.0):
    focals_table = pd.read_csv(focals_csv, sep=",",header=0)
    is_focal = []
    for i in range(0,len(focals_table['event_id'])):
        gene_log2 = focals_table['weighted_average_log2_copy_ratio'][i]
        if focals_table['stdevs'][i] > sd_cutoff and abs(gene_log2) > cnv_cutoff:
            focal = 1
        else:
            focal = 0
        is_focal.append(focal)
    focals_table['focal_status']=is_focal
    return focals_table

# Big function

def generate_chrom_plot(genes_table, cnvloh_dna_dfs, chrom_info_dict, chromosome, genes_list:list=None, outpath:str=None, output_type="json", cnvloh_dna_dfs_2=None, focals_dict:dict=None, male_adjustment:bool=False):

    genelist_field = TextInput(name = "genelist_field", syncable = True, visible = False, value_input = str(genes_list))

    genes_list = None

    if not outpath:
        outpath = f"./CNVLOH_Interactive_Bokeh_Chrom_{chromosome}.{output_type}"
    default_cnv_cutoff = 0.25
    default_sd_cutoff = 2
    default_baf_cutoff = 0.42
    
    num_points = len(cnvloh_dna_dfs['cnv_df']['MIDPOINT'])
    if cnvloh_dna_dfs_2:
        num_points += len(cnvloh_dna_dfs_2['cnv_df']['MIDPOINT'])
        
    if num_points > 50000:
        slider_step = 0.05
        num_stdev = 1
    else:
        slider_step = 0.01
        num_stdev = 1
    
    ### Define sliders ###
    cnv_slider, cnv_spinner = slider_spinner(start=0, end=3, value=default_cnv_cutoff, step=slider_step, title="CNV Log2 Cutoff", name='cnv_slider') #Determines the +/- Log2 that counts as CNV
    loh_slider, loh_spinner = slider_spinner(start=0, end=0.5, value=default_baf_cutoff, step=slider_step, title="LOH BAF Cutoff", name='loh_slider') #Determines the +/- VAF that counts as LOH
    cell_slider, cell_spinner = slider_spinner(start=0.10, end=1, value=1, step=slider_step, title="Cellularity", name='cell_slider') #Adjusts everything for samples with <100% cellularity
    offset_slider, offset_spinner = slider_spinner(start=-1, end=1, value=0, step=slider_step, title="Log2 Offset", name='offset_slider') #Changes zeroing point of CNV data
    if focals_dict:
        sd_slider, sd_spinner = slider_spinner(start=1,end=8, value=default_sd_cutoff,step=0.25,title='Focal SD Cutoff', name='sd_slider')
    if cnvloh_dna_dfs_2:
        offset_slider_2, offset_slider_spinner_2 = slider_spinner(start=-1, end=1, value=0, step=slider_step, title="Offset 2", name='offset_slider_2') #Changes zeroing point of overlay CNV data
        alpha_slider, alpha_spinner = slider_spinner(start=0, end=1, value=0.5, step=slider_step, title="Alpha", name='alpha_slider') #Changes transparency of overlay data)
 
    # Colors
    cnv_colors_list = []
    if (genes_table) is not None:
        for i in genes_table['seg_mean']:
            cnv_colors_list.append(calc_col_cnv(i))
    loh_colors_list = []
    if cnvloh_dna_dfs is not None:
        for i in cnvloh_dna_dfs['seg_df']['LOH']:
            loh_colors_list.append(calc_col_loh(i))

    # Gene data table
    # This is the table made from the 'CNV_LOH.breakpoints.csv' data
    
    text_positions = []
    for i in range(0,23):
        text_positions.append(i/50)
    
    print(genes_table)
    if (genes_table) is not None:
        genes_table = genes_table[genes_table['CHR'] == chromosome].copy()
    print(genes_table)

    gene_focals_datasource = None
    exon_focals_datasource = None
    stored_focals_datasource = None
    focal_effects_datasource = None

    if focals_dict:
        gene_focals_datasource, exon_focals_datasource, stored_focals_datasource, focal_effects_datasource = focals_json_to_datasources(focals_dict, genome_info={chromosome:chrom_info_dict}, filter_to_chrom=chromosome, genes_list=genes_list,stdev_cutoff=default_sd_cutoff,log2_cutoff=default_cnv_cutoff,male_adjustment=male_adjustment)
    
    gene_table_data = make_gene_table_data(genes_table,cnv_cutoff=default_cnv_cutoff,loh_cutoff=default_baf_cutoff,genes_list=genes_list,focals_dict=focals_dict,sd_cutoff=default_sd_cutoff,focal_effects_datasource=focal_effects_datasource)
    source_genetable = ColumnDataSource(gene_table_data, name="gene_table_data")

    if (genes_table) is not None and output_type != 'html':
        exclude_columns=['seg_start','seg_end','num_mark','seg_mean','LOH','loh_allele_fraction']
        genes_table_dict = table_to_dict(table=genes_table,key_1='seg_id',key_2='ENTREZ_ID',string_keys_only=True, remove_blanks=True, exclude_columns=exclude_columns)
        if genes_list:
            add_phenotype_gene_status(genes_table_dict, genes_list)
    else:
        genes_table_dict = None

    # Saves a stored copy of the data to refer to for calculations
    gene_table_data_stored = make_gene_table_data(genes_table,genes_list=genes_list,focals_dict=focals_dict,sd_cutoff=1,focal_effects_datasource=focal_effects_datasource)
    source_genetable_stored = ColumnDataSource(gene_table_data_stored, name="gene_table_data_stored")

    # Defines a table layout to display under the plots
    gene_table_columns = [
            TableColumn(field="SYMBOL", title="Gene"),
            TableColumn(field="CHR", title="Chrom"),
            TableColumn(field="START", title="Start"),
            TableColumn(field="END", title="End"),
            TableColumn(field="seg_mean", title="Segment Log2"),
            TableColumn(field="seg_num_copy", title="Segment Copies"),
            TableColumn(field="seg_baf",title='B-Allele Freq'),
            TableColumn(field="seg_loh",title='Loss-of-Heterozygosity')
        ]
    gene_data_table = DataTable(source=source_genetable, columns=gene_table_columns, sizing_mode="stretch_both")

    ### Data sources for plots. Changing these dynamically changes the plots###
    ### Values ending in _stored are to be referred back to when recalculating things that would overwrite data.
    
    print(cnvloh_dna_dfs)
    time.sleep(3)
    
    baseline_adj_list = []
    for i in list(cnvloh_dna_dfs['seg_df']['CHROM']):
        if male_adjustment and ('x' in i.lower() or 'y' in i.lower()):
            baseline_adj_list.append(1)
        else:
            baseline_adj_list.append(0)
    
    # Overall segments, including CNV and LOH
    source_seg = create_sourceseg_datasource(cnvloh_dna_dfs=cnvloh_dna_dfs, source_name='source_seg', num_stdev=num_stdev, genes_table_dict=genes_table_dict,genome_info=chrom_info_dict, baseline_adj_list=baseline_adj_list)
    
    if cnvloh_dna_dfs_2:
        source_seg_2 = create_sourceseg_datasource(cnvloh_dna_dfs=cnvloh_dna_dfs_2, source_name='source_seg_2', num_stdev=num_stdev)

    # CNV supporting sub-segments
    #print(datetime.datetime.now())
    #print('Filtering points by stdev')
    cnv_x_array,cnv_y_array,cnv_filtered_array = filter_cnv_points(cnvloh_dna_dfs,num_stdev)
    initial_points = len(cnvloh_dna_dfs['cnv_df']['LOG2'])
    kept_points = len(cnv_x_array)
    
    if cnvloh_dna_dfs_2:
        cnv_x_array_2,cnv_y_array_2 = filter_cnv_points(cnvloh_dna_dfs_2,num_stdev)
        initial_points += len(cnvloh_dna_dfs_2['cnv_df']['LOG2'])
        kept_points += len(cnv_x_array_2)
    
    #print(f"\tPlotted {kept_points} CNV points, out of {initial_points}.")
    #print(datetime.datetime.now())
    
    source_cnv = ColumnDataSource(data={'cnv_chrom':[chromosome for i in cnv_x_array],
                                        'cnv_x':np.array(cnv_x_array),
                                        'cnv_y_stored':np.array(cnv_y_array),
                                        'cnv_y':np.array(cnv_y_array),
                                        'cnv_alpha':[0.5]*len(cnv_y_array),
                                        'coordinates': [f"{chromosome}:{cnv_x_array[i]}" for i in range(0,len(cnv_x_array))],
                                        'cnv_reads': ['']*len(cnv_x_array),
                                        'cnv_stdev': ['']*len(cnv_x_array),
                                        'other_hovertext':['']*len(cnv_x_array),
                                        'type_label': ['Coverage Point']*len(cnv_x_array),
                                        'alpha': 1-np.array(cnv_filtered_array),
                                        'is_filtered': cnv_filtered_array
                                        }
                                 )
                                 
    source_bins = ColumnDataSource(data={'cnv_bin_chrom':[chromosome for i in list(cnvloh_dna_dfs['cnv_bins_df']['START'])],
                                         'cnv_bin_pos':( np.array(list(cnvloh_dna_dfs['cnv_bins_df']['START'])) + np.array(list(cnvloh_dna_dfs['cnv_bins_df']['END'])) )/2,
                                         'cnv_bin_start':list(cnvloh_dna_dfs['cnv_bins_df']['START']),
                                         'cnv_bin_end':list(cnvloh_dna_dfs['cnv_bins_df']['END']),
                                         'cnv_bin_width':( np.array(list(cnvloh_dna_dfs['cnv_bins_df']['END'])) - np.array(list(cnvloh_dna_dfs['cnv_bins_df']['START'])) ),
                                         'cnv_bin_sd': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV'])),
                                         'cnv_bin_height': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV']))*2,
                                         'cnv_y': list(cnvloh_dna_dfs['cnv_bins_df']['MEAN']),
                                         'cnv_bin_sd_stored': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV'])),
                                         'cnv_bin_sd': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV'])),
                                         'cnv_bin_height_stored': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV']))*2,
                                         'cnv_y_stored': list(cnvloh_dna_dfs['cnv_bins_df']['MEAN']),
                                         'cnv_bin_points':list(cnvloh_dna_dfs['cnv_bins_df']['POINTS']),
                                         'cnv_bin_reads': list(cnvloh_dna_dfs['cnv_bins_df']['READS']),
                                         'cnv_bin_depth': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['READS']))*150 / (np.array(list(cnvloh_dna_dfs['cnv_bins_df']['END'])) - np.array(list(cnvloh_dna_dfs['cnv_bins_df']['START']))),
                                         'alpha':[1 for i in list(cnvloh_dna_dfs['cnv_bins_df']['START'])]
                                        }
                                    )

    # LOH supporting points
    source_loh = ColumnDataSource(data={'loh_x':np.array(cnvloh_dna_dfs['loh_df']['POS']),
                                        'loh_y_stored':np.array(cnvloh_dna_dfs['loh_df']['VAF']),
                                       'loh_y':np.array(cnvloh_dna_dfs['loh_df']['VAF']),
                                       'color':['gray']*len(cnvloh_dna_dfs['loh_df']['POS']),
                                       'loh_alpha':[0.5]*len(cnvloh_dna_dfs['loh_df']['POS'])
                                       }
                                 )

    if cnvloh_dna_dfs_2:
        source_cnv_2 = ColumnDataSource(data={'cnv_x':np.array(cnv_x_array_2),
                                        'cnv_y_stored':np.array(cnv_y_array_2),
                                        'cnv_y':np.array(cnv_y_array_2),
                                        'cnv_alpha':[0.5]*len(cnv_y_array_2)}
                                 )
        source_loh_2 = ColumnDataSource(data={'loh_x':np.array(cnvloh_dna_dfs_2['loh_df']['POS']),
                                        'loh_y_stored':np.array(cnvloh_dna_dfs_2['loh_df']['VAF']),
                                       'loh_y':np.array(cnvloh_dna_dfs_2['loh_df']['VAF']),
                                       'color':['gray']*len(cnvloh_dna_dfs_2['loh_df']['POS']),
                                       'loh_alpha':[0.5]*len(cnvloh_dna_dfs_2['loh_df']['POS'])}
                                 )
    
    # Provides the lists of colors for JS to recalculate them when values change
    source_colors = ColumnDataSource(data={'cnv_colors':cnv_col_list,
                                           'loh_colors':loh_col_list},name='colors_source')

    ### Plotting ###
    # Get genome info for chromosome to plot centromeres, etc.
    chrom_info = chrom_info_dict

    ## CNV ##
    
    # Setup figure canvas
    chrom_padding = int(chrom_info['length']*0.01)
    cnv = figure(tools=TOOLS,y_axis_label="CNV Log2 Ratio", x_range=Range1d(0-chrom_padding,chrom_info['length']+chrom_padding), y_range=Range1d(-2,2), frame_height=300, name="cnv_plot")
    cnv.yaxis.formatter=NumeralTickFormatter(format="+0.0")
    cnv.xaxis.visible = False
    # Plot zero-line and centromeres
    cnv.add_layout(BoxAnnotation(left=int(chrom_info['centromere_start']), right=int(chrom_info['centromere_end']), fill_alpha=0.3, fill_color='grey',level="underlay"))
    cnv.line([0,chrom_info['length']],[0,0],color="black")

    # Plot confidence interval
    #cnv.rect(source=source_seg,x='seg_loh_x',
    #    width='seg_loh_width',
    #    y='seg_y',
    #    height='seg_y_height',
    #    fill_color="cornflowerblue",fill_alpha=1,line_color="cornflowerblue")
    #if cnvloh_dna_dfs_2:
    #    cnv.rect(source=source_seg_2,x='seg_loh_x',
    #        width='seg_loh_width',
    #        y='seg_y',
    #        height='seg_y_height',
    #        fill_color="goldenrod",line_color="goldenrod",fill_alpha="seg_alpha",line_alpha="seg_alpha")
    # Plot points
    hover_cnvpoints = HoverTool(tooltips="""
        <b>@coordinates</b><br>
        Log2: @cnv_y{0.00}<br>
        Reads: @cnv_reads
        """,
    names=['cnv_points1','cnv_points2'])
    hover_cnvbins = HoverTool(tooltips="""
        <b>@cnv_bin_chrom: @cnv_bin_start - @cnv_bin_end</b><br>
        Log2: @cnv_y{0.00}<br>
        St.Dev: @cnv_bin_sd{0.00}<br>
        Points: @cnv_bin_points{0}<br>
        Reads: @cnv_bin_reads{0}<br>
        Cov. Depth: @cnv_bin_depth{0.00}x
        """,
    names=['cnv_bins1','cnv_bins2'])
    hover_cnvsegs = HoverTool(tooltips="""
        <b><u>@cnvloh_status</u><br>
        @seg_chrom: @seg_x0 - @seg_x1</b><br>
        <i>Segment @seg_index_local of @seg_per_chrom</i><br>
        Length: @seg_length{0}<br>
        Log2: @cnv_y{0.00}<br>
        Copies: @cnv_num_copy{0.0}<br>
        St. Dev.: @cnv_stdev{0.00}<br>
        Points: @cnv_points{0}<br>
        Reads: @cnv_reads{0}<br>
        Cov. Depth: @cnv_cov_depth{0.00}x
        """,
    names=['cnv_segs1','cnv_segs2'])
    
    cnv.scatter(source=source_cnv,x="cnv_x",y="cnv_y",color="#88CCEE",alpha="alpha", name="cnv_points1")
    cnv.rect(source=source_bins, x='cnv_bin_pos', width='cnv_bin_width', height='cnv_bin_height', y='cnv_y', fill_color = "dodgerblue", alpha="alpha", name="cnv_bins1")
    
    if cnvloh_dna_dfs_2:
        cnv.scatter(source=source_cnv_2,x="cnv_x",y="cnv_y",color="#CCBB44",alpha="cnv_alpha", name="cnv_points2")
        cnv.rect(source=source_bins_2, x='cnv_bin_pos', width='cnv_bin_width', height='cnv_bin_height', y='cnv_y', fill_color = "gold",name="cnv_bins2")
    
    # Plot segments
    cnv.segment(source=source_seg,x0="seg_x0",x1="seg_x1",
                y0="cnv_y_display",y1="cnv_y_display",
               color="red",line_width=5,name="cnv_segs1")
    if cnvloh_dna_dfs_2:
        cnv.segment(source=source_seg_2,x0="seg_x0",x1="seg_x1",
                y0="cnv_y_display",y1="cnv_y_display",
               color="blue",line_width=5,alpha="seg_alpha",name="cnv_segs2")
               
    cnv_labels = LabelSet(x='seg_loh_x', y='cnv_label_pos', text='cnv_label_text',
                      source=source_seg, render_mode='canvas', text_align="center", text_alpha="cnv_label_alpha")
    cnv.add_layout(cnv_labels)
    
    # Plot focal events
    if focals_dict:
        # Filter to limit plotting of focal events to only show phenotype genes
        pheno_focal_filter = CustomJSFilter(args=dict(genelist_field=genelist_field), code='''
            const indices = [];

            const genelist_string = genelist_field.value_input;
            var genes_set = undefined;
        
            if (!(!genelist_string || genelist_string == '' || genelist_string == '[]'))
            {
                genes_set = new Set(genelist_string.replace('[','').replace(']','').split(','));
            }
            console.log(genes_set);
            
            // For each row in the genes table, see if the Entrez Id is in our genes_list.
            // If so, or if there is no genes_list, it passes the filter.
            for (let i = 0; i < source.get_length(); i++) {
                if (!genes_set || (source.data['entrez_id'][i] !== null && genes_set.has(source.data['entrez_id'][i].toString()))) {
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }

            return indices;
        ''')

        focal_sd_filter = CustomJSFilter(args=dict(sd_slider=sd_slider), code='''
            const indices = [];
            const focal_sd_cutoff = sd_slider.value;

            // For each row in the genes table, see if the focal sd absolute value is over our cutoff.
            // If so, or if there is no cutoff, it passes the filter.
            for (let i = 0; i < source.get_length(); i++) {
                if (!focal_sd_cutoff || Math.abs(source.data['delta_sd'][i]) >= focal_sd_cutoff) {
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }

            return indices;
        ''')

        gene_focals_view = CDSView(source=gene_focals_datasource, filters=[pheno_focal_filter, focal_sd_filter])
        exon_focals_view = CDSView(source=exon_focals_datasource, filters=[pheno_focal_filter, focal_sd_filter])
        
        cnv.scatter(source=gene_focals_datasource, view=gene_focals_view, x='event_midpoint',y='display_y',color='color',marker="triangle", size=12, name="focal_genes_scatter")
        cnv.scatter(source=exon_focals_datasource, view=exon_focals_view, x='event_midpoint',y='display_y',color='color',marker="inverted_triangle",size=8, name="focal_exons_scatter")
        hover_focal_gene = HoverTool(tooltips="""
            <b><u>@symbol (@entrez_id)</u><br>
            @text_status</b><br>
            @chrom: @start - @end<br>
            Log2: @{log2}{0.00} (@{delta_log2}{+0.0})<br>
            Copies: @{focal_copies}{0.0} (@{delta_copies}{+0.0})<br>
            Δ St.Dev.: @{delta_sd}{+0.00}<br>
            """,
            names=['focal_genes_scatter'])
        hover_focal_exon = HoverTool(tooltips="""
            <b><u>@symbol (@entrez_id)</u><br>
            @text_status</b><br>
            @chrom: @start - @end<br>
            Log2: @{log2}{0.00} (@{delta_log2}{+0.0})<br>
            Copies: @{focal_copies}{0.0} (@{delta_copies}{+0.0})<br>
            Δ St.Dev.: @{delta_sd}{+0.00}<br>
            @transcript_text<br>
            # Trans. Affected: @transcripts_affected<br>
            """,
            names=['focal_exons_scatter'])
        cnv.add_tools(hover_focal_gene)
        cnv.add_tools(hover_focal_exon)
    
    #cnv.add_tools(hover_cnvpoints)
    cnv.add_tools(hover_cnvbins)
    cnv.add_tools(hover_cnvsegs)
    
    ## Assist Plot ##
    
    # Setup figure canvas
    helper = figure(tools=TOOLS,y_axis_label="Affected Genes",x_range=cnv.x_range,y_range=Range1d(-0.6,0.6), frame_height=250,name='helper_plot')
    helper.yaxis.formatter=NumeralTickFormatter(format="+0.0")
    helper.xaxis.visible = False
    helper.yaxis.major_label_text_color="#FFFFFE"
    
    # Plot zero-line centromeres
    helper.add_layout(BoxAnnotation(left=int(chrom_info['centromere_start']), right=int(chrom_info['centromere_end']), fill_alpha=0.3, fill_color='grey',level="underlay"))
    helper.line([0,chrom_info['length']],[0,0],color="black")
    #helper.line([0,chrom_info['length']],[0.2,0.2],color="gray", line_dash="dashed")
    #helper.line([0,chrom_info['length']],[-0.2,-0.2],color="gray", line_dash="dashed")

    # Filter to limit helper plot to only show phenotype genes
    pheno_gene_filter = CustomJSFilter(args=dict(genes_list=genes_list), code='''
        const indices = [];

        // For each row in the genes table, see if the Entrez Id is in our genes_list.
        // If so, or if there is no genes_list, it passes the filter.
        for (let i = 0; i < source.get_length(); i++){
            if (!genes_list || genes_list.includes(source.data['ENTREZ_ID'][i].toString())) {
                indices.push(true);
            } else {
                indices.push(false);
            }
        }

        return indices;
    ''')
    custom_view = CDSView(source=source_genetable, filters=[pheno_gene_filter])

    # Plot affected genes in appropriate colors (blue = gain, red = loss)
    helper.rect(
        source=source_genetable,
        view=custom_view,
        x='gene_x',
        width='gene_width',
        y=0.3,
        height=0.6,
        fill_color="gene_cnv_color",
        fill_alpha=0.95,
        line_alpha=1,
        line_color="gene_cnv_color",
        name="gene_rect_cnv"
    )
             
    helper.rect(
        source=source_genetable,
        view=custom_view,
        x='gene_x',
        width='gene_width',
        y=-0.3,
        height=0.6,
        fill_color="gene_loh_color",
        fill_alpha=0.95,
        line_alpha=1,
        line_color="gene_loh_color",
        name="gene_rect_loh"
    )
            
    helper.rect(
        source=source_genetable,
        view=custom_view,
        x='gene_x',
        width='gene_width',
        y='gene_y',
        height='gene_height',
        fill_color="gene_cnv_color",
        fill_alpha=0,
        line_alpha=0,
        line_color="gene_loh_color",
        name="gene_rect_cnvloh"
    )

    hover_helper_rects = HoverTool(tooltips="""
            <b>@SYMBOL (@ENTREZ_ID)</b><br>
            @CHR: @START - @END<br>
            Log2: @gene_log2{0.00}<br>
            BAF: @seg_baf{0.00}<br>
            LOH: @seg_loh{0.0%}<br>
            """,
            names=['gene_rect_cnvloh'])

    # Plot labels of affected genes
    # This doesn't really work right yet, needs better positioning
    labels = LabelSet(
        x='gene_x',
        y='gene_y_label',
        text='SYMBOL',
        source=source_genetable,
        render_mode='canvas',
        text_align="center",
        text_alpha="label_alpha"
    )
    helper.add_layout(labels)
    helper.add_tools(hover_helper_rects)

    ## LOH ##
    
    # Setup figure canvass
    loh = figure(tools=TOOLS,y_axis_label="Variant Allele Frequency",x_range=cnv.x_range,y_range=Range1d(0, 1), frame_height=350)
    loh.yaxis.formatter=NumeralTickFormatter(format="0.00")
    # Plot zero-line and centromeres
    loh.add_layout(BoxAnnotation(left=int(chrom_info['centromere_start']), right=int(chrom_info['centromere_end']), fill_alpha=0.3, fill_color='grey',level="underlay"))
    loh.line([0,chrom_info['length']],[0,0],color="black")

    # Add box showing area of LOH based on center of 90% CI
    loh.rect(source=source_seg,x='seg_loh_x',
             width='seg_loh_width',
             y='seg_loh_y',
             height='seg_loh_height',
             fill_color="seg_loh_color",fill_alpha=1,line_color="seg_loh_line",name="seg_loh_rect1")
    if cnvloh_dna_dfs_2:
        loh.rect(source=source_seg_2,x='seg_loh_x',
             width='seg_loh_width',
             y='seg_loh_y',
             height='seg_loh_height',
             fill_color="seg_loh_color",fill_alpha="seg_alpha",line_color="seg_loh_line",line_alpha="seg_alpha",name="seg_loh_rect2")
    # Add variants
    loh.scatter(source=source_loh,x="loh_x", y="loh_y",color="black")
    if cnvloh_dna_dfs_2:
        loh.scatter(source=source_loh,x="loh_x", y="loh_y",color="black",alpha="loh_alpha")
        
    hover_lohrects = HoverTool(tooltips="""
        <b><u>@cnvloh_status</u><br>
        @seg_chrom: @seg_x0 - @seg_x1</b><br>
        <i>Segment @seg_index_local of @seg_per_chrom</i><br>
        Length: @seg_length<br>
        BAF: @seg_loh_baf{0.00}<br>
        LOH: @seg_loh_height{0.0%}<br>
        """,
        names=['seg_loh_rect1','seg_loh_rect2'])
    loh.add_tools(hover_lohrects)
    
    ### Source-updating JS ###
    
    # I don't know if I'm able to insert comments inside the JS code
    
    ## Seg+CNV+Loh Callback
    # Values that are entirely recalculated from other values don't need to be 'stored'
    # Values that are adjustments of an initial value need to be recalculated from a stored version
    callback = generate_callbackJS(source_seg=source_seg, source_cnv=source_cnv, source_cnv_bins=source_bins, source_loh=source_loh, source_colors=source_colors,
        cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider,offset_slider=offset_slider, alpha_slider=None)
    if cnvloh_dna_dfs_2:
        callback_2 = generate_callbackJS(source_seg=source_seg_2, source_cnv=source_cnv_2, source_loh=source_loh_2, source_colors=source_colors,
            cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider,offset_slider=offset_slider_2, alpha_slider=alpha_slider)

    point_bin_radio = RadioButtonGroup(labels=['Points','Bins','Both'], active=2)
    callback_cnv_radio = CustomJS(args=dict(source_points=source_cnv, source_bins=source_bins, radio_buttons = point_bin_radio, genelist_field = genelist_field),
                        code="""
    const button_state = radio_buttons.active;
    console.log('CNV Button State:' + button_state);
    
    const data_points = source_points.data;
    const data_bins = source_bins.data;
    
    for (var i = 0; i < data_points['alpha'].length; i++ ){
        if (button_state == 0){data_points['alpha'][i] = 1;}
        else if (button_state == 2){data_points['alpha'][i] = 1 - data_points['is_filtered'][i]}
        else {data_points['alpha'][i] = 0;}
    }
    for (var i = 0; i < data_bins['alpha'].length; i++ ){
        if (button_state == 1 || button_state == 2){data_bins['alpha'][i] = 1;}
        else {data_bins['alpha'][i] = 0;}
    }
    source_points.change.emit();
    source_bins.change.emit();
                        """)
    point_bin_radio.js_on_click(callback_cnv_radio)
    
    # Radio button to toggle gene labels in the helper plot on and off
    gene_label_radio = RadioButtonGroup(labels=['Hide Labels','Show Labels'], active=1)
    
    # Radio button to toggle between showing only phenotype genes vs all genes
    all_genes_radio = RadioButtonGroup(labels=['Phenotype Genes','All Genes'], active=0)
    callback_all_genes_radio = CustomJS(args=dict(
        source_genetable=source_genetable,
        source_genetable_stored=source_genetable_stored,
        all_genes_radio=all_genes_radio,
        custom_view=custom_view,
        pheno_gene_filter=pheno_gene_filter,
        gene_focals_view=gene_focals_view,
        exon_focals_view=exon_focals_view,
        pheno_focal_filter=pheno_focal_filter,
        focal_sd_filter=focal_sd_filter
    ), code="""
        const all_genes_selected_button = all_genes_radio.active;

        if (all_genes_selected_button == 0) {
            // Show phenotype genes
            custom_view.filters = [pheno_gene_filter]
            gene_focals_view.filters = [pheno_focal_filter, focal_sd_filter]
            exon_focals_view.filters = [pheno_focal_filter, focal_sd_filter]
        } else {
            // Show all genes
            custom_view.filters = []
            
            gene_focals_view.filters = [focal_sd_filter]
            exon_focals_view.filters = [focal_sd_filter]
        }
    """)
    all_genes_radio.js_on_click(callback_all_genes_radio)

    callback_update_label_transparency = CustomJS(args=dict(
            source_genetable=source_genetable, 
            source_genetable_stored=source_genetable_stored, 
            gene_label_radio=gene_label_radio, 
            all_genes_radio=all_genes_radio,
            genes_list=genes_list
        ), code="""
        const gene_label_radio_state = gene_label_radio.active;
        const all_genes_radio_state = all_genes_radio.active;
        
        const data_points = source_genetable.data;
        const stored_data_points = source_genetable_stored.data;
        
        // console.log('data_points')
        // console.log(data_points)
        // console.log('stored_data_points')
        // console.log(stored_data_points)
        // console.log('genes_list')
        // console.log(genes_list)

        for (var i = 0; i < data_points['label_alpha'].length; i++) {
            if (gene_label_radio_state == 0) {
                // Gene labels are off. Set alpha to 0 on all labels.
                data_points['label_alpha'][i] = 0
                stored_data_points['label_alpha'][i]=0
            } else {
                // Gene labels are on.
                if (all_genes_radio_state === 1) {
                    // We are displaying all genes, so set all the label alphas to 1
                    data_points['label_alpha'][i] = 1
                    stored_data_points['label_alpha'][i] = 1
                } else {
                    // We are displaying phenotype genes.  Label is displayed if either 
                    // - We have no genes_list ("All Genes" in Gene List dropdown)
                    // - Entrez Id is in the provided genes_list

                    const show_label = (!genes_list || genes_list.includes(data_points['ENTREZ_ID'][i].toString()))
                    const label_alpha = show_label ? 1 : 0

                    data_points['label_alpha'][i] = label_alpha
                    stored_data_points['label_alpha'][i] = label_alpha
                }
            }
        }
        source_genetable.change.emit();
    """)
    # Both the gene label on/off and all genes/pheno genes radio buttons need to update the gene label states
    gene_label_radio.js_on_click(callback_update_label_transparency)
    all_genes_radio.js_on_click(callback_update_label_transparency)

    ## Gene data table callback ##
    # Color-coding for Gain vs. Loss is calculated dynamically and obeys cellularity and centering changes
    # LOH plotting not yet implemented until I get an example of current data formatting
    callback_table = generate_callbackJS_gene_table(source_genetable=source_genetable, source_genetable_stored=source_genetable_stored, source_focal_effects=focal_effects_datasource,
        source_colors=source_colors, cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider,offset_slider=offset_slider, genelist_field=genelist_field)

    ## Gene list callback ##
    # Changes the 'alpha' value for labels in the stored table, which is used to determine if genes were initially listed or not.
    callback_genelist = generate_callbackJS_genelist(source_genetable_stored=source_genetable_stored, genelist_field=genelist_field)

    ### Slider callbacks ###
    # Note: You don't want to trigger recalculation for sliders that don't actually do anything.
    # I think each blip of movement recalculates the entire data source.
    # Can probably make things faster by reducing the number of steps in the sliders where reasonable.
    
    # Callback that alters Seg+CNV+LOH source
    if focals_dict:
        callback_focals = generate_callbackJS_focals(stored_focals_datasource, gene_focals_datasource, exon_focals_datasource, focal_effects_datasource, cnv_slider, cell_slider, offset_slider, sd_slider, genelist_field)

    ideogram = plot_chrom_ideogram(chrom_info_dict,x_range=cnv.x_range)

    ### Writing out figures ###
    
    ## Setup final layout ##
    # Grid for the three stacked plots, CNV+Helper+LOH
    #ideo = gridplot([[ideogram]],plot_width=800,plot_height=25)
    grid = gridplot([[ideogram],[cnv],[helper],[loh]], plot_height=900, sizing_mode="stretch_both",toolbar_options={"logo":None})

    # Sliders go above plots, table goes below
    
    if focals_dict:
        layout = column(
            row(
                cnv_slider, cnv_spinner,
                loh_slider, loh_spinner,
                cell_slider, cell_spinner,
                offset_slider, offset_spinner,
                sd_slider, sd_spinner
            ),
            row(point_bin_radio, gene_label_radio), #, all_genes_radio),
                grid,
            sizing_mode="stretch_both"
        )
    elif cnvloh_dna_dfs_2:
        layout = column(
            row(
                cnv_slider, cnv_spinner,
                loh_slider, loh_spinner,
                cell_slider, cell_spinner,
                offset_slider, offset_spinner,
                offset_slider_2, offset_spinner_2,
                alpha_slider, alpha_spinner
            ),
            row(point_bin_radio, gene_label_radio), #all_genes_radio),
                grid,
            sizing_mode="stretch_both"
        )
    else:
        layout = column(
            row(
                cnv_slider, cnv_spinner,
                loh_slider, loh_spinner,
                cell_slider, cell_spinner,
                offset_slider, offset_spinner,
            ),
            row(point_bin_radio, gene_label_radio), # all_genes_radio),
                grid,
            sizing_mode="stretch_both"
        )
    
    # Callback that alters Seg+CNV+LOH source
    if focals_dict:
        cell_slider.js_on_change('value',callback_focals)
        offset_slider.js_on_change('value',callback_focals)
        cnv_slider.js_on_change('value',callback_focals)
        sd_slider.js_on_change('value',callback_focals)
    cnv_slider.js_on_change('value', callback)
    loh_slider.js_on_change('value', callback)
    cell_slider.js_on_change('value',callback)
    offset_slider.js_on_change('value', callback)
    if cnvloh_dna_dfs_2:
        cell_slider.js_on_change('value',callback_2)
        offset_slider_2.js_on_change('value',callback_2)
        alpha_slider.js_on_change('value',callback_2)

    # Callback that alters Gene Table source
    genelist_field.js_on_change('value_input', callback_genelist)
    genelist_field.js_on_change('value_input', callback_focals)
    genelist_field.js_on_change('value_input', callback_table)
    
    cnv_slider.js_on_change('value', callback_table)
    loh_slider.js_on_change('value', callback_table)
    cell_slider.js_on_change('value',callback_table)
    offset_slider.js_on_change('value', callback_table)
    if focals_dict:
        sd_slider.js_on_change('value',callback_table)

    ### Writing out figures ###

    # Name and write out file
    if output_type == "html":
        output_file(os.path.join(outpath),title=f"Chromosome {chromosome}: Interactive CNV-LOH Display")
        save(layout)
    else:
        with open(outpath, 'w') as output_handle:
            output_handle.write(json.dumps(json_item(layout)))

    return None

# Big function -- genome ver.
def generate_genome_plot(genes_table, cnvloh_dna_dfs, genome_info, genes_list:list=None, outpath:str=None, output_type="json", cnvloh_dna_dfs_2=None, focals_dict:dict=None, male_adjustment:bool=False):
    
    genelist_field = TextInput(name = "genelist_field", syncable = True, visible = False, value_input = str(genes_list))
    
    genes_list = None
    
    plot_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    
    for i in list(genome_info.keys()):
        if i not in plot_chroms:
            genome_info.pop(i)
    
    if not outpath:
        outpath = f"./CNVLOH_Interactive_Bokeh_GenomeView.{output_type}"
    print("Generating genome plot")
    
    default_cnv_cutoff = 0.25
    default_sd_cutoff = 2
    default_baf_cutoff = 0.42
    
    num_points = len(cnvloh_dna_dfs['cnv_df']['MIDPOINT'])
    if cnvloh_dna_dfs_2:
        num_points += len(cnvloh_dna_dfs_2['cnv_df']['MIDPOINT'])
    
    if num_points > 50000:
        slider_step = 0.05
        num_stdev = 1
    else:
        slider_step = 0.01
        num_stdev = 1
    
    cnv_x_array = []
    cnv_y_array = []
    
    abs_position_dict = {}
    chrom_list = []
    length_list = []
    baseline_adj_list = []
    
    for i in genome_info:
        chrom_list.append(i)
        length_list.append(genome_info[i]['length'])
    chrom_list.append('end')
    length_list.append(0)
    
    for i in range(0,len(chrom_list)):
        chrom = chrom_list[i]
        abs_position_dict[chrom] = sum(list(length_list)[:i])
    
    abs_adj_list = []
    
    for i in range(0,len(cnvloh_dna_dfs['seg_df']['CHROM'])):
        current_chrom = cnvloh_dna_dfs['seg_df']['CHROM'][i]
        abs_adj_list.append(abs_position_dict[current_chrom])
        if male_adjustment is True and ('x' in current_chrom.lower() or 'y' in current_chrom.lower()):
            print(f'Male adjustment added to {current_chrom}')
            baseline_adj = 1
        else:
            baseline_adj = 0
        baseline_adj_list.append(baseline_adj)
    
    if cnvloh_dna_dfs_2:
        abs_adj_list_2 = []
        for i in range(0,len(cnvloh_dna_dfs_2['seg_df']['CHROM'])):
            abs_adj_list_2.append(abs_position_dict[cnvloh_dna_dfs_2['seg_df']['CHROM'][i]])
            

    print("\tCreating data tables")

    ### Data sources for plots. Changing these dynamically changes the plots###
    
    # Focals data
    
    if focals_dict:
        gene_focals_datasource, exon_focals_datasource, stored_focals_datasource, focal_effects_datasource = focals_json_to_datasources(focals_dict, genome_info=genome_info, genes_list=genes_list, abs_position_dict=abs_position_dict, stdev_cutoff=default_sd_cutoff, log2_cutoff=default_cnv_cutoff, male_adjustment=male_adjustment)
    else:
        gene_focals_datasource = None
        exon_focals_datasource = None
        stored_focals_datasource = None
        focal_effects_datasource = None

    # Gene data table
    # This is the table made from the 'CNV_LOH.breakpoints.csv' data
    gene_table_data = make_gene_table_data(genes_table,cnv_cutoff=default_cnv_cutoff,loh_cutoff=default_baf_cutoff,genes_list=genes_list,focal_effects_datasource=focal_effects_datasource)
    source_genetable = ColumnDataSource(gene_table_data, name="gene_table_data")

    if (genes_table) is not None and output_type != 'html':
        exclude_columns=['seg_start','seg_end','num_mark','seg_mean','LOH','loh_allele_fraction']
        genes_table_dict = table_to_dict(table=genes_table,key_1='seg_id',key_2='ENTREZ_ID',string_keys_only=True,remove_blanks=True, exclude_columns=exclude_columns)
        if genes_list:
            add_phenotype_gene_status(genes_table_dict, genes_list)
    else:
        genes_table_dict = None

    # Saves a stored copy of the data to refer to for calculations
    gene_table_data_stored = make_gene_table_data(genes_table)
    source_genetable_stored = ColumnDataSource(gene_table_data_stored, name="gene_table_data_stored")

    # Defines a table layout to display under the plots
    #gene_table_columns = [
    #        TableColumn(field="SYMBOL", title="Gene"),
    #        TableColumn(field="CHR", title="Chrom"),
    #        TableColumn(field="START", title="Start"),
    #        TableColumn(field="END", title="End"),
    #        TableColumn(field="seg_mean", title="Segment Log2")
    #    ]
    #gene_data_table = DataTable(source=source_genetable, columns=gene_table_columns, sizing_mode="stretch_width")
    
    ### Values ending in _stored are to be referred back to when recalculating things that would overwrite data.
    
    print("\tCreating segment source")

    # Overall segments, including CNV and LOH
    
    source_seg = create_sourceseg_datasource(cnvloh_dna_dfs=cnvloh_dna_dfs, source_name='source_seg', num_stdev=num_stdev, abs_adj_list=abs_adj_list, genes_table_dict=genes_table_dict, genome_info=genome_info, baseline_adj_list=baseline_adj_list)
    cnv_max = max(np.array(cnvloh_dna_dfs['seg_df']['LOG2'].astype(float)))
    cnv_min = min(np.array(cnvloh_dna_dfs['seg_df']['LOG2'].astype(float)))

    #print(cnvloh_dna_dfs['cnv_bins_df'])

    source_bins = ColumnDataSource(data={'cnv_bin_chrom':list(cnvloh_dna_dfs['cnv_bins_df']['CHROM']),
                                         'cnv_bin_pos':( np.array(list(cnvloh_dna_dfs['cnv_bins_df']['START'])) + np.array(list(cnvloh_dna_dfs['cnv_bins_df']['END'])) )/2 + np.array([abs_position_dict[c] for c in list(cnvloh_dna_dfs['cnv_bins_df']['CHROM'])]),
                                         'cnv_bin_start':list(cnvloh_dna_dfs['cnv_bins_df']['START']),
                                         'cnv_bin_end':list(cnvloh_dna_dfs['cnv_bins_df']['END']),
                                         'cnv_bin_width':( np.array(list(cnvloh_dna_dfs['cnv_bins_df']['END'])) - np.array(list(cnvloh_dna_dfs['cnv_bins_df']['START'])) ),
                                         'cnv_bin_sd': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV'])),
                                         'cnv_bin_height': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV']))*2,
                                         'cnv_y': list(cnvloh_dna_dfs['cnv_bins_df']['MEAN']),
                                         'cnv_bin_sd_stored': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV'])),
                                         'cnv_bin_sd': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV'])),
                                         'cnv_bin_height_stored': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['STDEV']))*2,
                                         'cnv_y_stored': list(cnvloh_dna_dfs['cnv_bins_df']['MEAN']),
                                         'cnv_bin_points':list(cnvloh_dna_dfs['cnv_bins_df']['POINTS']),
                                         'cnv_bin_reads': list(cnvloh_dna_dfs['cnv_bins_df']['READS']),
                                         'cnv_bin_depth': np.array(list(cnvloh_dna_dfs['cnv_bins_df']['READS']))*150 / (np.array(list(cnvloh_dna_dfs['cnv_bins_df']['END'])) - np.array(list(cnvloh_dna_dfs['cnv_bins_df']['START'])))
                                        }
                                    )

    if cnvloh_dna_dfs_2:
        source_seg_2 = create_sourceseg_datasource(cnvloh_dna_dfs_2, 'source_seg_2', num_stdev, abs_adj_list_2, fill_alpha=0.5)
        cnv_max_2 = max(np.array(cnvloh_dna_dfs_2['seg_df']['LOG2'].astype(float)))
        cnv_min_2 = min(np.array(cnvloh_dna_dfs_2['seg_df']['LOG2'].astype(float)))
        
        cnv_max = max(cnv_max,cnv_max_2)
        cnv_min = min(cnv_min,cnv_min_2)

    # Provides the lists of colors for JS to recalculate them when values change
    source_colors = ColumnDataSource(data={'cnv_colors':cnv_col_list,
                                           'loh_colors':loh_col_list})

    ### Define sliders ###
    cnv_slider, cnv_spinner = slider_spinner(start=0, end=3, value=default_cnv_cutoff, step=slider_step, title="CNV Log2 Cutoff", name='cnv_slider') #Determines the +/- Log2 that counts as CNV
    loh_slider, loh_spinner = slider_spinner(start=0, end=0.5, value=default_baf_cutoff, step=slider_step, title="LOH BAF Cutoff", name='loh_slider') #Determines the +/- VAF that counts as LOH
    cell_slider, cell_spinner = slider_spinner(start=0.10, end=1, value=1, step=slider_step, title="Cellularity", name='cell_slider') #Adjusts everything for samples with <100% cellularity
    offset_slider, offset_spinner = slider_spinner(start=-1, end=1, value=0, step=slider_step, title="Log2 Offset", name='offset_slider') #Changes zeroing point of CNV data
    if focals_dict:
        sd_slider, sd_spinner = slider_spinner(start=1,end=8, value=default_sd_cutoff,step=0.5,title='Focal SD Cutoff', name='sd_slider')
        callback_focals = generate_callbackJS_focals(stored_focals_datasource, gene_focals_datasource, exon_focals_datasource, focal_effects_datasource, cnv_slider, cell_slider, offset_slider, sd_slider, genelist_field)
    if cnvloh_dna_dfs_2:
        offset_slider, offset_spinner = slider_spinner(start=-1, end=1, value=0, step=slider_step, title="Offset 1", name='offset_slider') #Changes zeroing point of CNV data (auto-centering usually handles it though)
        offset_slider_2, offset_spinner_2 = slider_spinner(start=-1, end=1, value=0, step=slider_step, title="Offset 2", name='offset_slider_2') #Changes zeroing point of overlay CNV data
        alpha_slider, alpha_spinner = slider_spinner(start=0, end=1, value=0.5, step=slider_step, title="Alpha", name='alpha_slider') #Changes transparency of overlay data)
    
    ### Source-updating JS ###
    
    # I don't know if I'm able to insert comments inside the JS code
    
    ## Seg+CNV+Loh Callback
    # Values that are entirely recalculated from other values don't need to be 'stored'
    # Values that are adjustments of an initial value need to be recalculated from a stored version
    callback=generate_callbackJS(source_seg=source_seg, source_cnv=None, source_cnv_bins=source_bins, source_loh=None, source_colors=source_colors, cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider, offset_slider=offset_slider, alpha_slider=None, genome=True)
    if cnvloh_dna_dfs_2:
        callback_2=generate_callbackJS(source_seg_2, source_cnv=None, source_cnv_bins=None, source_loh=None, source_colors=source_colors, cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider, offset_slider=offset_slider_2, alpha_slider=alpha_slider, genome=True)

    ## Gene data table callback ##
    # Color-coding for Gain vs. Loss is calculated dynamically and obeys cellularity and centering changes
    callback_table = generate_callbackJS_gene_table(source_genetable=source_genetable, source_genetable_stored=source_genetable_stored, source_focal_effects=focal_effects_datasource, 
        source_colors=source_colors, cnv_slider=cnv_slider, loh_slider=loh_slider, cell_slider=cell_slider,offset_slider=offset_slider,genelist_field=genelist_field)

    ## Gene list callback ##
    # Changes the 'alpha' value for labels in the stored table, which is used to determine if genes were initially listed or not.
    callback_genelist = generate_callbackJS_genelist(source_genetable_stored=source_genetable_stored, genelist_field=genelist_field)

    ### Slider callbacks ###
    # Note: You don't want to trigger recalculation for sliders that don't actually do anything.
    # I think each blip of movement recalculates the entire data source.
    # Can probably make things faster by reducing the number of steps in the sliders where reasonable.
    
    # Callback that alters Seg+CNV+LOH source
    cnv_slider.js_on_change('value', callback)
    loh_slider.js_on_change('value', callback)
    cell_slider.js_on_change('value',callback)
    offset_slider.js_on_change('value', callback)
    if cnvloh_dna_dfs_2:
        cell_slider.js_on_change('value',callback_2)
        offset_slider_2.js_on_change('value',callback_2)
        alpha_slider.js_on_change('value',callback_2)
        
    if focals_dict:
        callback_focals = generate_callbackJS_focals(stored_focals_datasource, gene_focals_datasource, exon_focals_datasource, focal_effects_datasource, cnv_slider, cell_slider, offset_slider, sd_slider, genelist_field)
        cnv_slider.js_on_change('value',callback_focals)
        cell_slider.js_on_change('value',callback_focals)
        offset_slider.js_on_change('value',callback_focals)
        sd_slider.js_on_change('value',callback_focals)
    
    # Callback that alters Gene Table source
    genelist_field.js_on_change('value_input', callback_genelist)
    genelist_field.js_on_change('value_input', callback_focals)
    genelist_field.js_on_change('value_input', callback_table)
    
    cnv_slider.js_on_change('value', callback_table)
    loh_slider.js_on_change('value', callback_table)
    cell_slider.js_on_change('value',callback_table)
    offset_slider.js_on_change('value', callback_table)
    if focals_dict:
        sd_slider.js_on_change('value',callback_table)

    ### Plotting ###
    # Get genome info for chromosome to plot centromeres, etc.

    ## CNV ##
    
    print("\tCreating CNV canvas, plotting zero line and centromeres")

    # Setup figure canvas
    genome_padding = int(abs_position_dict['end']*0.01)
    cnv = figure(tools=TOOLS,y_axis_label="CNV Log2 Ratio",x_range=Range1d(0-genome_padding,abs_position_dict['end']+genome_padding),y_range=Range1d(-2,2),name='cnv_plot')
    cnv.toolbar.logo = None
    cnv.yaxis.formatter=NumeralTickFormatter(format="+0.0")
    cnv.xgrid.grid_line_color = None
    cnv.ygrid.grid_line_color = None
    cnv.xaxis.visible = False
    # Plot zero-line and centromeres
    for i in abs_position_dict:
        chromosome = i
        if chromosome in genome_info and chromosome in plot_chroms:
            chrom_info = genome_info[chromosome]

            if not math.isnan(chrom_info['centromere_start']) and not math.isnan(chrom_info['centromere_end']):
                cnv.add_layout(BoxAnnotation(left=int(chrom_info['centromere_start'])+abs_position_dict[chromosome], right=int(chrom_info['centromere_end']+abs_position_dict[chromosome]), fill_alpha=0.3, fill_color='grey',level="underlay"))
        cnv.line([abs_position_dict[chromosome],abs_position_dict[chromosome]],[-5,5],color="black",line_width=1)
    cnv.line([0,abs_position_dict['end']],[0,0],color="black",line_width=3)

    chrom_header = plot_chrom_headers(abs_position_dict, x_range = cnv.x_range)

    print("\tPlotting CNV bars")

    # Plot CNV bars
    #cnv.rect(source=source_seg,x='seg_loh_x_absolute',
    #         width='seg_loh_width',
    #         y='cnv_y',
    #         height='cnv_y_height',
    #         fill_color="cornflowerblue",fill_alpha=1,line_color="cornflowerblue")
    cnv.rect(source=source_bins, x='cnv_bin_pos', width='cnv_bin_width', height='cnv_bin_height', y='cnv_y', fill_color = "cornflowerblue", name="cnv_bins1")

    if cnvloh_dna_dfs_2:
        cnv.rect(source=source_seg_2,x='seg_loh_x_absolute',
                width='seg_loh_width',
                y='cnv_y',
                height='cnv_y_height',
                fill_color="goldenrod",fill_alpha="seg_alpha",line_color="goldenrod",line_alpha="seg_alpha",
                name="cnv_bins2")

    print("\tPlotting CNV segments")
    # Plot segments
    cnv.segment(source=source_seg,x0="seg_x0_absolute",x1="seg_x1_absolute",
                y0="cnv_y_display",y1="cnv_y_display",
               color="red",line_width=5, name="cnv_segs1")
    
    if cnvloh_dna_dfs_2:
        print("Plotting comparator CNV segments")
        # Plot segments
        cnv.segment(source=source_seg_2,x0="seg_x0_absolute",x1="seg_x1_absolute",
                y0="cnv_y",y1="cnv_y",
               color="blue",line_width=5,alpha="seg_alpha",name="cnv_segs2")
    
    cnv_labels = LabelSet(x='seg_loh_x_absolute', y='cnv_label_pos', text='cnv_label_text',
                      source=source_seg, render_mode='canvas', text_align="center", text_alpha="cnv_label_alpha")
    cnv.add_layout(cnv_labels)
    
    hover_cnvsegs = HoverTool(tooltips="""
        <b><u>@cnvloh_status</u><br>
        @seg_chrom: @seg_x0 - @seg_x1</b><br>
        <i>Segment @seg_index_local of @seg_per_chrom</i><br>
        Length: @seg_length{0}<br>
        Log2: @cnv_y{0.00}<br>
        Copies: @cnv_num_copy{0.0}<br>
        St. Dev.: @cnv_stdev{0.00}<br>
        Points: @cnv_points{0}<br>
        Reads: @cnv_reads{0}<br>
        Cov. Depth: @cnv_cov_depth{0.00}x
        """,
    names=['cnv_segs1','cnv_segs2'])
    hover_cnvbins = HoverTool(tooltips="""
        <b>@cnv_bin_chrom: @cnv_bin_start - @cnv_bin_end</b><br>
        Log2: @cnv_y{0.00}<br>
        St.Dev: @cnv_bin_sd{0.00}<br>
        Points: @cnv_bin_points{0}<br>
        Reads: @cnv_bin_reads{0}<br>
        Cov. Depth: @cnv_bin_depth{0.00}x
        """,
    names=['cnv_bins1','cnv_bins2'])
    cnv.add_tools(hover_cnvbins)
    cnv.add_tools(hover_cnvsegs)
    
    # Plot focal events
    if focals_dict:
        # Filter to limit plotting of focal events to only show phenotype genes
        pheno_focal_filter = CustomJSFilter(args=dict(genelist_field=genelist_field), code='''
            const indices = [];

            const genelist_string = genelist_field.value_input;
            var genes_set = undefined;
        
            if (!(!genelist_string || genelist_string == '' || genelist_string == '[]'))
            {
                genes_set = new Set(genelist_string.replace('[','').replace(']','').split(','));
            }
            console.log(genes_set);
            
            // For each row in the genes table, see if the Entrez Id is in our genes_list.
            // If so, or if there is no genes_list, it passes the filter.
            for (let i = 0; i < source.get_length(); i++) {
                if (!genes_set || (source.data['entrez_id'][i] !== null && genes_set.has(source.data['entrez_id'][i].toString()))) {
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }

            return indices;
        ''')

        focal_sd_filter = CustomJSFilter(args=dict(sd_slider=sd_slider), code='''
            const indices = [];
            const focal_sd_cutoff = sd_slider.value;

            // For each row in the genes table, see if the focal sd absolute value is over our cutoff.
            // If so, or if there is no cutoff, it passes the filter.
            for (let i = 0; i < source.get_length(); i++) {
                if (!focal_sd_cutoff || (source.data['delta_sd'][i] !== null && Math.abs(source.data['delta_sd'][i]) >= focal_sd_cutoff)) {
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }

            return indices;
        ''')

        gene_focals_view = CDSView(source=gene_focals_datasource, filters=[pheno_focal_filter, focal_sd_filter])
        exon_focals_view = CDSView(source=exon_focals_datasource, filters=[pheno_focal_filter, focal_sd_filter])

        cnv.scatter(source=gene_focals_datasource, view=gene_focals_view, x='event_midpoint',y='display_y',color='color',marker="triangle", size=12, name="focal_genes_scatter")
        cnv.scatter(source=exon_focals_datasource, view=exon_focals_view, x='event_midpoint',y='display_y',color='color',marker="inverted_triangle",size=6, name="focal_exons_scatter")
        hover_focal_gene = HoverTool(tooltips="""
            <b><u>@symbol (@entrez_id)</u><br>
            @text_status</b><br>
            @chrom: @start - @end<br>
            Log2: @{log2}{0.00} (@{delta_log2}{+0.0})<br>
            Copies: @{focal_copies}{0.0} (@{delta_copies}{+0.0})<br>
            Δ St.Dev.: @{delta_sd}{+0.00}<br>
            """,
            names=['focal_genes_scatter'])
        hover_focal_exon = HoverTool(tooltips="""
            <b><u>@symbol (@entrez_id)</u><br>
            @text_status</b><br>
            @chrom: @start - @end<br>
            Log2: @{log2}{0.00} (@{delta_log2}{+0.0})<br>
            Copies: @{focal_copies}{0.0} (@{delta_copies}{+0.0})<br>
            Δ St.Dev.: @{delta_sd}{+0.00}<br>
            @transcript_text<br>
            # Trans. Affected: @transcripts_affected<br>
            """,
            names=['focal_exons_scatter'])
        cnv.add_tools(hover_focal_gene)
        cnv.add_tools(hover_focal_exon)

    ## LOH ##
    
    print("\tCreating LOH canvas, plotting zero line and centromeres")

    # Setup figure canvass
    loh = figure(tools=TOOLS,y_axis_label="Variant Allele Frequency",x_range=cnv.x_range,y_range=Range1d(0, 1))
    loh.yaxis.formatter=NumeralTickFormatter(format="0.00")
    loh.xgrid.grid_line_color = None
    loh.ygrid.grid_line_color = None
    loh.xaxis.visible = False
    # Plot zero-line and centromeres
    for i in abs_position_dict:
        chromosome = i
        if chromosome in genome_info and chromosome in plot_chroms:
            chrom_info = genome_info[chromosome]

            if not math.isnan(chrom_info['centromere_start']) and not math.isnan(chrom_info['centromere_end']):
                loh.add_layout(BoxAnnotation(left=int(chrom_info['centromere_start'])+abs_position_dict[chromosome], right=int(chrom_info['centromere_end']+abs_position_dict[chromosome]), fill_alpha=0.3, fill_color='grey',level="underlay"))
        loh.line([abs_position_dict[chromosome],abs_position_dict[chromosome]],[-5,5],color="black",line_width=1)
    loh.line([0,abs_position_dict['end']],[0,0],color="black")

    print("\tPlotting LOH boxes")

    # Add box showing area of LOH based on center of 90% CI
    loh.rect(source=source_seg,x='seg_loh_x_absolute',
             width='seg_loh_width',
             y='seg_loh_y',
             height='seg_loh_height',
             fill_color="seg_loh_color",fill_alpha=0.85,line_color="seg_loh_line",name="seg_loh_rect1")
    
    if cnvloh_dna_dfs_2:
        print("Plotting comparator LOH boxes")
        loh.rect(source=source_seg_2,x='seg_loh_x_absolute',
             width='seg_loh_width',
             y='seg_loh_y',
             height='seg_loh_height',
             fill_color="seg_loh_color",fill_alpha="seg_alpha",line_color="seg_loh_line",line_alpha="seg_alpha",name="seg_loh_rect2")
             
    hover_lohrects = HoverTool(tooltips="""
        <b><u>@cnvloh_status</u><br>
        @seg_chrom: @seg_x0 - @seg_x1</b><br>
        <i>Segment @seg_index_local of @seg_per_chrom</i><br>
        Length: @seg_length<br>
        BAF: @seg_loh_baf{0.00}<br>
        LOH: @seg_loh_height{0.0%}<br>
        """,
        names=['seg_loh_rect1','seg_loh_rect2'])
    loh.add_tools(hover_lohrects)
    
    print("\tCreating controls and callbacks")
    
    ### Auto calculation mode select ###
    
    recalc_mode_radio = RadioButtonGroup(labels=['Cellularity','Baseline','Combined'], active=2)
    callback_mode_radio = CustomJS(args=dict(radio_buttons = recalc_mode_radio),
                        code="""
    const button_state = radio_buttons.active;
    console.log('Mode Button State:' + button_state);
    """)
    recalc_mode_radio.js_on_click(callback_mode_radio)
    
    ### Automatic cellularity calculation button
    
    cell_calc_button = Button(label="Auto-Adjustment", button_type="success")
    cell_calc = generate_cellcalcJS(source_seg, cnv_slider, loh_slider, cell_slider, offset_slider, recalc_mode_radio, cell_calc_button)
    cell_calc_button.js_on_click(cell_calc)

    ### Writing out figures ###
    
    ## Setup final layout ##
    # Grid for the three stacked plots, CNV+Helper+LOH
    grid = gridplot([[chrom_header],[cnv],[loh]], sizing_mode="stretch_width", plot_height=450, toolbar_options={"logo":None})

    # Sliders go above plots, table goes below
    if cnvloh_dna_dfs_2:
        layout = column(
            row(cnv_slider, cnv_spinner, 
            loh_slider, loh_spinner, 
            cell_slider, cell_spinner, 
            offset_slider, offset_spinner,
            offset_slider_2, offset_spinner_2,
            alpha_slider, alpha_spinner),
            row(grid),
            sizing_mode="stretch_width"
        )
    elif focals_dict:
        layout = column(
            row(cnv_slider, cnv_spinner,
            loh_slider, loh_spinner,
            cell_slider, cell_spinner,
            offset_slider, offset_spinner,
            sd_slider, sd_spinner),
            row(cell_calc_button, recalc_mode_radio),
            row(grid),
            sizing_mode="stretch_width"
        )
    else:
        layout = column(
            row(cnv_slider, cnv_spinner,
            loh_slider, loh_spinner,
            cell_slider, cell_spinner, 
            offset_slider, offset_spinner),
            row(cell_calc_button, recalc_mode_radio),
            row(grid),
            sizing_mode="stretch_width"
        )

    print("\tWriting plot file")

    # Name and write out file
    if output_type == "html":
        if chromosome == 'end':
            output_file(os.path.join(outpath),title=f"Genome View: Interactive CNV-LOH Display")
        else:
            output_file(os.path.join(outpath),title=f"Chromosome {chromosome}: Interactive CNV-LOH Display")
        save(layout)
    else:
        with open(outpath, 'w') as output_handle:
            output_handle.write(json.dumps(json_item(layout)))

    return None
