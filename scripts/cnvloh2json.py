import os
import json
import pandas as pd
import numpy as np
from statistics import median
from math import exp
from math import log
from math import sqrt
import argparse

def ci_from_data(data, confidence=0.90):
    from scipy import stats as st
    a = 1.0 * np.array(data)
    n = len(a)
    mean, se = np.mean(a), st.sem(a)
    h = se * st.t.ppf((1 + confidence) / 2., n-1)
    ci_lower = mean-h
    ci_upper = mean+h
    return mean, ci_lower, ci_upper

def pval_from_ci(ci_lower:float, ci_upper:float, ci_width:float=0.90):
    from scipy import stats as st
    ci_width = 1-((1-ci_width)/2)
    z_score = st.norm.ppf(ci_width)
    std_err = abs(max(ci_upper,ci_lower) - min(ci_upper,ci_lower))/(2*z_score)
    z_test = ((ci_upper + ci_lower)/2)/std_err
    p_val = exp( (-0.717*z_test) - (0.416 * (z_test ** 2)))
    return p_val

def ci_from_pval(mean_val:float, p_val:float, ci_width:float=0.90):
    from scipy import stats as st
    ci_width = 1-((1-ci_width)/2)
    z_score = st.norm.ppf(ci_width)
    z_from_p = -0.862 + ((0.743-(2.404*log(p_val))) ** 0.5)
    std_err = mean_val/z_from_p
    ci_lower = mean_val - z_score*std_err
    ci_upper = mean_val + z_score*std_err
    return ci_lower,ci_upper

def cnvloh2json_segments_gatk(modeled_seg_file, cnv_subseg_file, loh_calls_file):

    seg_count = 0
    subseg_count = 0
    var_count = 0
    
    # reading in the files
    with open(modeled_seg_file, 'r') as modeled_seg_data:
        modeled_seg_lines = modeled_seg_data.readlines()
    with open(cnv_subseg_file, 'r') as cnv_subseg_data:
        cnv_subseg_lines = cnv_subseg_data.readlines()
    with open(loh_calls_file, 'r') as loh_calls_data:
        loh_calls_lines = loh_calls_data.readlines()
        
    if os.path.isfile(cnv_subseg_file.replace('.denoisedCR.tsv','')):
        has_counts = True
        with open(cnv_subseg_file.replace('.denoisedCR.tsv',''), 'r') as cnv_subseg_reads_data:
            cnv_subseg_reads_lines = cnv_subseg_reads_data.readlines() 
    else:
        has_counts = False
    
    segment_array = []
    line_index = 0
    
    #process modeled segments
    for i in modeled_seg_lines:
        if i.startswith('@'):
            pass
        elif i.startswith('CONTIG\tSTART'):
            pass
        else:
            seg_count += 1
            seg_fields = i.strip().split()
            chrom = seg_fields[0]
            seg_start = int(seg_fields[1])
            seg_end = int(seg_fields[2])
            supporting_points_cnv = int(seg_fields[3])
            supporting_points_loh = int(seg_fields[4])
            log2_copy_ratio = float(seg_fields[6])
            log2_90ci_low = float(seg_fields[5])
            log2_90ci_high = float(seg_fields[7])
            log2_pval = pval_from_ci(log2_90ci_low,log2_90ci_high)
            loh_allele_fraction = float(seg_fields[9])
            loh_90ci_low  = float(seg_fields[8])
            loh_90ci_high = float(seg_fields[10])
            loh_pval = pval_from_ci(loh_90ci_low,loh_90ci_high)
            seg_key = f'{chrom}:{seg_start}-{seg_end}'
            segment_array.append({
                'position':{'seg_index':seg_count,'seg_id':seg_key,'chrom':chrom,'start':seg_start,'end':seg_end,'length':seg_end-seg_start+1}, 
                'cnv':{
                    'log2_copy_ratio': log2_copy_ratio,
                    'cnv_supporting_points':supporting_points_cnv,
                    'cnv_supporting_reads':0,
                    'log2_stdev':None, 
                    'log2_pval':log2_pval,
                    'log2_copy_ratio_90per_ci_low': log2_90ci_low,
                    'log2_copy_ratio_90per_ci_high': log2_90ci_high,
                    'cnv_supporting_subsegments':[]},
                'loh':{
                    'loh_allele_fraction':loh_allele_fraction,
                    'loh_supporting_points': supporting_points_loh,
                    'loh_supporting_reads':0,
                    'loh_allele_fraction_pval':loh_pval,
                    'loh_allele_fraction_90per_ci_low':loh_90ci_low,
                    'loh_allele_fraction_90per_ci_high':loh_90ci_high,
                    'loh_supporting_variants':[]
                    }
            })
    for i in cnv_subseg_lines:
        line_index += 1
        if i.startswith('@'):
            pass
        elif i.startswith('CONTIG\tSTART'):
            pass
        else:
            seg_fields = i.strip().split()
            chrom = seg_fields[0]
            seg_start = int(seg_fields[1])
            seg_end = int(seg_fields[2])
            log2_copy_ratio = float(seg_fields[3])
            seg_key = f'{chrom}:{seg_start}-{seg_end}'
            if has_counts:
                reads_count = int(cnv_subseg_reads_lines[line_index].strip().split()[3])
            else:
                reads_count = None
            for j in range(0,len(segment_array)):
                if segment_array[j]['position']['chrom'] == chrom and segment_array[j]['position']['start'] <= seg_start and segment_array[j]['position']['end'] > seg_end:
                    subseg_count += 1
                    segment_array[j]['cnv']['cnv_supporting_subsegments'].append({'subseg_index':subseg_count, 'subseg_id':seg_key,'chrom':chrom,'start':seg_start,'end':seg_end,'length':seg_end-seg_start+1,'reads':reads_count,'log2_copy_ratio':log2_copy_ratio})
                    segment_array[j]['cnv']['cnv_supporting_reads'] += reads_count
                    break
    for i in loh_calls_lines:
        if i.startswith('@'):
            pass
        elif i.startswith('CONTIG\tPOSITION'):
            pass
        else:
            var_count += 1
            seg_fields = i.strip().split()
            chrom = seg_fields[0]
            var_pos = int(seg_fields[1])
            ref_count = int(seg_fields[2])
            alt_count = int(seg_fields[3])
            ref_allele = seg_fields[4]
            alt_allele = seg_fields[5]
            var_id = f'{chrom}:{var_pos}_{ref_allele}>{alt_allele}'
            total_reads = ref_count + alt_count
            ref_vaf = ref_count / total_reads
            alt_vaf = alt_count / total_reads
            b_allele_freq = min([ref_vaf, alt_vaf])
            for j in range(0,len(segment_array)):
                if segment_array[j]['position']['chrom'] == chrom and segment_array[j]['position']['start'] <= var_pos and segment_array[j]['position']['end'] > var_pos:
                    segment_array[j]['loh']['loh_supporting_variants'].append({'var_index':var_count, 'var_id':var_id,'chrom':chrom, 'position':var_pos, 'total_reads':total_reads,'b_allele_freq':b_allele_freq,'ref_allele':ref_allele, 'alt_allele':alt_allele,'ref_count':ref_count,'alt_count':alt_count,'ref_vaf':ref_vaf,'alt_vaf':alt_vaf})
                    segment_array[j]['loh']['loh_supporting_reads'] += total_reads
                    break
    for i in range(0,len(segment_array)):
        segment_array[i]['cnv']['log2_stdev']=np.nan_to_num(np.std([j['log2_copy_ratio'] for j in segment_array[i]['cnv']['cnv_supporting_subsegments']]))
    return segment_array

def cnvloh2json_segments_varscan(modeled_seg_file, cnv_subseg_file, loh_calls_file):
    
    seg_count = 0
    bin_count = 0
    var_count = 0
        
    # reading in the files
    modeled_segs_data = pd.read_csv(modeled_seg_file, header=0, sep=",")
    modeled_segs_data.columns=["index","ID","chrom","start","end","num_mark","seg_mean","bstat","pval","lcl","ucl"]
    
    cnv_bins_data = pd.read_csv(cnv_subseg_file, header=0, sep="\t")
    cnv_bins_data.columns=["chrom","start","end","num_positions","n_depth","t_depth","seg_mean","gc","call","raw_ratio"]
    
    loh_vars_data = pd.read_csv(loh_calls_file, header=0, sep=",")
    loh_vars_data.columns=["index","chrom","pos","ref","alt","ref_reads","alt_reads","vaf"]
    
    segment_array = []
    
    segments_loh_allele_freq={}
    
    #process modeled segments
    for i in range(0,len(list(modeled_segs_data['index']))):
        chrom = list(modeled_segs_data['chrom'])[i-1]
        seg_start = int(list(modeled_segs_data['start'])[i-1])
        seg_end = int(list(modeled_segs_data['end'])[i-1])
        supporting_points_cnv = int(list(modeled_segs_data['num_mark'])[i-1])
        log2_copy_ratio = float(list(modeled_segs_data['seg_mean'])[i-1])
        seg_pval = float(list(modeled_segs_data['pval'])[i-1])
        if seg_pval == 0:
            log2_90ci_low,log2_90ci_high = log2_copy_ratio,log2_copy_ratio
        else:
            log2_90ci_low,log2_90ci_high = ci_from_pval(log2_copy_ratio,seg_pval)
        seg_key = f'{chrom}:{seg_start}-{seg_end}'
        segment_array.append({
            'position':{'seg_index':i,'seg_id':seg_key,'chrom':chrom,'start':seg_start,'end':seg_end,'length':seg_end-seg_start+1}, 
            'cnv':{
                'log2_copy_ratio': log2_copy_ratio,
                'cnv_supporting_points':supporting_points_cnv,
                'log2_pval':seg_pval,
                'log2_copy_ratio_90per_ci_low': log2_90ci_low,
                'log2_copy_ratio_90per_ci_high': log2_90ci_high,
                'cnv_supporting_subsegments':[]},
            'loh':{
                'loh_allele_fraction':None,
                'loh_supporting_points': 0,
                'loh_allele_fraction_pval':None,
                'loh_allele_fraction_90per_ci_low':None,
                'loh_allele_fraction_90per_ci_high':None,
                'loh_supporting_variants':[]
                }})
        seg_count += 1
    bin_chrom_list = list(cnv_bins_data['chrom'])
    bin_start_list = list(cnv_bins_data['start'])
    bin_end_list = list(cnv_bins_data['end'])
    bin_log2_list = list(cnv_bins_data['seg_mean'])
    current_segment = 0
    for i in range(0,len(list(cnv_bins_data['chrom']))):
        chrom = bin_chrom_list[i]
        seg_start = int(bin_start_list[i])
        seg_end = int(bin_end_list[i])
        log2_copy_ratio = float(bin_log2_list[i])
        seg_key = f'{chrom}:{seg_start}-{seg_end}'
        for j in range(current_segment,len(segment_array)):
            if segment_array[j]['position']['chrom'] == chrom and segment_array[j]['position']['start'] <= seg_start and segment_array[j]['position']['end'] > seg_end:
                bin_count += 1
                segment_array[j]['cnv']['cnv_supporting_subsegments'].append({'subseg_index':bin_count, 'subseg_id':seg_key,'chrom':chrom,'start':seg_start,'end':seg_end,'length':seg_end-seg_start+1,'log2_copy_ratio':log2_copy_ratio})
                current_segment = j
                break
    current_segment = 0
    var_chrom_list = list(loh_vars_data['chrom'])
    var_pos_list = list(loh_vars_data['pos'])
    var_ref_list = list(loh_vars_data['ref'])
    var_alt_list = list(loh_vars_data['alt'])
    var_refreads_list = list(loh_vars_data['ref_reads'])
    var_altreads_list = list(loh_vars_data['alt_reads'])
    for i in range(0,len(list(loh_vars_data['chrom']))):
        var_count = i
        chrom = var_chrom_list[i]
        var_pos = int(var_pos_list[i])
        ref_allele = var_ref_list[i]
        alt_allele = var_alt_list[i]
        ref_count = int(var_refreads_list[i])
        alt_count = int(var_altreads_list[i])
        var_id = f'{chrom}:{var_pos}_{ref_allele}>{alt_allele}'
        total_reads = ref_count + alt_count
        ref_vaf = ref_count / total_reads
        alt_vaf = alt_count / total_reads
        b_allele_freq = min([ref_vaf, alt_vaf])
        for j in range(current_segment,len(segment_array)):
            if segment_array[j]['position']['chrom'] == chrom and segment_array[j]['position']['start'] <= var_pos and segment_array[j]['position']['end'] > var_pos:
                segment_array[j]['loh']['loh_supporting_variants'].append({'var_index':var_count, 'var_id':var_id,'chrom':chrom, 'position':var_pos, 'total_reads':total_reads,'b_allele_freq':b_allele_freq,'ref_allele':ref_allele, 'alt_allele':alt_allele,'ref_count':ref_count,'alt_count':alt_count,'ref_vaf':ref_vaf,'alt_vaf':alt_vaf})
                if j not in segments_loh_allele_freq:
                    segments_loh_allele_freq[j] = []
                segments_loh_allele_freq[j].append(b_allele_freq)
                current_segment = j
                break
    for i in segments_loh_allele_freq:
        segment_array[i]['loh']['loh_supporting_points'] = len(segments_loh_allele_freq[i])
        segment_array[i]['loh']['loh_allele_fraction'],segment_array[i]['loh']['loh_allele_fraction_90per_ci_low'],segment_array[i]['loh']['loh_allele_fraction_90per_ci_high'] = ci_from_data(segments_loh_allele_freq[i])
        segment_array[i]['loh']['pval']=pval_from_ci(segment_array[i]['loh']['loh_allele_fraction_90per_ci_low'],segment_array[i]['loh']['loh_allele_fraction_90per_ci_high'])

    return segment_array

def cnvloh2json_segments_methyl(methyl_seg_file, methyl_bins_file):
    
    seg_count = 0
    bin_count = 0
    
    # reading in the files
    methyl_segs_data = pd.read_csv(methyl_seg_file, header=0, sep=",")
    methyl_segs_data.columns=["index","feature","chrom","start","end","num_mark","bstat","pval","seg_mean","seg_median"]
    
    methyl_bins_data = pd.read_csv(methyl_bins_file, header=0, sep=",")
    methyl_bins_data.columns=["index","chrom","start","end","feature","seg_mean"]
    
    segment_array = []
    
    #process modeled segments
    for i in list(methyl_segs_data['index']):
        chrom = list(methyl_segs_data['chrom'])[i-1]
        seg_start = int(list(methyl_segs_data['start'])[i-1])
        seg_end = int(list(methyl_segs_data['end'])[i-1])
        supporting_points_cnv = int(list(methyl_segs_data['num_mark'])[i-1])
        log2_copy_ratio = float(list(methyl_segs_data['seg_mean'])[i-1])
        seg_pval = float(list(methyl_segs_data['pval'])[i-1])
        if seg_pval == 0:
            log2_90ci_low,log2_90ci_high = log2_copy_ratio,log2_copy_ratio
        else:
            log2_90ci_low,log2_90ci_high = ci_from_pval(log2_copy_ratio,seg_pval)
        seg_key = f'{chrom}:{seg_start}-{seg_end}'
        segment_array.append({
            'position':{'seg_index':i,'seg_id':seg_key,'chrom':chrom,'start':seg_start,'end':seg_end,'length':seg_end-seg_start+1}, 
            'cnv':{
                'log2_copy_ratio': log2_copy_ratio,
                'cnv_supporting_points':supporting_points_cnv,
                'log2_pval':seg_pval,
                'log2_copy_ratio_90per_ci_low': log2_90ci_low,
                'log2_copy_ratio_90per_ci_high': log2_90ci_high,
                'cnv_supporting_subsegments':[]},
            })
    bin_chrom_list = list(methyl_bins_data['chrom'])
    bin_start_list = list(methyl_bins_data['start'])
    bin_end_list = list(methyl_bins_data['end'])
    bin_log2_list = list(methyl_bins_data['seg_mean'])
    for i in list(methyl_bins_data['index']):
        chrom = bin_chrom_list[i-1]
        seg_start = int(bin_start_list[i-1])
        seg_end = int(bin_end_list[i-1])
        log2_copy_ratio = float(bin_log2_list[i-1])
        seg_key = f'{chrom}:{seg_start}-{seg_end}'
        for j in range(0,len(segment_array)):
            if segment_array[j]['position']['chrom'] == chrom and segment_array[j]['position']['start'] <= seg_start and segment_array[j]['position']['end'] > seg_end:
                bin_count += 1
                segment_array[j]['cnv']['cnv_supporting_subsegments'].append({'subseg_index':subseg_count, 'subseg_id':seg_key,'chrom':chrom,'start':seg_start,'end':seg_end,'length':seg_end-seg_start+1,'log2_copy_ratio':log2_copy_ratio,'reads':reads_count})
                break
    return segment_array

def calc_recenter_value(segment_array:list,min_seg_length:int=25000000, max_loh=0.1,min_segments:int=1):
    log2_values = []
    for i in range(0,len(segment_array)):
        if segment_array[i]['position']['length'] >= min_seg_length and ('loh' not in segment_array[i] or (type(segment_array[i]['loh']['loh_allele_fraction']) is float and abs(0.5-segment_array[i]['loh']['loh_allele_fraction']) <= max_loh)):
            log2_values.append(segment_array[i]['cnv']['log2_copy_ratio'])
            #print(f"Length: {segment_array[i]['position']['length']} | LOH: {abs(0.5-segment_array[i]['loh']['loh_allele_fraction'])} | CNV: {segment_array[i]['cnv']['log2_copy_ratio']}")
    if len(log2_values) >= min_segments:
        recenter_value = -1 * median(log2_values)
    else:
        recenter_value = None
    print(f'From no-LOH segment values:\n{log2_values}\nCalculated recentering value: {recenter_value}.')
    return recenter_value
    
def apply_recenter_value(segment_array:list, recenter_value:float):
    for i in range(0,len(segment_array)):
        segment_array[i]['cnv']['log2_copy_ratio'] = segment_array[i]['cnv']['log2_copy_ratio'] + recenter_value
        segment_array[i]['cnv']['log2_copy_ratio_90per_ci_low'] = segment_array[i]['cnv']['log2_copy_ratio_90per_ci_low'] + recenter_value
        segment_array[i]['cnv']['log2_copy_ratio_90per_ci_high'] = segment_array[i]['cnv']['log2_copy_ratio_90per_ci_high'] + recenter_value
        for j in range(0,len(segment_array[i]['cnv']['cnv_supporting_subsegments'])):
            segment_array[i]['cnv']['cnv_supporting_subsegments'][j]['log2_copy_ratio'] = segment_array[i]['cnv']['cnv_supporting_subsegments'][j]['log2_copy_ratio'] + recenter_value
    return segment_array

def recenter_cnvloh_segments(segment_array,min_seg_length:int=25000000, max_loh=0.1,min_segments:int=1):
    recenter_value = calc_recenter_value(segment_array, min_seg_length, max_loh, min_segments)
    if recenter_value:
        segment_array = apply_recenter_value(segment_array, recenter_value)
    else:
        print('Skipping recentering. Insufficient valid segments.')
    return segment_array, recenter_value

def add_qc_dlrs(cnvloh_dict):
    prev_log2rat = 0
    prev_seg_index = ''
    dlrs_rows = []
    for segment in cnvloh_dict['segments']:
        segment_chr = segment['position']['chrom']
        seg_index = segment['position']['seg_index']
        for subsegment in segment['cnv']['cnv_supporting_subsegments']:
            log2rat = subsegment['log2_copy_ratio']
            if seg_index == prev_seg_index:
                dlrs_rows.append({
                    'spread': abs(log2rat - prev_log2rat),
                    'chr': segment_chr,
                    'log2r': log2rat,
                     },
                )
            prev_log2rat = log2rat
            prev_seg_index = seg_index
    dlrs = pd.DataFrame.from_records(dlrs_rows)
    quant_25, quant_75 = dlrs['spread'].quantile([0.25, 0.75])
    dlrsiq = dlrs[(dlrs['spread'] > quant_25) & (dlrs['spread'] < quant_75)]
    cnvloh_dict['metadata']["der_logratio_spread"] = dlrs['spread'].std()/sqrt(2)
    cnvloh_dict['metadata']["n_dlrs"] = len(dlrs['spread'])
    cnvloh_dict['metadata']["dlrs_interquartile"] = dlrsiq['spread'].std()/sqrt(2)
    return cnvloh_dict

def generate_cnvloh_json(sample_dir, sample_name, modeled_seg_file, cnv_subseg_file, loh_calls_file, data_type = "WES", case_name=None, cnv_params=None, software="GATK", software_version=None, sample_type=None):
    if (data_type.upper() == "WES" or data_type.upper() == "WGS") and software.upper() == "GATK":
        json_output_str = f'{sample_name}_cnvloh_output_{data_type}-GATK.json'
        min_seg_length=25000000
        max_loh = 0.1
        segment_array = cnvloh2json_segments_gatk(modeled_seg_file, cnv_subseg_file, loh_calls_file)
    elif (data_type.upper() == "WES" or data_type.upper() == "WGS") and software.upper() == "VARSCAN":
        json_output_str = f'{sample_name}_cnvloh_output_{data_type}-Varscan.json'
        segment_array = cnvloh2json_segments_varscan(modeled_seg_file, cnv_subseg_file, loh_calls_file)
        min_seg_length=5000000
        max_loh=0.2
    elif data_type.upper().startswith('METH'):
        json_output_str = f'{sample_name}_cnvloh_output_methylation.json'
        min_seg_length=25000000
        max_loh=0
        segment_array = cnvloh2json_segments_methyl(modeled_seg_file, cnv_subseg_file)
    else:
        raise ValueError(f'The specified data type and software combination is not supported: {data_type},{software}.')
    segment_array, recenter_value = recenter_cnvloh_segments(segment_array,min_seg_length=min_seg_length,max_loh=max_loh)
    metadata_dict={'data_type':data_type,'software':software,'software_version':software_version,'case_name':case_name,'sample_name':sample_name,'sample_type':sample_type,'analysis_parameters':cnv_params,'recentering_value':recenter_value}
    cnvloh_dict = {'metadata':metadata_dict, 'segments':segment_array}
    cnvloh_dict = add_qc_dlrs(cnvloh_dict)
    json_string = json.dumps(cnvloh_dict,indent=1)
    json_path = os.path.join(sample_dir, json_output_str)
    with open(json_path,'w') as cnvloh_json_out:
        cnvloh_json_out.write(json_string)
    return None

def run_cnv2json_argparser():
    parser = argparse.ArgumentParser(description="Creates CNV-LOH JSON file from raw outputs.")
    parser.add_argument('--sample-name', help="(Optional) Name of CNV-LOH sample being processed.", default='cnvloh_sample')
    parser.add_argument('--cnvloh-dir', help='Directory to output CNV-LOH JSON file.',required=True)
    parser.add_argument('--hets-file', help="Path to CNV-LOH *.hets.tsv File",default=None)
    parser.add_argument('--counts-file', help="Path to CNV-LOH *.counts.tsv.denoisedCR.tsv",required=True)
    parser.add_argument('--modelsegs-file', help="Path to CNV-LOH *.modelFinal.seg file",required=True)
    parser.add_argument('--sample-type', help='(Optional) Sample type -- typically Normal, Tumor, or Tumor_Normal', default=None)
    parser.add_argument('--case-name', help='(Optional) Name of the case to which the sample belongs.', default=None)
    parser.add_argument('--software', help='Software used to create data.', default='GATK')
    parser.add_argument('--version', help="Version of software used to create data.", default=None)
    parser.add_argument('--data-type', help='Type of input data', default='WES')
    parser.add_argument('--parameters', help='String of parameters used to run software.', default=None)
    return parser

#### Main Function ####

if __name__ == '__main__':
    args = run_cnv2json_argparser().parse_args()
    sample_dir = args.cnvloh_dir
    sample_name = args.sample_name
    case_name=args.case_name
    sample_type=args.sample_type
    software=args.software
    data_type=args.data_type
    cnv_params = args.parameters
    software_version = args.version
    modeled_seg_file = args.modelsegs_file
    cnv_subseg_file = args.counts_file
    loh_calls_file = args.hets_file
    data_type = args.data_type
    generate_cnvloh_json(sample_dir, sample_name, modeled_seg_file, cnv_subseg_file, loh_calls_file, data_type, case_name, cnv_params, software, software_version, sample_type)