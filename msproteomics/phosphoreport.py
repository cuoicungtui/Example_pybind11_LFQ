import sys
import argparse
import pandas as pd

import msproteomics

def main():
    
    parser = argparse.ArgumentParser(description='Create phosphosite and phosphopeptide report.', usage = 'phosphoreport inputfile [...]')
    
    parser.add_argument('file')
    parser.add_argument('-tool', '--tool', default = 'generic', choices=['generic', 'diann', 'sn'], type = str, help = 'Processing tool, diann for DIA-NN and sn for Spectronaut')
    parser.add_argument('-output_site', '--output_site', default = 'phosphosite-report.txt', type = str, help = 'Filename for phosphosite report')
    parser.add_argument('-output_peptide', '--output_peptide', default = 'phosphopeptide-report.txt', type = str, help = 'Filename for phosphopeptide report')
    parser.add_argument('-sample_id_col', '--sample_id_col', default = None, type = str, help = 'The column containing sample id')
    parser.add_argument('-intensity_col', '--intensity_col', default = None, type = str, help = 'The column containing intensity')
    parser.add_argument('-secondary_id_cols', '--secondary_id_cols', default = None, nargs='*', type = str, help = 'Columns forming secondary ids')
    parser.add_argument('-annotation_cols', '--annotation_cols', default = None, nargs='*', type = str, help = 'Annotation columns')

    parser.add_argument('-protein_id_col', '--protein_id_col', default = None, type = str, help = 'The column containing protein id')
    parser.add_argument('-site_id_col', '--site_id_col', default = None, type = str, help = 'The column containing site id')

    parser.add_argument('-site_filter_double_less', '--site_filter_double_less', default = None, nargs=2, action='append', type = str, help = 'Site filtering double less than or equal to')
    parser.add_argument('-site_filter_double_greater', '--site_filter_double_greater', default = None, nargs=2, action='append', type = str, help = 'Site filtering double greater than or equal to')
    parser.add_argument('-site_filter_string_equal', '--site_filter_string_equal', default = None, nargs=2, action='append', type = str, help = 'Site filtering string equal')
    parser.add_argument('-site_filter_string_not_equal', '--site_filter_string_not_equal', default = None, nargs=2, action='append', type = str, help = 'Site filtering string not equal')

    parser.add_argument('-modified_sequence_col', '--modified_sequence_col', default = None, type = str, help = 'The column containing modified sequences')
    parser.add_argument('-regex_str', '--regex_str', default = None, type = str, help = 'Regular expression to extract sites')
    parser.add_argument('-target_modification', '--target_modification', default = None, type = str, help = 'Target modification')
    parser.add_argument('-modifications', '--modifications', default = None, nargs='*', type = str, help = 'List of modifications')

    parser.add_argument('-peptide_filter_double_less', '--peptide_filter_double_less', default = None, nargs=2, action='append', type = str, help = 'Peptide filtering double less than or equal to')
    parser.add_argument('-peptide_filter_double_greater', '--peptide_filter_double_greater', default = None, nargs=2, action='append', type = str, help = 'Peptide filtering double greater than or equal to')
    parser.add_argument('-peptide_filter_string_equal', '--peptide_filter_string_equal', default = None, nargs=2, action='append', type = str, help = 'Peptide filtering string equal')
    parser.add_argument('-peptide_filter_string_not_equal', '--peptide_filter_string_not_equal', default = None, nargs=2, action='append', type = str, help = 'Peptide filtering string not equal')

    parser.add_argument('-normalize', '--normalize', default = 'none', type = str, choices=['none', 'median'], help = 'Normalization method')
    parser.add_argument('-quant_method', '--quant_method', default = 'maxLFQ', type = str, help = 'Quantification method')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        print('Missing input filename. Try phosphoreport -h.')
        sys.exit(0)
    
    inputfile = sys.argv[1]

    if args.tool == "generic" or args.tool == "diann":
        if args.sample_id_col is None:
            args.sample_id_col = 'Run'
        if args.intensity_col is None:
            args.intensity_col = 'Fragment.Intensity'
        if args.secondary_id_cols is None:
            args.secondary_id_cols = ['Precursor.Id', 'Fragment.Rel.Id']
        if args.annotation_cols is None:
            args.annotation_cols = ['Fasta.Files']

        # site paramters
        if args.protein_id_col is None:
            args.protein_id_col = 'Protein.Ids'
        if args.site_id_col is None:
            args.site_id_col = 'Phospho.Site.Specs'

        if args.site_filter_double_less is None:
            args.site_filter_double_less = [['Global.Q.Value', '0.01']]
        if args.site_filter_double_greater is None:
            args.site_filter_double_greater = [['PTM.Site.Confidence', '0.01']]
        if args.site_filter_string_not_equal is None:
            args.site_filter_string_not_equal = [['Fragment.Intensity', '0.0']]

        # peptide paramters
        if args.modified_sequence_col is None:
            args.modified_sequence_col = 'Modified.Sequence'
        if args.regex_str is None:
            args.regex_str = '\\(UniMod:[0-9]*\\)'
        if args.target_modification is None:
            args.target_modification = '(UniMod:21)'                

        if args.peptide_filter_double_less is None:
            args.peptide_filter_double_less = [['Global.Q.Value', '0.01']]
        if args.peptide_filter_string_not_equal is None:
            args.peptide_filter_string_not_equal = [['Fragment.Intensity', '0.0']]
    
    elif args.tool == 'sn':

        if args.sample_id_col is None:
            args.sample_id_col = 'R.FileName'
        if args.intensity_col is None:
            args.intensity_col = 'F.PeakArea'
        if args.secondary_id_cols is None:
            args.secondary_id_cols = ['EG.PrecursorId', 'EG.Library', 'FG.Charge', 'F.FrgIon', 'F.Charge', 'F.FrgLossType']
        if args.annotation_cols is None:
            args.annotation_cols = ['PG.Organisms']

        # site paramters
        if args.protein_id_col is None:
            args.protein_id_col = 'PG.ProteinAccessions'
        if args.site_id_col is None:
            args.site_id_col = 'EG.ProteinPTMLocations'
        
        if args.site_filter_double_greater is None:
            args.site_filter_double_greater = [['EG.PTMAssayProbability', '0.75']]

        # peptide paramters
        if args.modified_sequence_col is None:
            args.modified_sequence_col = 'EG.ModifiedSequence'
        if args.regex_str is None:
            args.regex_str = '\[[^\\[]+\]'
        if args.target_modification is None:
            args.target_modification = '[Phospho (STY)]'
    
    else:
        raise Exception('Unsupported tool.')

    print('\nGeneral setting:')
    print('  input data = ', inputfile)
    print('  processing tool = ', args.tool)
    print('  sample_id_col = ', args.sample_id_col)
    print('  intensity_col = ', args.intensity_col)
    print('  secondary_id_cols = ', args.secondary_id_cols)        
    print('  annotation_cols = ', args.annotation_cols)        
    
    print('\nPhosphosite setting:')
    print('  protein_id_col = ', args.protein_id_col) 
    print('  site_id_col = ', args.site_id_col) 
    print('  site_filter_double_less = ', args.site_filter_double_less) 
    print('  site_filter_double_greater = ', args.site_filter_double_greater) 
    print('  site_filter_string_not_equal = ', args.site_filter_string_not_equal) 
    print('  output_site = ', args.output_site) 

    print('\nPhosphopeptide setting:')
    print('  modified_sequence_col = ', args.modified_sequence_col) 
    print('  regex_str = ', args.regex_str) 
    if args.modifications is None:
        print('  target_modification = ', args.target_modification) 
    else:
        print('  modifications = ', args.modifications) 

    print('  peptide_filter_double_less = ', args.peptide_filter_double_less) 
    print('  peptide_filter_string_not_equal = ', args.peptide_filter_string_not_equal) 
    print('  output_peptide = ', args.output_peptide) 
    
    #TODO: check filter datatype
    
    #TODO: check or read required columns
    print('\nLoading data file:', inputfile)
    d = pd.read_csv(inputfile, sep = '\t')                
    print(d.shape[0], 'rows x', d.shape[1], 'columns read')

    #TODO: check filter datatype

    import numpy as np
    if args.normalize == 'median':
        log2_intensity = msproteomics.normalize(d, args.sample_id_col, args.intensity_col)
    else:    
        log2_intensity = np.log2(d[args.intensity_col],  out = np.zeros_like(d[args.intensity_col]), where = (d[args.intensity_col] != 0))
    
    if sum(np.isnan(log2_intensity)) > 0:
        print("Warning: missing values in the intensity column.")

    d = d.assign(log2_intensity = log2_intensity)    

    #-- site filtering
    index = pd.Series([True] * d.shape[0])

    if args.site_filter_double_less is not None:
        for f in args.site_filter_double_less:
            col = f[0]
            val = float(f[1])
            index = index & (~d[col].isna()) & (d[col] <= val)

    if args.site_filter_double_greater is not None:
        for f in args.site_filter_double_greater:
            col = f[0]
            val = float(f[1])
            index = index & (~d[col].isna()) & (d[col] >= val)

    if args.site_filter_string_not_equal is not None:
        for f in args.site_filter_string_not_equal:
            col = f[0]
            val = f[1]
            index = index & (~d[col].isna()) & (d[col].astype(str) != val)

    if args.site_filter_string_equal is not None:
        for f in args.site_filter_string_equal:
            col = f[0]
            val = f[1]
            index = index & (~d[col].isna()) & (d[col].astype(str) == val)

    d_filtered = d[index]
    
    print('\n{0} rows after site filtering'.format(d_filtered.shape[0]))

    protein_ids = d_filtered[args.protein_id_col].fillna('').tolist()
    
    ptm_locations = d_filtered[args.site_id_col].fillna('').tolist()

    print('Creating site identifiers')
    (all_sites, unique_sites, first_sites, unique_first_sites) = msproteomics.create_site_key(protein_ids, ptm_locations)
    d_filtered = msproteomics.create_site_report_longformat(d_filtered.assign(all_sites = all_sites), first_sites)
    d_filtered = d_filtered.loc[d_filtered['site'] != '']

    print('Creating site report')
    result = msproteomics.create_report_wideformat(d_filtered, method = args.quant_method,
                                            sample_id_col = args.sample_id_col,
                                            intensity_col = "log2_intensity", 
                                            secondary_id_cols = args.secondary_id_cols, annotation_cols = args.annotation_cols + ['all_sites'])

    result.to_csv(args.output_site, sep = '\t', index = False)


    #-- peptide filtering
    index = pd.Series([True] * d.shape[0])

    if args.peptide_filter_double_less is not None:
        for f in args.peptide_filter_double_less:
            col = f[0]
            val = float(f[1])
            index = index & (~d[col].isna()) & (d[col] <= val)

    if args.peptide_filter_double_greater is not None:
        for f in args.peptide_filter_double_greater:
            col = f[0]
            val = float(f[1])
            index = index & (~d[col].isna()) & (d[col] >= val)

    if args.peptide_filter_string_not_equal is not None:
        for f in args.peptide_filter_string_not_equal:
            col = f[0]
            val = f[1]
            index = index & (~d[col].isna()) & (d[col].astype(str) != val)

    if args.peptide_filter_string_equal is not None:
        for f in args.peptide_filter_string_equal:
            col = f[0]
            val = f[1]
            index = index & (~d[col].isna()) & (d[col].astype(str) == val)

    d_filtered = d[index]        
    print('\n{0} rows after peptide filtering'.format(d_filtered.shape[0]))

    print('Creating peptide identifiers')
    if args.modifications is None:
        (ptm_key, STY_count) = msproteomics.create_peptide_key(d_filtered[args.modified_sequence_col].tolist(),
                                                               regex_str = args.regex_str,
                                                               target_modification = args.target_modification)
    else:                                                               
        (ptm_key, STY_count) = msproteomics.create_peptide_key(d_filtered[args.modified_sequence_col].tolist(),
                                                               regex_str = args.regex_str,
                                                               modifications = args.modifications)
    print('Creating peptide report')
    result = msproteomics.create_report_wideformat(d_filtered.assign(ptm_key = ptm_key, STY_count = STY_count),
                                                    sample_id_col = args.sample_id_col,
                                                    intensity_col = "log2_intensity",
                                                    primary_id = 'ptm_key', 
                                                    secondary_id_cols = args.secondary_id_cols,
                                                    annotation_cols = ['STY_count'] + args.annotation_cols,
                                                    method = args.quant_method)
    result.to_csv(args.output_peptide, sep = '\t', index = False)
