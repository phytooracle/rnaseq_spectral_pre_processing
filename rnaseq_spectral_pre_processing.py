#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez
Date   : 2022-03-30
Purpose: RNA-Seq & Spectral preprocessing
"""

import argparse
import os
import sys
import pandas as pd
import multiprocessing
import warnings
warnings.filterwarnings("ignore")


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='RNA-Seq & Spectral preprocessing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-r',
                        '--rnaseq_csv',
                        help='RNA-Seq input CSV.',
                        metavar='str',
                        type=str,
                        default='https://data.cyverse.org/dav-anon/iplant/home/emmanuelgonzalez/rnaseq_spectral/rnaseq/Cotton_TPM_TOP25.csv')

    parser.add_argument('-s',
                        '--spectral_csv',
                        help='Spectral input CSV.',
                        metavar='str',
                        type=str,
                        default='https://data.cyverse.org/dav-anon/iplant/home/emmanuelgonzalez/rnaseq_spectral/spectral/2019averagescans.csv')

    parser.add_argument('-f',
                        '--fieldbook_csv',
                        help='Fieldbook input CSV.',
                        metavar='str',
                        type=str,
                        default='https://data.cyverse.org/dav-anon/iplant/home/emmanuelgonzalez/rnaseq_spectral/Duke_Fieldbook_2019.csv')

    parser.add_argument('-o',
                        '--out_dir',
                        help='Output directory',
                        metavar='str',
                        type=str,
                        default='2019_cotton_rnaseq_spectra')
    

    parser.add_argument('-t',
                        '--treatment',
                        help='Treatment zones to output.',
                        metavar='list',
                        nargs='+',
                        default=['WW', 'WL'])

    return parser.parse_args()


# --------------------------------------------------
def get_rnaseq_data(csv_path, fb_path):
    '''
    Reads the RNA-Seq file and adds fieldbook data to it using the join function within Pandas.

    Input:
        - csv_path: path to CSV file containing RNA-Seq data
        - fb_path: path to CSV file containing fieldbook information
    Output: 
        - rna_fb: Dataframe containing RNA-Seq and fieldbook information
    '''

    fb = get_fieldbook_data(csv_path=fb_path, index='entry')
    rna_df = pd.read_csv(csv_path)
    rna_df['Gene'] = rna_df['Gene'].str.replace('Gohir.', 'Gh_')
    rna_df = rna_df.set_index('Gene').T.reset_index().rename(columns={'index': 'entry'}).set_index('entry', drop=True)
    
    if 'TPM' in csv_path:
        rna_df = rna_df.reset_index()
        rna_df['treatment'] = ['_'.join(entry.split('_')[-1:]) for entry in rna_df['entry']]
        rna_df['entry'] = ['_'.join(entry.split('_')[:-1]) for entry in rna_df['entry']]
        rna_fb = fb.reset_index().set_index(['entry', 'treatment']).join(rna_df.set_index(['entry', 'treatment'])).reset_index().set_index('plot')
    
    else:
        rna_fb = join_dataframes(df1=fb, df2=rna_df, convert_dtypes=True).set_index('plot')
        
    return rna_fb


# --------------------------------------------------
def get_spectral_data(csv_path, fb_path):
    '''
    Reads hyperspectral file and adds fieldbook data to it using the join function within Pandas.

    Input: 
        - csv_path: path to CSV file containing hyperspectral data
        - fb_path: path to CSV file containing fieldbook information
    Output: 
        - spectra_df: Dataframe containing only hyperspectral data
        - spectra_fb: Dataframe containing both hyperspectral and fieldbook information
        - fb: Dataframe containing only fieldbook data
    '''

    fb = get_fieldbook_data(csv_path=fb_path, index='plot')
    spectra_df = pd.read_csv(csv_path)

    spectra_df = fix_plot_column(df=spectra_df)

    spectra_df = spectra_df.set_index('Lambda').T.reset_index().rename(columns={'index': 'plot'}).set_index('plot', drop=True)

    spectra_df.columns = spectra_df.columns.astype(str)
    spectra_fb = join_dataframes(df1=fb, df2=spectra_df, convert_dtypes=True)
    
    return spectra_df, spectra_fb, fb


# --------------------------------------------------
def get_fieldbook_data(csv_path, index):
    '''
    Reads fieldbook file.

    Input:
        - csv_path: path to CSV file containing the fieldbook information
        - index: column name for column to be set as index
    Output: 
        - fb: Dataframe containing fieldbook data
    '''

    fb = pd.read_csv(csv_path).dropna(subset=['plot'])
    fb['plot'] = pd.to_numeric(fb['plot']).astype(int).astype(str)
    fb = fb.set_index(index)
    

    return fb


# --------------------------------------------------
def fix_plot_column(df):
    '''
    Extracts plot number from hyperspectral collection ID (column)

    Input: 
        - df: dataframe which needs to be fixed
    Output: 
        - df: fixed dataframe with plot number
    '''

    plot_number_columns = [column.split('_')[1].replace('00000.asd', '').replace('p', '') for column in df.columns if 'Cotton' in column]
    plot_number_columns.insert(0, 'Lambda')

    df.columns = plot_number_columns
    
    return df


# --------------------------------------------------
def join_dataframes(df1, df2, convert_dtypes=False, groupby_list=None, stat=None):
    '''
    Joins two dataframes with options for converting data types and group by.

    Input:
        - df1: first dataframe with shared index
        - df2: second dataframe with shared index
        - convert_dtypes: bool to indicate of converting data types is necessary
        - groupby_list: list containing column names by which to group by
        - stat: statistic to use for groupby, either "mean" or "median"
    Output: 
        - joined_df: Joined dataframe with selected operations (convert_dtypes and/or groupby_list)
    '''

    joined_df = df1.join(df2).reset_index()

    if convert_dtypes:
        print('Converting data types.')
        joined_df = joined_df.convert_dtypes()
    
    if groupby_list:
        
        if stat=='mean':
            joined_df = joined_df.groupby(by=groupby_list).mean()
        else:
            joined_df = joined_df.groupby(by=groupby_list).median()
    
    return joined_df


# --------------------------------------------------
def get_transcript_list(df, substring):
    '''
    Returns a list of transcript names. 

    Input: 
        - df: Dataframe containing transcript names as columns 
        - substring: shared substring within transcripts, for example "Gh_"
    Output: 
        - transcript_list: list of all transcripts contained within the dataframe
    '''

    transcript_list = [transcript for transcript in df.columns if substring in transcript]
    
    return transcript_list 


# --------------------------------------------------
def generate_save_transcript_csv(transcript):
    '''
    Generates RNA-Seq merged file for subsequent model training.

    Input: 
        - transcript: a single transcript name to process
    Output: 
        - CSV file saved to disk, within the specified output directory command line argument 
    '''

    args = get_args()
    # Collect response (RNAseq TPM) and explanatory variables (spectra) for a single transcript
    temp_data = rnaseq_spectra.dropna(subset=[transcript])
    temp_data = temp_data.dropna(subset=['350', '2500'])
    
    
    column_list = [str(i) for i in range(350, 2501)]
    column_list.insert(0, transcript)
    result = temp_data[column_list]
    
    if not os.path.isdir(args.out_dir):
        try:
            os.makedirs(args.out_dir)
        except OSError:
            pass

    if 'TPM' in args.rnaseq_csv:
        out_file = os.path.join(args.out_dir, '_'.join([transcript, 'tpm_spectra.csv']))

    else:
        out_file = os.path.join(args.out_dir, '_'.join([transcript, 'logfc_spectra.csv']))

    result.to_csv(out_file)


# --------------------------------------------------
def main():
    """
    Run the RNA-Seq and Hyperspectral merging here.
    """

    args = get_args()
    # Get RNA-seq data
    rna_fb = get_rnaseq_data(csv_path=args.rnaseq_csv,
                            fb_path=args.fieldbook_csv)
    
    rna_fb = rna_fb[rna_fb['treatment'].isin(args.treatment)]

    # Get spectral data
    spectra_df, spectra_fb, fb = get_spectral_data(csv_path=args.spectral_csv,
                                fb_path=args.fieldbook_csv)

    # Join RNA-seq and spectral data into single dataframe
    global rnaseq_spectra
    rnaseq_spectra = join_dataframes(df1=rna_fb, df2=spectra_df, groupby_list=['entry', 'treatment'], stat='mean')


    # Generate and save transcript CSV
    transcript_list = get_transcript_list(df=rnaseq_spectra, substring='Gh_')

    with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
        p.map(generate_save_transcript_csv, transcript_list)
        
    print('Processing complete.')


# --------------------------------------------------
if __name__ == '__main__':
    main()
