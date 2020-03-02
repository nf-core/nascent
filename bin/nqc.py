from __future__ import print_function
import argparse
import numpy as np
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
import os
import re
import concurrent.futures
import warnings
import subprocess
from collections import defaultdict
from operator import itemgetter
from argparse import RawTextHelpFormatter
from io import StringIO

# Plotting imports
import bokeh.io as bk_io
from bokeh.plotting import figure, output_file, show, save
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, CDSView, BooleanFilter
from bokeh.transform import factor_cmap, factor_mark
from bokeh.palettes import Spectral6, Dark2, cividis, Viridis, viridis

def get_bedgraph_data(params):
    
    (bedgraph, genome) = params
    subsample_blacklist_vals =  ['_001', '_002', '_003', '_004', '_005', '_006', '_007', '_008', '_009', \
                                 '_01', '_02', '_03', '_04', '_05', '_06', '_07', '_08', '_09', \
                                 '_1', '_2', '_3', '_4', '_5', '_6', '_7', '_8', '_9']     
    sample_name = os.path.splitext(os.path.basename(bedgraph))[0]
    splt = sample_name.split('_')
    sample_id = '_'.join(splt[:-1]) if '_{}'.format(splt[-1:][0]) in subsample_blacklist_vals else '_'.join(splt) 
    
    coverage_file = pd.read_csv(bedgraph, comment='#', names=['chr','start','stop','coverage'], sep='\t', header=None \
                                 ).assign(sample=sample_name)
    coverage_variance = coverage_file['coverage'].var(ddof=1)
    
# Made a mistake previously where I did not take the abs -- resulted in linear relationship -- should explore this    
    log_coverage_variance = (np.log2(abs(coverage_file['coverage']))).var(ddof=1)
    
    complement = subprocess.Popen(['complementBed', '-i', bedgraph, '-g', genome], 
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)
    compout,comperr = complement.communicate()
    compdata = StringIO(str(compout,'utf-8')) 
    gaps=pd.read_csv(compdata, sep='\t', header=None, \
                  names = ['chr','start','end'])
    
    gaps['gap_size'] = gaps.end - gaps.start
    
    gap_std = gaps['gap_size'].std(ddof=1)
    log_gap_std = np.log2(gaps['gap_size']).std(ddof=1)
    gap_var = gaps['gap_size'].var(ddof=1)
    log_gap_var = np.log2(gaps['gap_size']).var(ddof=1)
    gap_median_size = gaps['gap_size'].median()    
    log_gap_median_size = np.log2((gaps['gap_size']).median())
    number_of_1bp_gaps = (gaps['gap_size'] == 0).sum()
    log_number_of_1bp_gaps = np.log2((gaps['gap_size'] == 0).sum()) 
    
    coverage_stats = { 'sample': sample_name, \
                       'sample_id': sample_id, \
                       'gap_std': float(gap_std), \
                       'log_gap_std': float(log_gap_std), \
                       'gap_var': float(gap_var), \
                       'log_gap_var': float(log_gap_var), \
                       'gap_median_size': int(gap_median_size), \
                       'log_gap_median_size': float(log_gap_median_size), \
                       'coverage_variance': float(coverage_variance), \
                       'log_coverage_variance': float(log_coverage_variance), \
                       'number_of_1bp_gaps': int(number_of_1bp_gaps), \
                       'log_number_of_1bp_gaps': float(log_number_of_1bp_gaps) }  
    
    print("Done with sample %s bedGraph." % sample_name)
    
    return coverage_stats

def get_mapped_complexity_data(dup_file):
    
    sample_name = os.path.splitext(os.path.basename(dup_file))[0]
    if sample_name.endswith('.marked_dup_metrics'):
        sample_name = sample_name[:-19]
    dup_file_df = pd.read_csv(dup_file, comment='#', sep='\t', usecols = range(1,9), skiprows=[6], header=None, \
                                names=['unpaired_reads','paired_reads','secondary_reads', \
                                'unmapped_reads', 'unpaired_read_dups', 'paired_read_dups', \
                                'read_pair_optical_dups', 'percent_dup'] \
                                 ).assign(sample=sample_name)
    
    print("Done with sample %s mapped duplication file." % sample_name)
    
    dup_stats = { 'sample': sample_name, \
                         'unpaired_reads': int(dup_file_df.iloc[0]['unpaired_reads']), \
                         'paired_reads': int(dup_file_df.iloc[0]['paired_reads']), \
                         'secondary_reads': int(dup_file_df.iloc[0]['secondary_reads']), \
                         'unmapped_reads': int(dup_file_df.iloc[0]['unmapped_reads']), \
                         'unpaired_read_dups': int(dup_file_df.iloc[0]['unpaired_read_dups']), \
                         'paired_read_dups': int(dup_file_df.iloc[0]['paired_read_dups']), \
                         'read_pair_optical_dups': int(dup_file_df.iloc[0]['read_pair_optical_dups']), \
                         'percent_dup': float(dup_file_df.iloc[0]['percent_dup']) }
    return dup_stats

def get_read_complexity_data(dup_file):
    
    sample_name = os.path.splitext(os.path.basename(dup_file))[0]
    if sample_name.endswith('.fastqc_stats'):
        sample_name = sample_name[:-13]    
    read_dup_file_df = pd.read_csv(dup_file, comment='#', sep='\t')
    
    print("Done with sample %s read duplication file." % sample_name)
    
    read_dup_stats = { 'sample': read_dup_file_df.iloc[0]['SRR'], \
                         'percent_gc': int(read_dup_file_df.iloc[0]['%GC']), \
                         'total_sequences': int(read_dup_file_df.iloc[0]['Total_Sequences']), \
                         'percent_deduplicated': float(read_dup_file_df.iloc[0]['%Total_Deduplicated']) }
    return read_dup_stats


def main():
    parser = argparse.ArgumentParser(description='Nascent QC\n\n===============================================================================================================\n\nAn algorithm which provides a metric for nascent data quality.', epilog='@Dowell Lab, Margaret Gruca, margaret.gruca@colorado.edu\nFor questions and issues, see https://github.com/Dowell-Lab/nqc', usage='%(prog)s --bedgraph /my/sample/dir/*.bedGraph --output /my/out/dir/', formatter_class=RawTextHelpFormatter)
    
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    required.add_argument('-b', '--bedgraph', dest='bedgraphdir', \
                        help='Path to bedgraph file(s) of interest. Wildcard may be used to process multiple bedGraphs simultaneously with one output report.', required=True)
    
    required.add_argument('-d', '--mapped-duplicates', dest='mappeddupdir', \
                        help='Path to picard duplication summary file(s). Wildcard may be used to process multiple duplication files simultaneously.', required=True)
    
    required.add_argument('-r', '--read-duplicates', dest='readdupdir', \
                        help='Path to fastqc duplication summary file(s). Wildcard may be used to process multiple duplication files simultaneously.', required=True)        
    
    required.add_argument('-o', '--output', dest='output', \
                        help='Path to where output stats file and plots will be saved.', required=False, default='./')
    
    required.add_argument('-g', '--genome', dest='genome', \
                        help='Chromosome sizes file -- bedGraph file sort must match genome file. See README for details on generating this file for your specific genome.', \
                        default='', required=False, type=str)
    
    optional.add_argument('-png', '--save-png', dest='save_png', action='store_true', \
                        help='Save png in addition to the default html output.', default=False, required=False)
    
    optional.add_argument('-t', '--threads', type=int, dest='threads', metavar='<THREADS>', \
                        help='Number of threads for multi-processing. Default=1', default=1, required=False)  
    
    args = parser.parse_args()
    
    genome = args.genome
    warnings.simplefilter(action='ignore', category=(FutureWarning,RuntimeWarning))
    
    threads = args.threads
    bedgraph_list = glob.glob(args.bedgraphdir + '*.bedGraph')
    mapped_duplication_file_list = glob.glob(args.mappeddupdir + '*.marked_dup_metrics.txt')
    read_duplication_file_list = glob.glob(args.readdupdir + '*.fastqc_stats.txt')    
    
    #subsample_blacklist_vals =  ['_001', '_002', '_003', '_004', '_005', '_006', '_007', '_008', '_009', \
    #                             '_01', '_02', '_03', '_04', '_05', '_06', '_07', '_08', '_09', \
    #                             '_1', '_2', '_3', '_4', '_5', '_6', '_7', '_8', '_9']    
    
    print("Starting to process bedGraphs...\n" + str(datetime.datetime.now()))
    
    with concurrent.futures.ProcessPoolExecutor(threads) as executor:
        jobs = [executor.submit(get_bedgraph_data, [file,genome])
                for file in bedgraph_list]
        coverage_stats = [r.result() for r in jobs]
        
    with concurrent.futures.ProcessPoolExecutor(threads) as executor:
        jobs = [executor.submit(get_mapped_complexity_data, file)
                for file in mapped_duplication_file_list]
        mapped_duplication_stats = [r.result() for r in jobs]
        
    with concurrent.futures.ProcessPoolExecutor(threads) as executor:
        jobs = [executor.submit(get_read_complexity_data, file)
                for file in read_duplication_file_list]
        read_duplication_stats = [r.result() for r in jobs]        
    
    d = defaultdict(dict)
    for l in (coverage_stats, mapped_duplication_stats):
        for elem in l:
            d[elem['sample']].update(elem)
    mapdup_bg_stats = d.values()
    
    e = defaultdict(dict)
    for l in (mapdup_bg_stats, read_duplication_stats):
        for elem in l:
            e[elem['sample']].update(elem)
    nqc_stats = e.values()
    
    sorted_nqc_stats = sorted(nqc_stats, key=itemgetter('sample'))

    nqc_stats_file = open(args.output + '/nqc_stats.txt', 'w')
    for stat in sorted_nqc_stats:
        nqc_stats_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                        (stat['sample'], \
                         stat['sample_id'], \
                         stat['gap_std'], stat['log_gap_std'], \
                         stat['gap_var'], stat['log_gap_var'], \
                         stat['gap_median_size'], stat['log_gap_median_size'], \
                         stat['coverage_variance'], stat['log_coverage_variance'], \
                         stat['number_of_1bp_gaps'], stat['log_number_of_1bp_gaps'], \
                         stat['unpaired_reads'], stat['paired_reads'], stat['secondary_reads'], \
                         stat['unmapped_reads'], stat['unpaired_read_dups'], \
                         stat['paired_read_dups'], stat['read_pair_optical_dups'], \
                         stat['percent_dup'], stat['percent_gc'], stat['total_sequences'], stat['percent_deduplicated']))
    nqc_stats_file.close()
    
    #sample = [d['sample'] for d in sorted_nqc_stats]
    
    #identifier = []
    #for i in sample:
    #    splt = i.split('_')
    #    value = '_'.join(splt[:-1]) if '_{}'.format(splt[-1:][0]) in subsample_blacklist_vals else '_'.join(splt)
    #    identifier.append(value)
    #    
    #df["id"] = identifier
    
    df = pd.DataFrame(sorted_nqc_stats).replace([np.inf, -np.inf], np.nan).fillna(0)
    print(list(df.columns.values))
    df['log2_unique_mapped_reads'] = np.log2(df['unpaired_reads'] - df['unpaired_read_dups'])
    df['percent_raw_read_dup'] = 1 - df['percent_deduplicated']
    df['log2_unique_raw_reads'] = np.log2(df['total_sequences'] * df['percent_deduplicated'])
    df['log2_unmapped_reads'] = np.log2(df['unmapped_reads'])
    
    output_file(args.output + '/nqc_plot.html')
    source = ColumnDataSource(df)
    SAMPLE = df['sample_id'].unique()
    pal=viridis(len(SAMPLE))

    # create a view of the source for one plot to use
    #view = CDSView(source=source, filters=[BooleanFilter([True if y > 250 or y < 100 else False for y in y1])])
    
    TOOLS = "box_select,lasso_select,hover,help,save"
    #colors = factor_cmap('sample_id', palette=Category20b_20, factors=sample_id.unique()) 
    
    # create a new plot and add a renderer
    s1 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s1.xaxis.axis_label = 'Log\u2082 Number of Unique Reads'
    s1.yaxis.axis_label = 'Log\u2082 Gap Standard Deviation'       
    s1.circle('log2_unique_mapped_reads', 'log_gap_std', size=6, hover_color="firebrick", legend='sample_id',
                color=factor_cmap('sample_id', pal, SAMPLE),
                source=source)
    
    s2 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s2.xaxis.axis_label = 'Log\u2082 Number of Unique Reads'
    s2.yaxis.axis_label = "Log\u2082 Gap Variance"       
    s2.circle('log2_unique_mapped_reads', 'log_gap_var', size=6, hover_color="firebrick", legend='sample_id',
                  color=factor_cmap('sample_id', pal, SAMPLE),
                  source=source)
    
    s3 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s3.xaxis.axis_label = 'Log\u2082 Number of Unique Reads'
    s3.yaxis.axis_label = "Log\u2082 Gap Median Size" 
    s3.circle('log2_unique_mapped_reads', 'log_gap_median_size', size=6, hover_color="firebrick", legend='sample_id',
                    color=factor_cmap('sample_id', pal, SAMPLE),
                    source=source)
    
    s4 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s4.xaxis.axis_label = 'Log\u2082 Unmapped Reads'
    s4.yaxis.axis_label = "Log\u2082 Coverage Variance"
    s4.circle('log2_unmapped_reads', 'log_coverage_variance', size=6, hover_color="firebrick", legend='sample_id',
                    color=factor_cmap('sample_id', pal, SAMPLE),
                    source=source)
    
    s5 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s5.xaxis.axis_label = 'Log\u2082 Number of Unique Reads'
    s5.yaxis.axis_label = "Log\u2082 Coverage Variance"
    s5.circle('log2_unique_mapped_reads', 'log_coverage_variance', size=6, hover_color="firebrick", legend='sample_id',
                    color=factor_cmap('sample_id', pal, SAMPLE),
                    source=source)
    
    s6 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s6.xaxis.axis_label = 'Log\u2082 Number of Unique Reads'
    s6.yaxis.axis_label = "Log\u2082 Number of 1bp Gaps"      
    s6.circle('log2_unique_raw_reads', 'log_number_of_1bp_gaps', size=6, hover_color="firebrick", legend='sample_id',
                    color=factor_cmap('sample_id', pal, SAMPLE),
                    source=source)
    s6.legend.location = "bottom_right"
    
    s7 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s7.xaxis.axis_label = 'Log\u2082 Number of Unique Raw Reads'
    s7.yaxis.axis_label = "Log\u2082 Gap Median Size"      
    s7.circle('log2_unique_raw_reads', 'log_gap_median_size', size=6, hover_color="firebrick", legend='sample_id',
                    color=factor_cmap('sample_id', pal, SAMPLE),
                    source=source)
    
    s8 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s8.xaxis.axis_label = 'Log\u2082 Number of Unique Raw Reads'
    s8.yaxis.axis_label = "Log\u2082 Coverage Variance" 
    s8.circle('log2_unique_raw_reads', 'log_coverage_variance', size=6, hover_color="firebrick", legend='sample_id',
                    color=factor_cmap('sample_id', pal, SAMPLE),
                    source=source)
    
    s9 = figure(tools=TOOLS, plot_width=600, plot_height=600, title=None)
    s9.xaxis.axis_label = 'Log\u2082 Number of Unique Raw Reads'
    s9.yaxis.axis_label = "Log\u2082 Number of Unique Mapped Reads"
    s9.circle('log2_unique_raw_reads', 'log2_unique_mapped_reads', size=6, hover_color="firebrick", legend='sample_id',
                    color=factor_cmap('sample_id', pal, SAMPLE),
                    source=source)
    s9.legend.location = "bottom_right"
    
    # make a grid
    grid = gridplot([s1, s2, s3, s4, s5, s6, s7, s8, s9], ncols=3, plot_width=600, plot_height=600)
    save(grid, browser=args.output + 'nqc_plot.html')
    if args.save_png:
        bk_io.export_png(grid, filename='./nqc_plot.png')
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    ##ax2 = fig.add_subplot(3, 3, 2)
    ## We use ax parameter to tell seaborn which subplot to use for this plot
    #sns.scatterplot(y='log_gap_std', x='log2_unique_reads', hue='id', data=df)
    #sns.set_style("white")
    #sns.set_style("ticks")    
    #plt.tight_layout()    
    ##ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));
    #plt.savefig('./plot_test' + '.png')
    
    print("Done.\n" + str(datetime.datetime.now()))
    
    sys.exit(0)
    
if __name__ == '__main__':
    main()   
