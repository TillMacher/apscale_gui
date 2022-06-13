import subprocess, gzip, datetime, pickle, glob, os, openpyxl, shutil, math
import pandas as pd
from plotly.subplots import make_subplots
from pathlib import Path
from joblib import Parallel, delayed
import plotly.graph_objects as go
import plotly.express as px
from statistics import mean
from statistics import median
from statistics import stdev
import numpy as np
import PySimpleGUI as sg

############################################################################
## general function

def calculate_read_stats(lst):
    minimum = min(lst)
    maximum = max(lst)
    average = round(mean(lst),2)
    med = round(median(lst),2)
    deviation = round(stdev(lst),2)

    return([minimum, maximum, average, med, deviation])

## Plot a heatmap of the ESV or OTU table
def plot_heatmap(file, unit, project):
    """ Function to plot a heatmap of the OTU or ESV table """

    ## load OTU table and extract read numbers
    df = pd.read_excel(file)
    ID_list = df['ID'].values.tolist()
    df = df.drop(columns=['ID', 'Seq'])
    Sample_list = df.columns.tolist()

    z = []
    for i in df.values.tolist():
        sublist=[]
        for num in i:
            if num == 0:
                sublist.append(0)
            else:
                sublist.append(math.log(num))
        z.append(sublist)

    ## create heatmap
    fig = px.imshow(z, y=ID_list, x=Sample_list, aspect="auto")

    ## calculate optimal height and trim maximum height and adjust layout
    h = len(ID_list)*15
    if h >= 3000:
        h = 3000
        fig.update_layout(template='simple_white', width=1500, height=h, title='log('+unit+')', coloraxis_showscale=False)
        fig.update_yaxes(showticklabels=False, title=unit+'s')
    else:
        fig.update_layout(template='simple_white', width=1500, height=h, title='log('+unit+')', coloraxis_showscale=False)
        fig.update_yaxes(tickmode='linear')
        fig.update_xaxes(tickmode='linear')

    ## write image
    out_pdf = Path(project).joinpath('0_statistics', 'Summary_statistics', unit + '_heatmap.pdf')
    fig.write_image(str(out_pdf))
    out_html = Path(project).joinpath('0_statistics', 'Summary_statistics', unit + '_heatmap.html')
    fig.write_html(str(out_html))

def plot_reads(Project_report, project, ESV_table_lulu, OTU_table_lulu):

    df_3 = pd.read_excel(Project_report, sheet_name='3_PE merging')
    df_4 = pd.read_excel(Project_report, sheet_name='4_primer_trimming')
    df_5 = pd.read_excel(Project_report, sheet_name='5_quality_filtering')
    df_7 = pd.read_excel(Project_report, sheet_name='7_otu_clustering')
    df_9_OTUs = pd.read_excel(OTU_table_lulu)
    df_9_ESVs = pd.read_excel(ESV_table_lulu)

    samples = df_7['File'].values.tolist()
    raw_reads = df_3['processed reads'].values.tolist()
    merged_reads = df_3['merged reads'].values.tolist()
    trimmed_reads = df_4['trimmed reads'].values.tolist()
    filtered_reads = df_5['passed reads'].values.tolist()
    mapped_reads_OTUs = [sum(df_9_OTUs[sample].values.tolist()) for sample in samples]
    mapped_reads_ESVs = [sum(df_9_ESVs[sample].values.tolist()) for sample in samples]

    stats_raw_reads = calculate_read_stats(raw_reads)
    stats_merged_reads = calculate_read_stats(merged_reads)
    stats_trimmed_reads = calculate_read_stats(trimmed_reads)
    stats_filtered_reads = calculate_read_stats(filtered_reads)
    stats_mapped_reads_OTUs = calculate_read_stats(mapped_reads_OTUs)
    stats_mapped_reads_ESVs = calculate_read_stats(mapped_reads_ESVs)

    ## dataframe
    df_stats = pd.DataFrame()
    df_stats['Raw reads'] = raw_reads + stats_raw_reads
    df_stats['Merged reads'] = merged_reads + stats_merged_reads
    df_stats['Trimmed reads'] = trimmed_reads + stats_trimmed_reads
    df_stats['Filtered reads'] = filtered_reads + stats_filtered_reads
    df_stats['Mapped reads (OTUs)'] = mapped_reads_OTUs + stats_mapped_reads_OTUs
    df_stats['Mapped reads (ESVs)'] = mapped_reads_ESVs + stats_mapped_reads_ESVs
    df_stats.index = samples + ['_minimum', '_maximum', '_average', '_median', '_deviation']
    out_xlsx = Path(project).joinpath('0_statistics', 'Summary_statistics/summary_stats.xlsx')
    df_stats.to_excel(out_xlsx)

    ## plot OTUs
    fig = go.Figure()
    for sample in samples:
        y_values = df_stats.loc[sample].values.tolist()[:-1]
        x_values = df_stats.columns.tolist()[:-1]
        fig.add_trace(go.Scatter(x=x_values, marker_color='navy', y=y_values, name=sample))
    fig.update_layout(template='simple_white', width=1000, height=800, title='Reads per sample for each module')
    fig.update_yaxes(title='Reads')
    out_html = Path(project).joinpath('0_statistics', 'Summary_statistics/{}_scatter.html'.format('OTUs'))
    fig.write_html(str(out_html))

    ## plot ESVs
    fig = go.Figure()
    for sample in samples:
        y_values = df_stats.loc[sample].values.tolist()[:-2] + [df_stats.loc[sample].values.tolist()[-1]]
        x_values = df_stats.columns.tolist()[:-2] + [df_stats.columns.tolist()[-1]]
        fig.add_trace(go.Scatter(x=x_values, marker_color='navy', y=y_values, name=sample))
    fig.update_layout(template='simple_white', width=1000, height=800, title='Reads per sample for each module')
    fig.update_yaxes(title='Reads')
    out_html = Path(project).joinpath('0_statistics', 'Summary_statistics/{}_scatter.html'.format('ESVs'))
    fig.write_html(str(out_html))

    ## plot stats
    fig = go.Figure()
    for category in df_stats.columns.tolist():
        y_values = df_stats.loc[samples][category].values.tolist()
        fig.add_trace(go.Box(y=y_values, name=category, marker_color = 'navy'))
    fig.update_layout(template='simple_white', width=1000, height=800, title='Reads per module')
    fig.update_yaxes(title='Reads')
    out_html = Path(project).joinpath('0_statistics', 'Summary_statistics/boxplot.html')
    fig.write_html(str(out_html))

def lulu_stats(ESV_table, OTU_table, ESV_table_lulu, OTU_table_lulu):
    ESVs = len(pd.read_excel(ESV_table))
    ESVs_filtered = len(pd.read_excel(ESV_table_lulu))
    OTUs = len(pd.read_excel(OTU_table))
    OTUs_filtered = len(pd.read_excel(OTU_table_lulu))

    fig = go.Figure()
    x_values = ['ESVs', 'ESVs filtered', 'OTUs', 'OTUs filtered']
    y_values = [ESVs, ESVs_filtered, OTUs, OTUs_filtered]
    text = [ESVs, ESVs_filtered, OTUs, OTUs_filtered]
    fig.add_trace(go.Bar(x=x_values, y=y_values, marker_color='navy', text=text))
    fig.update_layout(template='simple_white', width=1000, height=800, title='LULU filtering')
    fig.update_traces(textposition='outside')
    fig.update_yaxes(title='OTUs/ESVs')
    out_html = Path(project).joinpath('0_statistics', 'Summary_statistics/LULU_filtering_stats.html')
    fig.write_html(str(out_html))

## main function to call the script
def main(project = Path.cwd()):
    """Main function of the script. Default working directory is the current working directory."""

    print('{}: Starting to collect summary statistics.'.format(datetime.datetime.now().strftime("%H:%M:%S")))

    ## create output folders
    try:
        os.mkdir(Path(project).joinpath('0_statistics'))
    except FileExistsError:
        pass
    try:
        os.mkdir(Path(project).joinpath('0_statistics', 'Summary_statistics'))
    except FileExistsError:
        pass

    ## get name of the project
    project_name = Path(project).stem

    ## collect required files
    Project_report = Path(project).joinpath('Project_report.xlsx')
    OTU_table = Path(project).joinpath('7_otu_clustering/', project_name + '_OTU_table.xlsx')
    ESV_table = Path(project).joinpath('8_denoising/', project_name + '_ESV_table.xlsx')
    OTU_table_lulu = Path(project).joinpath('9_lulu_filtering/otu_clustering', project_name + '_OTU_table_filtered.xlsx')
    ESV_table_lulu = Path(project).joinpath('9_lulu_filtering/denoising', project_name + '_ESV_table_filtered.xlsx')

    ############################################################################
    ## OTUs and ESVs
    plot_heatmap(OTU_table_lulu, 'OTU', project)
    plot_heatmap(ESV_table_lulu, 'ESV', project)

    ############################################################################
    ## read progress
    plot_reads(Project_report, project, ESV_table_lulu, OTU_table_lulu)

if __name__ == "__main__":
    main()
