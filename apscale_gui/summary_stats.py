import subprocess, gzip, datetime, pickle, glob, os, openpyxl, shutil, math
import pandas as pd
from plotly.subplots import make_subplots
from pathlib import Path
from joblib import Parallel, delayed
import plotly.graph_objects as go
import plotly.express as px
from statistics import mean
import numpy as np
import PySimpleGUI as sg

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

## Plot the number of reads of the whole dataset across the processing
def pipeline_stats_bar(Project_report, OTU_table, ESV_table, project):
    """ Function to plot a the number of processed reads along the processing """

    # Project_report = '/Users/tillmacher/Dropbox/Insta_Upload/Project_report.xlsx'
    # OTU_table = '/Users/tillmacher/Dropbox/Insta_Upload/Test_OTU_table.xlsx'
    # ESV_table = '/Users/tillmacher/Dropbox/Insta_Upload/Test_ESV_table.xlsx'

    ## raw reads and PE merging results
    pe_merging = pd.read_excel(Project_report, sheet_name='3_PE merging')
    reads_1 = pe_merging['processed reads'].values.tolist()
    reads_total = sum(reads_1)
    reads_2 = pe_merging['merged reads'].values.tolist()
    ## primer trimming results
    primer_trimming = pd.read_excel(Project_report, sheet_name='4_primer_trimming')
    reads_3 = primer_trimming['trimmed reads'].values.tolist()
    ## quality filtering results
    quality_filtering = pd.read_excel(Project_report, sheet_name='5_quality_filtering')
    reads_4 = quality_filtering['passed reads']
    ## Denoising results
    denoising = pd.read_excel(ESV_table).drop(columns=['ID', 'Seq'])
    reads_5 = [sum(denoising[sample].values.tolist()) for sample in denoising.columns.tolist()]
    ## OTU clustering results
    otu_clustering = pd.read_excel(OTU_table).drop(columns=['ID', 'Seq'])
    reads_6 = [sum(otu_clustering[sample].values.tolist()) for sample in otu_clustering.columns.tolist()]

    ## bar chart categories
    ## Bar no 1: Sum of reads left after the processing step
    ## Bar no 2: Difference to number of total reads
    # PE merging
    bar_1_1 = sum(reads_2)
    bar_1_2 = reads_total - bar_1_1
    # Primer trimming
    bar_2_1 = sum(reads_3)
    bar_2_2 = reads_total - bar_2_1
    # Quality filtering
    bar_3_1 = sum(reads_4)
    bar_3_2 = reads_total - bar_3_1
    # Denoising
    bar_4_1 = sum(reads_5)
    bar_4_2 = reads_total - bar_4_1
    # OTU clustering
    bar_5_1 = sum(reads_6)
    bar_5_2 = reads_total - bar_5_1

    fig = go.Figure()
    fig.add_trace(go.Bar(x=['PE merging'], y=[bar_1_1], marker_color='Teal'))
    fig.add_trace(go.Bar(x=['PE merging'], y=[bar_1_2], marker_color='Orange'))

    fig.add_trace(go.Bar(x=['Primer trimming'], y=[bar_2_1], marker_color='Teal'))
    fig.add_trace(go.Bar(x=['Primer trimming'], y=[bar_2_2], marker_color='Orange'))

    fig.add_trace(go.Bar(x=['Quality filtering'], y=[bar_3_1], marker_color='Teal'))
    fig.add_trace(go.Bar(x=['Quality filtering'], y=[bar_3_2], marker_color='Orange'))

    fig.add_trace(go.Bar(x=['Denoising'], y=[bar_4_1], marker_color='Teal'))
    fig.add_trace(go.Bar(x=['Denoising'], y=[bar_4_2], marker_color='Orange'))

    fig.add_trace(go.Bar(x=['OTU clustering'], y=[bar_5_1], marker_color='Teal'))
    fig.add_trace(go.Bar(x=['OTU clustering'], y=[bar_5_2], marker_color='Orange'))

    fig.update_layout(template='simple_white', width=800, height=500, showlegend=False, barmode='stack')
    fig.update_yaxes(title='Reads')

    ## write image
    out_pdf = Path(project).joinpath('0_statistics', 'Summary_statistics', 'Summary_statistics_bar.pdf')
    fig.write_image(str(out_pdf))
    out_html = Path(project).joinpath('0_statistics', 'Summary_statistics', 'Summary_statistics_bar.html')
    fig.write_html(str(out_html))

## Plot the number of reads per sample (as boxplot) across the processing
def pipeline_stats_box(Project_report, OTU_table, ESV_table, project):
    """ Function to plot a the number of processed reads along the processing """

    ## raw reads and PE merging results
    pe_merging = pd.read_excel(Project_report, sheet_name='3_PE merging')
    reads_1 = pe_merging['processed reads'].values.tolist()
    reads_total = sum(reads_1)
    reads_2 = pe_merging['merged reads'].values.tolist()
    ## primer trimming results
    primer_trimming = pd.read_excel(Project_report, sheet_name='4_primer_trimming')
    reads_3 = primer_trimming['trimmed reads'].values.tolist()
    ## quality filtering results
    quality_filtering = pd.read_excel(Project_report, sheet_name='5_quality_filtering')
    reads_4 = quality_filtering['passed reads']
    ## Denoising results
    denoising = pd.read_excel(ESV_table).drop(columns=['ID', 'Seq'])
    reads_5 = [sum(denoising[sample].values.tolist()) for sample in denoising.columns.tolist()]
    ## OTU clustering results
    otu_clustering = pd.read_excel(OTU_table).drop(columns=['ID', 'Seq'])
    reads_6 = [sum(otu_clustering[sample].values.tolist()) for sample in otu_clustering.columns.tolist()]

    fig = go.Figure()
    fig.add_trace(go.Box(y=reads_1, name='Raw reads', marker_color='Teal'))
    fig.add_trace(go.Box(y=reads_2, name='PE merging', marker_color='Teal'))
    fig.add_trace(go.Box(y=reads_3, name='Primer trimming', marker_color='Teal'))
    fig.add_trace(go.Box(y=reads_4, name='Quality filtering', marker_color='Teal'))
    fig.add_trace(go.Box(y=reads_5, name='Denoising', marker_color='Teal'))
    fig.add_trace(go.Box(y=reads_6, name='OTU clustering', marker_color='Teal'))

    fig.update_layout(template='simple_white', width=800, height=500, showlegend=False)
    fig.update_yaxes(title='Reads')

    ## write image
    out_pdf = Path(project).joinpath('0_statistics', 'Summary_statistics', 'Summary_statistics_box.pdf')
    fig.write_image(str(out_pdf))
    out_html = Path(project).joinpath('0_statistics', 'Summary_statistics', 'Summary_statistics_box.html')
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
    OTU_table = Path(project).joinpath('7_otu_clustering', project_name + '_OTU_table.xlsx')
    ESV_table = Path(project).joinpath('8_denoising', project_name + '_ESV_table.xlsx')

    ## check if all files exist
    if os.path.isfile(Project_report) and os.path.isfile(OTU_table) and os.path.isfile(ESV_table):
        ## Plot ESV and OTUs as a heatmap
        print('{}: Creating ESV heatmap.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
        plot_heatmap(ESV_table, 'ESV', project)
        print('{}: Creating OTU heatmap.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
        plot_heatmap(OTU_table, 'OTU', project)

        ## Plot the number of reads of the whole dataset across the processing
        print('{}: Creating bar charts.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
        pipeline_stats_bar(Project_report, OTU_table, ESV_table, project)

        ## Plot the number of reads per sample (as boxplot) across the processing
        print('{}: Creating box plots.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
        pipeline_stats_box(Project_report, OTU_table, ESV_table, project)

    else:
        print('{}: Warning: Please first run the whole pipeline and make sure that following files exist:\n -Project report\n -OTU table\n -ESV table'.format(datetime.datetime.now().strftime("%H:%M:%S")))
        sg.Popup('Warning: Please first run the whole pipeline and make sure that following files exist:\n -Project report\n -OTU table\n -ESV table', title='Error')

if __name__ == "__main__":
    main()
