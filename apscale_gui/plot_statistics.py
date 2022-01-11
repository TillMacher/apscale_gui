import pandas as pd
import PySimpleGUI as sg
import plotly.graph_objects as go
import subprocess, gzip, glob, pickle, datetime, psutil, os, shutil
import pandas as pd
from pathlib import Path
from demultiplexer import file_pairs
from joblib import Parallel, delayed

def plot_reads_processing(file):
    " Plot the report stats for all different modules "

    ## open file
    sheet_names = pd.ExcelFile(file).sheet_names

    ## open file
    if '3_PE merging' in sheet_names:
        df = pd.read_excel(file, sheet_name='3_PE merging')
        x1_values = [i[1]/i[0]*100 for i in df[['processed reads', 'merged reads']].values.tolist()]
        x2_values = [100-i for i in x1_values]
        y_values = [i.replace('_PE.fastq.gz', '') for i in df['File'].values.tolist()]

        fig = go.Figure(data=[
            go.Bar(name='Merged', x=x1_values, y=y_values, orientation='h', marker_color='Teal'),
            go.Bar(name='Discarded', x=x2_values, y=y_values, orientation='h', marker_color='Orange')])
        fig.update_layout(barmode='stack', template='plotly_white', width=1000, height=1000)
        fig.update_xaxes(title='reads (%)')
        fig.show()

    if '4_primer_trimming' in sheet_names:
        df = pd.read_excel(file, sheet_name='4_primer_trimming')
        x1_values = [i[1]/i[0]*100 for i in df[['processed reads', 'trimmed reads']].values.tolist()]
        x2_values = [100-i for i in x1_values]
        y_values = [i.replace('_PE.fastq.gz', '') for i in df['File'].values.tolist()]

        fig = go.Figure(data=[
            go.Bar(name='Trimmed', x=x1_values, y=y_values, orientation='h', marker_color='Teal'),
            go.Bar(name='Discarded', x=x2_values, y=y_values, orientation='h', marker_color='Orange')])
        fig.update_layout(barmode='stack', template='plotly_white', width=1000, height=1000)
        fig.update_xaxes(title='reads (%)')
        fig.show()

    if '5_quality_filtering' in sheet_names:
        df = pd.read_excel(file, sheet_name='5_quality_filtering')
        x1_values = [i[1]/i[0]*100 for i in df[['processed reads', 'passed reads']].values.tolist()]
        x2_values = [100-i for i in x1_values]
        y_values = [i.replace('_PE.fastq.gz', '') for i in df['File'].values.tolist()]

        fig = go.Figure(data=[
            go.Bar(name='Passed', x=x1_values, y=y_values, orientation='h', marker_color='Teal'),
            go.Bar(name='Discarded', x=x2_values, y=y_values, orientation='h', marker_color='Orange')])
        fig.update_layout(barmode='stack', template='plotly_white', width=1000, height=1000)
        fig.update_xaxes(title='reads (%)')
        fig.show()

    if '6_dereplication' in sheet_names:
        df = pd.read_excel(file, sheet_name='6_dereplication')
        x1_values = [i[1]/i[0]*100 for i in df[['processed sequences', 'unique sequences']].values.tolist()]
        x2_values = [100-i for i in x1_values]
        y_values = [i.replace('_PE.fastq.gz', '') for i in df['File'].values.tolist()]

        fig = go.Figure(data=[
            go.Bar(name='Unique', x=x1_values, y=y_values, orientation='h', marker_color='Teal'),
            go.Bar(name='Non-unique', x=x2_values, y=y_values, orientation='h', marker_color='Orange')])
        fig.update_layout(barmode='stack', template='plotly_white', width=1000, height=1000)
        fig.update_xaxes(title='sequences (%)')
        fig.show()
