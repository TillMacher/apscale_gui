import subprocess, gzip, datetime, pickle, glob, os, openpyxl, shutil
import pandas as pd
from plotly.subplots import make_subplots
from pathlib import Path
from joblib import Parallel, delayed
import plotly.graph_objects as go
import statistics
import numpy as np

## calculate weighted average
def weighted_avg_per_file(df):
    rate = df['Mean_EE'].values.tolist()
    amount = df['Recs'].values.tolist()
    return np.average(rate, weights=amount)

## quality statistics function for each file
def fastq_eestats(file, folder, project = None):
    """ Function to calculate fastq statistics """

    ## extract the filename from the sample path / name and convert to output name
    sample_name_out = Path(str(project) + '/0_statistics/' + folder + '/log/' + str(Path(file).with_suffix('').with_suffix('').name) + '.csv')

    ## run vsearch --fastq_eestats2 to get length distributions
    f = subprocess.run(['vsearch',
                    '--fastq_eestats', file,
                    '--output', sample_name_out,
                    '--threads', str(1)], capture_output = True)

    ## Give user output, if 0 reads are the output handle Zero division exception
    print(datetime.datetime.now().strftime("%H:%M:%S") + ': ' + str(Path(file).with_suffix('').with_suffix('').name) + ' - done')

    ################################################################################################
    ## figure 1
    ## max EE and PctRecs
    df = pd.read_csv(sample_name_out, delimiter='\t')
    x_values = df['Pos'].values.tolist()
    y1_values = df['Mean_Q'].values.tolist()
    y1_bad = len(x_values)*[18]
    y1_medium = len(x_values)*[9]
    y1_good = len(x_values)*[12]

    ## create new figure
    fig = go.Figure()

    # Add traces
    fig.add_trace(go.Bar(x=x_values, y=y1_bad, name="Bad", marker_color='rgb(229, 195, 196)', marker_line_color='rgb(229, 195, 196)', base=2))
    fig.add_trace(go.Bar(x=x_values, y=y1_medium, name="Medium", marker_color='rgb(230, 220, 196)', marker_line_color='rgb(230, 220, 196)', base=20))
    fig.add_trace(go.Bar(x=x_values, y=y1_good, name="Good", marker_color='rgb(196, 230, 196)', marker_line_color='rgb(196, 230, 196)', base=29))
    fig.add_trace(go.Scatter(x=x_values, y=y1_values, name='Mean quality score', line=dict(color='black', width=1)))

    # Set axis titles
    fig.update_xaxes(title_text="Position (bp)")
    fig.update_yaxes(title_text="Quality score", range=(2,44))
    fig.update_layout(template='simple_white', width=800, height=500, showlegend=True, barmode='stack')

    ## write image
    out_pdf = Path(str(project) + '/0_statistics/' + folder + '/plot_files/qscore_' + str(Path(file).with_suffix('').with_suffix('').name) + '.pdf')
    fig.write_image(str(out_pdf))

    ################################################################################################
    ## figure 2
    ## calculate the number of reads per base (i.e. length distriubution)
    n_recs = df['Recs'].values.tolist()[0]
    len_dist_dict = {}
    for i,j in df[['Pos', 'Recs']].values.tolist():
        len_dist_dict[i] = n_recs - j
        n_recs -= n_recs - j
    mean_ee = df['Mean_EE'].values.tolist()
    recs = df['Recs'].values.tolist()
    ee_1 = len(x_values)*[1]

    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # Add traces
    fig.add_trace(go.Scatter(x=x_values, y=mean_ee, name="Mean EE", marker_color="Grey"), secondary_y=True)
    fig.add_trace(go.Scatter(x=x_values, y=ee_1, name="EE = 1", line=dict(color='black', dash='dot')), secondary_y=True)
    fig.add_trace(go.Scatter(x=x_values, y=list(len_dist_dict.values()), name="Records (%)", marker_color="Blue"), secondary_y=False)
    fig.update_layout(template='simple_white', width=800, height=500)

    fig.update_yaxes(title_text="Mean EE", secondary_y=True)
    fig.update_yaxes(title_text="Reads", secondary_y=False)
    fig.update_xaxes(title_text="Position (bp)")

    ## write image
    out_pdf = Path(str(project) + '/0_statistics/' + folder + '/plot_files/maxee_' + str(Path(file).with_suffix('').with_suffix('').name) + '.pdf')
    fig.write_image(str(out_pdf))

## quality statistics function for the whole dataset
def dataset_stats(folder, input, project = None):
    ## extract the filename from the sample path / name and convert to output name
    files = [Path(str(project) + '/0_statistics/' + folder + '/log/' + str(Path(file).with_suffix('').with_suffix('').name) + '.csv') for file in input]

    number_of_reads = [pd.read_csv(file, delimiter='\t')['Recs'][0] for file in files]
    file_names = [Path(file).with_suffix('').name for file in files]
    weighted_avg = [weighted_avg_per_file(pd.read_csv(file, delimiter='\t')) for file in files]
    ee_1 = len(file_names)*[1]

    weighted_avg, number_of_reads, file_names = zip(*sorted(zip(weighted_avg, number_of_reads, file_names)))

    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # Add traces
    fig.add_trace(go.Bar(x=file_names, y=number_of_reads, name="Reads", marker_color="Grey"), secondary_y=False)
    fig.add_trace(go.Scatter(x=file_names, y=ee_1, name="EE = 1", line=dict(color='black', dash='dot')), secondary_y=True)
    fig.add_trace(go.Scatter(x=file_names, y=weighted_avg, name="Mean EE", marker_color="Blue"), secondary_y=True)
    fig.update_layout(template='simple_white', width=800, height=500)

    fig.update_yaxes(title_text="Mean EE", secondary_y=True)
    fig.update_yaxes(title_text="Reads", secondary_y=False)
    fig.update_xaxes(title_text="Samples", tickmode='linear', showticklabels=False)

    ## write image
    out_pdf = Path(str(project) + '/0_statistics/' + folder + '/plot_dataset/mean_ee_reads.pdf')
    fig.write_image(str(out_pdf))
    out_html = Path(str(project) + '/0_statistics/' + folder + '/plot_dataset/mean_ee_reads.html')
    fig.write_html(str(out_html))

## main function to call the script
def main(folder, project = Path.cwd()):
    """Main function of the script.Default working directory is the current working directory."""

    ## collect variables from the settings file
    gen_settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '0_general_settings')
    cores, comp_lvl = gen_settings['cores to use'].item(), gen_settings['compression level'].item()

    ## collect the input files from primer trimming step
    input = glob.glob(str(Path(project).joinpath(folder, 'data', '*.fastq.gz')))

    print('{}: Starting to collect statistics for {} input files.'.format(datetime.datetime.now().strftime("%H:%M:%S"), len(input)))

    ## create output folders
    try:
        os.mkdir(Path(project).joinpath('0_statistics'))
    except FileExistsError:
        pass
    try:
        os.mkdir(Path(project).joinpath('0_statistics', folder))
    except FileExistsError:
        pass
    try:
        os.mkdir(Path(project).joinpath('0_statistics', folder, 'log'))
    except FileExistsError:
        pass
    try:
        os.mkdir(Path(project).joinpath('0_statistics', folder, 'plot_files'))
    except FileExistsError:
        pass
    try:
        os.mkdir(Path(project).joinpath('0_statistics', folder, 'plot_dataset'))
    except FileExistsError:
        pass

    ## parallelize the fastq_eestats for each file
    Parallel(n_jobs = cores)(delayed(fastq_eestats)(file, folder, project = project) for file in input)

    ## calculate whole data set statistics
    dataset_stats(folder, input, project = project)

if __name__ == "__main__":
    main()
