import subprocess, gzip, datetime, pickle, glob, os, openpyxl, shutil
import pandas as pd
from plotly.subplots import make_subplots
from pathlib import Path
from joblib import Parallel, delayed
import plotly.graph_objects as go


## quality filtering function to quality filter the specified file
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
    out_pdf = Path(str(project) + '/0_statistics/' + folder + '/plot/qscore_' + str(Path(file).with_suffix('').with_suffix('').name) + '.pdf')
    fig.write_image(str(out_pdf))

    ################################################################################################
    ## figure 2
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
    out_pdf = Path(str(project) + '/0_statistics/' + folder + '/plot/maxee_' + str(Path(file).with_suffix('').with_suffix('').name) + '.pdf')
    fig.write_image(str(out_pdf))



## main function to call the script
def main(folder, project = Path.cwd()):
    """Main function of the script.Default working directory is the current working directory."""

    ## collect variables from the settings file
    gen_settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '0_general_settings')
    cores, comp_lvl = gen_settings['cores to use'].item(), gen_settings['compression level'].item()

    ## collect the input files from primer trimming step
    input = glob.glob(str(Path(project).joinpath(folder, 'data', '*.fastq.gz')))

    print('{}: Starting to collect statistics for {} input files.'.format(datetime.datetime.now().strftime("%H:%M:%S"), len(input)))

    ## create temporal output folder
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
        os.mkdir(Path(project).joinpath('0_statistics', folder, 'plot'))
    except FileExistsError:
        pass

    ## parallelize the fastq_eestats
    Parallel(n_jobs = cores)(delayed(fastq_eestats)(file, folder, project = project) for file in input)

if __name__ == "__main__":
    main()
