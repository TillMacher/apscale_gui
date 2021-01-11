def get_length_distribution(length_trimming_folder, path_to_outdirs, num_cores_to_use):

    import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import pandas as pd
    import webbrowser
    from Bio import SeqIO
    from Bio.Seq import Seq
    import plotly.graph_objects as go

    length_trimming_folder = Path(length_trimming_folder)
    input_files = sorted(glob.glob(str(length_trimming_folder) + "/*.fastq.gz"))
    n_files = len(input_files)
    stats_file = str(path_to_outdirs) + "/5_Length_trimming/_stats/length_distribution.html"

    def get_length_distribution_command(i, file, path_to_outdirs):

        output_file = str(path_to_outdirs) + "/5_Length_trimming/_log/" + Path(file).stem + "_length_dist.log"
        stderr_file = str(path_to_outdirs) + "/5_Length_trimming/_log/" + Path(file).stem + ".log"
        f = open(stderr_file, "w")
        subprocess.call(["vsearch", "--fastq_eestats2", file, "--output", output_file], stderr=f)
        f.close()

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(output_file, "r")
        all_lines = log_file_to_read.readlines()
        max_reads = int(all_lines[0].split()[4].replace(",", ""))
        avg_reads = round(float(all_lines[0].split()[6]))
        log_file_to_read.close()

        print("(" + str(i+1) + "/" + str(n_files) + ")\tMax length (" + str(max_reads) + ") Avg length (" + str(avg_reads) + ") -", Path(file).stem)

    results = Parallel(n_jobs=num_cores_to_use)(delayed(get_length_distribution_command)(i, file, path_to_outdirs) for i, file in enumerate(input_files))

    x_values, y_max_reads, y_avg_reads, y_size_1, y_size_2 = [], [], [], [], []

    for file in input_files:

        output_file = str(path_to_outdirs) + "/5_Length_trimming/_log/" + Path(file).stem + "_length_dist.log"

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(output_file, "r")
        all_lines = log_file_to_read.readlines()
        max_reads = int(all_lines[0].split()[4].replace(",", ""))
        avg_reads = round(float(all_lines[0].split()[6]))
        log_file_to_read.close()

        x_values.append(str(Path(file).stem))
        y_max_reads.append(max_reads)
        y_avg_reads.append(avg_reads)
        y_size_1.append(9)
        y_size_2.append(13)

    # draw the stats plot
    fig = go.Figure(data=[
        go.Scatter(name='max read length', x=x_values, y=y_max_reads, mode='markers',  marker=dict(size=y_size_1)),
        go.Bar(name='avg read length', x=x_values, y=y_avg_reads)
        ])
    fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
    fig.write_html(str(stats_file))

    sg.Popup("Length distribution plots is found under:", Path(str(path_to_outdirs) + "/5_Length_trimming/_stats/length_distribution.html"), title="Finished", keep_on_top=True)
