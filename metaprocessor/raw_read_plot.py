def raw_read_plot(raw_data_folder, path_to_outdirs, num_cores_to_use):

    import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import pandas as pd
    import webbrowser
    from Bio import SeqIO
    from Bio.Seq import Seq
    import plotly.graph_objects as go

    paired_end = True

    raw_data_folder = Path(raw_data_folder)
    input_files = sorted(glob.glob(str(raw_data_folder) + "/*.fastq.gz"))
    n_files = len(input_files)
    output_file = Path(str(path_to_outdirs) + "/0_raw_data/_stats/raw_reads.html")
    output_xlsx = Path(str(path_to_outdirs) + "/0_raw_data/_stats/raw_reads.xlsx")

    print("Information:\nOnly the read numbers for the forward reads are shown.\nThe reverse reads have the equal number of reads by default.\n")

    def vsearch_fastq_stats_command(i, file):

        stderr_file = Path(str(path_to_outdirs) + "/0_raw_data/_log/" + Path(file).stem + ".log")
        f = open(stderr_file, "w")
        subprocess.call(["vsearch", "--fastq_stats", file ],  stderr=f)
        f.close()

        stderr_file = Path(str(path_to_outdirs) + "/0_raw_data/_log/" + Path(file).stem + ".log")
        f = open(stderr_file, "r")
        all_lines = f.readlines()
        n_reads = all_lines[4].strip().split()[1]
        f.close()

        print("(" + str(i+1) + "/" + str(n_files) + ")\traw reads (" + str(n_reads) + ")", "-", Path(file).stem)

    results = Parallel(n_jobs=num_cores_to_use)(delayed(vsearch_fastq_stats_command)(i, file) for i, file in enumerate(input_files))

    reads_dict = {}

    for file in input_files:

        stderr_file = Path(str(path_to_outdirs) + "/0_raw_data/_log/" + Path(file).stem + ".log")
        f = open(stderr_file, "r")
        all_lines = f.readlines()
        n_reads = all_lines[4].strip().split()[1]
        f.close()

        reads_dict[Path(file).stem] = n_reads

    x_values = list(reads_dict.keys())
    y_values = list(reads_dict.values())
    fig = go.Figure([go.Bar(x=x_values, y=y_values)])
    fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
    fig.write_html(str(output_file))

    reads_df = pd.DataFrame([[key, int(value)] for key, value in reads_dict.items()], columns=(["Sample", "Reads"]))
    if paired_end == True:
        total_reads = sum(reads_df["Reads"].values.tolist()) / 2
        nc_reads = sum([entry[1] for entry in reads_df.values.tolist() if "nc_" in entry[0]]) / 2
        sample_reads = sum([entry[1] for entry in reads_df.values.tolist() if "nc_" not in entry[0]]) / 2
    else:
        total_reads = sum(reads_df["Reads"].values.tolist())
        nc_reads = sum([entry[1] for entry in reads_df.values.tolist() if "nc_" in entry[0]])
        sample_reads = sum([entry[1] for entry in reads_df.values.tolist() if "nc_" not in entry[0]])

    reads_df = pd.DataFrame(reads_df.values.tolist() + [["Sum", total_reads]] + [["Samples", sample_reads]] + [["Negative controls", nc_reads]], columns=(["Sample", "Reads"]))
    reads_df.to_excel(output_xlsx, index=False)

    sg.Popup("The raw read plot is found under:", Path(str(output_file)), title="Finished")









#
