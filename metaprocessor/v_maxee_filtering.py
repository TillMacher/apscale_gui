def v_maxee_filtering(input_folder, path_to_outdirs, maxee_value, num_cores_to_use, print_handle, window):

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    import pandas as pd

    ## define variables
    input_folder = Path(input_folder)
    input_files = sorted(glob.glob(str(input_folder) + "/*.fastq.gz"))
    n_files = len(input_files)

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting quality filtering for " + str(n_files) + " files.")
    window.refresh()

    ## start main command
    x_values, y_passed_reads, y_discarded_reads = [], [], []
    stats_file = str(path_to_outdirs) + "/6_Quality_filtering/_stats/maxee_trimmed_reads.html"
    stats_file_xlsx = str(path_to_outdirs) + "/6_Quality_filtering/_stats/maxee_trimmed_reads.xlsx"

    def vsearch_maxee_filtering(i, file, path_to_outdirs, min_length):

        output_file = str(path_to_outdirs) + "/6_Quality_filtering/_data/" + Path(file).stem
        stderr_file = str(path_to_outdirs) + "/6_Quality_filtering/_log/" + Path(file).stem + ".log"

        f = open(stderr_file, "w")
        subprocess.call(["vsearch", "-fastq_filter", file, "-fastaout", output_file, "--fastq_maxee" , maxee_value, "-fastq_qmax", "64", "-sizeout"], stderr=f)
        f.close()

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(stderr_file, "r")
        all_lines = log_file_to_read.readlines()

        try:
            # read the log file to collect the stats to draw a plot later
            log_file_to_read = open(stderr_file, "r")
            all_lines = log_file_to_read.readlines()
            passed_reads = all_lines[4].split()[0]
            discarded_reads = all_lines[4].split()[7]
            written_reads_pct = "(" + str(int(int(passed_reads) / (int(passed_reads) + int(discarded_reads)) *100)) + "%)"
            log_file_to_read.close()
        except:
            written_reads_pct = "(0%)"

        log_file_to_read.close()

        print("(" + str(i+1) + "/" + str(n_files) + ")\tPassed reads", written_reads_pct, "-", Path(file).stem)

        # compress the file using gzip
        input = open(output_file, 'rb')
        s = input.read()
        input.close()
        output = gzip.GzipFile(Path(str(output_file) + ".gz"), 'wb')
        output.write(s)
        output.close()
        os.remove(output_file)

    results = Parallel(n_jobs=num_cores_to_use)(delayed(vsearch_maxee_filtering)(i, file, path_to_outdirs, maxee_value) for i, file in enumerate(input_files))

    for file in input_files:

        output_file = str(path_to_outdirs) + "/6_Quality_filtering/_data/" + Path(file).stem + ".gz"
        stderr_file = str(path_to_outdirs) + "/6_Quality_filtering/_log/" + Path(file).stem + ".log"

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(stderr_file, "r")
        all_lines = log_file_to_read.readlines()
        passed_reads = all_lines[4].split()[0]
        discarded_reads = all_lines[4].split()[7]
        x_values.append(str(Path(file).stem))
        y_passed_reads.append(passed_reads)
        y_discarded_reads.append(discarded_reads)
        log_file_to_read.close()

        # check if the outputfile is empty
        # if True delete it
        with gzip.open(output_file, 'rb') as f:
            file_content = f.read()
            f.close()
        if file_content == b'':
            print("Removed:", Path(file).stem)
            os.remove(Path(str(output_file)))
            discarded_file = str(path_to_outdirs) + "/6_Quality_filtering/_data_discarded/" + Path(file).stem
            f = open(discarded_file, "w")
            f.write("")
            f.close()

    # draw the stats plot
    fig = go.Figure(data=[
        go.Bar(name='passed reads', x=x_values, y=y_passed_reads),
        go.Bar(name='discarded reads', x=x_values, y=y_discarded_reads)
        ])
    fig.update_layout(barmode='stack')
    fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
    fig.write_html(str(stats_file))

    # stats_df = pd.DataFrame([[sample, int(passed_reads), int(discarded_reads)] for sample, passed_reads, discarded_reads in zip(x_values, y_passed_reads, y_discarded_reads)], columns=(["Sample", "Passed", "Discarded"]))
    # total_passed_reads = sum(stats_df["Passed"].values.tolist())
    # total_discarded_reads = sum(stats_df["Discarded"].values.tolist())
    # nc_reads = sum([entry[1] for entry in stats_df.values.tolist() if "nc_" in entry[0]])
    # sample_reads = sum([entry[1] for entry in stats_df.values.tolist() if "nc_" not in entry[0]])
    # stats_df = pd.DataFrame(stats_df.values.tolist() + [["Sum passed", total_passed_reads]] + [["Sum discarded", total_discarded_reads]] + [["Samples", sample_reads]] + [["Negative controls", nc_reads]], columns=(["Sample", "Passed", "Discarded"]))
    # stats_df.to_excel(stats_file_xlsx, index=False)

    ## print standard closing message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished quality filtering for " + str(n_files) + " files.")
    window.refresh()
    sg.Popup("Quality trimmed files are found under:", Path(str(path_to_outdirs) + "/5_Length_filtering/_data/"), title="Finished", keep_on_top=True)












#
