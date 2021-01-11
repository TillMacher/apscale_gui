def v_merge(input_folder, path_to_outdirs, num_cores_to_use, minovlen, maxdiffs, maxdiffpct, print_handle, window):

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    from metaprocessor import file_pairs

    # collect the files from the inout folder
    input_folder = Path(input_folder)
    input_files = sorted(glob.glob(str(input_folder) + "/*.fastq.gz"))

    ## collect the pairs
    pairs = [i for i in file_pairs.main(input_files) if len(i) == 2]
    n_files = len(input_files)
    n_pairs = len(pairs)

    ## print standard closing message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting paired-end merging for " + str(n_files) + " files.")
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Found " + str(n_pairs) + " pairs.")
    window.refresh()

    # check if there is a PE-file for each sample
    if n_pairs * 2 != n_files:
        print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Error: Could not find all paired-end read files. Please check your data.")
        sg.PopupError("Error: Could not find all paired-end read files. Please check your data.", keep_on_top=True)

    else:
        x_values = []
        y_pairs_values = []
        y_merged_values = []
        y_not_merged_values = []
        stats_file = str(path_to_outdirs) + "/3_Paired-end_merging/_stats/merged_reads.html"

        def vsearch_merging_command(i, pe_files, path_to_outdirs):
            # define file variables
            sample_name = Path(pe_files[0]).stem.replace(".fastq", "")
            r1_read, r2_read = pe_files[0], pe_files[1]
            merged_file_fastq = Path(str(path_to_outdirs) + "/3_Paired-end_merging/_data/" + sample_name + ".fastq")
            log_file = Path(str(path_to_outdirs) + "/3_Paired-end_merging/_log/" + sample_name + ".log")

            # call VSEARCH for the paired-end merging
            # write the output to a temporary file
            f = open(log_file, "w")
            subprocess.call(["vsearch", "--fastq_mergepairs", str(r1_read), "--reverse", str(r2_read), "--fastqout", str(merged_file_fastq), "--fastq_maxdiffpct", maxdiffpct, "--fastq_maxdiffs", maxdiffs, "--fastq_minovlen", minovlen, "--quiet"], stderr=f)
            f.close()

            # read the log file to print progress
            f = open(log_file, "r")
            pairs = f.readline().strip().split()
            merged_pairs = f.readline().strip().split()
            not_merged_pairs = f.readline().strip().split()
            f.close()

            # compress the file using gzip
            input = open(merged_file_fastq, 'rb')
            s = input.read()
            input.close()
            output = gzip.GzipFile(Path(str(merged_file_fastq) + ".gz"), 'wb')
            output.write(s)
            output.close()
            os.remove(merged_file_fastq)

            print("(" + str(i+1) + "/" + str(n_files) + ")\tMerged", merged_pairs[2], "Not merged", not_merged_pairs[3], "-", sample_name)

        results = Parallel(n_jobs=num_cores_to_use)(delayed(vsearch_merging_command)(i, pe_files, path_to_outdirs) for i, pe_files in enumerate(pairs))

        for pe_files in pairs:
            sample_name = Path(pe_files[0]).stem.replace(".fastq", "")
            r1_read, r2_read = pe_files[0], pe_files[1]
            merged_file_fastq = Path(str(path_to_outdirs) + "/3_Paired-end_merging/_data/" + str(Path(pe_files[0]).stem) + ".fastq")
            log_file = Path(str(path_to_outdirs) + "/3_Paired-end_merging/_log/" + sample_name + ".log")

            # read the log file to collect the stats to draw a plot later
            f = open(log_file, "r")
            pairs = f.readline().strip().split()
            merged_pairs = f.readline().strip().split()
            not_merged_pairs = f.readline().strip().split()
            f.close()

            x_values.append(sample_name)
            y_pairs_values.append(pairs[0])
            y_merged_values.append(merged_pairs[0])
            y_not_merged_values.append(not_merged_pairs[0])

            if not_merged_pairs[3] == "(100.0%)":
                print("Removed:", sample_name)
                os.remove(Path(str(merged_file_fastq) + ".gz"))
                not_merged_file_fastq = Path(str(path_to_outdirs) + "/3_Paired-end_merging/_data_discarded/" + sample_name + ".fastq")
                f = open(not_merged_file_fastq, "w")
                f.write("")
                f.close()

        # draw the stats plot
        fig = go.Figure(data=[
            go.Bar(name='merged pairs', x=x_values, y=y_merged_values),
            go.Bar(name='not merged pairs', x=x_values, y=y_not_merged_values)
            ])
        fig.update_layout(barmode='stack')
        fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
        fig.write_html(str(stats_file))

        ## print standard closing message
        print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished paired-end merging for " + str(n_files) + " files.")
        window.refresh()
        sg.Popup("Merged files are found under:", Path(str(path_to_outdirs) + "/3_Paired-end_merging/_data/"), title="Finished", keep_on_top=True)
