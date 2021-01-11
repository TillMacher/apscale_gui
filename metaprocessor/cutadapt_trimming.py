def cutadapt_trimming(input_folder, path_to_outdirs, p5_forward_primer, p7_reverse_primer, anchoring, num_cores_to_use, print_handle, window):

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go

    ## define variables
    input_folder = Path(input_folder)
    input_files = sorted(glob.glob(str(input_folder) + "/*.fastq.gz"))
    n_files = len(input_files)

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting primer trimming for " + str(n_files) + " files.")
    window.refresh()

    ## start command
    x_values, y_passed_reads, y_discarded_reads = [], [], []
    stats_file = str(path_to_outdirs) + "/4_Primer_trimming/_stats/primer_trimmed_reads.html"

    def cutadapt_command(i, file, path_to_outdirs, p5_forward_primer, p7_reverse_primer, anchoring):
        output_file_g = str(path_to_outdirs) + "/4_Primer_trimming/_data/" + Path(file).stem + "_g.gz"
        output_file_g_a = str(path_to_outdirs) + "/4_Primer_trimming/_data/" + Path(file).stem + ".gz"
        stdout_file = str(path_to_outdirs) + "/4_Primer_trimming/_log/" + Path(file).stem + ".log"
        stderr_file = str(path_to_outdirs) + "/4_Primer_trimming/_log/" + Path(file).stem + "_stderr.log"

        if anchoring == False:
            f = open(stdout_file, "w")
            f2 = open(stderr_file, "w")
            # first cut the forward primer and write it to a new file
            subprocess.call(["cutadapt", "-g", p5_forward_primer, "-o", output_file_g, file, "--discard-untrimmed", "--core", "1"], stdout=f, stderr=f2)
            # then cut the reverse primer and write it to the final file
            subprocess.call(["cutadapt", "-a", p7_reverse_primer, "-o", output_file_g_a, output_file_g, "--discard-untrimmed", "--core", "1"], stdout=f, stderr=f2)
            # remove the forwad only trimmed file
            os.remove(output_file_g)
            f.close()
            f2.close()
        elif anchoring == True:
            f = open(stdout_file, "w")
            f2 = open(stderr_file, "w")
            # first cut the forward primer and write it to a new file
            subprocess.call(["cutadapt", "-g", "^" + p5_forward_primer, "-o", output_file_g, file, "--discard-untrimmed", "--core", "1"], stdout=f, stderr=f2)
            # then cut the reverse primer and write it to the final file
            subprocess.call(["cutadapt", "-a", p7_reverse_primer + "$", "-o", output_file_g_a, output_file_g, "--discard-untrimmed", "--core", "1"], stdout=f, stderr=f2)
            # remove the forwad only trimmed file
            os.remove(output_file_g)
            f.close()
            f2.close()

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(stdout_file, "r")
        all_lines = log_file_to_read.readlines()
        processed_reads = [line.split(" ")[-1].replace("\n", "").replace(",", "") for line in all_lines if "Total reads processed:" in line][0]
        written_reads = [line.split(" ")[-2].replace(",", "") for line in all_lines if "Reads written (passing filters):" in line][1]
        written_reads_pct = round(int(written_reads) / int(processed_reads) * 100, 1)
        print("(" + str(i+1) + "/" + str(n_files) + ")\tTrimmed reads (" + str(written_reads_pct) + "%) -", Path(file).stem)
        log_file_to_read.close()

    results = Parallel(n_jobs=num_cores_to_use)(delayed(cutadapt_command)(i, file, path_to_outdirs, p5_forward_primer, p7_reverse_primer, anchoring) for i, file in enumerate(input_files))

    for file in input_files:
        output_file = str(path_to_outdirs) + "/4_Primer_trimming/_data/" + Path(file).stem + ".gz"
        stdout_file = str(path_to_outdirs) + "/4_Primer_trimming/_log/" + Path(file).stem + ".log"

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(stdout_file, "r")
        all_lines = log_file_to_read.readlines()

        total_reads = all_lines[7].split()[3].replace(',','')
        passed_reads = all_lines[9].split()[4].replace(',','')
        discarded_reads = float(total_reads) - float(passed_reads)
        written_reads_pct = all_lines[9].split()[5]

        log_file_to_read.close()

        x_values.append(str(Path(file).stem))
        y_passed_reads.append(passed_reads)
        y_discarded_reads.append(discarded_reads)

        # check if the outputfile is empty
        # if True delete it
        with gzip.open(output_file, 'rb') as f:
            file_content = f.read()
            f.close()
        if file_content == b'':
            print("Removed:", Path(file).stem)
            os.remove(Path(str(output_file)))
            discarded_file = Path(str(path_to_outdirs) + "/4_Primer_trimming/_data_discarded/" + str(Path(file).stem))
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

    ## print standard closing message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished primer trimming for " + str(n_files) + " files.")
    window.refresh()
    sg.Popup("Primer trimmed files are found under:", Path(str(path_to_outdirs) + "/4_Primer_trimming/_data/"), title="Finished", keep_on_top=True)
