def v_derep_singletons_uchime(input_folder, path_to_outdirs, num_cores_to_use, large_file_option, print_handle, window):

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
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting dereplication and chimera removal for " + str(n_files) + " files.")
    window.refresh()

    ##################################
    ##################################
    ## dereplication and chimera removal command

    x_values, y_passed_reads, y_discarded_reads, y_chimeras, y_non_chimeras = [], [], [], [], []
    stats_file_dereplication = str(path_to_outdirs) + "/7_Clustering/_stats/dereplication.html"
    stats_file_chimeras = str(path_to_outdirs) + "/7_Clustering/_stats/chimeras.html"

    def derep_uchime_command(i, file, path_to_outdirs):

        dirname = str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/')
        if not os.path.exists(dirname):
            os.mkdir(dirname)

        file_derep = str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/') + Path(file).stem.replace('.fastq', "_derep.fasta")
        file_derep_singletons = str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/') + Path(file).stem.replace('.fastq', "_derep_singletons.fasta")
        file_derep_singletons_nochimeras =  str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/') + Path(file).stem.replace('.fastq', "_derep_singletons_nochimeras.fasta")
        stderr_file = str(path_to_outdirs) + "/7_Clustering/_log/" + Path(file).stem.replace('.fastq', ".log1")

        f = open(stderr_file, "w")
        # DEREPLICATION
        subprocess.call(["vsearch", "--derep_fulllength", file, "-output", file_derep, "-sizeout"], stderr=f)
        # DEREPLICATION, DISCARD SINGLETONS
        subprocess.call(["vsearch", "--derep_fulllength", file_derep, "-output", file_derep_singletons, "-sizein", "-sizeout", "-minuniquesize", "2"], stderr=f)
        # CHIMERA REMOVAL
        subprocess.call(["vsearch", "--uchime_denovo", file_derep_singletons, "--nonchimeras", file_derep_singletons_nochimeras], stderr=f)
        f.close()

        try:
            stderr_file = str(path_to_outdirs) + "/7_Clustering/_log/" + Path(file).stem.replace('.fastq', ".log1")

            # read the log file to collect the stats to draw a plot later
            log_file_to_read = open(stderr_file, "r")
            all_lines = log_file_to_read.readlines()

            unique_reads = all_lines[16].split()[0]
            discarded_reads = all_lines[16].split()[3]
            discarded_reads_pct = all_lines[16].split()[6]
            chimeras = all_lines[26].split()[1]
            non_chimeras = all_lines[26].split()[4]
            chimeras_pct = all_lines[26].split()[2]

            print("(" + str(i+1) + "/" + str(n_files) + ")\tDereplicated", discarded_reads_pct, "Chimeras", chimeras_pct, "-", Path(file).stem.replace('.fastq', ""))

        except:
            print("(" + str(i+1) + "/" + str(n_files) + ")\tDiscarded", "(100%) - ", Path(file).stem.replace('.fastq', ""))

    results = Parallel(n_jobs=num_cores_to_use)(delayed(derep_uchime_command)(i, file, path_to_outdirs) for i, file in enumerate(input_files))

    for file in input_files:

        output_file =  str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/') + Path(file).stem.replace('.fastq', "_derep_singletons_nochimeras.fasta")
        stderr_file = str(path_to_outdirs) + "/7_Clustering/_log/" + Path(file).stem.replace('.fastq', ".log1")

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(stderr_file, "r")
        all_lines = log_file_to_read.readlines()

        try:

            unique_reads = all_lines[16].split()[0]
            discarded_reads = all_lines[16].split()[3]
            discarded_reads_pct = all_lines[16].split()[6]
            chimeras = all_lines[26].split()[1]
            non_chimeras = all_lines[26].split()[4]
            chimeras_pct = all_lines[26].split()[2]

            x_values.append(str(Path(file).stem))
            y_passed_reads.append(unique_reads)
            y_discarded_reads.append(discarded_reads)
            y_chimeras.append(chimeras)
            y_non_chimeras.append(non_chimeras)

        except:
            x_values.append(str(Path(file).stem))
            y_passed_reads.append(0)
            y_discarded_reads.append(0)
            y_chimeras.append(0)
            y_non_chimeras.append(0)

        # check if the outputfile is empty
        # if True delete it
        with open(output_file, 'r') as f:
            file_content = f.read()
            f.close()
        if file_content == '':
            print("Removed:", Path(file).stem)

            # remove all files
            file_derep = str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/') + Path(file).stem.replace('.fastq', "_derep.fasta")
            os.remove(Path(str(file_derep)))
            file_derep_singletons = str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/') + Path(file).stem.replace('.fastq', "_derep_singletons.fasta")
            os.remove(Path(str(file_derep_singletons)))
            file_derep_singletons_nochimeras =  str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/') + Path(file).stem.replace('.fastq', "_derep_singletons_nochimeras.fasta")
            os.remove(Path(str(file_derep_singletons_nochimeras)))

            # create an empty file
            discarded_file = Path(str(path_to_outdirs) + "/7_Clustering/_data_discarded/" + str(Path(file).stem))
            f = open(discarded_file, "w")
            f.write("")
            f.close()

            # remove the dir
            dirname = str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace('.fastq', '/')
            os.rmdir(dirname)

    # draw the stats plot
    fig = go.Figure(data=[
        go.Bar(name='passed reads', x=x_values, y=y_passed_reads),
        go.Bar(name='singletons', x=x_values, y=y_discarded_reads)
        ])
    fig.update_layout(barmode='stack')
    fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
    fig.write_html(str(stats_file_dereplication))

    # draw the stats plot
    fig = go.Figure(data=[
        go.Bar(name='non-chimeras', x=x_values, y=y_non_chimeras),
        go.Bar(name='chimeras', x=x_values, y=y_chimeras)
        ])
    fig.update_layout(barmode='stack')
    fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
    fig.write_html(str(stats_file_chimeras))

    ##################################
    ##################################
    ## pooling command

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting pooling of " + str(n_files) + " files.")
    window.refresh()

    input_folder = Path(input_folder)
    input_files = sorted(glob.glob(str(input_folder) + "/*.fastq.gz"))
    pooled_files_derep_folder = str(path_to_outdirs) + "/7_Clustering/_data/_clusters"
    pooled_files_fasta = str(path_to_outdirs) + "/7_Clustering/_data/_clusters/pooled_files.fasta"
    pooled_files_fasta_derep = str(path_to_outdirs) + "/7_Clustering/_data/_clusters/pooled_files_derep.fasta"
    pooled_files_fasta_derep_nochimeras = str(path_to_outdirs) + "/7_Clustering/_data/_clusters/pooled_files_derep_nochimeras.fasta"
    stats_file_pooled_files_fasta_derep_nochimeras = str(path_to_outdirs) + "/7_Clustering/_stats/pooled_files_dereplication_chimeras.html"
    stderr_file = str(path_to_outdirs) + "/7_Clustering/_log/pooled_files.log2"

    if not os.path.exists(pooled_files_derep_folder):
        os.mkdir(pooled_files_derep_folder)

    ######################################################
    # pool, dereplicate and uchime all files

    derep_files = sorted(glob.glob(str(path_to_outdirs) + "/7_Clustering/_data/*/*nochimeras.fasta"))

    with open(pooled_files_fasta, "wb") as outfile:
        for f in derep_files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())

    f = open(stderr_file, "w")


    ##################################
    ##################################
    ## dereplication of pooled files

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting dereplication of pooled files.")
    window.refresh()

    subprocess.call(["vsearch", "--derep_fulllength", pooled_files_fasta, "-output", pooled_files_fasta_derep, "-sizein", "-sizeout", "-minuniquesize", "2"], stderr=f)

    if large_file_option != True:
        # chimera REMOVAL of the pooled files only if user wants this option
        ## print standard starting message
        print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting chimera removal of pooled files.")
        window.refresh()

        subprocess.call(["vsearch", "--uchime_denovo", pooled_files_fasta_derep, "--nonchimeras", pooled_files_fasta_derep_nochimeras], stderr=f)
        f.close()

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(stderr_file, "r")
        all_lines = log_file_to_read.readlines()

        unique_reads = all_lines[6].split()[0]
        discarded_reads = str(int(all_lines[4].split()[3]) - int(unique_reads))
        unique_reads_pct = int(unique_reads) / (int(unique_reads) + int(discarded_reads)) * 100
        discarded_reads_pct = int(discarded_reads) / (int(unique_reads) + int(discarded_reads)) * 100

        chimeras = all_lines[17].split()[1]
        non_chimeras = all_lines[17].split()[4]
        chimeras_pct = int(chimeras) / (int(chimeras) + int(non_chimeras)) * 100
        non_chimeras_pct = int(non_chimeras) / (int(chimeras) + int(non_chimeras)) * 100

        y_passed_reads = [unique_reads_pct]
        y_discarded_reads = [discarded_reads_pct]
        y_chimeras = [chimeras_pct]
        y_non_chimeras = [non_chimeras_pct]

        log_file_to_read.close()

    else:
        ## else just pool the files and rename it
        ## print standard starting message
        print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Skipping dereplication of pooled files.")
        window.refresh()

        os.rename(pooled_files_fasta_derep, pooled_files_fasta_derep_nochimeras)

        # read the log file to collect the stats to draw a plot later
        log_file_to_read = open(stderr_file, "r")
        all_lines = log_file_to_read.readlines()

        unique_reads = all_lines[6].split()[0]
        discarded_reads = str(int(all_lines[4].split()[3]) - int(unique_reads))
        unique_reads_pct = int(unique_reads) / (int(unique_reads) + int(discarded_reads)) * 100
        discarded_reads_pct = int(discarded_reads) / (int(unique_reads) + int(discarded_reads)) * 100

        y_passed_reads = [unique_reads_pct]
        y_discarded_reads = [discarded_reads_pct]
        y_chimeras = [0]
        y_non_chimeras = [0]

        log_file_to_read.close()

    # draw the stats plot
    fig = go.Figure(data=[
        go.Bar(name='passed reads', x=["dereplication"], y=y_passed_reads),
        go.Bar(name='discarded reads', x=["dereplication"], y=y_discarded_reads),
        go.Bar(name='chimeras', x=["uchime"], y=y_chimeras),
        go.Bar(name='non-chimeras', x=["uchime"], y=y_non_chimeras)
        ])
    fig.update_layout(barmode='stack')
    fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
    fig.write_html(str(stats_file_pooled_files_fasta_derep_nochimeras))

    ## print standard closing message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished dereplication, chimera removal and pooling.")
    window.refresh()
    sg.Popup("Processed files are found under:", Path(str(path_to_outdirs) + "/7_Clustering/_data/"), title="Finished", keep_on_top=True)
