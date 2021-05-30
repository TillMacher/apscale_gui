def cutadapt_trimming(path_to_outdirs, p5_forward_primer, p7_reverse_primer, anchoring, num_cores_to_use, main_log_file):

    # main_log_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/log.xlsx"
    # path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST"

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    import pandas as pd

    ## define variables
    input_folder = Path(str(path_to_outdirs) + "/2_pe_merging/_data/")
    input_files = sorted(glob.glob(str(input_folder) + "/*.fastq.gz"))
    n_files = len(input_files)

    ## collect the number of raw reads
    log_df = pd.read_excel(main_log_file)

    dirName = Path(str(path_to_outdirs) + "/3_primer_trimming/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    data_folder = Path(str(path_to_outdirs) + "/3_primer_trimming/_data")
    if not os.path.exists(data_folder):
        os.mkdir(Path(data_folder))

    log_folder = Path(str(path_to_outdirs) + "/3_primer_trimming/_log")
    if not os.path.exists(log_folder):
        os.mkdir(Path(log_folder))

    if n_files == 0:
        sg.PopupError("Warning: No files found! Please check your the merged reads!")

    else:
        ## print standard starting message
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting primer trimming for " + str(n_files) + " files.")

        ## start command
        x_values, y_passed_reads, y_discarded_reads = [], [], []
        stats_file = str(path_to_outdirs) + "/4_Primer_trimming/_stats/primer_trimmed_reads.html"

        def cutadapt_command(i, file, data_folder, log_folder, p5_forward_primer, p7_reverse_primer, anchoring, log_df):
            # file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/3_primer_trimming/_data/neg_2nd_4.fastq.gz"
            # path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST"
            # p5_forward_primer = "AAACTCGTGCCAGCCACC"
            # p7_reverse_primer = "CAAACTGGGATTAGATACCC"
            # anchoring = False

            sample_name = str(Path(file).stem).replace(".fastq", "")
            output_file = Path(str(data_folder) + "/" + sample_name + ".fastq.gz")
            stdout_file = Path(str(log_folder) + "/" + sample_name + ".log")
            stderr_file = Path(str(log_folder) + "/" + sample_name + "_stderr.log")

            ## count the merged reads from the log file
            try:
                reads = log_df.loc[log_df['Sample ID'] == sample_name]["Merged reads"].values.tolist()[0]
            except:
                reads = 0

            if reads == 0:
                print("(" + str(i+1) + "/" + str(n_files) + ")\tNo reads to process:", Path(file).stem)
            else:

                if anchoring == False:
                    primer_pair = p5_forward_primer + "..." + p7_reverse_primer
                    f = open(stdout_file, "w")
                    f2 = open(stderr_file, "w")
                    # first cut the forward primer and write it to a new file
                    subprocess.call(["cutadapt", "-a", primer_pair, "-o", str(output_file), str(file), "--discard-untrimmed", "--core", "1"], stdout=f, stderr=f2)
                    f.close()
                    f2.close()
                elif anchoring == True:
                    primer_pair = "^" + p5_forward_primer + "..." + p7_reverse_primer
                    f = open(stdout_file, "w")
                    f2 = open(stderr_file, "w")
                    # first cut the forward primer and write it to a new file
                    subprocess.call(["cutadapt", "-a", primer_pair, "-o", str(output_file), str(file), "--discard-untrimmed", "--core", "1"], stdout=f, stderr=f2)
                    f.close()
                    f2.close()

                ## count the written reads
                try:
                    with gzip.open(output_file, 'rb') as f:
                        for j, l in enumerate(f):
                            pass
                    written_reads = int((j+1) / 4)
                except:
                    written_reads = 0

                ## write the reads to the log file
                f = open(stderr_file, "a")
                f.write(str(written_reads))
                f.close()

                ## calculate percentage
                try:
                    written_reads = round(written_reads / reads * 100, 1)
                except:
                    written_reads = 0

                print("(" + str(i+1) + "/" + str(n_files) + ")\tWritten:", written_reads , "% -", Path(file).stem)

        results = Parallel(n_jobs=int(num_cores_to_use))(delayed(cutadapt_command)(i, file, str(data_folder), str(log_folder), p5_forward_primer, p7_reverse_primer, anchoring, log_df) for i, file in enumerate(input_files))

        ## store the merged reads for the log file
        trimmed_reads_list = []

        for file in log_df["Sample ID"].values.tolist():
            sample_name = str(Path(file).stem).replace(".fastq", "")
            log_file = Path(str(log_folder) + "/" + sample_name + "_stderr.log")

            try:
                f = open(log_file, "r")
                for line in f:
                    pass
                trimmed_reads_list.append(int(line))
            except:
                trimmed_reads_list.append(0)

        ## write to log sheet
        log_df["Primer trimming"] = trimmed_reads_list
        log_df.to_excel(main_log_file, index=False)

        ## print standard closing message
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished primer trimming for " + str(n_files) + " files.")
