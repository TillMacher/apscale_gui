def v_merge(raw_data_folder, path_to_outdirs, main_log_file, num_cores_to_use, minovlen, maxdiffs, maxdiffpct):

    raw_data_folder = "/Users/tillmacher/Desktop/MP_Projects/Projects/EGLV_diatoms/1_raw_data/_demultiplexing"
    path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/EGLV_diatoms/"
    num_cores_to_use = "4"
    minovlen = "5"
    maxdiffs = "199"
    maxdiffpct = "25"
    main_log_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/EGLV_diatoms/log.xlsx"

    import datetime, glob, subprocess, gzip, os, re
    import PySimpleGUI as sg
    from pathlib import Path
    import pandas as pd
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    from metaprocessor import file_pairs

    # collect the files from the inout folder
    raw_data_folder = Path(raw_data_folder)
    input_files = sorted(glob.glob(str(raw_data_folder) + "/*.fastq.gz"))

    ## collect the pairs
    pairs = [i for i in file_pairs.main(input_files) if len(i) == 2]
    n_files = len(input_files)
    n_pairs = len(pairs)

    ## collect the number of raw reads
    log_df = pd.read_excel(main_log_file)

    ## print standard closing message
    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting paired-end merging for " + str(n_files) + " files.")
    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Found " + str(n_pairs) + " pairs.")

    dirName = Path(str(path_to_outdirs) + "/2_pe_merging/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    data_folder = Path(str(path_to_outdirs) + "/2_pe_merging/_data")
    if not os.path.exists(data_folder):
        os.mkdir(Path(data_folder))

    log_folder = Path(str(path_to_outdirs) + "/2_pe_merging/_log")
    if not os.path.exists(log_folder):
        os.mkdir(Path(log_folder))

    # check if there is a PE-file for each sample
    if n_pairs * 2 != n_files:
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Error: Could not find all paired-end read files. Please check your data.")
        sg.PopupError("Error: Could not find all paired-end read files. Please check your data.", keep_on_top=True)

        # flat_pairs = [item for sublist in pairs for item in sublist]
        # [file for file in input_files if file not in flat_pairs]

    else:
        def vsearch_merging_command(i, pe_files, data_folder, log_folder, log_df):
            ## define file variables
            ## pe_files = pairs[0]
            r1_read, r2_read = pe_files[0], pe_files[1]
            sample_name = "_".join(str(r1_read).replace(".fastq.gz", "").split("/")[-1].split("_")[:-1])
            merged_file_fastq = Path(str(data_folder) + "/" + sample_name + ".fastq")
            log_file = Path(str(log_folder) + "/" + sample_name + ".log")

            ## call VSEARCH for the paired-end merging
            ## write the output to a temporary file
            f = open(log_file, "w")
            subprocess.call(["vsearch", "--fastq_mergepairs", str(r1_read), "--reverse", str(r2_read), "--fastqout", str(merged_file_fastq), "--fastq_maxdiffpct", str(maxdiffpct), "--fastq_maxdiffs", str(maxdiffs), "--fastq_minovlen", str(minovlen), "--quiet"], stderr=f)
            f.close()

            try:
                ## compress the file using gzip
                input = open(merged_file_fastq, 'rb')
                s = input.read()
                input.close()
                gzip_file = Path(str(merged_file_fastq) + ".gz")
                output = gzip.GzipFile(gzip_file, 'wb')
                output.write(s)
                output.close()

                ## count the raw reads from the log file
                try:
                    raw_reads = log_df.loc[log_df['Sample ID'] == sample_name]["Raw reads"].values.tolist()[0]
                except:
                    raw_reads = 0

                ## count the merged reads
                with gzip.open(gzip_file, 'rb') as f:
                    j = 0
                    for j, l in enumerate(f):
                        pass
                merged_reads = int((j+1) / 4)

                ## calculate percentage
                try:
                    written_reads = round(merged_reads / raw_reads * 100, 1)
                except:
                    written_reads = 0

                ## write to sample log file
                f = open(log_file, "a")
                f.write("\n" + str(merged_reads))
                f.close()

                ## remove fastq file and keep fastq.gz file
                os.remove(merged_file_fastq)

                ## print update
                print("(" + str(i+1) + "/" + str(n_pairs) + ")\tWritten:", written_reads , "% -", sample_name)

            except:
                ## write to sample log file
                f = open(log_file, "a")
                f.write("\n" + str(0))
                f.close()
                print("(" + str(i+1) + "/" + str(n_pairs) + ")\tWritten:", "0 % -", sample_name)



        results = Parallel(n_jobs=int(num_cores_to_use))(delayed(vsearch_merging_command)(i, pe_files, data_folder, log_folder, log_df) for i, pe_files in enumerate(pairs))

        ## store the merged reads for the log file
        merged_reads_list = []

        for pe_files in pairs:
            sample_name = "_".join(str(pe_files[0]).replace(".fastq.gz", "").split("/")[-1].split("_")[:-1])
            log_file = Path(str(log_folder) + "/" + sample_name + ".log")

            f = open(log_file, "r")
            for line in f:
                pass
            merged_reads_list.append(int(line))

        ## write to log sheet
        log_df["Merged reads"] = merged_reads_list
        log_df.to_excel(main_log_file, index=False)

        ## print standard closing message
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished paired-end merging for " + str(n_files) + " files.")
