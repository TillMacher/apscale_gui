def v_read_filtering(path_to_outdirs, min_length, max_length, maxee_value, num_cores_to_use, main_log_file):

    # main_log_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/log.xlsx"
    # path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST"
    # min_length = 130
    # max_length = 210
    # maxee_value = 1

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    import pandas as pd

    ## define functions
    def count_fasta_reads(file):
        def blocks(files, size=65536):
            while True:
                b = files.read(size)
                if not b: break
                yield b
        with open(file, "r",encoding="utf-8",errors='ignore') as f:
            return sum(bl.count(">") for bl in blocks(f))

    ## define variables
    input_folder = Path(str(path_to_outdirs) + "/3_primer_trimming/_data/")
    input_files = sorted(glob.glob(str(input_folder) + "/*.fastq.gz"))
    n_files = len(input_files)

    ## collect the number of raw reads
    log_df = pd.read_excel(main_log_file)

    dirName = Path(str(path_to_outdirs) + "/4_read_filtering/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    data_folder = Path(str(path_to_outdirs) + "/4_read_filtering/_data")
    if not os.path.exists(data_folder):
        os.mkdir(Path(data_folder))

    log_folder = Path(str(path_to_outdirs) + "/4_read_filtering/_log")
    if not os.path.exists(log_folder):
        os.mkdir(Path(log_folder))

    if n_files == 0:
        sg.PopupError("Warning: No files found! Please check your the merged reads!")

    ## print standard starting message
    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting read filtering for " + str(n_files) + " files.")

    else:
        def vsearch_read_filtering_command(i, file, data_folder, log_folder, min_length, max_length, maxee_value, log_df):
            # file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/3_primer_trimming/_data/Dessau_MuldeFL_5A_180419_a.fastq.gz"
            sample_name = str(Path(file).stem).replace(".fastq", "")
            output_file = Path(str(data_folder) + "/" + sample_name + ".fasta")
            stderr_file = Path(str(log_folder) + "/" + sample_name + "_stderr.log")

            ## count the merged reads from the log file
            try:
                reads = log_df.loc[log_df['Sample ID'] == sample_name]["Primer trimming"].values.tolist()[0]
            except:
                reads = 0

            if reads == 0:
                print("(" + str(i+1) + "/" + str(n_files) + ")\tNo reads to process:", Path(file).stem)

            else:
                ## run vsearch command to filter reads according to maxee and length
                f = open(stderr_file, "w")
                subprocess.call(["vsearch", "-fastq_filter", file, "-fastaout", output_file, "-fastq_maxee" , str(maxee_value), "-fastq_minlen", str(max_length), "-fastq_maxlen", str(min_length), "-fastq_qmax", "64", "-sizeout"], stderr=f)
                f.close()

                ## count the written reads
                written_reads = count_fasta_reads(output_file)

                ## write the reads to the log file
                f = open(stderr_file, "a")
                f.write(str(written_reads))
                f.close()

                ## calculate percentage
                try:
                    written_reads = round(written_reads / reads * 100, 1)
                except:
                    written_reads = 0

                print("(" + str(i+1) + "/" + str(n_files) + ")\tPassed reads:", written_reads, "% -", sample_name)

        results = Parallel(n_jobs=int(num_cores_to_use))(delayed(vsearch_read_filtering_command)(i, file, data_folder, log_folder, min_length, max_length, maxee_value, log_df) for i, file in enumerate(input_files))

        ## store the merged reads for the log file
        filtered_reads_list = []

        for sample_name in log_df["Sample ID"].values.tolist():
            log_file = Path(str(log_folder) + "/" + sample_name + "_stderr.log")
            try:
                f = open(log_file, "r")
                for line in f:
                    pass
                filtered_reads_list.append(int(line))
            except:
                filtered_reads_list.append(0)

        ## write to log sheet
        log_df["Read filtering"] = filtered_reads_list
        log_df.to_excel(main_log_file, index=False)

        ## print standard closing message
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished read filtering for " + str(n_files) + " files.")
