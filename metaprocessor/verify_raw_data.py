def verify_raw_data(raw_data_folder, path_to_outdirs, main_log_file):

    # raw_data_folder = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/1_raw_data/_data"
    # path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST"
    # main_log_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/log.xlsx"

    from pathlib import Path
    import datetime, glob, subprocess, gzip, os, re
    import PySimpleGUI as sg
    from metaprocessor import file_pairs
    import pandas as pd

    # collect the files from the inout folder
    raw_data_folder = Path(raw_data_folder)
    input_files = sorted(glob.glob(str(raw_data_folder) + "/*.fastq.gz"))
    pairs = [i for i in file_pairs.main(input_files) if len(i) == 2]
    n_files = len(input_files)
    n_pairs = len(pairs)

    ## print standard closing message
    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Counting raw reads for " + str(n_files) + " files.")
    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Found " + str(n_pairs) + " pairs.")

    ## store the number of raw reads in a dict
    reads_dict = {}

    if input_files != []:
        for i, pair in enumerate(pairs):
            file = Path(pair[0])
            sample_name = "_".join(str(pair).replace(".fastq.gz", "").split("/")[-1].split("_")[:-1])

            ## write the number of raw reads to the log file
            ## count the reads
            with gzip.open(file, 'rb') as f:
                for j, l in enumerate(f):
                    pass
            reads = int((j+1) / 4)
            reads_dict[sample_name] = reads
            print("(" + str(i+1) + "/" + str(n_pairs) + ")\tRaw reads:", reads, "-", sample_name)
    else:
        sg.PopupError("Warning: No .fastq.gz files were found!")

    log_df = pd.DataFrame(reads_dict.keys(), columns=["Sample ID"])
    log_df["Raw reads"] = reads_dict.values()
    log_df.to_excel(main_log_file, index=False)
