def function_rename_samples(rename_data_folder, sample_renaming_sheet, path_to_outdirs, suffix_list):

    import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import pandas as pd
    import webbrowser
    from Bio import SeqIO
    from Bio.Seq import Seq
    import plotly.graph_objects as go

    suffix_list = [".1", ".2"]
    rename_data_folder = "/Users/tillmacher/Desktop/MP_Projects/Projects/Robot_LFU_data/0_raw_data/_data"
    sample_renaming_sheet = "/Users/tillmacher/Desktop/MP_Projects/Projects/Robot_LFU_data/rename_sheet.xlsx"
    path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/Robot_LFU_data/"

    sample_renaming_sheet = Path(sample_renaming_sheet)
    rename_sheet_df = pd.read_excel(sample_renaming_sheet)
    input_files = sorted(glob.glob(str(rename_data_folder) + "/*.fastq.gz"))

    files_list = []

    #\xa0

    for i, file in enumerate(rename_sheet_df.values.tolist()):
        # check if r1 reads file exists
        # then rename it
        file_r1 = Path(rename_data_folder + "/" + file[0] + suffix_list[0] + ".fastq.gz")
        if file_r1.is_file():
            new_name = Path(rename_data_folder + "/" + file[1] + suffix_list[0] + ".fastq.gz")
            os.rename(file_r1, new_name)
            print(file_r1.name + "\t>>\t" + new_name.name)
        else:
            print("Warning: r1 reads are missing -", file[0])

        # check if r2 reads file exists
        # then rename it
        file_r2 = Path(rename_data_folder + "/" + file[0] + suffix_list[1] + ".fastq.gz")
        if file_r2.is_file():
            new_name = Path(rename_data_folder + "/" + file[1] + suffix_list[1] + ".fastq.gz")
            os.rename(file_r2, new_name)
            print(file_r2.name + "\t>>\t" + new_name.name)
        else:
            print("Warning: r2 reads are missing -", file[0])

    sg.Popup("Renamed all files in:", Path(str(rename_data_folder)), title="Finished")














#
