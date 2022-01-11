def create_rename_sheet(rename_data_folder, path_to_outdirs):

    try:
        import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
        import PySimpleGUI as sg
        from pathlib import Path
        from joblib import Parallel, delayed
        import pandas as pd
        import webbrowser
        from Bio import SeqIO
        from Bio.Seq import Seq
        import plotly.graph_objects as go
    except:
        print("\nError: You need to install the required packages first!\n")
        exit()

    suffix_list = [".1", ".2"]
    file_ending = ".fastq.gz"
    rename_data_folder = "/Users/tillmacher/Desktop/MP_Projects/Projects/Robot_LFU_data/0_raw_data/_data"
    sample_renaming_sheet = "/Users/tillmacher/Desktop/MP_Projects/Projects/Robot_LFU_data/rename_sheet.xlsx"
    path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/Robot_LFU_data/"

    rename_data_folder = Path(rename_data_folder)
    input_files = sorted(glob.glob(str(rename_data_folder) + "/*.fastq.gz"))

    files_list = list(set([Path(file).name.replace(suffix_list[0] + file_ending, "").replace(suffix_list[1] + file_ending, "") for file in input_files]))
    rename_sheet_df = pd.DataFrame(files_list, columns=["Original name"])
    rename_sheet_df.insert(1, "New name", [""] * int(len(input_files)/2))
    rename_sheet_df.to_excel(Path(str(path_to_outdirs) + "/rename_sheet.xlsx"), index=False)

    sg.Popup("Please add the new names to the following rename sheet:", Path(str(path_to_outdirs) + "/rename_sheet.xlsx"), title="Finished")
