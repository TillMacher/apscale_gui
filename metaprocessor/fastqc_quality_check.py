def fastqc_quality_check(raw_data_folder, path_to_outdirs, num_cores_to_use):

    import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import pandas as pd
    import webbrowser
    from Bio import SeqIO
    from Bio.Seq import Seq
    import plotly.graph_objects as go

    sg.PopupAutoClose("Running FastQC - this may take a while", auto_close_duration=5, title="Notification")

    raw_data_folder = Path(raw_data_folder)
    input_files = sorted(glob.glob(str(raw_data_folder) + "/*.fastq.gz"))
    n_files = len(input_files)

    def fastqc_command(i, file, path_to_outdirs):
        output_dir = Path(str(path_to_outdirs) + "/1_Quality_check/_data/")
        subprocess.call(["fastqc", "--quiet", "-t", str(num_cores_to_use), "--outdir", output_dir, file])
        print("(" + str(i+1) + "/" + str(n_files) + ")   ", Path(file).stem)

    results = Parallel(n_jobs=num_cores_to_use)(delayed(fastqc_command)(i, file, path_to_outdirs) for i, file in enumerate(input_files))

    sg.Popup("Finished FastQC", title="Finished")
