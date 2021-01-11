def md5sum(raw_data_folder):

    import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import pandas as pd
    import webbrowser
    from Bio import SeqIO
    from Bio.Seq import Seq
    import plotly.graph_objects as go

    raw_data_folder = Path(raw_data_folder)
    input_files = sorted(glob.glob(str(raw_data_folder) + "/*.fastq.gz"))

    md5sum_list = []

    for file in input_files:
        with open(file, 'rb') as fh:
            m = hashlib.md5()
            while True:
                data = fh.read(8192)
                if not data:
                    break
                m.update(data)
        md5sum_list.append(Path(file).stem + " = " + m.hexdigest() + "\n")

    sg.PopupScrolled(' '.join(md5sum_list), title = "md5sums", size=(80, 30))
