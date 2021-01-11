def create_demultiplexing_sheet(path_to_outdirs):

    import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import pandas as pd
    import webbrowser
    from Bio import SeqIO
    from Bio.Seq import Seq
    import plotly.graph_objects as go

    df_1 = pd.DataFrame([["index1", "AGTC"]], columns=["index", "sequences"])
    df_2 = pd.DataFrame([["sample1", "index1", "index1"]], columns=["name", "i5", "i7"])
    writer = pd.ExcelWriter(str(path_to_outdirs) + "/2_Demultiplexing/demultiplexing_sheet_empty.xlsx", engine='xlsxwriter')
    df_1.to_excel(writer, sheet_name='index', index=False)
    df_2.to_excel(writer, sheet_name='combinations', index=False)
    writer.save()

    sg.Popup("The new demultiplexing sheet is found under: /2_Demultiplexing/demultiplexing_sheet_empty.xlsx/", title="Finished")
