def check_installations():

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

    missing = False
    try:
        subprocess.call(["swarm", "--version"])
    except:
        missing = True
        sg.Popup("Oops! Swarm seems to be missing.")

    print("\n##############################\n")

    try:
        subprocess.call(["vsearch", "--version"])
    except:
        missing = True
        sg.Popup("Oops! Vsearch seems to be missing.")

    print("\n##############################\n")

    try:
        subprocess.call(["cutadapt", "--version"])
    except:
        missing = True
        sg.Popup("Oops! Cutadapt seems to be missing.")

    print("\n##############################\n")

    try:
        subprocess.call(["fastqc", "--version"])
    except:
        missing = True
        sg.Popup("Oops! FastQC seems to be missing.")

    if missing == False:
        sg.Popup("Awesome! Everything is installed and you are ready to go!")
