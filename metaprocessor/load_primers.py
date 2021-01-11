def load_primers(primer_sheet_xlsx):

    import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import pandas as pd
    import webbrowser
    from Bio import SeqIO
    from Bio.Seq import Seq
    import plotly.graph_objects as go

    # collect the available primers from the primer sheet
    # store them in a fw primer and a rv primer dict

    df = pd.read_excel(Path(primer_sheet_xlsx))
    available_fw_primers_dict = {}
    available_rv_primers_dict = {}
    for primer in df.values.tolist():
        direction = primer[0]
        ID = primer[1]
        sequence = primer[2]

        if direction == 'forward':
            available_fw_primers_dict[ID] = sequence
        elif direction == 'reverse':
            available_rv_primers_dict[ID] = sequence
        else:
            print('Error: no direction given.')

    def slices(list, slice):
        for i in range(0, len(list), slice):
            yield list[i : i + slice]

    available_fw_primers_list = list(slices([sg.Radio(name, "fw", key=name, default=True, size=(15,1)) for name in sorted(list(available_fw_primers_dict.keys()))], 10))
    available_rv_primers_list = list(slices([sg.Radio(name, "rv", key=name, default=True, size=(15,1)) for name in sorted(list(available_rv_primers_dict.keys()))], 10))

    win2_active = True
    #window.Hide()
    layout2 = [[sg.Text("Available primers", size=(20,1))],
    [sg.Frame(layout = available_fw_primers_list, title = 'Choose your forward primer')],
    [sg.Frame(layout = available_rv_primers_list, title = 'Choose your reverse primer')],
    [sg.Text("")],
    [sg.Button('Trim primers'), sg.Text(""), sg.Button('Back')]]

    win2 = sg.Window('Available primers', layout2, keep_on_top=True)
    while True:

        event2, values2 = win2.Read()

        ##################################################
        # break event
        # exit the function

        if event2 is None or event2 == 'Back':
            win2.Close()
            win2_active = False
            #window.UnHide()
            return "Back"
            break

        ##################################################
        # filter event
        # start the actual filter script

        if event2 == 'Trim primers':

            for key, value in values2.items():
                if value == True:
                    if key in available_fw_primers_dict.keys():
                        fw_primer = available_fw_primers_dict[key]
                    if key in available_rv_primers_dict.keys():
                        rv_primer = available_rv_primers_dict[key]

            win2.Close()
            win2_active = False
            return [fw_primer, rv_primer]
            #window.UnHide()
            break
