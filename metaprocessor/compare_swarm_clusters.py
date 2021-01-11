def compare_swarm_clusters(path_to_outdirs, print_handle, window):

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    from Bio import SeqIO
    import pandas as pd

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting SWARM cluster comparison.")
    window.refresh()

    input_files = sorted(glob.glob(str(path_to_outdirs) + "/8_Read_tables/_data/*/*SWARM*.xlsx"))
    input_files_short = []

    for file in input_files:
        input_files_short.append(Path(file).stem + ".xlsx")

    ##################################################
    # start a second window to ask for the read tables to process

    win2_active = True
    #window.Hide()
    input_files_list = list(slices([sg.CB(name, default=False) for name in sorted(input_files_short)], 2))
    layout2 = [[sg.Text("OTU table filering", size=(20,1))],
    [sg.Frame(layout = input_files_list, title = 'Check read tables to filter')],
    [sg.Button('Compare clusters')],
    [sg.Button('Back')]]

    win2 = sg.Window('Taxon table filtering', layout2, keep_on_top=True)
    while True:
        event2, values2 = win2.Read()

        ##################################################
        # filter event
        # start the actual filter script
        files_to_process = []

        if event2 == 'Compare clusters':
            for file, key in zip(input_files, list(values2.values())):
                if key == True:
                    files_to_process.append(file)
            if files_to_process == []:
                sg.Popup("Error: Please choose a file", title="Error", keep_on_top=True)
            else:

                ##################################################
                # start the script

                n_OTUs_dict = {}

                for file in files_to_process:
                    read_table_df = pd.read_excel(file)
                    if type(int(Path(file).stem.split("_")[2])) == int:
                        n_OTUs_dict[int(Path(file).stem.split("_")[2])] = (len(read_table_df["OTUs"].values.tolist()))
                    else:
                        n_OTUs_dict[Path(file).stem] = (len(read_table_df["OTUs"].values.tolist()))

                x_values = list(n_OTUs_dict.keys())
                y_values = list(n_OTUs_dict.values())

                fig = go.Figure([go.Bar(x=x_values, y=y_values)])
                fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
                fig.update_layout(xaxis_title="swarm d", yaxis_title="# OTUs")
                output_html = str(path_to_outdirs) + "/8_Read_tables/_stats/swarm_cluster_comparison.html"
                fig.write_html(output_html)

                ## print standard closing message
                print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished SWARM cluster comparison.")
                window.refresh()
                sg.Popup("Processed files are found under:", Path(str(path_to_outdirs) + "/8_Read_tables/_stats/"), title="Finished", keep_on_top=True)

                win2.Close()
                win2_active = False
                #window.UnHide()
                break

        ##################################################
        # break event
        # exit the function

        if event2 is None or event2 == 'Back':
            win2.Close()
            win2_active = False
            #window.UnHide()
            break
