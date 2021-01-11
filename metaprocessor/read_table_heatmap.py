def read_table_heatmap(path_to_outdirs, print_handle, window):

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    from Bio import SeqIO
    import pandas as pd

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Creating read table heatmap.")
    window.refresh()

    # slice function for lists to split up lists
    def slices(list, slice):
        for i in range(0, len(list), slice):
            yield list[i : i + slice]

    input_files = sorted(glob.glob(str(path_to_outdirs) + "/8_Read_tables/_data/*/*.xlsx"))
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
    [sg.Button('Create Heatmap')],
    [sg.Button('Back')]]

    win2 = sg.Window('Taxon table filtering', layout2, keep_on_top=True)
    while True:
        event2, values2 = win2.Read()
        print(values2.values())

        ##################################################
        # filter event
        # start the actual filter script
        files_to_process = []

        if event2 == 'Create Heatmap':
            for file, key in zip(input_files, list(values2.values())):
                if key == True:
                    files_to_process.append(file)
            if files_to_process == []:
                sg.Popup("Error: Please choose a file", title="Error")
            else:

                ##################################################
                # start the script

                for file in files_to_process:
                    print(file)
                    input_df = pd.read_excel(file)

                    x_values = input_df.columns.tolist()[1:-1]
                    z_values = input_df[x_values].values.tolist()
                    y_values = input_df["ID"].values.tolist()

                    fig = go.Figure(data=go.Heatmap(
                                       z=z_values,
                                       x=x_values,
                                       y=y_values,
                                       hoverongaps = False, colorscale="YlGnBu"))

                    output_html = str(path_to_outdirs) + "/8_Read_tables/_stats/" + str(Path(file).stem) + "_heatmap.html"
                    fig.write_html(output_html)

                ## print standard closing message
                print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished read table heatmap(s).")
                window.refresh()
                sg.Popup("Processed files are found under:", Path(str(path_to_outdirs) + "/8_Read_tables/_stats/"), title="Finished")

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
