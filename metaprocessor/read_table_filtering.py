def read_table_filtering(filter_treshold, path_to_outdirs, print_handle, window):

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    from Bio import SeqIO
    import pandas as pd

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting read table filtering.")
    window.refresh()

    # slice function for lists to split up lists
    def slices(list, slice):
        for i in range(0, len(list), slice):
            yield list[i : i + slice]

    # SWARM_raw_read_table = Path("/home/till/Desktop/Projects/Projects_Development/dev_Meta_Tools/Projects/Default_project/8_Read_tables/_data/Swarm_cluster_d_1/SWARM_output_d_1_raw_read_table.xlsx")
    # path_to_outdirs = Path("/home/till/Desktop/Projects/Projects_Development/dev_Meta_Tools/Projects/Default_project")
    input_files = sorted(glob.glob(str(path_to_outdirs) + "/8_Read_tables/_data/*/*.xlsx"))
    input_files_short = []
    files_to_process = []
    filter_treshold_p = filter_treshold
    filter_treshold = float(filter_treshold) / 100

    for file in input_files:
        input_files_short.append(Path(file).stem + ".xlsx")

    ##################################################
    # start a second window to ask for the read tables to process

    win2_active = True
    #window.Hide()
    input_files_list = list(slices([sg.CB(name, default=True) for name in sorted(input_files_short)], 2))
    layout2 = [[sg.Text("OTU table filering", size=(20,1))],
    [sg.Frame(layout = input_files_list, title = 'Check read tables to filter')],
    [sg.Button('Filter')],
    [sg.Button('Back')]]

    win2 = sg.Window('Taxon table filtering', layout2, keep_on_top=False)
    while True:
        event2, values2 = win2.Read()

        ##################################################
        # filter event
        # start the actual filter script

        if event2 == 'Filter':
            for file, key in zip(input_files, list(values2.values())):
                if key == True:
                    files_to_process.append(file)
            if files_to_process == []:
                sg.Popup("Error: Please choose a file", title="Error")
            else:

                for file in files_to_process:

                    ## print standard closing message
                    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Filtering " + Path(file).stem + " with " + str(filter_treshold*100) + "%." )
                    window.refresh()

                    read_table_xlsx = Path(file)

                    read_table_df = pd.read_excel(read_table_xlsx)
                    read_table_df_new = read_table_df["ID"]
                    samples = read_table_df.columns.tolist()[1:-1]

                    for sample in samples:
                        sum_of_reads = sum(read_table_df[sample].values.tolist())
                        cutoff = round(sum_of_reads * float(filter_treshold))
                        reads_list = []
                        for n_reads in read_table_df[sample].values.tolist():
                            if n_reads <= cutoff:
                                reads_list.append(0)
                            else:
                                reads_list.append(n_reads)
                        sample_df = pd.DataFrame(reads_list, columns=[sample])
                        read_table_df_new = pd.concat([read_table_df_new, sample_df], axis=1)

                    read_table_df_new = pd.concat([read_table_df_new, read_table_df["Sequences"]], axis=1)

                    for row in read_table_df_new.values.tolist():
                        if sum(row[1:-1]) == 0:
                            read_table_df_new = read_table_df_new[read_table_df_new.ID != str(row[0])]

                    output_xlsx = str(Path(file).parent) + "/" + str(Path(file).stem) + "_" + str(filter_treshold_p) + "_p.xlsx"
                    read_table_df_new.to_excel(output_xlsx, index=False)

                    output_fasta = str(Path(file).parent) + "/" + str(Path(file).stem) + "_" + str(filter_treshold_p) + "_p.fasta"
                    f = open(output_fasta, "w")
                    for row in read_table_df_new[["ID", "Sequences"]].values.tolist():
                        f.write(">" + row[0] + "\n" + row[1] + "\n")
                    f.close()

                ## print standard closing message
                print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished read table filtering.")
                window.refresh()
                sg.Popup("Processed files are found under:", Path(str(path_to_outdirs) + "/8_Read_tables/_data/"), title="Finished")

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
