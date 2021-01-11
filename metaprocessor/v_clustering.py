def v_clustering(clustering_folder, cluster_threshold, path_to_outdirs, representative, print_handle, window):

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    from Bio import SeqIO
    import pandas as pd

    ## define variables
    pooled_files_fasta_derep_nochimeras = str(path_to_outdirs) + "/7_Clustering/_data/_clusters/pooled_files_derep_nochimeras.fasta"
    pooled_files_derep_folder = str(path_to_outdirs) + "/7_Clustering/_data/_clusters"
    v_cluster_OTUs_fasta = Path(pooled_files_derep_folder + "/v_cluster/v_cluster_OTUs_" + cluster_threshold.replace(".", "") + ".fasta")
    v_cluster_OTUs_relabeled_fasta = Path(pooled_files_derep_folder + "/v_cluster/v_cluster_OTUs_" + cluster_threshold.replace(".", "") + "_relabeled.fasta")
    stats_file_OTUs = str(path_to_outdirs) + "/7_Clustering/_stats/v_cluster_" + cluster_threshold.replace(".", "") + "_OTUs.html"

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting " + cluster_threshold + " clustering")
    window.refresh()

    dirName = Path(pooled_files_derep_folder)
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    dirName = Path(pooled_files_derep_folder + "/v_cluster/")
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    dirName = Path(str(path_to_outdirs) + "/8_Read_tables/_data" + "/v_cluster")
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    if representative == "Consensus":
        representative = "--consout"
    else:
        representative = "--centroids"

    subprocess.call(["vsearch", "--cluster_size", pooled_files_fasta_derep_nochimeras, "--id", cluster_threshold, "--sizein", "--relabel", "OTU_", representative, v_cluster_OTUs_fasta])

    OTU = 1
    x_values =[]
    y_values = []

    outfile =  open(v_cluster_OTUs_relabeled_fasta, "w")
    for record in SeqIO.parse (v_cluster_OTUs_fasta, "fasta"):
        x_values.append(">OTU_" + str(OTU))
        y_values.append(record.id.split(":")[-1].split("=")[-1].replace(";", ""))
        outfile.write(">OTU_" + str(OTU) + "\n")
        outfile.write(str(record.seq) + "\n")
        OTU += 1

    ##################################################
    # draw the stats plot
    fig = go.Figure(data=[
        go.Bar(name='OTUs', x=x_values, y=y_values),
        ])
    fig.update_layout(barmode='stack')
    fig.update_xaxes(tickangle=90, tickfont=dict(size=10))
    fig.write_html(str(stats_file_OTUs))

    ##################################################
    # re-mapping

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Remapping OTUs to files")
    window.refresh()

    input_files = sorted(glob.glob(str(path_to_outdirs) + "/7_Clustering/_data/*/*_derep_singletons_nochimeras.fasta"))

    all_OTUs_list = []
    all_sequences_list = []
    for record in SeqIO.parse(v_cluster_OTUs_relabeled_fasta, "fasta"):
        all_OTUs_list.append(record.id)
        all_sequences_list.append(str(record.seq))

    currentDirectory = os.getcwd()

    for file in input_files:
        hits = Path(str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace("_derep_singletons_nochimeras", "/") + Path(file).stem.replace("_derep_singletons_nochimeras", "_v_clustering_hits") + ".fasta")
        OTU_stats = Path(str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace("_derep_singletons_nochimeras", "/") + Path(file).stem.replace("_derep_singletons_nochimeras", "_v_clustering_stats") + ".txt")
        stderr_file = Path(str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace("_derep_singletons_nochimeras", "/") + Path(file).stem.replace("_derep_singletons_nochimeras", "_v_clustering") + ".log")

        ### Correct? ID = 97% is not a problem?
        # MAP EACH FILE AGAINST OTUs
        f = open(stderr_file, "w")
        subprocess.call(["vsearch", "--usearch_global", file, "-db", v_cluster_OTUs_relabeled_fasta, "--id", "0.97", "--blast6out", str(hits), "--maxhits", "1", "--otutabout", str(OTU_stats), "--output_no_hits"], stderr=f)
        f.close()

    ##################################################
    # create raw OTU table

    read_table_dict = {}

    for file in input_files:
        filename = Path(file).stem.replace("_derep_singletons_nochimeras", "")
        OTU_stats = Path(str(path_to_outdirs) + "/7_Clustering/_data/" + Path(file).stem.replace("_derep_singletons_nochimeras", "/") + Path(file).stem.replace("_derep_singletons_nochimeras", "_v_clustering_stats") + ".txt")
        v_clustering_raw_read_table = Path(str(path_to_outdirs) + "/8_Read_tables/_data/v_cluster/v_clustering_" + cluster_threshold.replace(".", "") + ".xlsx")

        if os.stat(OTU_stats).st_size != 0:
            OTU_stats_df = pd.read_csv(OTU_stats, delimiter="\t").transpose()
            OTU_stats_df = OTU_stats_df.rename(columns=OTU_stats_df.iloc[0])
            OTU_stats_df = OTU_stats_df.drop(OTU_stats_df.index[0])
            OTU_stats_df_headers = OTU_stats_df.columns.tolist()
            read_table_list = []

            for OTU in all_OTUs_list:
                if OTU in OTU_stats_df_headers:
                    read_table_list.append(int(''.join(map(str, OTU_stats_df[OTU].values.tolist()))))
                else:
                    read_table_list.append(0)

            read_table_dict[filename] = read_table_list

            # # lastly append the number of unassigned reads!
            # # those are collected from the merge log file!!
            # merge_log = Path(str(path_to_outdirs) + "/3_Paired-end_merging/_log/" + filename + ".log")
            # f = open(merge_log, "r")
            # sum_of_reads = f.readline()
            # f.close()
            # sum_of_reads = int(sum_of_reads.split(" ")[-3])
            # sum_of_assigned_reads = sum(read_table_list)
            # sum_of_unassinged_reads = sum_of_reads - sum_of_assigned_reads
            # read_table_list.append(sum_of_unassinged_reads)

        else:
            print("Discarded:", filename, "- no matches after remapping")

    # now collect all previously discarded files
    # add them to the table
    discarded_files = sorted(glob.glob(str(path_to_outdirs) + "/*/_data_discarded/*"))
    # add 0 only
    for file in discarded_files:
        filename = Path(file).name.replace(".fastq", "")
        read_table_list = []
        for OTU in all_OTUs_list:
            read_table_list.append(0)

        read_table_dict[filename] = read_table_list

        # lastly append the number of unassigned reads!
        # those are collected from the merge log file!!
        # merge_log = Path(str(path_to_outdirs) + "/3_Paired-end_merging/_log/" + filename + ".log")
        # f = open(merge_log, "r")
        # sum_of_reads = f.readline()
        # f.close()
        # sum_of_reads = int(sum_of_reads.split(" ")[-3])
        # sum_of_assigned_reads = sum(read_table_list)
        # sum_of_unassinged_reads = sum_of_reads - sum_of_assigned_reads
        # read_table_list.append(sum_of_unassinged_reads)
        # read_table_dict[filename] = read_table_list

    # create the read table and finish the script
    read_table_df = pd.DataFrame.from_dict(read_table_dict, orient='columns')
    read_table_df.insert(0, "ID", all_OTUs_list, True)
    read_table_df['Sequences'] = all_sequences_list
    read_table_df.to_excel(v_clustering_raw_read_table, index=False)

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Writing read table.")
    window.refresh()

    v_clustering_raw_read_table_fasta = Path(str(path_to_outdirs) + "/8_Read_tables/_data/v_cluster/v_clustering_" + cluster_threshold.replace(".", "") + ".fasta")
    f = open(v_clustering_raw_read_table_fasta, "w")
    for row in read_table_df[["ID", "Sequences"]].values.tolist():
        f.write(">" + row[0] + "\n" + row[1] + "\n")
    f.close()

    ## print standard starting message
    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished " + cluster_threshold + " clustering.")
    window.refresh()
    sg.Popup("Read table is found under:", Path(str(path_to_outdirs) + "/8_Read_tables/_data/v_cluster"), title="Finished", keep_on_top=True)
