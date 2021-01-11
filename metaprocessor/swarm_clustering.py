def swarm_clustering(clustering_folder, d_min, d_max, path_to_outdirs, num_cores_to_use):

    import glob, sys, subprocess, os, hashlib, gzip, multiprocessing, shutil, re
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import pandas as pd
    import webbrowser
    from Bio import SeqIO
    from Bio.Seq import Seq
    import plotly.graph_objects as go

    print("Swarm: Clustering")

    for d_val in range(int(d_min), int(d_max) + 1):

        pooled_files_fasta_derep_nochimeras = str(path_to_outdirs) + "/7_Clustering/_data/_clusters/pooled_files_derep_nochimeras.fasta"
        pooled_files_derep_folder = str(path_to_outdirs) + "/7_Clustering/_data/_clusters"
        SWARM_OTUs_fasta = Path(pooled_files_derep_folder + "/Swarm_cluster_d_" + str(d_val) + "/SWARM_OTUs.fasta")
        SWARM_OTUs_fasta_relabel = Path(pooled_files_derep_folder + "/Swarm_cluster_d_" + str(d_val) + "/SWARM_OTUs_relabel_d_" + str(d_val) + ".fasta")
        SWARM_output_txt = Path(pooled_files_derep_folder + "/Swarm_cluster_d_" + str(d_val) + "/SWARM_output_d_" + str(d_val) + ".txt")
        swarm_stats_file = Path(pooled_files_derep_folder + "/Swarm_cluster_d_" + str(d_val) + "/SWARM_output_d_" + str(d_val) + ".stats")
        stats_file_OTUs = str(path_to_outdirs) + "/7_Clustering/_stats/" + "d_" + str(d_val) + "_OTUs.html"
        input_files = glob.glob(str(path_to_outdirs) + "/7_Clustering/_data/*/*_derep_singletons_nochimeras.fasta")

        dirName = Path(pooled_files_derep_folder + "/Swarm_cluster_d_" + str(d_val))
        if not os.path.exists(dirName):
            os.mkdir(dirName)

        dirName = Path(str(path_to_outdirs) + "/8_Read_tables/_data" + "/Swarm_cluster_d_" + str(d_val))
        if not os.path.exists(dirName):
            os.mkdir(dirName)

        dirName = Path(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "/OTU_fastas")
        if not os.path.exists(dirName):
            os.mkdir(dirName)

        dirName = Path(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "/OTU_mapping")
        if not os.path.exists(dirName):
            os.mkdir(dirName)

        # CLUSTERING     subprocess.call(["swarm", "-z", "-t", cores, "--statistics-file", swarm_stats_file, "-o", pooled_files_OTUs, pooled_files_derep])
        subprocess.call(["swarm", "-z", "-d", str(d_val), "-t", str(num_cores_to_use), "-w", SWARM_OTUs_fasta, "-o", SWARM_output_txt, pooled_files_fasta_derep_nochimeras])

        ##################################################
        # relabel the SWARM_OTUs_fasta

        OTU = 1
        x_values =[]
        y_values = []

        outfile =  open(SWARM_OTUs_fasta_relabel, "w")
        for record in SeqIO.parse (SWARM_OTUs_fasta, "fasta"):
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

        # first create a dictionary that contains all OTUs and their respective underlying hits
        OTU_hits_dict = {}
        for i, hit in enumerate(open(SWARM_output_txt, "r")):
            OTU = "OTU_" + str(i + 1)
            # split the OTU swarm into a list
            # that list contains the [ID] and the [size]
            # get the ID by collecting every second element from 0
            # get the size by collecting every second element starting from 1
            hit_list = re.split(r'[ ;]', hit)
            hit_list_ID = hit_list[0::2]
            OTU_hits_dict[OTU] = hit_list_ID

        # then create a fasta file for each OTU
        def create_OTU_fastas(OTU, hits):

            ## these files can be used to adress the intraspecific distance within OTUs

            OTU_list = []
            # open a new fasta file for each OTU
            output_fasta = Path(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "/OTU_fastas/" + OTU + ".fasta")
            f = open(output_fasta, "w")
            # iterate through the pooled_files fasta (the input for SWARM)
            input_seq_iterator = SeqIO.parse(pooled_files_fasta_derep_nochimeras, "fasta")
            for record in input_seq_iterator:
                # exctract the record (and drop the size)
                record_id = (str(record.id).split(";")[0])
                # collect the record if found in hits
                if record_id in hits:
                    # write to the OTU fasta file
                    f.write(">" + OTU + ":" + str(record.id) + "\n")
                    f.write(str(record.seq) + "\n")
            f.close()

        # execute the create_OTU_fastas command
        results = Parallel(n_jobs=num_cores_to_use)(delayed(create_OTU_fastas)(OTU, hits) for OTU, hits in OTU_hits_dict.items())

        # create a main_OTU file by concatenating all previously created files
        main_OTU_fasta = Path(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "/OTU_fastas/main_OTUs.fasta")
        with open(main_OTU_fasta, 'w') as outfile:
            for file in glob.glob(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "/OTU_fastas/*.fasta"):
                if file != main_OTU_fasta:
                    with open(file, 'r') as readfile:
                        outfile.write(readfile.read())

        def vsearch_mapping(file):

            main_dict = {}
            abundance_dict = {}

            # get filename
            filename = str(Path(file).stem).replace("_derep_singletons_nochimeras", "")

            # define output files
            # overwrite each file for every OTU to reduce file abundance (and disk space)
            hits = Path(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "/OTU_mapping/" + filename + "_hits.txt")
            OTU_stats = Path(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "//OTU_mapping/" + filename + "_stats.txt")
            stderr_file = Path(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "/OTU_mapping/" + filename + "_log.txt")

            # map each OTU fasta
            f = open(stderr_file, "w")
            subprocess.call(["vsearch", "--search_exact", file, "-db", str(main_OTU_fasta),"--blast6out", str(hits), "--maxhits", "1", "--otutabout", str(OTU_stats), "--output_no_hits"], stderr=f)
            f.close()

        # execute the vsearch_mapping command
        results = Parallel(n_jobs=num_cores_to_use)(delayed(vsearch_mapping)(file) for file in input_files)

        all_OTUs_list = []
        all_sequences_list = []
        for record in SeqIO.parse(SWARM_OTUs_fasta_relabel, "fasta"):
            all_OTUs_list.append(record.id)
            all_sequences_list.append(str(record.seq))

        read_table_dict = {}

        for file in input_files:
            filename = str(Path(file).stem).replace("_derep_singletons_nochimeras", "")
            OTU_stats = Path(str(path_to_outdirs) + "/7_Clustering/_data/_clusters/Swarm_cluster_d_" + str(d_val) + "/OTU_mapping/" + filename + "_stats.txt")

            try:
                OTU_stats_df = pd.read_csv(OTU_stats, delimiter="\t").transpose()
                OTU_stats_df = OTU_stats_df.rename(columns=OTU_stats_df.iloc[0])
                OTU_stats_df = OTU_stats_df.drop(OTU_stats_df.index[0])
                OTU_stats_df_headers = OTU_stats_df.columns.tolist()
                renamed_headers_list = []
                for header in OTU_stats_df_headers:
                    renamed_headers_list.append(header.split(":")[0])
                OTU_stats_df.columns = renamed_headers_list

                read_table_list = []

                for OTU in all_OTUs_list:
                    if OTU in renamed_headers_list:
                        read_numbers = OTU_stats_df[OTU].values.tolist()
                        test = type(read_numbers[0])
                        # check if there are more than one hit
                        if test == list:
                            # flatten the list
                            flattened_list = [val for sublist in read_numbers for val in sublist]
                            # calculate the sum of all hits, which results in the number of reads for this OTU
                            OTU_abundance = sum(flattened_list)
                            # add the result to the output list
                            read_table_list.append(OTU_abundance)
                        else:
                            # if there is only one hit (this is not a singleton per se, just for this sample. There are more reads in another file):
                            # add the result to the output list
                            OTU_abundance = int(''.join(map(str, read_numbers)))
                            read_table_list.append(OTU_abundance)
                    else:
                        read_table_list.append(0)

                read_table_dict[filename] = read_table_list

            except:
                print("Discarded:", filename, "- no matches after remapping")

        read_table_df = pd.DataFrame.from_dict(read_table_dict, orient='columns')
        read_table_df.insert(0, "ID", all_OTUs_list, True)
        read_table_df['Sequences'] = all_sequences_list
        SWARM_raw_read_table = Path(str(path_to_outdirs) + "/8_Read_tables/_data" + "/Swarm_cluster_d_" + str(d_val) + "/SWARM_d_" + str(d_val) + "_raw_read_table.xlsx")
        read_table_df.to_excel(SWARM_raw_read_table, index=False)

        SWARM_raw_read_table_fasta = Path(str(path_to_outdirs) + "/8_Read_tables/_data" + "/Swarm_cluster_d_" + str(d_val) + "/SWARM_d_" + str(d_val) + "_raw_read_table.fasta")
        f = open(SWARM_raw_read_table_fasta, "w")
        for row in read_table_df[["ID", "Sequences"]].values.tolist():
            f.write(">" + row[0] + "\n" + row[1] + "\n")
        f.close()
