def v_clustering_raw(path_to_outdirs, clustering_threshold, main_log_file):

    # main_log_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/log.xlsx"
    # path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST"
    # clustering_threshold = 0.97
    # representative = "Centroid"

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    from Bio import SeqIO
    import pandas as pd

    ## define functions
    def count_fasta_reads(file):
        def blocks(files, size=65536):
            while True:
                b = files.read(size)
                if not b: break
                yield b
        with open(file, "r",encoding="utf-8",errors='ignore') as f:
            return sum(bl.count(">") for bl in blocks(f))

    ## define variables
    v_unoise_ESVs_fasta_no_chimera = Path(str(path_to_outdirs) + "/5_pre_processing/_data/pooled_files.fasta")
    OTU_folder = Path(str(path_to_outdirs) + "/7_clustering/OTUs_a" + str(clustering_threshold))
    v_OTUs_fasta = Path(str(OTU_folder) + "/OTU_a" + str(clustering_threshold) + ".fasta")
    stderr_file = Path(str(path_to_outdirs) + "/7_clustering/_log/OTU_a" + str(clustering_threshold) + ".log")

    ## collect the number of raw reads
    log_df = pd.read_excel(main_log_file)

    dirName = Path(str(path_to_outdirs) + "/7_clustering/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    if not os.path.exists(OTU_folder):
        os.mkdir(Path(OTU_folder))

    log_folder = Path(str(path_to_outdirs) + "/7_clustering/_log")
    if not os.path.exists(log_folder):
        os.mkdir(Path(log_folder))

    if not v_unoise_ESVs_fasta_no_chimera.exists():
        sg.PopupError("Warning: No files found! Please check your the merged reads!")

    else:
        f = open(stderr_file, "w")

        ## CLUSTERING
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting to cluster the pooled data.")

        subprocess.call(["vsearch", "--cluster_size", str(v_unoise_ESVs_fasta_no_chimera), "--id", str(clustering_threshold), "--sizein", "--relabel", "OTU_", "--centroids", str(v_OTUs_fasta)], stderr=f)

        n_OTUs = count_fasta_reads(v_OTUs_fasta)

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Data was cluster into", str(n_OTUs), "OTUs.")
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished to cluster the pooled data.")

        f.close()

        # RE-MAPPING and READ TABLE GENERATION
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting to re-map OTUs to files.")

        ## collect sequences from the fasta file
        all_OTUs_dict = {}
        for record in SeqIO.parse(v_OTUs_fasta, "fasta"):
            all_OTUs_dict[str(record.id)] = str(record.seq)

        input_files = log_df["Sample ID"].values.tolist()

        read_table_df = pd.DataFrame(all_OTUs_dict.keys(), columns=["ID"])

        for sample_name in input_files:
            remapping_folder = Path(str(OTU_folder) + "/" + sample_name)
            stderr_file = Path(str(OTU_folder) + "/" + sample_name + "/" + sample_name + ".log")
            blastout_file = Path(str(OTU_folder) + "/" + sample_name + "/" + sample_name + "_blast.txt")
            otutabout_file = Path(str(OTU_folder) + "/" + sample_name + "/" + sample_name + "_otu.txt")
            sample_fasta = Path(str(path_to_outdirs) + "/5_pre_processing/_data/" + sample_name + "/" + sample_name + "_derep_singletons_nochimeras.fasta")

            if not os.path.exists(remapping_folder):
                os.mkdir(Path(remapping_folder))

            # MAP EACH FILE AGAINST OTUs
            f = open(stderr_file, "w")
            subprocess.call(["vsearch", "--usearch_global", sample_fasta, "-db", str(v_OTUs_fasta), "--id", str(clustering_threshold), "--blast6out", str(blastout_file), "--output_no_hits", "--maxhits", "1", "--otutabout", str(otutabout_file)], stderr=f)
            f.close()

            try:
                ## open the mapping results
                otutabout_df = pd.read_csv(otutabout_file, sep='\t')
                ## store the number of reads for each OTU
                reads_list = []
                for OTU in all_OTUs_dict.keys():
                    ## test if the OTU exists
                    try:
                        ## collect the reads
                        reads = otutabout_df.loc[otutabout_df['#OTU ID'] == OTU].values.tolist()[0][1]
                    except:
                        ## if not store 0 reads
                        reads = 0
                    ## append the reads to the list
                    reads_list.append(reads)
                ## add a new column for each sample to the df, containing all read numbers
                read_table_df[sample_name] = reads_list
            except:
                ## if file does not exist, add 0 reads for all OTUs
                read_table_df[sample_name] = [0] * len(all_OTUs_dict.keys())

        ## sort the read table
        read_table_df["sort"] = [int(OTU.split("_")[1]) for OTU in read_table_df["ID"]]
        read_table_df = read_table_df.sort_values("sort")
        read_table_df = read_table_df.drop(columns=['sort'])
        read_table_df["Sequences"] = [all_OTUs_dict[OTU] for OTU in read_table_df["ID"].values.tolist()]

        ## write the read table
        read_table_xlsx = Path(str(OTU_folder) + "/OTU_a" + str(clustering_threshold) + ".xlsx")
        read_table_df.to_excel(read_table_xlsx, index=False)

        ## write to log file
        OTU_list = []
        for sample in log_df["Sample ID"]:
            OTU_list.append(len([n for n in read_table_df[sample].values.tolist() if n != 0]))
        id = "OTUs (id=" + str(clustering_threshold) + ")"
        log_df[id] = OTU_list
        log_df.to_excel(main_log_file, index=False)

        ## Finish script
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished to re-map OTUs to files.")
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Wrote OTU read table.")

def v_clustering_denoised(path_to_outdirs, clustering_threshold, unoise_alpha, main_log_file):

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
    from Bio import SeqIO
    import pandas as pd

    ## define functions
    def count_fasta_reads(file):
        def blocks(files, size=65536):
            while True:
                b = files.read(size)
                if not b: break
                yield b
        with open(file, "r",encoding="utf-8",errors='ignore') as f:
            return sum(bl.count(">") for bl in blocks(f))

    ## define variables
    ESV_ID = "ESV_a" + str(unoise_alpha)
    v_unoise_ESVs_fasta_no_chimera = Path(str(path_to_outdirs) + "/6_denoising/" + ESV_ID + "/" + ESV_ID + "_no_chimera.fasta")
    OTU_folder = Path(str(path_to_outdirs) + "/7_clustering/OTU_id" + str(clustering_threshold) + "_a" + str(unoise_alpha))
    v_OTUs_fasta = Path(str(OTU_folder) + "/OTU_id" + str(clustering_threshold)  + "_a" + str(unoise_alpha) + ".fasta")
    stderr_file = Path(str(path_to_outdirs) + "/7_clustering/_log/OTU_id" + str(clustering_threshold) + ".log")

    ## collect the number of raw reads
    log_df = pd.read_excel(main_log_file)

    dirName = Path(str(path_to_outdirs) + "/7_clustering/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    if not os.path.exists(OTU_folder):
        os.mkdir(Path(OTU_folder))

    log_folder = Path(str(path_to_outdirs) + "/7_clustering/_log")
    if not os.path.exists(log_folder):
        os.mkdir(Path(log_folder))

    if not v_unoise_ESVs_fasta_no_chimera.exists():
        sg.PopupError("Warning: No files found! Please check your the merged reads!")

    else:
        f = open(stderr_file, "w")

        ## CLUSTERING
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting to cluster the pooled data.")

        subprocess.call(["vsearch", "--cluster_size", str(v_unoise_ESVs_fasta_no_chimera), "--id", str(clustering_threshold), "--sizein", "--relabel", "OTU_", "--centroids", str(v_OTUs_fasta)], stderr=f)

        n_OTUs = count_fasta_reads(v_OTUs_fasta)

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Data was cluster into", str(n_OTUs), "OTUs.")
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished to cluster the pooled data.")

        f.close()

        # RE-MAPPING and READ TABLE GENERATION
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting to re-map OTUs to files.")

        ## collect sequences from the fasta file
        all_OTUs_dict = {}
        for record in SeqIO.parse(v_OTUs_fasta, "fasta"):
            all_OTUs_dict[str(record.id)] = str(record.seq)

        input_files = log_df["Sample ID"].values.tolist()

        read_table_df = pd.DataFrame(all_OTUs_dict.keys(), columns=["ID"])

        for sample_name in input_files:
            remapping_folder = Path(str(OTU_folder) + "/" + sample_name)
            stderr_file = Path(str(OTU_folder) + "/" + sample_name + "/" + sample_name + ".log")
            blastout_file = Path(str(OTU_folder) + "/" + sample_name + "/" + sample_name + "_blast.txt")
            otutabout_file = Path(str(OTU_folder) + "/" + sample_name + "/" + sample_name + "_otu.txt")
            sample_fasta = Path(str(path_to_outdirs) + "/5_pre_processing/_data/" + sample_name + "/" + sample_name + "_derep_singletons_nochimeras.fasta")

            if not os.path.exists(remapping_folder):
                os.mkdir(Path(remapping_folder))

            # MAP EACH FILE AGAINST OTUs
            f = open(stderr_file, "w")
            subprocess.call(["vsearch", "--usearch_global", sample_fasta, "-db", str(v_OTUs_fasta), "--id", str(clustering_threshold), "--blast6out", str(blastout_file), "--output_no_hits", "--maxhits", "1", "--otutabout", str(otutabout_file)], stderr=f)
            f.close()

            try:
                ## open the mapping results
                otutabout_df = pd.read_csv(otutabout_file, sep='\t')
                ## store the number of reads for each OTU
                reads_list = []
                for OTU in all_OTUs_dict.keys():
                    ## test if the OTU exists
                    try:
                        ## collect the reads
                        reads = otutabout_df.loc[otutabout_df['#OTU ID'] == OTU].values.tolist()[0][1]
                    except:
                        ## if not store 0 reads
                        reads = 0
                    ## append the reads to the list
                    reads_list.append(reads)
                ## add a new column for each sample to the df, containing all read numbers
                read_table_df[sample_name] = reads_list
            except:
                ## if file does not exist, add 0 reads for all OTUs
                read_table_df[sample_name] = [0] * len(all_OTUs_dict.keys())

        ## sort the read table
        read_table_df["sort"] = [int(OTU.split("_")[1]) for OTU in read_table_df["ID"]]
        read_table_df = read_table_df.sort_values("sort")
        read_table_df = read_table_df.drop(columns=['sort'])
        read_table_df["Sequences"] = [all_OTUs_dict[OTU] for OTU in read_table_df["ID"].values.tolist()]

        ## write the read table
        read_table_xlsx = Path(str(OTU_folder) + "/OTU_id" + str(clustering_threshold)  + "_a" + str(unoise_alpha) + ".xlsx")
        read_table_df.to_excel(read_table_xlsx, index=False)

        ## write to log file
        OTU_list = []
        for sample in log_df["Sample ID"]:
            OTU_list.append(len([n for n in read_table_df[sample].values.tolist() if n != 0]))
        id = "OTUs (id=" + str(clustering_threshold) + ";a=" + str(unoise_alpha) + ")"
        log_df[id] = OTU_list
        log_df.to_excel(main_log_file, index=False)

        ## Finish script
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished to re-map OTUs to files.")
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Wrote OTU read table.")
