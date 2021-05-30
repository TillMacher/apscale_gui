def v_unoise(path_to_outdirs, unoise_minsize, unoise_alpha, main_log_file):

    # main_log_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/log.xlsx"
    # path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST"
    # unoise_minsize = 8
    # unoise_alpha = 2
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
    pooled_files_fasta = Path(str(path_to_outdirs) + "/5_pre_processing/_data/pooled_files.fasta")
    ESV_folder = Path(str(path_to_outdirs) + "/6_denoising/ESV_a" + str(unoise_alpha))
    v_unoise_ESVs_fasta = Path(str(ESV_folder) + "/ESV_a" + str(unoise_alpha) + ".fasta")
    v_unoise_ESVs_fasta_no_chimera = Path(str(ESV_folder) + "/ESV_a" + str(unoise_alpha) + "_no_chimera.fasta")
    stderr_file = Path(str(path_to_outdirs) + "/6_denoising/_log/ESV_a" + str(unoise_alpha) + ".log")

    ## collect the number of raw reads
    log_df = pd.read_excel(main_log_file)

    dirName = Path(str(path_to_outdirs) + "/6_denoising/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    if not os.path.exists(ESV_folder):
        os.mkdir(Path(ESV_folder))

    log_folder = Path(str(path_to_outdirs) + "/6_denoising/_log")
    if not os.path.exists(log_folder):
        os.mkdir(Path(log_folder))

    if not pooled_files_fasta.exists():
        sg.PopupError("Warning: No files found! Please check your the merged reads!")

    else:
        f = open(stderr_file, "w")

        ## DENOISING
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting to denoise the pooled data.")
        subprocess.call(["vsearch", "--cluster_unoise", str(pooled_files_fasta), "--unoise_alpha", str(unoise_alpha), "--minsize", str(unoise_minsize), "--sizeout", "--centroids", str(v_unoise_ESVs_fasta)], stderr=f)

        n_ESVs = count_fasta_reads(v_unoise_ESVs_fasta)
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Data was denoised into", str(n_ESVs), "ESVs.")
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished to denoise the pooled data.")

        # CHIMERA REMOVAL
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting chimera removal.")
        subprocess.call(["vsearch", "--uchime3_denovo", str(v_unoise_ESVs_fasta), "--relabel", "ESV_", "--nonchimeras", str(v_unoise_ESVs_fasta_no_chimera)], stderr=f)
        n_ESVs = count_fasta_reads(v_unoise_ESVs_fasta_no_chimera)
        print(datetime.datetime.now().strftime("%H:%M:%S") + ":", str(n_ESVs), "ESVs remained after chimera removal.")
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished chimera removal.")
        f.close()

        # RE-MAPPING and READ TABLE GENERATION
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting to re-map ESVs to files.")

        ## collect sequences from the fasta file
        all_ESVs_dict = {}
        for record in SeqIO.parse(v_unoise_ESVs_fasta_no_chimera, "fasta"):
            all_ESVs_dict[str(record.id)] = str(record.seq)

        input_files = log_df["Sample ID"].values.tolist()

        read_table_df = pd.DataFrame(all_ESVs_dict.keys(), columns=["ID"])

        for sample_name in input_files:
            remapping_folder = Path(str(ESV_folder) + "/" + sample_name)
            stderr_file = Path(str(ESV_folder) + "/" + sample_name + "/" + sample_name + ".log")
            blastout_file = Path(str(ESV_folder) + "/" + sample_name + "/" + sample_name + "_blast.txt")
            otutabout_file = Path(str(ESV_folder) + "/" + sample_name + "/" + sample_name + "_otu.txt")
            sample_fasta = Path(str(path_to_outdirs) + "/5_pre_processing/_data/" + sample_name + "/" + sample_name + "_derep_singletons_nochimeras.fasta")

            if not os.path.exists(remapping_folder):
                os.mkdir(Path(remapping_folder))

            # MAP EACH FILE AGAINST ESVs
            f = open(stderr_file, "w")
            subprocess.call(["vsearch", "--search_exact", sample_fasta, "-db", str(v_unoise_ESVs_fasta_no_chimera), "--blast6out", str(blastout_file), "--output_no_hits", "--maxhits", "1", "--otutabout", str(otutabout_file)], stderr=f)
            f.close()

            try:
                ## open the mapping results
                otutabout_df = pd.read_csv(otutabout_file, sep='\t')
                ## store the number of reads for each ESV
                reads_list = []
                for ESV in all_ESVs_dict.keys():
                    ## test if the ESV exists
                    try:
                        ## collect the reads
                        reads = otutabout_df.loc[otutabout_df['#OTU ID'] == ESV].values.tolist()[0][1]
                    except:
                        ## if not store 0 reads
                        reads = 0
                    ## append the reads to the list
                    reads_list.append(reads)
                ## add a new column for each sample to the df, containing all read numbers
                read_table_df[sample_name] = reads_list
            except:
                ## if file does not exist, add 0 reads for all ESVs
                read_table_df[sample_name] = [0] * len(all_ESVs_dict.keys())

        ## sort the read table
        read_table_df["sort"] = [int(OTU.split("_")[1]) for OTU in read_table_df["ID"]]
        read_table_df = read_table_df.sort_values("sort")
        read_table_df = read_table_df.drop(columns=['sort'])
        read_table_df["Sequences"] = [all_ESVs_dict[ESV] for ESV in read_table_df["ID"].values.tolist()]

        ## write the read table
        read_table_xlsx = Path(str(ESV_folder) + "/ESV_a" + str(unoise_alpha) + ".xlsx")
        read_table_df.to_excel(read_table_xlsx, index=False)

        ## write to log file
        ESV_list = []
        for sample in log_df["Sample ID"]:
            ESV_list.append(len([n for n in read_table_df[sample].values.tolist() if n != 0]))
        id = "ESVs (a=" + str(unoise_alpha) + ")"
        log_df[id] = ESV_list
        log_df.to_excel(main_log_file, index=False)

        ## Finish script
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished to re-map ESVs to files.")
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Wrote ESV read table.")



                #
