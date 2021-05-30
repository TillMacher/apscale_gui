def v_derep_singletons_uchime(path_to_outdirs, num_cores_to_use, main_log_file):

    # main_log_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/log.xlsx"
    # path_to_outdirs = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST"
    # num_cores_to_use = 4

    import datetime, glob, subprocess, gzip, os
    import PySimpleGUI as sg
    from pathlib import Path
    from joblib import Parallel, delayed
    import plotly.graph_objects as go
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
    input_folder = Path(str(path_to_outdirs) + "/4_read_filtering/_data/")
    input_files = sorted(glob.glob(str(input_folder) + "/*.fasta"))
    n_files = len(input_files)

    ## collect the number of raw reads
    log_df = pd.read_excel(main_log_file)

    dirName = Path(str(path_to_outdirs) + "/5_pre_processing/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    data_folder = Path(str(path_to_outdirs) + "/5_pre_processing/_data")
    if not os.path.exists(data_folder):
        os.mkdir(Path(data_folder))

    log_folder = Path(str(path_to_outdirs) + "/5_pre_processing/_log")
    if not os.path.exists(log_folder):
        os.mkdir(Path(log_folder))

    ## print standard starting message
    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting dereplication and chimera removal for " + str(n_files) + " files.")

    if n_files == 0:
        sg.PopupError("Warning: No files found! Please check your the merged reads!")

    else:
        ## dereplication and chimera removal command
        def derep_uchime_command(i, file, data_folder, path_to_outdirs, log_df):
            # file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/4_read_filtering/_data/Dessau_MuldeFL_5A_180419_a.fasta"

            sample_name = str(Path(file).stem).replace(".fasta", "")

            ## create a new folder for each file
            dirname = Path(str(data_folder) + "/" + sample_name)
            if not os.path.exists(dirname):
                os.mkdir(dirname)

            ## create new files
            file_derep_singletons = Path(str(dirname) + "/" + sample_name + "_derep_singletons.fasta")
            file_derep_singletons_nochimeras = Path(str(dirname) + "/" + sample_name + "_derep_singletons_nochimeras.fasta")
            stderr_file = Path(str(log_folder) + "/" + sample_name + "_stderr.log")

            ## count the filtered reads from the log file
            try:
                reads = log_df.loc[log_df['Sample ID'] == sample_name]["Read filtering"].values.tolist()[0]
            except:
                reads = 0

            if reads == 0:
                print("(" + str(i+1) + "/" + str(n_files) + ")\tNo reads to process:", Path(file).stem)

            else:
                f = open(stderr_file, "w")
                # DEREPLICATION, DISCARD SINGLETONS
                subprocess.call(["vsearch", "--derep_fulllength", file, "-output", file_derep_singletons, "-sizein", "-sizeout", "-minuniquesize", "2"], stderr=f)
                # CHIMERA REMOVAL
                subprocess.call(["vsearch", "--uchime_denovo", file_derep_singletons, "--nonchimeras", file_derep_singletons_nochimeras], stderr=f)
                f.close()

                ## count the written reads
                dereplicated_reads = count_fasta_reads(file_derep_singletons)
                non_chimera_reads = count_fasta_reads(file_derep_singletons_nochimeras)
                try:
                    non_chimera_reads_perc = round(100 - round(non_chimera_reads / dereplicated_reads *100, 1),1)
                except:
                    non_chimera_reads_perc = 0

                print("(" + str(i+1) + "/" + str(n_files) + ")\tUnique reads:", dereplicated_reads, "- Chimeras removed:", non_chimera_reads_perc, "% -", sample_name)

        results = Parallel(n_jobs=int(num_cores_to_use))(delayed(derep_uchime_command)(i, file, data_folder, path_to_outdirs, log_df) for i, file in enumerate(input_files))


        ## store the merged reads for the log file
        singletons_removed_list, chimeras_removed_list = [], []

        for sample_name in log_df["Sample ID"].values.tolist():
            dirname = Path(str(data_folder) + "/" + sample_name)
            file_derep_singletons = Path(str(dirname) + "/" + sample_name + "_derep_singletons.fasta")
            file_derep_singletons_nochimeras = Path(str(dirname) + "/" + sample_name + "_derep_singletons_nochimeras.fasta")

            try:
                f = open(file_derep_singletons, "r")
                reads = sum([int(line.split("size=")[-1].strip()) for line in f if "size=" in line])
                singletons_removed_list.append(reads)
            except:
                singletons_removed_list.append(0)

            try:
                f = open(file_derep_singletons_nochimeras, "r")
                reads = sum([int(line.split("size=")[-1].strip()) for line in f if "size=" in line])
                chimeras_removed_list.append(reads)
            except:
                chimeras_removed_list.append(0)

        ## write to log sheet
        log_df["Singleton removal"] = singletons_removed_list
        log_df["Chimera removal"] = chimeras_removed_list
        log_df.to_excel(main_log_file, index=False)

        ## print standard closing message
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Dereplication, singleton removal and chimera removal for " + str(n_files) + " files.")


        #########################
        ## pooling
        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Starting pooling reads of " + str(n_files) + " files.")

        files = [file for file, reads in log_df[["Sample ID", "Chimera removal"]].values.tolist() if reads != 0]
        pooled_files_fasta = Path(str(data_folder) + "/pooled_files.fasta")

        with open(pooled_files_fasta, 'w') as outfile:
            for fname in files:
                fasta = Path(str(data_folder) + "/" + fname + "/" + fname + "_derep_singletons_nochimeras.fasta")
                with open(fasta) as infile:
                    for line in infile:
                        outfile.write(line)

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished pooling reads of " + str(n_files) + " files.")

##
