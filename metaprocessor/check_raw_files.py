def check_raw_files(nc_ID, path_to_outdirs, print_handle, window):

    from pathlib import Path
    import glob, gzip, os, datetime
    from metaprocessor import file_pairs
    import PySimpleGUI as sg

    ## collect the files from the inout folder
    raw_data_folder = Path(str(path_to_outdirs) + "/0_raw_data/_data")
    input_files = sorted(glob.glob(str(raw_data_folder) + "/*"))

    ## count the number of files
    n_files = len(input_files)

    ## check if files are present in the raw data folder
    if n_files == 0:
        print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": No files were found in the raw data folder!")
        sg.PopupError("No files were found in the raw data folder!")

    else:
        ## check if all are in .fastq.gz format
        check_suffix = [file for file in input_files if Path(file).suffixes == ['.fastq', '.gz']]

        if len(check_suffix) == n_files:
            print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": All files are in .fastq.gz format.")

        else:
            fastq_format = [file for file in input_files if Path(file).suffix == ".fastq"]
            if len(fastq_format + check_suffix) == n_files:
                answer = sg.PopupOKCancel("All files must be compressed.\nUncompressed files will be converted to .fastq.gz format.\nContinue?")
                if answer == "OK":
                    for file in fastq_format:
                        file_gzip = Path(file + ".gz")
                        with open(file, 'rb') as f_in, gzip.open(file_gzip, 'wb') as f_out:
                            f_out.writelines(f_in)
                        os.remove(file)
                    print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": All files are in .fastq.gz format.")
            else:
                print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": All files must be in .fastq.gz format.")
                sg.PopupError("All files must be in .fastq.gz format!")

        ## check for the suffixes
        pairs = [i for i in file_pairs.main(input_files) if len(i) == 2]
        n_pairs = len(pairs)
        if n_pairs * 2 != n_files:
            print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Could not find pairs for all files!")
            sg.PopupError("Could not find pairs for all files!")
        else:
            print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Found " + str(n_files) + " files and " + str(n_pairs) + " pairs.")

        ## count the number of negative controls
        n_ncs = len([file for file in input_files if nc_ID in file])
        print_handle.print(datetime.datetime.now().strftime("%H:%M:%S") + ": Found " + str(n_ncs) + " labeled negative controls.")
