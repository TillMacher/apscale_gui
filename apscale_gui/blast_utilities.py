fasta_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/7_clustering/OTUs_a0.97/OTU_a0.97.fasta"
read_table = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/7_clustering/OTUs_a0.97/OTU_a0.97.xlsx"
blast_xml_files = ["/Users/tillmacher/Downloads/BJ1H0WUB01R-Alignment.xml"]
taxonomy_table = "/Users/tillmacher/Desktop/MP_Projects/Projects/z_TEST/7_clustering/OTUs_a0.97/OTU_a0.97_taxonomy.xlsx"
limit = 5
batch_size = 10


def subset_fasta(fasta_file, batch_size):

    from pathlib import Path
    from Bio import SeqIO
    import datetime

    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Create subsets from the fasta file.")

    batch_size = int(batch_size)

    fasta_file = Path(fasta_file)
    ## create batches from the main fasta file
    with open(fasta_file) as handle:
        i = 1
        n = 1
        chunk_fasta_files = []
        for record in SeqIO.parse(handle, "fasta"):
            ## create a new fasta file for each chunk
            chunk_fasta = Path(str(fasta_file.parent) + "/" + str(i) + ".fasta")
            ## save the name of all batches
            if chunk_fasta not in chunk_fasta_files:
                chunk_fasta_files.append(chunk_fasta)
            ## write the record to the respective file
            with open(chunk_fasta, "a") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
            ## proceed to next chunk
            if n == batch_size:
                n = 1
                i += 1
            else:
                n += 1

def blast_xml_to_taxonomy(fasta_file, blast_xml_files, read_table, limit):

    import datetime, sys, re, subprocess, itertools
    import pandas as pd
    from ete3 import NCBITaxa
    import requests
    import xmltodict
    import time
    import PySimpleGUI as sg
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
    from Bio import SeqIO
    from pathlib import Path
    from Bio.Blast import NCBIXML

    ########################################################################################
    ## define function

    def sort_df(df):
        sort_col = [int(n.split("_")[1]) for n in df["ID"]]
        df["sort"] = sort_col
        df = df.sort_values("sort")
        return df.drop(["sort"], axis=1)

    def open_file(file):
        if sys.platform == "win32":
            os.startfile(file)
        else:
            opener = "open" if sys.platform == 'darwin' else 'xdg-open'
            subprocess.call([opener, file])

    def get_desired_ranks(taxid, desired_ranks):
        ncbi = NCBITaxa()
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

    def ncbi_taxid_request(hit):
        desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        taxonomy_list = []
        taxid = hit[1]
        try:
            results = get_desired_ranks(taxid, desired_ranks)
            taxids = [str(taxid) for taxid in list(results.values())]

            # if the taxonomy is not present
            # DO THIS
            if '<not present>' in taxids:
                for taxid in taxids:
                    if taxid != '<not present>':
                        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=" + str(taxid)
                        response = requests.get(url)
                        data = xmltodict.parse(response.content)
                        for entry in data['eSummaryResult']['DocSum']['Item']:
                            if entry['@Name'] == 'ScientificName':
                                name = entry['#text']
                                taxonomy_list.append(name)
                        time.sleep(0.2)
                    else:
                        taxonomy_list.append("")

            # if all taxonomy information is present
            # DO THIS
            else:
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=" + ','.join(taxids)
                response = requests.get(url)
                data = xmltodict.parse(response.content)
                for entry in data['eSummaryResult']['DocSum']:
                    for item in entry['Item']:
                        if item['@Name'] == 'ScientificName':
                            name = item['#text']
                            taxonomy_list.append(name)
            return taxonomy_list
        except:
            return "taxid_not_found"

    ########################################################################################
    ## start script

    limit = int(limit)
    taxonomy_table = str(Path(fasta_file)).replace(".fasta", "_taxonomy_table.xlsx")

    with pd.ExcelWriter(taxonomy_table) as writer:
        hit_list = []

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Collecting Taxids from the blast result file(s).")

        if type(blast_xml_files) == "str":
            blast_xml_files = [blast_xml_files]

        for results_xml in blast_xml_files:
            ################################################################################
            # Convert the xml to a hit table
            f = open(results_xml)
            for line in f:
                # collect the query name
                if "<query-title>" in line:
                    query = re.split('>|<', line)[2]
                # collect the query sequence length
                if "<query-len>" in line:
                    query_len = re.split('>|<', line)[2]
                # collect the taxonomy id
                if "<taxid>" in line:
                    taxid = re.split('>|<', line)[2]
                # calculate the similarity
                if "<identity>" in line:
                    identity = re.split('>|<', line)[2]
                    diff = int(query_len) - int(identity)
                    perc_id = round(100 - diff / int(query_len) * 100, 2)
                    hit_list.append([query, taxid, perc_id, "NCBI"])
            f.close()

        ## filter list according to limit set by user
        hit_list_filtered = []
        for hit in hit_list:
            ID = hit[0]
            n_occurences = len([i for i in hit_list_filtered if ID in i])
            if n_occurences < limit:
                hit_list_filtered.append(hit)

        hit_table_df = pd.DataFrame(hit_list_filtered, columns=["ID", "Taxonomy", "Similarity", "Status"])
        hit_table_df.to_excel(writer, sheet_name='Raw hits', index=False)

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Pre-sorting hits by similarity.")

        query_list = list(set(hit_table_df["ID"].tolist()))
        sorted_hits_dict = {}

        for query in sorted(query_list):
            prev_hit = False
            hit_list = hit_table_df.loc[hit_table_df["ID"] == query][["Taxonomy", "Similarity"]].values.tolist()
            similar_hit_list = [hit_list[0]]

            if len([list(x) for x in set(tuple(x) for x in hit_list)]) == 1:
                sorted_hits_dict[query] = similar_hit_list

            else:
                for hit in hit_list:
                    if prev_hit == False:
                        prev_hit = hit
                    else:
                        if hit[1] < prev_hit[1]:
                            sorted_hits_dict[query] = similar_hit_list
                            break
                        else:
                            similar_hit_list.append(hit)
                        prev_hit = hit

            if query not in sorted_hits_dict.keys():
                similar_hit_list.sort()
                similar_hit_list = [list(x) for x in set(tuple(x) for x in similar_hit_list)]
                sorted_hits_dict[query] = similar_hit_list

        similarity_filtered_list = []
        for key, values in sorted_hits_dict.items():
            unique_values = [list(x) for x in set(tuple(x) for x in values)]
            if len(unique_values) == 1:
                similarity_filtered_list.append([key] + values[0])
            else:
                similarity_filtered_list = similarity_filtered_list + [[key] + value for value in unique_values]

        ################################################################################
        # Download the NCBI taxonomy for all hits

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Downloading taxonomy from NCBI.")

        hit_list_2 = []

        for hit in similarity_filtered_list:
            print(hit[0])
            try:
                taxonomy_list = ncbi_taxid_request(hit)
            except:
                time.sleep(3)
                taxonomy_list = ncbi_taxid_request(hit)

            if taxonomy_list != "taxid_not_found":
                hit_list_2.append([hit[0]] + taxonomy_list + [hit[2]] + ["NCBI"])

        hit_table_2_df = pd.DataFrame(hit_list_2, columns=["ID","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
        hit_table_2_df = sort_df(hit_table_2_df)
        hit_table_2_df.to_excel(writer, sheet_name='Taxonomy added', index=False)

        ################################################################################
        # Filter the table according to the JAMP filtering method

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Filtering hits by similarity thresholds.")

        hit_list_3 = []
        for hit in hit_list_2:
            identity = hit[-2]
            # remove empty hits first
            if hit[1] != '':
                if 100 >= identity >= 98:
                    if hit[1:7] != ['']*6:
                        hit_list_3.append(hit)
                elif 98 >= identity >= 95:
                    hit[6] = ''
                    if hit[1:7] != ['']*6:
                        hit_list_3.append(hit)
                elif 95 >= identity >= 90:
                    hit[5], hit[6] = '', ''
                    if hit[1:7] != ['']*6:
                        hit_list_3.append(hit)
                elif 90 >= identity >= 85:
                    hit[4], hit[5], hit[6] = '', '', ''
                    if hit[1:7] != ['']*6:
                        hit_list_3.append(hit)
                else:
                    hit[3], hit[4], hit[5], hit[6] = '', '', '', ''
                    if hit[1:7] != ['']*6:
                        hit_list_3.append(hit)

        hit_list_3 = [list(x) for x in set(tuple(x) for x in hit_list_3)]

        hit_table_3_df = pd.DataFrame(hit_list_3, columns=["ID","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
        hit_table_3_df = sort_df(hit_table_3_df)
        hit_table_3_df.to_excel(writer, sheet_name='JAMP filtering', index=False)

        ################################################################################
        # Remove remaining duplicates and ambigious taxonomies

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Removing ambigious taxonomy.")

        IDs = hit_table_3_df["ID"].values.tolist()
        duplicates = list(set([x for n, x in enumerate(IDs) if x in IDs[:n]]))
        duplicates_dict = {}

        for hit in sorted(hit_list_3):
            if (hit[0] in duplicates and hit[0] not in duplicates_dict.keys()):
                duplicates_dict[hit[0]] = [hit]
            elif hit[0] in duplicates:
                duplicates_dict[hit[0]] = duplicates_dict[hit[0]] + [hit]

        for key in duplicates_dict.keys():
            # remove: species
            for entry in duplicates_dict[key]:
                entry[6] = ""
            if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
            else:
                # remove: genus
                for entry in duplicates_dict[key]:
                    entry[5] = ""
                if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                    duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                else:
                    # remove: family
                    for entry in duplicates_dict[key]:
                        entry[4] = ""
                    if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                        duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                    else:
                        # remove: order
                        for entry in duplicates_dict[key]:
                            entry[3] = ""
                        if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                            duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                        else:
                            # remove: class
                            for entry in duplicates_dict[key]:
                                entry[2] = ""
                                duplicates_dict[key] = entry

        hit_list_4 = []
        for hit in hit_list_3:
            if hit[0] in duplicates_dict.keys():
                hit_list_4.append(duplicates_dict[hit[0]])
            else:
                hit_list_4.append(hit)
        hit_list_4 = [list(x) for x in set(tuple(x) for x in hit_list_4)]
        present_OTUs = [OTU[0] for OTU in hit_list_4]

        # lastly add the OTU that had ch against the NCBI database
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id not in present_OTUs:
                    hit_list_4.append([record.id] + ['No Match'] * 7 + ['NCBI'])

        hit_table_4_df = pd.DataFrame(hit_list_4, columns=["ID","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
        #hit_table_4_df = sort_df(hit_table_4_df)
        hit_table_4_df.to_excel(writer, sheet_name='JAMP hit', index=False)

        print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished writing taxonomy table.")

        answer = sg.PopupYesNo("Open taxonomy table?")
        if answer == "Yes":
            open_file(Path(taxonomy_table))

def ESV_reference_taxonomy(reference_fasta, blast_xml_files, limit):
    import datetime, sys, re, subprocess, itertools
    import pandas as pd
    from ete3 import NCBITaxa
    import requests
    import xmltodict
    import time
    import PySimpleGUI as sg
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
    from Bio import SeqIO
    from pathlib import Path
    from Bio.Blast import NCBIXML

    # reference_fasta = "/Users/tillmacher/Documents/GitHub/MetaProcessor/metaprocessor/user_data/ESV_references/tele02_reference.fasta"
    # reference_xml = "/Users/tillmacher/Downloads/reference.xml"
    # limit = 10


    ########################################################################################
    ## define function

    def sort_df(df):
        sort_col = [int(n.split("_")[1]) for n in df["ID"]]
        df["sort"] = sort_col
        df = df.sort_values("sort")
        return df.drop(["sort"], axis=1)

    def open_file(file):
        if sys.platform == "win32":
            os.startfile(file)
        else:
            opener = "open" if sys.platform == 'darwin' else 'xdg-open'
            subprocess.call([opener, file])

    def get_desired_ranks(taxid, desired_ranks):
        ncbi = NCBITaxa()
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

    def ncbi_taxid_request(hit):
        desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        taxonomy_list = []
        taxid = hit[1]
        try:
            results = get_desired_ranks(taxid, desired_ranks)
            taxids = [str(taxid) for taxid in list(results.values())]

            # if the taxonomy is not present
            # DO THIS
            if '<not present>' in taxids:
                for taxid in taxids:
                    if taxid != '<not present>':
                        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=" + str(taxid)
                        response = requests.get(url)
                        data = xmltodict.parse(response.content)
                        for entry in data['eSummaryResult']['DocSum']['Item']:
                            if entry['@Name'] == 'ScientificName':
                                name = entry['#text']
                                taxonomy_list.append(name)
                        time.sleep(0.2)
                    else:
                        taxonomy_list.append("")

            # if all taxonomy information is present
            # DO THIS
            else:
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=" + ','.join(taxids)
                response = requests.get(url)
                data = xmltodict.parse(response.content)
                for entry in data['eSummaryResult']['DocSum']:
                    for item in entry['Item']:
                        if item['@Name'] == 'ScientificName':
                            name = item['#text']
                            taxonomy_list.append(name)
            return taxonomy_list
        except:
            return "taxid_not_found"

    ########################################################################################
    ## start script

    limit = int(limit)
    hit_list = []

    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Collecting Taxids from the blast result file(s).")

    ################################################################################
    # Convert the xml to a hit table
    f = open(reference_xml)
    for line in f:
        # collect the query name
        if "<query-title>" in line:
            query = re.split('>|<', line)[2]
        # collect the query sequence length
        if "<query-len>" in line:
            query_len = re.split('>|<', line)[2]
        # collect the taxonomy id
        if "<taxid>" in line:
            taxid = re.split('>|<', line)[2]
        # calculate the similarity
        if "<identity>" in line:
            identity = re.split('>|<', line)[2]
            diff = int(query_len) - int(identity)
            perc_id = round(100 - diff / int(query_len) * 100, 2)
            hit_list.append([query, taxid, perc_id, "NCBI"])
    f.close()

    ## filter list according to limit set by user
    hit_list_filtered = []
    for hit in hit_list:
        ID = hit[0]
        n_occurences = len([i for i in hit_list_filtered if ID in i])
        if n_occurences < limit:
            hit_list_filtered.append(hit)

    hit_table_df = pd.DataFrame(hit_list_filtered, columns=["ID", "Taxonomy", "Similarity", "Status"])

    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Pre-sorting hits by similarity.")

    query_list = list(set(hit_table_df["ID"].tolist()))
    sorted_hits_dict = {}

    for query in sorted(query_list):
        prev_hit = False
        hit_list = hit_table_df.loc[hit_table_df["ID"] == query][["Taxonomy", "Similarity"]].values.tolist()
        similar_hit_list = [hit_list[0]]

        if len([list(x) for x in set(tuple(x) for x in hit_list)]) == 1:
            sorted_hits_dict[query] = similar_hit_list

        else:
            for hit in hit_list:
                if prev_hit == False:
                    prev_hit = hit
                else:
                    if hit[1] < prev_hit[1]:
                        sorted_hits_dict[query] = similar_hit_list
                        break
                    else:
                        similar_hit_list.append(hit)
                    prev_hit = hit

        if query not in sorted_hits_dict.keys():
            similar_hit_list.sort()
            similar_hit_list = [list(x) for x in set(tuple(x) for x in similar_hit_list)]
            sorted_hits_dict[query] = similar_hit_list

    similarity_filtered_list = []
    for key, values in sorted_hits_dict.items():
        unique_values = [list(x) for x in set(tuple(x) for x in values)]
        if len(unique_values) == 1:
            similarity_filtered_list.append([key] + values[0])
        else:
            similarity_filtered_list = similarity_filtered_list + [[key] + value for value in unique_values]

    ################################################################################
    # Download the NCBI taxonomy for all hits

    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Downloading taxonomy from NCBI.")

    hit_list_2 = []

    for hit in similarity_filtered_list:
        print(hit[0])
        try:
            taxonomy_list = ncbi_taxid_request(hit)
        except:
            time.sleep(3)
            taxonomy_list = ncbi_taxid_request(hit)

        if taxonomy_list != "taxid_not_found":
            hit_list_2.append([hit[0]] + taxonomy_list + [hit[2]] + ["NCBI"])

    hit_table_2_df = pd.DataFrame(hit_list_2, columns=["ID","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
    hit_table_2_df = sort_df(hit_table_2_df)

    ################################################################################
    # Filter the table according to the JAMP filtering method

    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Filtering hits by similarity thresholds.")

    hit_list_3 = []
    for hit in hit_list_2:
        identity = hit[-2]
        # remove empty hits first
        if hit[1] != '':
            if 100 >= identity >= 98:
                if hit[1:7] != ['']*6:
                    hit_list_3.append(hit)
            elif 98 >= identity >= 95:
                hit[6] = ''
                if hit[1:7] != ['']*6:
                    hit_list_3.append(hit)
            elif 95 >= identity >= 90:
                hit[5], hit[6] = '', ''
                if hit[1:7] != ['']*6:
                    hit_list_3.append(hit)
            elif 90 >= identity >= 85:
                hit[4], hit[5], hit[6] = '', '', ''
                if hit[1:7] != ['']*6:
                    hit_list_3.append(hit)
            else:
                hit[3], hit[4], hit[5], hit[6] = '', '', '', ''
                if hit[1:7] != ['']*6:
                    hit_list_3.append(hit)

    hit_list_3 = [list(x) for x in set(tuple(x) for x in hit_list_3)]

    hit_table_3_df = pd.DataFrame(hit_list_3, columns=["ID","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
    hit_table_3_df = sort_df(hit_table_3_df)

    ################################################################################
    # Remove remaining duplicates and ambigious taxonomies

    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Removing ambigious taxonomy.")

    IDs = hit_table_3_df["ID"].values.tolist()
    duplicates = list(set([x for n, x in enumerate(IDs) if x in IDs[:n]]))
    duplicates_dict = {}

    for hit in sorted(hit_list_3):
        if (hit[0] in duplicates and hit[0] not in duplicates_dict.keys()):
            duplicates_dict[hit[0]] = [hit]
        elif hit[0] in duplicates:
            duplicates_dict[hit[0]] = duplicates_dict[hit[0]] + [hit]

    for key in duplicates_dict.keys():
        # remove: species
        for entry in duplicates_dict[key]:
            entry[6] = ""
        if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
            duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
        else:
            # remove: genus
            for entry in duplicates_dict[key]:
                entry[5] = ""
            if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
            else:
                # remove: family
                for entry in duplicates_dict[key]:
                    entry[4] = ""
                if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                    duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                else:
                    # remove: order
                    for entry in duplicates_dict[key]:
                        entry[3] = ""
                    if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                        duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                    else:
                        # remove: class
                        for entry in duplicates_dict[key]:
                            entry[2] = ""
                            duplicates_dict[key] = entry

    hit_list_4 = []
    for hit in hit_list_3:
        if hit[0] in duplicates_dict.keys():
            hit_list_4.append(duplicates_dict[hit[0]])
        else:
            hit_list_4.append(hit)
    hit_list_4 = [list(x) for x in set(tuple(x) for x in hit_list_4)]

    hit_table_4_df = pd.DataFrame(hit_list_4, columns=["ID","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
    hit_table_4_df = sort_df(hit_table_4_df)

    ## write reference fasta to a dict
    reference_fasta_dict = {}
    for record in SeqIO.parse(reference_fasta, "fasta"):
        reference_fasta_dict[record.id] = str(record.seq)

    ## write new fasta file
    reference_fasta_taxonomy = Path(str(reference_fasta).replace(".fasta", "_taxonomy.fasta"))
    f = open(reference_fasta_taxonomy, "w")
    ## write the new reference fasta
    for hit in hit_table_4_df.values.tolist():
        for tax_level in hit[1:7][::-1]:
            if tax_level != '':
                taxonomy = tax_level
                break
        new_header = taxonomy.replace(" ", "_") + "_" + hit[0]
        f.write(">" + new_header + "\n")
        sequence = reference_fasta_dict[hit[0]]
        f.write(sequence + "\n")
    f.close()

    print(datetime.datetime.now().strftime("%H:%M:%S") + ": Finished adding taxonomy to reference fasta.")
