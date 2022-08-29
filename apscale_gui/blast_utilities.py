from pathlib import Path
from Bio import SeqIO
import datetime, sys, re, subprocess, itertools, os, time, shutil, glob, pkg_resources
import pandas as pd
from ete3 import NCBITaxa
import requests
import xmltodict
import PySimpleGUI as sg
from pathlib import Path
from tqdm import tqdm
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from requests_html import HTMLSession
import numpy as np

################################################################################
## general blast functions

def sort_df(df):
    sort_col = [int(n.split('_')[1]) for n in df['ID']]
    df['sort'] = sort_col
    df = df.sort_values('sort')
    return df.drop(['sort'], axis=1)

def open_file(file):
    if sys.platform == 'win32':
        os.startfile(file)
    else:
        opener = 'open' if sys.platform == 'darwin' else 'xdg-open'
        subprocess.call([opener, file])

def get_desired_ranks(taxid, desired_ranks):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

def ncbi_taxid_request(taxid):
    desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    taxonomy_list = []
    try:
        results = get_desired_ranks(taxid, desired_ranks)
        taxids = [str(taxid) for taxid in list(results.values())]

        # if the taxonomy is not present
        # DO THIS
        if '<not present>' in taxids:
            for taxid in taxids:
                if taxid != '<not present>':
                    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=' + str(taxid)
                    response = requests.get(url)
                    data = xmltodict.parse(response.content)
                    for entry in data['eSummaryResult']['DocSum']['Item']:
                        if entry['@Name'] == 'ScientificName':
                            name = entry['#text']
                            taxonomy_list.append(name)
                    time.sleep(0.2)
                else:
                    taxonomy_list.append('')

        # if all taxonomy information is present
        # DO THIS
        else:
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=' + ','.join(taxids)
            response = requests.get(url)
            data = xmltodict.parse(response.content)
            for entry in data['eSummaryResult']['DocSum']:
                for item in entry['Item']:
                    if item['@Name'] == 'ScientificName':
                        name = item['#text']
                        taxonomy_list.append(name)
        return taxonomy_list
    except:
        return ['No Match'] * 6

################################################################################
## online blast functions

def subset_fasta(fasta_file, batch_size):

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Create subsets from the fasta file.')

    batch_size = int(batch_size)

    fasta_file = Path(fasta_file)

    ## create a new folder for the database
    subset_folder = Path(fasta_file.parent).joinpath('fasta_subsets')
    try:
        os.mkdir(subset_folder)
    except FileExistsError:
        pass

    ## create batches from the main fasta file
    with open(fasta_file) as handle:
        i = 1
        n = 1
        chunk_fasta_files = []
        for record in SeqIO.parse(handle, 'fasta'):
            ## create a new fasta file for each chunk
            chunk_fasta = '{}/subset_{}.fasta'.format(subset_folder, i)
            ## save the name of all batches
            if chunk_fasta not in chunk_fasta_files:
                chunk_fasta_files.append(chunk_fasta)
            ## write the record to the respective file
            with open(chunk_fasta, 'a') as output_handle:
                SeqIO.write(record, output_handle, 'fasta')
            ## proceed to next chunk
            if n == batch_size:
                n = 1
                i += 1
            else:
                n += 1

    sg.Popup('{} subsets were created from the original .fasta file.'.format(i))

def blast_xml_to_taxonomy(fasta_file, blast_xml_files, read_table, limit):

    ########################################################################################
    ## start script

    limit = int(limit)
    taxonomy_table = Path(str(Path(fasta_file)).replace('.fasta', '_taxonomy_table.xlsx'))

    writer = pd.ExcelWriter(taxonomy_table, engine='xlsxwriter')

    hit_list = []

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Loading result file(s).')

    if isinstance(blast_xml_files, str) == True:
        blast_xml_files = [blast_xml_files]

    for results_xml in tqdm(blast_xml_files):
        ################################################################################
        # Convert the xml to a hit table
        f = open(results_xml)
        for line in f:
            # collect the query name
            if '<query-title>' in line:
                query = re.split('>|<', line)[2]
            # collect the taxonomy id
            if '<taxid>' in line:
                taxid = re.split('>|<', line)[2]
            # calculate the similarity
            if '<identity>' in line:
                identity = re.split('>|<', line)[2]
            # collect the query sequence length
            if '<align-len>' in line:
                align_len = re.split('>|<', line)[2]
                perc_id = round(int(identity) / int(align_len) * 100, 2)
                hit_list.append([query, taxid, perc_id, 'NCBI'])
        f.close()

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Collecting Taxids from the blast result file(s).')

    ## filter list according to limit set by user
    hit_list_filtered = []
    for hit in tqdm(hit_list):
        ID = hit[0]
        n_occurences = len([i for i in hit_list_filtered if ID in i])
        if n_occurences < limit:
            hit_list_filtered.append(hit)

    hit_table_df = pd.DataFrame(hit_list_filtered, columns=['ID', 'Taxonomy', 'Similarity', 'Status'])
    hit_table_df.to_excel(writer, sheet_name='Raw hits', index=False)

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Pre-sorting hits by similarity.')

    query_list = list(set(hit_table_df['ID'].tolist()))
    sorted_hits_dict = {}

    for query in tqdm(sorted(query_list)):
        prev_hit = False
        hit_list = hit_table_df.loc[hit_table_df['ID'] == query][['Taxonomy', 'Similarity']].values.tolist()
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

    hit_table_df = pd.DataFrame(similarity_filtered_list, columns=['ID', 'taxid', 'Similarity'])

    ################################################################################
    # Download the NCBI taxonomy for all hits

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Downloading taxonomy from NCBI.')

    hit_list_2 = []
    taxid_dict = {}

    for hit in tqdm(similarity_filtered_list):
        taxid = hit[1]

        if taxid in taxid_dict.keys():
            taxonomy_list = taxid_dict[taxid]
        else:
            try:
                taxonomy_list = ncbi_taxid_request(hit[1])
                taxid_dict[taxid] = taxonomy_list
            except ValueError:
                taxonomy_list = 'taxid_not_found'
                taxid_dict[taxid] = taxonomy_list

        if taxonomy_list != 'taxid_not_found':
            hit_list_2.append([hit[0]] + taxonomy_list + [hit[2]] + ['NCBI'])

    ################################################################################
    ## remove 'uncultured' hits
    hit_list_2_filtered =[]
    for hit in hit_list_2:
        if 'uncultured' in hit[6]:
            hit[1] = hit[6]
            hit[6] = ''
            hit_list_2_filtered.append(hit)
        else:
            hit_list_2_filtered.append(hit)

    hit_table_2_df = pd.DataFrame(hit_list_2_filtered, columns=['ID','Phylum','Class','Order','Family','Genus','Species', 'Similarity', 'Status'])
    hit_table_2_df = sort_df(hit_table_2_df)

    ################################################################################
    # Filter the table according to the JAMP filtering method

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Filtering hits by similarity thresholds.')

    hit_list_3 = []
    for hit in tqdm(hit_list_2_filtered):
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

    hit_table_3_df = pd.DataFrame(hit_list_3, columns=['ID','Phylum','Class','Order','Family','Genus','Species', 'Similarity', 'Status'])
    hit_table_3_df = sort_df(hit_table_3_df)

    # hit_table_3_df.to_excel(writer, sheet_name='JAMP filtering', index=False)

    ################################################################################
    # Remove remaining duplicates and ambigious taxonomies

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Removing ambigious taxonomy.')

    IDs = hit_table_3_df['ID'].values.tolist()
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
            entry[6] = ''
        if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
            duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
        else:
            # remove: genus
            for entry in duplicates_dict[key]:
                entry[5] = ''
            if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
            else:
                # remove: family
                for entry in duplicates_dict[key]:
                    entry[4] = ''
                if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                    duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                else:
                    # remove: order
                    for entry in duplicates_dict[key]:
                        entry[3] = ''
                    if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                        duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                    else:
                        # remove: class
                        for entry in duplicates_dict[key]:
                            entry[2] = ''
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
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id not in present_OTUs:
                hit_list_4.append([record.id] + ['No Match'] * 7 + ['NCBI'])

    hit_table_4_df = pd.DataFrame(hit_list_4, columns=['ID','Phylum','Class','Order','Family','Genus','Species', 'Similarity', 'Status'])
    #hit_table_4_df = sort_df(hit_table_4_df)
    hit_table_4_df.to_excel(writer, sheet_name='JAMP hit', index=False)

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Finished writing taxonomy table.')

    answer = sg.PopupYesNo('Open taxonomy table?')
    if answer == 'Yes':
        open_file(Path(taxonomy_table))

    sg.Popup('Taxonomy table is found under {}.'.format(taxonomy_table))

    writer.save()

################################################################################
## general blastn function
def blastn(query_fasta, blastn_database, project_folder, n_threads, task):

    if task == 'Highly similar sequences (megablast)':
        task = 'megablast'
    elif task == 'More dissimilar sequences (discontiguous megablast)':
        task = 'dc-megablast'
    elif task == 'Somewhat similar sequences (blastn)':
        task = 'blastn'

    ## replace spaces and dots in the filename
    filename = Path(query_fasta).stem
    filename = filename.replace('.', '_').replace(' ', '_')

    print('{}: Starting {} for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), task, filename))
    print('{}: Your database: {}'.format(datetime.datetime.now().strftime('%H:%M:%S'),  Path(blastn_database).stem))

    ## collect files and folders
    query_fasta = Path(query_fasta)
    db_folder = Path(blastn_database).joinpath('db')

    ## create a new folder for each blast search
    blast_search = 'BLAST_' + query_fasta.stem + '_(' + datetime.datetime.now().strftime('%D_%H_%M').replace('/', '_') + ')'
    blast_folder = Path(project_folder).joinpath('10_local_BLAST', blast_search)
    try:
        os.mkdir(blast_folder)
    except FileExistsError:
        pass

    ## create the output file
    blast_csv = blast_folder.joinpath(query_fasta.stem + '_' + task + '.csv')

    blastn_exe = "blastn"

    ## run blast search
    subprocess.call([blastn_exe, '-task', task, '-db', str(db_folder), '-query', str(query_fasta), '-num_threads', str(n_threads),  '-outfmt', '6 delim=;; qseqid sseqid pident evalue', '-out', str(blast_csv)])

    ## write the name of the database
    blast_log = blast_folder.joinpath('log.txt')
    f = open(blast_log, 'w')
    f.write('Your database: {}\n'.format(Path(blastn_database).stem))
    f.write('Your task: {}\n'.format(task))
    f.close()

    print('{}: Finished {} for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), task, filename))
    print()

################################################################################
## databases
## diat.barcode
def create_database_diat_barcode(diatbarcode_xlsx, project_folder, databases_path):
    ' Convert the diat.barcode file to fasta and create an new database using makeblastdb '

    ## collect files
    diatbarcode_xlsx = Path(diatbarcode_xlsx)
    filename = Path(diatbarcode_xlsx).stem
    print('{}: Creating fasta file from \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), filename))

    ## replace spaces and dots in the filename
    filename = filename.replace('.', '_').replace(' ', '_')

    ## load dataframe
    diat_barcode_df = pd.read_excel(diatbarcode_xlsx, sheet_name='data')

    ## create a new folder for the database
    db_folder = Path(databases_path).joinpath('DB_' + filename)
    try:
        os.mkdir(db_folder)
    except FileExistsError:
        pass

    ## create a fasta file
    fasta_file = Path(db_folder).joinpath('db.fasta')
    f = open(fasta_file, 'w')
    for sequence in diat_barcode_df[['Sequence ID', 'Sequence']].values.tolist():
        header = '>{}\n'.format(sequence[0])
        sequence = '{}\n'.format(sequence[1])
        f.write(header)
        f.write(sequence)
    f.close()

    ## build a new database
    print('{}: Starting to build a new database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    makeblastdb_exe = "makeblastdb"

    db_name = Path(db_folder).joinpath('db')
    subprocess.call([makeblastdb_exe, '-in', str(fasta_file), '-dbtype', 'nucl', '-out', str(db_name)])

    ## create a copy of the database xlsx
    diatbarcode_xlsx_db = db_folder.joinpath('db.xlsx')
    shutil.copyfile(diatbarcode_xlsx, diatbarcode_xlsx_db)

    print('{}: Finished building database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    sg.Popup('Finished building database.', title='Finished')

## midori2 database
def create_database_midori2(midori2_fasta, project_folder, databases_path):
    ' Create new midori2 database using makeblastdb '

    ## collect files
    midori2_fasta = Path(midori2_fasta)
    filename = Path(midori2_fasta).stem

    ## replace spaces and dots in the filename
    filename = filename.replace('.', '_').replace(' ', '_')

    ## create a new folder for the database
    db_folder = Path(databases_path).joinpath('DB_' + filename)
    try:
        os.mkdir(db_folder)
    except FileExistsError:
        pass

    ## build a new database
    print('{}: Starting to build a new database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    makeblastdb_exe = "makeblastdb"

    db_name = Path(db_folder).joinpath('db')
    subprocess.call([makeblastdb_exe, '-in', str(midori2_fasta), '-dbtype', 'nucl', '-out', str(db_name)])

    print('{}: Finished building database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    sg.Popup('Finished building database.', title='Finished')

## trnl P6 loop EMBL EBI
def create_database_trnl_P6(trnl_p6_xlsx, project_folder, databases_path):
    ' Convert the diat.barcode file to fasta and create an new database using makeblastdb '

    ## collect files
    trnl_p6_xlsx = Path(trnl_p6_xlsx)
    filename = Path(trnl_p6_xlsx).stem
    print('{}: Creating fasta file from \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), filename))

    ## replace spaces and dots in the filename
    filename = filename.replace('.', '_').replace(' ', '_')

    ## load dataframe
    trnl_p6_df = pd.read_excel(trnl_p6_xlsx)

    ## create a new folder for the database
    db_folder = Path(databases_path).joinpath('DB_' + filename)
    try:
        os.mkdir(db_folder)
    except FileExistsError:
        pass

    ## create a fasta file
    fasta_file = Path(db_folder).joinpath('db.fasta')
    f = open(fasta_file, 'w')
    for line in trnl_p6_df.values.tolist():
        if '>' in line[0]:
            taxid = [i.replace(';', '') for i in line[0].split() if 'taxid=' in i and 'merged' not in i]
            writer = '>{}\n'.format(taxid[0])
        else:
            writer = '{}\n'.format(line[0])
        f.write(writer)
    f.close()

    ## build a new database
    print('{}: Starting to build a new database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    makeblastdb_exe = "makeblastdb"

    db_name = Path(db_folder).joinpath('db')
    subprocess.call([makeblastdb_exe, '-in', str(fasta_file), '-dbtype', 'nucl', '-out', str(db_name)])

    print('{}: Finished building database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    sg.Popup('Finished building database.', title='Finished')

## local NCBI-mined DB

def create_database_NCBI(genbank_file, project_folder, databases_path):
    ' Create a database from the NCBI fasta using makeblastdb'

    ## build a new database
    print('{}: Collecting records from genbank file to build a new database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect files
    genbank_file = Path(genbank_file)
    filename = Path(genbank_file).stem

    ## create a new folder for the database
    db_folder = Path(databases_path).joinpath('DB_' + filename)
    try:
        os.mkdir(db_folder)
    except FileExistsError:
        pass

    ## collect taxids -> should be possible with seqio but isnt!
    fasta_file = Path(db_folder).joinpath('db.fasta')
    f_fasta = open(fasta_file, 'w')
    f_gb = open(genbank_file, 'r')
    skipped = []

    sequences = False
    for line in f_gb.readlines():

        if line.startswith('ACCESSION'):
            accession = line.split()[1]

        if '/db_xref="taxon' in line:
            taxid = line.strip().split(':')[-1].replace('\"', '')

        if line.startswith('//'):

            if accession == 'NW_019713564':
                break

            seq = ''.join(record)
            allowed_s = "acgtuwsmkrybdhvn"
            if all(ch in allowed_s for ch in seq) == True and seq != '':
                id = '>{}_{}\n'.format(accession, taxid)
                seq = '{}\n'.format(seq)
                f_fasta.write(id)
                f_fasta.write(seq)
                sequences = False
                record = []
            else:
                skipped.append(accession)
                sequences = False
                record = []

        if sequences == True:
            record.append(''.join(line.split()[1::]))

        if line.startswith('ORIGIN'):
            sequences = True
            record = []

    f_gb.close()
    f_fasta.close()

    print('{}: Skipped {} malformated record(s).'.format(datetime.datetime.now().strftime('%H:%M:%S'), len(skipped)))

    makeblastdb_exe = "makeblastdb"

    ## create db
    print('{}: Starting to build a new database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    db_name = Path(db_folder).joinpath('db')
    subprocess.call([makeblastdb_exe, '-in', str(fasta_file), '-dbtype', 'nucl', '-out', str(db_name)])

    print('{}: Finished building database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    sg.Popup('Finished building database.', title='Finished')

def accession2taxid(accession):
    url = 'https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(accession)
    as_session = HTMLSession()
    as_session.headers.update({'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/98.9.4758.82 Safari/537.36'})
    retry_strategy = Retry(total = 10, status_forcelist = [400, 401, 403, 404, 429, 500, 502, 503, 504], backoff_factor = 1)
    adapter = HTTPAdapter(max_retries = retry_strategy)
    as_session.mount('https://', adapter)
    as_session.mount('http://', adapter)
    r = as_session.get(url, timeout = 300)
    data = r.text.split(';')
    taxid = [i for i in data if '?ORGANISM' in i][0].split('?')[-1].replace('&amp', '').replace('ORGANISM=', '')

    return taxid

def filter_blastn_results_NCBI(blast_csv, read_table, project_folder):

    " Filter the BLAST results and download the according taxonomy "

    ## collect files and folders
    blast_csv = Path(blast_csv)
    read_table = Path(read_table)

    ## replace spaces and dots in the filename
    filename = Path(blast_csv).stem

    ## create an output file
    taxonomy_table_xlsx = blast_csv.parent.joinpath(filename + '_taxonomy_table.xlsx')

    print('{}: Starting to filter BLAST results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), filename))

    ## load blast results
    blast_df = pd.read_csv(blast_csv, header=None)
    blast_df.columns = ['ID', 'Sequence ID', 'Similarity']

    print('{}: Building reference dictionary.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    print('{}: Filtering hits by similarity.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## filter by similarity
    sorted_hits_dict = {}
    for query in tqdm(set(blast_df['ID'].values.tolist())):
        prev_hit = False
        hit_list = blast_df.loc[blast_df['ID'] == query][['Sequence ID','Similarity']].values.tolist()
        hit_list_sorted = pd.DataFrame(hit_list).sort_values(1).values.tolist()
        similar_hit_list = [hit_list[0]]

        ## collect all hits that have the same similarity (beginning from the best hit)
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

        ## if there are only identical hits, add them all!
        if query not in sorted_hits_dict.keys():
            similar_hit_list.sort()
            similar_hit_list = [list(x) for x in set(tuple(x) for x in similar_hit_list)]
            sorted_hits_dict[query] = similar_hit_list

    print('{}: Fetching taxonomy from NCBI.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## replace accession number with taxids
    for OTU, hit in tqdm(sorted_hits_dict.items()):
        sublist = []
        for i in hit:
            if 'gb|' in i[0]:
                accession = i[0].split('|')[1]
            else:
                accession = i
            try:
                time.sleep(1)
                taxid = accession2taxid(accession)
                taxonomy = ncbi_taxid_request(taxid)
                sublist.append([taxonomy, i[1]])
            except IndexError:
                time.sleep(60)
                taxid = accession2taxid(accession)
                taxonomy = ncbi_taxid_request(taxid)
                sublist.append([taxonomy, i[1]])

        sorted_hits_dict[OTU] = sublist


    print('{}: Filtering hits by taxonomy filter.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## reduce multiple hits
    hit_dict_1 = {}
    for query, hits in tqdm(sorted_hits_dict.items()):

        ## store the similarity
        similarity = hits[0][1]

        ## collect the taxonomy from the dict
        taxonomy_list = [i[0] for i in hits]

        ## remove duplicates
        taxonomy_list_set = [list(x) for x in set(tuple(x) for x in taxonomy_list)]

        filtered_list = []
        for taxonomy in taxonomy_list_set:
            # reduce taxonomy according to JAMP
            if 100 >= similarity >= 98:
                if taxonomy != ['']*6:
                    filtered_list.append(taxonomy)
            elif 98 >= similarity >= 95:
                taxonomy[5] = ''
                if taxonomy != ['']*6:
                    filtered_list.append(taxonomy)
            elif 95 >= similarity >= 90:
                taxonomy[4], taxonomy[5] = '', ''
                if taxonomy != ['']*6:
                    filtered_list.append(taxonomy)
            elif 90 >= similarity >= 85:
                taxonomy[3], taxonomy[4], taxonomy[5] = '', '', ''
                if taxonomy[1:7] != ['']*6:
                    filtered_list.append(taxonomy)
            else:
                taxonomy[2], taxonomy[3], taxonomy[4], taxonomy[5] = '', '', '', ''
                if taxonomy != ['']*6:
                    filtered_list.append(taxonomy)
            ## apend the filtered taxonomy to the main table
            hit_dict_1[query] = filtered_list + [similarity]

    print('{}: Downgrading conflicting hits.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## downgrading of remaining multiple hits
    hit_dict_2 = {}
    for query, taxonomy in tqdm(hit_dict_1.items()):
        ## check if there are multiple hits remaining
        if len(taxonomy) == 2:
            ## remove species entries that were identified to Genus
            if 'sp.' in taxonomy[0][5]:
                taxonomy[0][5] = ''
                hit_dict_2[query] = taxonomy[0] + [taxonomy[1]]
            else:
                hit_dict_2[query] = taxonomy[0] + [taxonomy[1]]
        else:
            ## collect information
            taxonomy_list = taxonomy[:-1]
            identity = taxonomy[-1]

            # remove: species
            for entry in taxonomy_list:
                entry[5] = ''
            if len([list(x) for x in set(tuple(x) for x in taxonomy_list)]) == 1:
                taxonomy_list = [list(x) for x in set(tuple(x) for x in taxonomy_list)][0]
            else:
                # remove: genus
                for entry in taxonomy_list:
                    entry[4] = ''
                if len([list(x) for x in set(tuple(x) for x in taxonomy_list)]) == 1:
                    taxonomy_list = [list(x) for x in set(tuple(x) for x in taxonomy_list)][0]
                else:
                    # remove: family
                    for entry in taxonomy_list:
                        entry[3] = ''
                    if len([list(x) for x in set(tuple(x) for x in taxonomy_list)]) == 1:
                        taxonomy_list = [list(x) for x in set(tuple(x) for x in taxonomy_list)][0]
                    else:
                        # remove: order
                        for entry in taxonomy_list:
                            entry[2] = ''
                        if len([list(x) for x in set(tuple(x) for x in taxonomy_list)]) == 1:
                            taxonomy_list = [list(x) for x in set(tuple(x) for x in taxonomy_list)][0]
                        else:
                            # remove: class
                            for entry in taxonomy_list:
                                entry[1] = ''
                                taxonomy_list = entry

            ## append the back-ranked hit to the dict
            hit_dict_2[query] = taxonomy_list + [identity]

    print('{}: Collecting OTUs/ESVs without a hit.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## check which OTUs did not produce a blast hit
    read_table_df = pd.read_excel(read_table)
    read_table_otus = read_table_df['ID'].values.tolist()

    blast_df_out_list = []
    for OTU in tqdm(read_table_otus):
        if OTU in hit_dict_2.keys():
            blast_df_out_list.append([OTU] + hit_dict_2[OTU] + ['diat.barcode'])
        else:
            blast_df_out_list.append([OTU] + ['No Match']*7 + ['diat.barcode'])

    print('{}: Writing taxonomy table.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## create the final dataframe
    blast_df_out = pd.DataFrame(blast_df_out_list)
    blast_df_out.columns = ['ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'Status']
    blast_df_out.to_excel(taxonomy_table_xlsx, index=False, sheet_name='Taxonomy table')

    print('{}: Finished filtering BLAST results.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

################################################################################
## result filtering

def filter_blast_results(blast_csv, read_table, project_folder, blastn_database, filter_tresholds, database_type, threshold_e_value):

    ## collect files and folders
    blast_csv = Path(blast_csv)
    read_table = Path(read_table)
    project_folder = Path(project_folder)

    ## replace spaces and dots in the filename
    filename = Path(blast_csv).stem

    ## create an output file
    taxonomy_table_xlsx = blast_csv.parent.joinpath(filename + '_taxonomy_table.xlsx')

    print('{}: Starting to filter BLAST results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), filename))

    ## load blast results
    col_names = ['ID', 'Sequence ID', 'Similarity', 'evalue']
    blast_df = pd.read_csv(blast_csv, header=None, sep=';;', names=col_names, engine='python')
    blast_df['Similarity'] = [float(i) for i in blast_df['Similarity'].values.tolist()]

    #############################################################################
    ## LOAD TAXIDS
    #############################################################################

    ## taxids
    ## collect stored taxids
    taxids_xlsx = Path(blastn_database).joinpath('taxids.xlsx')
    print(taxids_xlsx)

    ## print a message
    if database_type == 'Midori2' or database_type == 'Custom NCBI':
        print('{}: Converting taxids to taxonomy.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        print('{}: This can take a while if taxids need to be downloaded from NCBI.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        print('{}: Taxids will be stored to improve runtimes for future datasets.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

        if taxids_xlsx.is_file() == True:
            taxid_df = pd.read_excel(taxids_xlsx, index_col=0)
            print('{}: Found {} available taxids.'.format(datetime.datetime.now().strftime('%H:%M:%S'), len(taxid_df)))

            ## retry no matches?
            indexNames = taxid_df.loc[taxid_df[0] == 'No Match'].index.tolist()
            n_no_matches = len(indexNames)
            if n_no_matches != 0:
                answer = sg.PopupYesNo('Found {} taxids without a matching taxonomy!\nDo you want to retry those?'.format(n_no_matches), title='No Matches')
                if answer == 'Yes':
                    taxid_df.drop(indexNames , inplace=True)

    #############################################################################
    ## PRE FILTER HITS
    #############################################################################

    #############################################################################
    ## filter by evalue
    print('{}: Filtering hits by e-value.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    sorted_hits_dict = {}
    for query in tqdm(set(blast_df['ID'].values.tolist())):

        ## extract all hits for the respective OTU
        df = blast_df.loc[blast_df['ID'] == query][['Sequence ID','Similarity', 'evalue']]
        best_evalue = df['evalue'].values.tolist()[0]

        if best_evalue <= threshold_e_value:
            ## collect all hits that have the same similarity (beginning from the best hit)
            sorted_hits_dict[query] = df.loc[df['evalue'] == best_evalue].values.tolist()

    #############################################################################
    ## DATABASE SPECIFIC PROCESSING
    #############################################################################

    #############################################################################
    ## assign taxonomy to hits
    print('{}: Assigning taxonomy to hits.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## process hits based on the database
    if database_type == 'Midori2':
        ## replace IDs with taxonomy
        ## reduce server connections by storing taxids
        ## write taxids to table for faster processing

        if taxids_xlsx.is_file() == True:
            taxid_df['INDEX'] = taxid_df.index.tolist()
            taxid_storage = {}
            for i in taxid_df.values.tolist():
                taxid_storage[str(i[-1])] = i[:-1]
        else:
            taxid_storage = {}

        # blast_df['Sequence ID'] = [i.split('_')[-1] for i in blast_df['Sequence ID'].values.tolist()]

        for key, values in tqdm(sorted_hits_dict.items()):
            for value in values:

                taxid = value[0].split('_')[-1]

                if taxid not in taxid_storage.keys():
                    res = ncbi_taxid_request(taxid)
                    if len(res) == 6:
                        value[0] = res
                        taxid_storage[taxid] = res
                    else:
                        n = 6 - len(res)
                        res + [''] * n
                        value[0] = res
                        taxid_storage[taxid] = res
                else:
                    value[0] = taxid_storage[taxid]

        taxid_df = pd.DataFrame.from_dict(taxid_storage, orient='index')
        taxid_df.to_excel(taxids_xlsx)

    elif database_type == 'Diat.barcode':
        ## load diat.barcode xlsx
        diat_barcode_xlsx = Path(blastn_database).joinpath('db.xlsx')
        diat_barcode_df = pd.read_excel(diat_barcode_xlsx, sheet_name='data')

        ## collect taxonomy from database file
        df = diat_barcode_df[['Sequence ID', 'Phylum (following Algaebase 2018)', 'Class (following Round, Crawford & Mann 1990)', 'Order (following Round, Crawford & Mann 1990)', 'Family (following Round, Crawford & Mann 1990)', 'Genus', 'Species']]
        reference_taxonomy_dict = {}
        for i in df.values.tolist():
            reference_taxonomy_dict[i[0]] = i[1:]

        ## replace IDs with taxonomy
        for key, values in sorted_hits_dict.items():
            for value in values:
                value[0] = reference_taxonomy_dict[value[0]]

    elif database_type == 'Custom NCBI':
        ## replace IDs with taxonomy
        ## reduce server connections by storing taxids
        ## write taxids to table for faster processing

        if taxids_xlsx.is_file() == True:
            taxid_df['INDEX'] = taxid_df.index.tolist()
            taxid_storage = {}
            for i in taxid_df.values.tolist():
                taxid = int(i[-1])
                taxonomy = i[:-1]
                taxid_storage[taxid] = taxonomy
        else:
            taxid_storage = {}

        for key, values in tqdm(sorted_hits_dict.items()):
            for value in values:

                taxid = int(value[0].split('_')[-1])

                if taxid not in taxid_storage.keys():
                    res = ncbi_taxid_request(taxid)
                    if len(res) == 6:
                        value[0] = res
                        taxid_storage[taxid] = res
                    else:
                        n = 6 - len(res)
                        res + [''] * n
                        value[0] = res
                        taxid_storage[taxid] = res

        taxid_df = pd.DataFrame.from_dict(taxid_storage, orient='index')
        taxid_df.to_excel(taxids_xlsx)


    #############################################################################
    ## FILTER HITS BY TAXONOMY
    #############################################################################

    #############################################################################
    ## filter by taxonomy and  reduce multiple hits
    print('{}: Filtering hits by taxonomy.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    hit_dict_1 = {}

    for query, hits in tqdm(sorted_hits_dict.items()):
        ## store the similarity
        try:
            similarity = float(hits[0][1])
        except:
            similarity = 0

        ## collect taxonomies
        taxonomy_list = [i[0] for i in hits]
        taxonomy_list_set = [list(x) for x in set(tuple(x) for x in taxonomy_list)]

        filtered_list = []
        for taxonomy in taxonomy_list_set:
            # reduce taxonomy according to JAMP
            if 100 >= similarity >= filter_tresholds[0]:
                if taxonomy != ['']*6:
                    filtered_list.append(taxonomy)
            elif filter_tresholds[0] >= similarity >= filter_tresholds[1]:
                taxonomy[5] = ''
                if taxonomy != ['']*6:
                    filtered_list.append(taxonomy)
            elif filter_tresholds[1] >= similarity >= filter_tresholds[2]:
                taxonomy[4], taxonomy[5] = '', ''
                if taxonomy != ['']*6:
                    filtered_list.append(taxonomy)
            elif filter_tresholds[2] >= similarity >= filter_tresholds[3]:
                taxonomy[3], taxonomy[4], taxonomy[5] = '', '', ''
                if taxonomy[1:7] != ['']*6:
                    filtered_list.append(taxonomy)
            else:
                taxonomy[2], taxonomy[3], taxonomy[4], taxonomy[5] = '', '', '', ''
                if taxonomy != ['']*6:
                    filtered_list.append(taxonomy)
            ## apend the filtered taxonomy to the main table
            hit_dict_1[query] = filtered_list + [similarity]

    #############################################################################
    ## FILTER HITS BY TAXONOMY
    #############################################################################

    #############################################################################
    ## Downgrade conflicting hits
    print('{}: Downgrading conflicting hits.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## downgrading of remaining multiple hits
    hit_dict_2 = {}
    for query, taxonomy in tqdm(hit_dict_1.items()):
        ## check if there are multiple hits remaining
        if len(taxonomy) == 2:
            taxonomy[0][5] = ' '.join(taxonomy[0][5].split('_'))
            ## remove species entries that were identified to Genus
            if 'sp.' in taxonomy[0][5]:
                taxonomy[0][5] = ''
                hit_dict_2[query] = taxonomy[0] + [taxonomy[1]]
            else:
                hit_dict_2[query] = taxonomy[0] + [taxonomy[1]]
        else:
            ## collect information
            taxonomy_list = taxonomy[:-1]
            identity = taxonomy[-1]

            # remove: species
            for entry in taxonomy_list:
                entry[5] = ''
            if len([list(x) for x in set(tuple(x) for x in taxonomy_list)]) == 1:
                taxonomy_list = [list(x) for x in set(tuple(x) for x in taxonomy_list)][0]
            else:
                # remove: genus
                for entry in taxonomy_list:
                    entry[4] = ''
                if len([list(x) for x in set(tuple(x) for x in taxonomy_list)]) == 1:
                    taxonomy_list = [list(x) for x in set(tuple(x) for x in taxonomy_list)][0]
                else:
                    # remove: family
                    for entry in taxonomy_list:
                        entry[3] = ''
                    if len([list(x) for x in set(tuple(x) for x in taxonomy_list)]) == 1:
                        taxonomy_list = [list(x) for x in set(tuple(x) for x in taxonomy_list)][0]
                    else:
                        # remove: order
                        for entry in taxonomy_list:
                            entry[2] = ''
                        if len([list(x) for x in set(tuple(x) for x in taxonomy_list)]) == 1:
                            taxonomy_list = [list(x) for x in set(tuple(x) for x in taxonomy_list)][0]
                        else:
                            # remove: class
                            for entry in taxonomy_list:
                                entry[1] = ''
                                taxonomy_list = entry

            ## append the back-ranked hit to the dict
            hit_dict_2[query] = taxonomy_list + [identity]

    #############################################################################
    ## COLLECT NO MATCHES
    #############################################################################

    #############################################################################
    ## Collect remaining OTUs
    print('{}: Collecting OTUs/ESVs without a hit.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## check which OTUs did not produce a blast hit
    read_table_df = pd.read_excel(read_table)
    read_table_otus = read_table_df['ID'].values.tolist()

    blast_df_out_list = []
    for OTU in tqdm(read_table_otus):
        if OTU in hit_dict_2.keys():
            blast_df_out_list.append([OTU] + hit_dict_2[OTU] + [database_type])
        else:
            blast_df_out_list.append([OTU] + ['No Match']*7 + [database_type])

    #############################################################################
    ## WRITE OTU TABLE
    #############################################################################

    #############################################################################
    ## write OTU table
    print('{}: Writing taxonomy table.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## create the final dataframe
    blast_df_out = pd.DataFrame(blast_df_out_list)
    blast_df_out.columns = ['ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'Status']
    blast_df_out.to_excel(taxonomy_table_xlsx, index=False, sheet_name='Taxonomy table')

    print('{}: Finished filtering BLAST results.'.format(datetime.datetime.now().strftime('%H:%M:%S')))











#
