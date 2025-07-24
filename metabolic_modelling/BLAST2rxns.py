__author__ = "Emilian Paulitz"
__license__ = "MIT"

import argparse
from unidecode import unidecode
import ast
import pandas as pd
import numpy as np
import requests
import io
from time import sleep
import datetime

"""
Design choices:
Exclude symbols with less than three letters, as they lead to spurious search results (1)

UniProt:
First three letters of symbol have to be contained in the line (because sometimes symbol NAC001 corresponds to uniprot string NAC1) (2)
Only accept reviewed entries (3)
Only take first entry of a search (4)

KEGG:
Exclude general EC-numbers (e.g. 1.1.1.-) (5)
KEGG reactions after this procedure with either no substrates or no products are discarded. (8)
Only parse the net amount of a metabolite in a reaction. 2 H2O + ... -> 1 H2O + ... would be parsed as H2O + ... -> ... (10)
Reaction stoichiometries including n as a factor are unique per case and need to be dealt with manually. Currently, unknown factors are always parsed as 1.
If reactions contain the same compound on both reaction sides, print out the reaction and have the compound be consumed (11)

fixed by always using the KEGG reactions denoted in EC API response:
Sometimes, multiple reaction IDs of the same enzyme are mentioned. Use all (9) (because e.g. for 2.1.1.4 both reactions are mentioned under https://www.genome.jp/kegg-bin/get_htext#D2008)
metabolites without KEGG ID in KEGG reactions are not parsed. These are often times either polymers with n, or protein-residues. (6)
    Their stoichiometry might be incomplete, and metabolites with spaces often times lead to incomplete stoichiometries. (7)
    We try to fix by looking up reaction ID, if not possible, it is simply denoted in comment (7.1)
"""

def dissect_KEGG_unused(s, cpds):

    ls = s.split('\n')

    # extract substrates/products from response
    prods, subs = [], []
    sub, prod = False, False
    compounds_incomplete = False
    rxn_name = ''
    rxn_string = ''
    name_to_cpd = dict()
    S = dict()
    for line in ls:
        if not line.startswith(' '):
            sub, prod = False, False
            if line.startswith('SUBSTRATE'):
                sub = True
                line = line[len('SUBSTRATE'):]
            elif line.startswith('PRODUCT'):
                prod = True
                line = line[len('PRODUCT'):]
            elif line.startswith('NAME'):
                if ';' in line:
                    rxn_name = line[len('NAME'):line.rfind(';')]
                else:
                    rxn_name = line[len('NAME'):]
            elif line.startswith('REACTION'):
                rxn_id = ''
                if line.rfind('[RN:') >= 0 and line.rfind(']') >= 0:
                    rxn_id = line[line.rfind('[RN:') + len('[RN:'): line.rfind(']')]
                    if " " in rxn_id:  # (9)
                        rxn_id = rxn_id.split(" ")[0]
                rxn_string = line
        if sub:
            try:
                name, cpd = line.strip().split(' [CPD:')
                if ' ' in name:
                    compounds_incomplete = True  # (7)
            except ValueError:
                compounds_incomplete = True  # (7)
                continue  # (6)
            cpd = cpd[:cpd.rfind(']')]
            name_to_cpd[name] = cpd
            if cpd not in cpds:
                cpds[cpd] = name
            subs.append(cpd)
            
        elif prod:
            try:
                name, cpd = line.strip().split(' [CPD:')
                if ' ' in name:
                    compounds_incomplete = True  # (7)
            except ValueError:
                compounds_incomplete = True  # (7)
                continue  # (6)
            cpd = cpd[:cpd.rfind(']')]
            name_to_cpd[name] = cpd
            if cpd not in cpds:
                cpds[cpd] = name
            prods.append(cpd)
    if subs and prods:  # (8)
        # Parse stoichiometry
        ls = rxn_string.split()
        factor = 1
        for word in ls[1:]:
            if word.isnumeric():
                factor = int(word)
            elif word in name_to_cpd:
                if name_to_cpd[word] in subs:
                    S[name_to_cpd[word]] = S.get(name_to_cpd[word], 0) + factor  # (10)
                elif name_to_cpd[word] in prods:
                    S[name_to_cpd[word]] = S.get(name_to_cpd[word], 0) - factor  # (10)
                factor = 1
            # Other cases could be +, =, [RN:R00006]

        return rxn_name, rxn_id, S, compounds_incomplete
    
    return [None] * 4

def is_KEGG_rid(rid):
    if len(rid) != 6:
        return False
    if rid[0] != 'R':
        return False
    for i in range(1, len(rid)):
        if not rid[i].isnumeric():
            return False
    return True

# parses response text from a request for an EC number to retrieve all rxn IDs
def get_all_reac(ec_text):
    res = set()
    parse_line = False
    for line in ec_text.split('\n'):
        curr_rids = None
        if line.startswith('ALL_REAC'):
            curr_rids = line[len('ALL_REAC'):].strip().strip(';').split()
        elif parse_line:
            parse_line = False
            curr_rids = line.strip().strip(';').split()

        if curr_rids is not None:
            for rid in curr_rids:
                if rid.endswith('(G)') and is_KEGG_rid(rid[:-len('(G)')]):
                    rid = rid[:-len('(G)')]

                if is_KEGG_rid(rid):
                    res.add(rid)
                elif rid not in ['(other)', '>']:
                    print(f'Found unexpected string {rid} in ALL_REAC field')

            if line.strip().endswith(';'): # there is another line of reaction IDs
                parse_line = True
            else:
                break
    return res

# performs API request to retrieve reaction information from reaction ID
def get_rxn(rid, prev_rxn_name):
    success_fix = False
    rxn_name = ''
    new_entry = dict()
    if rid:
        print(f'Attempting {rid}...', end=' ')

        perform_request = True
        retries = 0
        while perform_request:
            if retries > 5:
                print(f'Got error 403 for {retries} attempts for {rid}')
                break
            if retries > 0:
                sleep(31)

            retries += 1
            response = requests.get(f'https://rest.kegg.jp/get/{rid}')
            perform_request = response.status_code == 403               

        if response.status_code == 200:
            for line in response.text.split('\n'):
                if line.startswith('NAME'):
                    rxn_name = line[len('NAME'):].strip()
                if line.startswith('EQUATION'):
                    S = dict()
                    ls = line[len('EQUATION'):].strip().split()
                    factor = 1
                    subs_or_prods = -1  # start with substrates
                    for ele in ls:
                        if ele.isnumeric():
                            factor = int(ele)
                        elif '=' in ele:
                            subs_or_prods = 1  # now we are in prods
                        elif ele.startswith('C') or ele.startswith('G'):
                            # (11)
                            if ele not in S:
                                S[ele] = factor * subs_or_prods
                            else:
                                print(f'Error for reaction {rid}: The same compound {ele} occured multiple times.')
                            factor = 1
                        elif ele != '+':
                            print(f'Unexpected word encountered in equation of KEGG API request for reaction {rid}: {ele}')
                    success_fix = True
                    break

        print('Success!' if success_fix else 'Failed!')

    return rxn_name.strip(';') if rxn_name else prev_rxn_name, rid, S, not success_fix


def build_model(ec_numbers, db_dir, working_dir, accession, reasoning = None):
    if reasoning is None:
        reasoning = {ec: '' for ec in ec_numbers}

    # Retrieve reaction from EC-number from KEGG
    rxns = []
    cpds = dict()

    # read ec2rxn db
    try:
        ec2rxn = pd.read_csv(db_dir + 'ec2rxn.csv').set_index('EC')
        ec2rxn['comment'] = ec2rxn['comment'].apply(lambda x: x if x == x else '')
    except FileNotFoundError:
        ec2rxn = pd.DataFrame()

    # add ec number data to rxns for all ecs in the db, except if it is an invalid reaction
    def add_row_if_valid(row):  # axis keyword is ignored; purpose: the function can be called via Series.apply and DatatFrame.apply
        
        if (not (type(row['stoichiometry']) == float and np.isnan(row['stoichiometry']))
            and row['stoichiometry'] != '{}'):  # empty dictionary
            
            rxns.append([row.name,  # EC number
                        row['name'],
                        row['id'],
                        row['stoichiometry'],
                        row['comment'] + reasoning[row.name]])
        return None
    
    ec_done = ec2rxn.index.intersection(ec_numbers)
    ec2rxn.loc[ec_done].apply(add_row_if_valid, axis=1)
    ec_done = set(ec_done)
    print(f'Found {len(rxns)} reactions for {len(ec_numbers)} level-4 enzyme numbers in the database.')

    # proceed with API requests if ECs were not in the db
    # You can submit up to 10 ec numbers at a time
    ec_numbers = ec_numbers.difference(ec2rxn.index)
    print(f'Performing API requests for remaining {len(ec_numbers)} enzyme numbers')
    
    # some EC numbers might be obsolete, in which case their new number needs to be searched
    try:
        obsolete2ec = pd.read_csv(f'{db_dir}obsolete2ec.csv').set_index('obsolete_ec')['ec'].to_dict()
        new_ec_numbers = set()
        for ec in ec_numbers:
            if ec not in obsolete2ec:
                new_ec_numbers.add(ec)
            elif obsolete2ec[ec] != '-':
                new_ec_numbers.add(obsolete2ec[ec])
                reasoning[obsolete2ec[ec]] = f'new name of obsolete EC {ec}\n' + reasoning[ec]
        ec_numbers = new_ec_numbers
    except:
        obsolete2ec = dict()
    num_new_obsolete = 0
    ec2rxns_to_add = list()
    rids_to_request = dict()  # rid: ec
    while ec_numbers:

        batch = []
        while len(batch) < 10 and ec_numbers:
            curr_ec = ec_numbers.pop()
            if curr_ec in ec_done:
                continue

            # check whether information for ec number is already on disk; else add to API query
            if curr_ec in ec2rxn.index:
                # ec2rxn.loc[curr_ec] can be either a Series or a DataFrame
                tmp = ec2rxn.loc[curr_ec]
                if len(tmp.shape) == 1:
                    add_row_if_valid(tmp)
                else:
                    tmp.apply(add_row_if_valid, axis=1)
            else:
                batch.append(curr_ec)
                ec_done.add(curr_ec)
        
        # first get reaction IDs from EC numbers
        retries = 0
        perform_request = True
        while perform_request:
            retries += 1
            response = requests.get(f'https://rest.kegg.jp/get/{"+".join(batch)}')
            perform_request = response.status_code == 403
            if perform_request:
                if retries > 5:
                    print(f'Got error 403 for {retries} attempts for {batch}')
                    break
                sleep(31)

        if response.status_code != 200:
            print('error:', response.status_code, batch)
            continue

        responses = response.text.split('///')
        for i, s in enumerate(responses):
            # find out EC number of response text. necessary because responses of unknown 
            # EC numbers are empty and do not even produce a separator symbol, therefore 
            # not increasing i. Still use batch mode because it is faster (~2x)
            found_entry = False
            obsolete = False
            for line in s.split('\n'):
                if obsolete and line.startswith('NAME'):
                    ls = line.split()
                    if len(ls) <= 3 or ls[0] != 'NAME' or ls[1] != 'Transferred' or ls[2] != 'to':
                        if ls == ['NAME', 'Deleted', 'entry']:
                            print(f'Warning: Got deleted obsolete enzyme {obsolete_ec}')
                            obsolete2ec[obsolete_ec] = '-'
                            num_new_obsolete += 1
                        else:
                            print('Warning: Obsolete enzyme but NAME line was unexpected:')
                            print(line) 
                    else:
                        ec_numbers.add(ls[3])
                        obsolete2ec[obsolete_ec] = ls[3]
                        num_new_obsolete += 1
                        reasoning[ls[3]] = f'new name of obsolete EC {obsolete_ec}\n' + reasoning[obsolete_ec]
                        print(f'Found obsolete enzyme {obsolete_ec} and added new number {ls[3]} to query')
                    found_entry = True
                    break

                if line.startswith('ENTRY'):
                    ls = line.split()
                    if ls[0] != 'ENTRY' or ls[1] != 'EC' or ls[3] != 'Enzyme':
                        if ls[3] == 'Obsolete' and ls[4] == 'Enzyme':
                            obsolete_ec = ls[2]
                            obsolete = True
                            continue
                        else:
                            print('Warning: ENTRY line was unexpected:')
                            print(line)
                            found_entry = False
                            break
                    else:
                        ec = ls[2]
                        found_entry = True
                    break
            
            if not found_entry:
                # the last element is usually empty, so not problem
                if i + 1 != len(responses):
                    print('Warning: ENTRY line was not found or unexpected:')
                    print(s)
                    print('Aborted enzyme')
                continue
            if obsolete:
                continue

            curr_rids = get_all_reac(s)

            if not curr_rids:
                print(f'Warning: Did not find ALL_REAC for EC {ec}')
                rxn_to_add = dict()
                rxn_to_add['EC'] = ec
                rxn_to_add['name'] = ''
                rxn_to_add['id'] = ''
                rxn_to_add['stoichiometry'] = ''
                rxn_to_add['comment'] = 'No proper reaction\n'
                ec2rxns_to_add.append(rxn_to_add)

            for rid in curr_rids:
                rids_to_request[rid] = rids_to_request.get(rid, set()).union([ec])

    if num_new_obsolete:
        print(f'Found {num_new_obsolete} new obsolete EC numbers, adding them to obsolete2ec.csv')
        df = pd.DataFrame(obsolete2ec, index=['ec']).T
        df.index.name = 'obsolete_ec'
        df.to_csv(f'{db_dir}obsolete2ec.csv')

    # Then call get_rxn to obtain the following variables
    print(f'Found {len(rids_to_request)} reactions to request')
    for rid in rids_to_request:
        
        rxn_name, rxn_id, S, compounds_incomplete = get_rxn(rid, rid)

        for ec in rids_to_request[rid]:
            if compounds_incomplete:  # (7.1)
                print(f"Reaction {rid} for EC {ec} is likely incomplete.")
                comment = "The reaction might be incomplete\n" + reasoning[ec]
            else:
                comment = reasoning[ec]

            rxns.append([ec, rxn_name, rxn_id, str(S), comment])
            
            rxn_to_add = dict()
            rxn_to_add['EC'] = ec
            rxn_to_add['name'] = rxn_name.strip() if type(rxn_name) == str else rxn_name
            rxn_to_add['id'] = rxn_id.strip() if type(rxn_id) == str else rxn_id
            rxn_to_add['stoichiometry'] = str(S).strip()
            rxn_to_add['comment'] = 'The reaction might be incomplete\n' if comment.startswith('The reaction might be incomplete') else ''
            ec2rxns_to_add.append(rxn_to_add)

    print("Number of rxns retrieved from ec_numbers:", len(rxns))

    # update ec2rxn
    if ec2rxns_to_add:
        df = pd.DataFrame(ec2rxns_to_add)
        df.set_index('EC', inplace=True)
        pd.concat([ec2rxn, df]).to_csv(db_dir + 'ec2rxn.csv')
        print(f'{datetime.datetime.now()}: Updated ec2rxn with {df.shape[0]} entries')

    # create final df
    df = pd.DataFrame(rxns, columns=['EC', 'name', 'id', 'stoichiometry', 'comment']).sort_values('EC')

    # The following merges multiple enzymes with the same reaction ID. I decided that it might be better to keep the long format however
    if False:
        # For the following step every reaction need an ID, so fill in the EC number if ID is missing
        df.loc[df['id'].isna(), 'id'] = df.loc[df['id'].isna(), 'EC']

        # Merging enzymes with the same reaction ID 
        # A function to merge the comments including the GPR information
        def merge_comment(s, sep = '$$$'):

            res = ''
            # For each string, skip to table and append to the result
            for item in s:
                ec, actual_comment = item.split(sep)
                after_header = False
                for line in actual_comment.split('\n'):
                    if after_header:
                        res += f'{line},{ec}\n'
                    if not after_header and 'Query-proteins' in line:
                        after_header = True

            # get header from the last element used
            header = ''
            for line in actual_comment.split('\n'):
                header += line + '\n'
                if 'Query-proteins' in line:
                    # remove trailing newline and add new column
                    header = header[:-1] + ', EC\n'
                    break
            return header + res[:-1]  # remove trailing newline

        # Solve the problem with an involved groupby, that chooses the first element for every column, except for...
        groupby_dict = {c: lambda x: x.iloc[0] for c in df.columns}
        groupby_dict.pop('id')

        # names: merge the names if different (by joining by ', ')
        groupby_dict['name'] = lambda x: ', '.join(set(x))

        # comment: merge the reasoning by: adding a column "EC" and concatenate all of the reasoning tables 
        groupby_dict['comment'] = merge_comment

        # temporarily add the EC number with a separator to the comment to be able to access that information
        df['comment'] = df['EC'] + '$$$' + df['comment']
        df = df.groupby('id').agg(groupby_dict).reset_index()    

    df.to_csv(f'{working_dir}/{accession}.rxns.csv', index=False)
    print('Model reactions written to', f'{working_dir}/{accession}.rxns.csv')


def main(args):
    db_dir = args.db_dir + ('' if args.db_dir.endswith('/') else '/')

    # Process Araport11, gather symbols, and find out corresponding uniprot entries
    if args.input_format == 'araport':
        # Process Araport11 and collect information about headers
        headers = dict()
        with open(f'{args.wd}/Araport11.faa', 'r') as inp:
            with open(f'{args.wd}/Araport11.processed.faa', 'w') as outp:
                for line in inp:
                    # keep newlines around header lines, but remove them for all other lines
                    if line.startswith('>'):
                        replaced = line.replace('&#64257;', 'fi').replace("&#8208;", "-").replace("&#946;", "beta").replace("&#916;", "delta").replace("&#948;", "delta").replace("&#945;", "alpha").replace("&#8242;", "-prime").replace("&#64258;", "fl").replace("&#949;", "epsilon").replace("&#8260;", "/")
                        outp.write('\n' + unidecode(replaced))
                        header = unidecode(replaced[1:-1])
                        headers[header.split(' | ')[0]] = header.split(' | ')[1][len('Symbols: '):].split(', ')
                    else:
                        outp.write(unidecode(line[:-1]))

        # cd ~/Desktop/PhD-Synch/Arabidopsis_panGEM/model_reconstruction/data/blast
        # makeblastdb -in Araport11.processed.faa -dbtype prot
        # blastp -query query.faa -db ./Araport11.processed.faa -evalue 1e-2 -num_threads 8 -mt_mode 1 -max_target_seqs 5 -outfmt "6 qseqid sseqid pident evalue" > result.tsv

        # Extract symbols from blast result
        df = pd.read_csv(args.input_path, sep='\t', header=None)
        df.columns =['query', 'subject', 'pid', 'e']
        df[['pid', 'e']] = df[['pid', 'e']].astype(float)

        # Filter df for valid hits
        valid_hits = df.loc[df['pid'] >= args.BLAST_pid_threshold, ].groupby(by='query').first().reset_index()

        # Extract dictionary with aggregated pid and query entries
        genes = valid_hits.groupby('subject')[['query', 'pid']].agg({'query': list, 'pid': 'max'}).to_dict(orient='index')
        genes.pop('no symbol available', None)

        # Build dict symbol: gene
        symbols = dict()
        for gene in genes:
            for symbol in headers[gene]:
                if len(symbol.replace("#", " ").replace(" ", "")) >= 3:  # (1) filter out ambiguous symbols
                    symbols[symbol.replace("#", " ")] = gene
        print("Number of symbols from BLAST search:", len(symbols))

        # Search symbol(s) on Uniprot
        # TODO integrate writing results into a database
        uniprot_entries = dict()
        for symbol in symbols:

            retries = 0
            perform_request = True
            while perform_request:
                retries += 1
                response = requests.get(f'https://rest.uniprot.org/uniprotkb/search?query="{symbol}"+AND+organism_id:3702&format=tsv')
                perform_request = response.status_code == 403
                if perform_request:
                    if retries > 5:
                        print(f'Got error 403 for {retries} attempts for {symbol}')
                        break
                    sleep(31)

            if response.status_code != 200:
                print('error:', response.status_code, symbol)
                continue
            output = io.BytesIO(response.content)

            # Parse tsv and select entrie(s)
            df = pd.read_csv(output, sep='\t')

            # filter for entries that exactly contain the start of the search string (2)
            contain_start = df[['Entry Name', 'Protein names', 'Gene Names']].apply(lambda x: ''.join(map(str, x.to_numpy())), axis=1).apply(lambda x: symbol[:3] in x)
            # Another filtering idea: df['Entry Name'] == f'{symbol}_ARATH'
            # Filter for reviewed entries, containing the start of the symbol, then take the first result (3,4)
            valid_entries = df.loc[(df['Reviewed'] == 'reviewed') & contain_start, 'Entry'].to_numpy()

            if valid_entries.size > 0:
                if valid_entries[0] in uniprot_entries:
                    uniprot_entries[valid_entries[0]].append(symbol)
                else:
                    uniprot_entries[valid_entries[0]] = [symbol]

        print("Number of uniprot_entries from collected symbols:", len(uniprot_entries))

    # Process input DIAMOND .tsv
    if args.input_format == 'uniprot':
        # parse diamond tsv
        df = pd.read_csv(args.input_path, sep='\t', header=None)
        df.columns =['query', 'subject', 'pid', 'e']
        df[['pid', 'e']] = df[['pid', 'e']].astype(float)

        # extract uniprot ID / gene ID 
        df['subject'] = df['subject'].apply(lambda x: x.split('|')[1])
        df['query'] = df['query'].apply(lambda x: x.split('|')[0])

        # Filter df for valid hits
        valid_hits = df.loc[df['pid'] >= args.BLAST_pid_threshold, ].groupby(by='query').first().reset_index()

        # Extract dictionary with aggregated pid and query entries: uniprotID:{'query':[genes*],'pid':max_pid}
        uniprot_entries = valid_hits.groupby('subject')[['query', 'pid']].agg({'query': list, 'pid': 'max'}).to_dict(orient='index')

        print("Number of uniprot entries from DIAMOND run:", len(uniprot_entries))

    # Map uniprot entry IDs to EC numbers
    if args.input_format in ['araport', 'uniprot']:
        ec_numbers = dict()

        # read upe2ec
        try:
            upe2ec = pd.read_csv(db_dir + 'upe2ec.csv').set_index('UniprotID')['EC'].apply(lambda x: ast.literal_eval(x)).T.to_dict()
        except FileNotFoundError:
            upe2ec = dict()
        added_entries = 0

        for upe in uniprot_entries:
            ecs_to_add = set()
            
            # check if upe is in upe2ec, and add the result if any is not nan
            if upe in upe2ec:
                ecs_to_add = ecs_to_add.union(upe2ec[upe])
            else:  # otherwise perform an API request
                retries = 0
                perform_request = True
                while perform_request:
                    retries += 1
                    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{upe}')
                    perform_request = response.status_code == 403
                    if perform_request:
                        if retries > 5:
                            print(f'Got error 403 for {retries} attempts for {upe}')
                            break
                        sleep(31)
                
                if response.status_code != 200:
                    print('error:', response.status_code, upe)
                    continue
                
                # Search for EC-number in response
                if 'ecNumber' in response.text:
                    try:
                        ecs_to_add.add(response.json()['proteinDescription']['recommendedName']['ecNumbers'][0]['value'])
                    except KeyError:
                        success = False
                        for comment in response.json()['comments']:
                            if comment['commentType'] == "CATALYTIC ACTIVITY":
                                try:
                                    ecs_to_add.add(comment['reaction']['ecNumber'])
                                    success = True
                                except KeyError:
                                    continue
                        if not success:
                            print('Extracting EC number from JSON failed for uniprot entry', upe)
                            
                # And add the result to the table
                added_entries += 1
                upe2ec[upe] = list(ecs_to_add)

            for ec in ecs_to_add:
                if ec in ec_numbers:
                    ec_numbers[ec].append(upe)
                else:
                    ec_numbers[ec] = [upe]
        print("Number of ec_numbers found in uniprot_entries:", len(ec_numbers))

        # write the updated table 
        if added_entries:
            df = pd.DataFrame({upe: str(upe2ec[upe]) for upe in upe2ec}, index=['EC']).T
            df.index.name = 'UniprotID'
            df.to_csv(db_dir + 'upe2ec.csv')
            print(f'Updated upe2ec with {added_entries} new entries')

        # output functional annotation:
        if args.input_format == 'uniprot':
            with open(f'{args.wd}/{args.accession}.functional.csv', 'w') as file:
                file.write('gene,EC,BLAST_pid\n')
                for ec in ec_numbers:
                    for upe in ec_numbers[ec]:
                        file.write(f'"{uniprot_entries[upe]["query"]}",{ec},{uniprot_entries[upe]["pid"]}\n')
                print('Functional annotations file written to', f'{args.wd}/{args.accession}.functional.csv')

        # synthesize the chain of reasoning as a comment
        reasoning = dict()
        for ec in ec_numbers:
            if args.input_format == 'uniprot':
                s = "Chains of reasoning:\nUniProt-entry, Query-proteins, maximum BLAST-pid"
                for upe in ec_numbers[ec]:
                    s += f'\n{upe},"{uniprot_entries[upe]["query"]}",{uniprot_entries[upe]["pid"]}'
                reasoning[ec] = s
            elif args.input_format == 'araport':
                s = "Chains of reasoning:\nUniProt-entry, Araport11-symbol, Araport11-gene, Query-proteins, maximum BLAST-pid"
                for upe in ec_numbers[ec]:
                    for symbol in uniprot_entries[upe]:
                        gene = symbols[symbol]
                        s += f'\n{upe},{symbol},{gene},"[{",".join(genes[gene]["query"])}]",{genes[gene]["pid"]}'
                reasoning[ec] = s

        # These general ec numbers will create trouble for KEGG (5)
        for ec in list(ec_numbers.keys()):  # this is necessary to not change size of iterated object during loop execution
            if '-' in ec:
                ec_numbers.pop(ec)
        ec_numbers = set(ec_numbers.keys())

    # Read deepec .tsv to generate ec_numbers list and reasoning for GPR rules
    if args.input_format == 'deepec':
        # read deepec input table
        df = pd.read_csv(args.input_path, sep='\t', header=0)
        df['Predicted EC number'] = df['Predicted EC number'].apply(lambda s: s[len('EC:'):])

        # Gather list of level-4 EC numbers
        ec_numbers = set([ec for ec in df['Predicted EC number'] if not '-' in ec])

        # synthesize reasoning as a dict {ec: 'csv table containing the header "Query-proteins"'}
        reasoning = df.groupby('Predicted EC number').agg(list)['Query ID'].apply(lambda s: f"Query-proteins\n\"{s}\"").to_dict()

    # Read deepEcTransformer .tsv to generate ec_numbers list and reasoning for GPR rules
    if args.input_format == 'deepectransformer':
        # read deepec input table
        df = pd.read_csv(args.input_path, sep='\t', header=0)
        df = df.loc[~df['prediction'].isna()]
        df['prediction'] = df['prediction'].apply(lambda s: [ec.removeprefix('EC:') for ec in s.split(';')])
        df = df.explode('prediction')

        # Gather list of level-4 EC numbers
        ec_numbers = set([ec for ec in df['prediction'] if not '-' in ec])

        # synthesize reasoning as a dict {ec: 'csv table containing the header "Query-proteins"'}
        reasoning = df.groupby('prediction').agg(list)['sequence_ID'].apply(lambda s: f"Query-proteins\n\"{s}\"").to_dict()

    build_model(ec_numbers, db_dir, args.wd, args.accession, reasoning)


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("accession", help="Name of accession to process. Example: Tibet-0")
    parser.add_argument("input_path", help="Path to tsv to process")
    parser.add_argument("wd", help="Path to working directory. Araport faa is expected "
                        "there (if -f araport), and the output will be stored there.")
    parser.add_argument("db_dir", help="Path to directory where a database of uniprot/kegg "
                        "already exists or should be created, so the same API requests do "
                        "not have to be made multiple times ")

    # Optional argument flag which defaults to False
    parser.add_argument("-f", "--input-format", action="store", default='uniprot',
                        help="Format of input. Either a mapping against Araport11 ('araport'), "
                        "a mapping against UniProtKB (default, 'uniprot'), or a tsv by DeepEC "
                        "already containing EC numbers ('deepec').")
    parser.add_argument("-b", "--blast_pid_threshold", action="store", type=int, default=90, dest="BLAST_pid_threshold",
                         help='BLAST pid threshold to filter results by (default 90)')

    args = parser.parse_args()
    main(args)
