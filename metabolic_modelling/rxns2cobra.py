# Define function to read in result from BLAST pipeline as model

def parse_genes(s):
    return ast.literal_eval(s.replace('[', '[\'').replace(']', '\']').replace(',', '\',\''))

# parses genes from a list of string lines as in the .rxns.csv
def parse_gpr(ls):
    
    # Skip to the headers
    for i, line in enumerate(ls):
        if 'Query-proteins' in line:
            break

    # Parse reasoning as csv to get a df
    output = io.StringIO('\n'.join([line.replace('""', '"') for line in ls[i:]]))
    reasoning = pd.read_csv(output)
    reasoning.columns = [c.strip() for c in reasoning.columns]

    # Read every entry in reasoning['Query-proteins'] as a list, then flatten the list and build GPR
    genes = [gene for sub_list in reasoning["Query-proteins"].apply(parse_genes).to_list() for gene in sub_list]
    return genes

# This function takes a long time, the majority is spent adding (any) gene_reaction_rule to reactions for some reason
def build_cobra(acc, tool='diamond'):

    model = cobra.Model(f"{acc}_BLAST_KEGG")
    df = pd.read_csv(f'{wd}{tool}_models/{acc}.rxns.csv')

    rxns = []
    rids = set()
    incomplete_rxn = []
    for idx, row in df.iterrows():

        rid = row['id'] if row['id'] == row['id'] else row['EC']

        # Do not involve incomplete reactions
        if 'The reaction might be incomplete' in row['comment']:
            incomplete_rxn.append(rid)
            continue

        # The same reaction ID can turn up multiple times because of different EC numbers. Concat GPR and names.
        if rid in rids:
            for rxn in rxns:
                if rxn.id == rid:
                    rxn.name += ', ' + row['name']
                    genes = parse_gpr(row['comment'].split('\n'))
                    if genes:
                        genes_prev = rxn.gene_reaction_rule.split(' or ')
                        rxn.gene_reaction_rule = f'( {" or ".join(set(genes_prev + genes))} )'
                    break
            continue

        # Build reaction
        rxn = cobra.Reaction(rid)
        rxn.name = row['name']
        curr_mets = ast.literal_eval(row['stoichiometry'])
        rxn.add_metabolites({cobra.Metabolite(m, compartment='s'): curr_mets[m] for m in curr_mets})
        genes = parse_gpr(row['comment'].split('\n'))
        rxn.gene_reaction_rule = f'( {" or ".join(genes)} )'

        rxns.append(rxn)
        rids.add(rid)

    model.add_reactions(rxns)
    return model, incomplete_rxn
