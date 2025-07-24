__author__ = "Emilian Paulitz"
__license__ = "MIT"

import argparse
import pandas as pd
from os.path import exists

def main(args):
    wd = args.wd
    out = args.out
    if exists(f'{out}/rxn_translation.csv'):
        if not args.override:
            exit('File to write into exists. Please specify --override to override it or rename the existing file')

    # read table and drop unnecessary first row
    df = pd.read_csv(f'{wd}/reac_xref.tsv', sep = '\t', comment='#', header=None)
    df.columns = 'source ID      description'.split()
    df = df.loc[df['ID'] != 'EMPTY', :]

    # infer the type of identifier and extract only the raw identifier (seed.reaction:R12345 -> seed, R12345)
    tmp_d = {'MNXR01':'mnx', 'MNXR02': 'mnx', 'MNXR03':'mnx', 'bigg.reaction':'bigg', 'biggR':'bigg', 'kegg.reaction':'kegg', 'keggR':'kegg',
     'metacyc.reaction':'metacyc', 'metacycR':'metacyc', 'mnx':'mnx', 'rhea':'rhea', 'rheaR':'rhea', 'sabiork.reaction':'sabiork',
     'sabiorkR':'sabiork', 'seed.reaction':'seed', 'seedR':'seed'}
    df['db'] = df['source'].apply(lambda x: tmp_d[x.split(':')[0]])
    df['xrefid'] = df['source'].apply(lambda x: x.split(':')[1] if ':' in x else x)

    # drop unnecessary data
    df = df.loc[df['db'] != 'mnx', :]
    df = df.drop('source', axis=1).drop_duplicates()

    # read EC number information
    ecs = pd.read_csv(f'{wd}/reac_prop.tsv', sep = '\t', comment='#', header=None)
    ecs.columns='ID     mnx_equation    reference       ec        is_balanced     is_transport'.split()
    ecs = ecs.loc[~ecs['ec'].isna(), ]

    # merge with df and write to file
    df = df.merge(ecs[['ID', 'ec', 'mnx_equation']], on='ID', how='left')
    df.to_csv(f'{out}/rxn_translation.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("wd", help="Path to folder where reac_*.tsv files are stored. Required positional argument.")
    parser.add_argument("out", help="Path to output folder. Will generate the translation table rxn_translation_new.csv there. Required positional argument.")

    # Optional argument flag which defaults to False
    parser.add_argument("-O", "--override", action="store_true", default=False, help="Override existing files without asking.")

    args = parser.parse_args()
    main(args)
