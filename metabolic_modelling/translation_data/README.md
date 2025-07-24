Tables translating metabolites and reactions of different widely used namespaces like KEGG, Bigg, MetaCyc can be downloaded from https://www.metanetx.org/mnxdoc/mnxref.html (e.g. `reac_xref.tsv` and `reac_prop.tsv` for reactions/enzymes)
The script `format_mnx_tables.py` can be used to format the files `reac_prop.tsv` and `reac_xref.tsv` into a translation table. It requires pandas to run.
Also, Philipp Wendering has the data in a nice format in his COMMIT project: https://github.com/pwendering/COMMIT/tree/master/data/tables/MNXref

The following tables were not uploaded due to copyright, but can be obtained from the databases directly (or by asking Emil):

- metacyc2ec.csv can be obtained from PGDB local install. The file protein-seq-ids-reduced-70.dat is in the pathway-tools folder and can be parsed with a custom python script. Ask Emil about it.

**KEGG API is a great resource for translations**
For documentation, see https://www.kegg.jp/kegg/rest/keggapi.html. Especially the link and list command are very useful:
- ec2pw.tsv, mapping EC numbers to KEGG pathway identifiers can be downloaded from https://rest.kegg.jp/link/pathway/ec
- kegg2ec.tsv can be obtained from https://rest.kegg.jp/link/enzyme/reaction
- pw2names.tsv to obtain human-readable names for pathways, e.g. for plotting, can be obtained from https://rest.kegg.jp/list/pathway

- ec2names.dict.txt maps EC numbers to human understandable names. It was compiled from various sources:
    - https://www.sigmaaldrich.com/DE/de/technical-documents/technical-article/protein-biology/enzyme-activity-assays/enzyme-commission-numbers
    - https://www-archiv.fdm.uni-hamburg.de/b-online/e18_1/ec4.htm
    - KEGG API requests (but probably https://rest.kegg.jp/list/enzyme would have been easier)
    It is a text-file that can be read into a python dictionary with
```
import ast
with open('ec2names.dict.txt', 'r') as infile:
    ec2names = infile.read()
ec2names = ast.literal_eval(ec2names)
```
