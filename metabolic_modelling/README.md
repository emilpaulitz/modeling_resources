# Resources and scripts for metabolic modelling


In the folder translation_data, data translating between different namespaces, EC numbers, human readable names etc. is collected.

The script **BLAST2rxns.py** by Emilian Paulitz automatically annotates a set of protein sequences with EC numbers. Input:
- The script expects as input a .tsv table obtained from running DIAMOND (or BLAST) with your proteins against Araport11 or Uniprot, or a .tsv file coming from deepec (adapt `--format` in any case).
	- The DIAMOND runs could be done like this (`--outfmt` creates tsv forma as expected by BLAST2rxns.py): `diamond blastp --query query.faa --db ./diamond_db.dmnd --evalue <e-val, e.g. 1e-2> --threads <threads> --max-target-seqs <e.g. 5> --outfmt 6 qseqid sseqid pident evalue --out ./diamond.tsv` 
	- `diamond_db.dmnd` is created by `diamond makedb --in db.faa --db ./diamond_db`
		- `db.faa` can be the protein fasta file of Araport11 or Uniprot database. Remember to adapt `--format` of `BLAST2rxns.py` accordingly
	- Output:
		- A file containing a gene to EC number mapping
		- A file containing information needed to build a model from this mapping, mostly the reaction stoichiometries. This output can be used to build a cobra model using the code in rxns2cobra.py (please adapt paths etc. to your needs)

- The script **rxns2cobra.py** by Emilian Paulitz builds a cobra model from the output of BLAST2rxns.py
- The script **cobrapy_convenience.py** by Emilian Paulitz contains a bunch of convenience functions for working with cobrapy models. Import by copying the file to your script's location and `from cobrapy_convenience import *`. For some functions, certain databases might be needed, please ask Emil about it. Most notable functions include: 
	- `get_rid` and `get_mid` for getting a reaction or metabolite by ID without having to type `model.reactions.get_by_id()` every time
	- `match_mname` and `match_rname` are useful for working with models from other namespaces. It searches for metabolites/reactions by their name.
	- `find_rxn_from_mids` will return all reactions that contain the given metabolites. Useful for working with models in other namespaces and manual curation
	- `check_production` is extremely useful for curating models: check if a specific metabolite can be produced. You can specify metabolites can should be imported for the check to see if specific connections are in the model. It can return either 1) the maximum amount of the metabolite that can be produced, 2) the solutions object to check which reactions are involved in the production, or 3) the whole model with temporaryly added reactions for deeper analysis. It has an option to perform pFBA, so that only relevant reactions will be included in the solution. 
	- `analyze_uptake` is useful for making sure your model does not generate any metabolites in infinite loops by checking which of the import reactions are essential for growth
	- `dissect_bio` checks which of the metabolites in the biomass reaction can be produced and which can not
	- `silence` class: use `with silence():` to suppress all standard output. Reading in older models can produce enormous amounts of warnings and errors that significantly slows read-in and your whole IDE (if you are working in notebooks)
	- `Formula` class: useful for comparing formulas. `Formula(formula_string)` will create a Formula object. These can be added, subtracted, and compared without the need for the same order of elements etc.

