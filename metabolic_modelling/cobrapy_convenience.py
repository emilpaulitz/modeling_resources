import os
from cobra import Reaction
import cobra
import pandas as pd
import ast
import io
import sys
import numpy as np
from itertools import chain, combinations

class silence():
    def __enter__(self):
        self._stdout = sys.stdout
        self._stderr = sys.stderr

        self.text_trap = io.StringIO()
        sys.stdout = self.text_trap
        sys.stderr = self.text_trap
        return self

    def __exit__(self, *args):
        sys.stdout = self._stdout
        sys.stderr = self._stderr

# Holds a string to represent the formula, but:
#   implements equivalence more leniently
#   sorts input on initialization to be able to compare different Formulas
#   also implements addition and subtraction
class Formula:
    @staticmethod
    def _display_formula(elem_composition, sep = ''):

        res = ''
        if 'H' in elem_composition:
            res += 'H' + str(elem_composition['H'])
            res += sep

        for elem in sorted(elem_composition.keys()):
            if elem in ['H', 'R']:
                continue
            res += elem + str(elem_composition[elem])
            res += sep
        
        if 'R' in elem_composition:
            res += 'R' + str(elem_composition['R'])
            res += sep

        return res

    @staticmethod
    def _get_dict(s: str):
        if not s:
            return ''

        curr_elem = ''
        curr_factor = None
        elem_composition = dict()
        for c in s:
            if c.isalpha():
                if c.isupper():
                    if curr_elem:
                        elem_composition[curr_elem] = 1 if curr_factor is None else int(curr_factor)
                        curr_factor = None
                    curr_elem = c
                elif c.islower():
                    curr_elem += c
            elif c.isnumeric():
                if curr_factor is None:
                    curr_factor = c
                else:
                    curr_factor += c

        # add last element
        elem_composition[curr_elem] = 1 if curr_factor is None else int(curr_factor)

        return elem_composition

    # Allows H count to be +-1 equal
    @staticmethod
    def formulae_equal(f1, f2, H_err = 1):
        if f1 == f2:
            return True
        
        # H are in the error bounds
        if abs((f1['H'] if 'H' in f1 else 0) - (f2['H'] if 'H' in f2 else 0)) <= H_err:

            for k1, v1 in f1.items():
                if k1 != 'H' and (k1 not in f2 or f2[k1] != v1):
                    return False
                
            for k2, v2 in f2.items():
                if k2 != 'H' and k2 not in f1:
                    return False
                
            return True
        return False
        
    def __init__(self, f=None):
        if f is None:
            self.d = dict()
        elif type(f) == str:
            self.d = Formula._get_dict(f)
        elif type(f) == dict:
            self.d = f
        else:
            raise TypeError(f'Only dict or str are allowed for constructor, provided {type(f)}')
        self.s = self.__str__()

    @classmethod
    def from_met_id(cls, mid, model):
        return cls(model.metabolites.get_by_id(mid).formula)

    def __eq__(self, other):
        if type(other) == Formula:
            return Formula.formulae_equal(self.d, other.d, 1)
        return False
    
    def equals(self, other, H_err = 1):
        if type(other) == Formula:
            return Formula.formulae_equal(self.d, other.d, H_err)
        return False
    
    def __str__(self):
        return Formula._display_formula(self.d)
    
    def pretty(self):
        return Formula._display_formula(self.d, ' ')
    
    def __repr__(self):
        return self.s
    
    def __hash__(self):
        return hash(self.s)
    
    def __lt__(self, other):
        return self.s > other.s
    
    def __add__(self, other):
        if type(other) != Formula:
            raise TypeError()
        
        res_d = self.d.copy()

        for k in other.d:
            res_d[k] = res_d.get(k, 0) + other.d[k]
        
        return Formula(res_d)
    
    def __sub__(self, other):
        if type(other) != Formula:
            raise TypeError()
        
        res_d = self.d.copy()

        for k in other.d:
            res_d[k] = res_d.get(k, 0) - other.d[k]
        
        return Formula(res_d)
    
    def __mul__(self, other):
        if type(other) not in [int, float]:
            raise TypeError()
        
        res_d = {k: int(val * other) for k, val in self.d.items()}
        
        return Formula(res_d)

def get_rname(model, rname, verbose = True):
    rs = []
    for r in model.reactions:
        if rname == r.name:
            if verbose:
                print(r)
            rs.append(r)
    return rs
    
def match_rname(model, rname, verbose = True):
    rs = []
    for r in model.reactions:
        if rname.lower() in r.name.lower():
            rs.append(r)

    rs.sort(key=lambda r: len(r.name))

    if verbose:
        for r in rs:
            print(r, r.name)
    return rs

def get_rid(model, rid):
    return model.reactions.get_by_id(rid)

def get_mname(model, mname, verbose = True):
    ms = []
    for m in model.metabolites:
        if mname == m.name:
            if verbose:
                print(m)
            ms.append(m)
    return ms

def match_mname(model, mname, verbose = True):
    ms = []
    search_terms = mname.split(', ')
    for m in model.metabolites:
        if all([search_term.lower() in m.name.lower() for search_term in search_terms]):
            ms.append(m)
    ms.sort(key=lambda m: len(m.name))
    
    if verbose:
        for m in ms:
            print(m, m.name)

    return ms

def get_mid(model, mid):
    return model.metabolites.get_by_id(mid)

def print_rxns_mid(model, mid, verbose=False):
    for r in get_mid(model, mid).reactions:
        print(r)
        if verbose:
            print(r.build_reaction_string(use_metabolite_names=True))
            print()

def print_rxn(r, print_genes = False):
    print(r.id)
    print(r.build_reaction_string(use_metabolite_names=True))
    print(r.build_reaction_string(use_metabolite_names=False))
    if print_genes:
        print(r.gene_reaction_rule)
    print()

# if rs is a list of reaction IDs, model has to be a model
def print_rxns(rs, model = None, print_genes=False):
    for r in rs:
        if type(r) == str:
            r = get_rid(model, r)
        print_rxn(r, print_genes=print_genes)

# return all reactions that all of the given metabolites are involved in
def find_rxn_from_mids(model, mids, allow_others = True):
    if type(mids) in [set, list]:
        found_rxns = list()
        for r in model.reactions:
            curr_mids = {m.id for m in r.metabolites}
            if all([mid in curr_mids for mid in mids]) and (
                allow_others or all([mid in mids for mid in curr_mids])):
                found_rxns.append(r)
        return found_rxns

def find_rxn_from_subs_prods(model, subs, prods, _second_round = False):
    if type(subs) in [set, list] and type(prods) in [set, list]:
        found_rxns = list()
        for r in model.reactions:
            curr_subs = {m.id for m in r.metabolites if r.metabolites[m] < 0}
            if all([mid in curr_subs for mid in subs]):
                curr_prods = {m.id for m in r.metabolites if r.metabolites[m] > 0}
                if all([mid in curr_prods for mid in prods]):
                    found_rxns.append(r)
        if _second_round:
            return found_rxns
        return found_rxns + find_rxn_from_subs_prods(model, prods, subs, _second_round = True)

def merge_metabolites(model, obsolete_id, kept_id):
    kept_m = get_mid(model, kept_id)
    obsolete_m = get_mid(model, obsolete_id)
    for r in obsolete_m.reactions:
        r.add_metabolites({obsolete_m: -r.metabolites[obsolete_m], kept_m: r.metabolites[obsolete_m]})
    model.remove_metabolites([obsolete_m])
    return model

def _parse_rxn_str(s):
    res = dict()
    next_is_factor = True
    product = False
    for word in s.split():
        if word in ['+', '=']:
            res[curr_met] = curr_factor if product else -curr_factor
            next_is_factor = True
            if word == '=':
                product = True
        elif next_is_factor:
            curr_factor = float(word)
            next_is_factor = False
        else:
            curr_met = word.split('@')[0]
    res[curr_met] = curr_factor
    return res

def add_rxns_from_mnx_table(model, table_path, rids, lower_bound=0, upper_bound=1000, gpr=''):
    # read reac_prop, a collection of all mnxref reactions
    mid2mets = {m.id: m for m in model.metabolites}
    all_rxns = list()
    new_mids = list()
    for comp in model.compartments.copy():
        break
    with open(table_path, 'r') as outf:
        for i, line in enumerate(outf):
            if line.startswith('#'):
                continue
            ls = line.split('\t')
            if ls[0] in rids:
                curr_rxn = cobra.Reaction(ls[0], lower_bound=lower_bound, upper_bound=upper_bound)
                mid2fact = _parse_rxn_str(ls[1])
                d = dict()
                for mid in mid2fact:
                    if mid not in mid2mets:
                        mid2mets[mid] = cobra.Metabolite(mid, compartment=comp)
                        new_mids.append(mid)
                    d[mid2mets[mid]] = mid2fact[mid]
                curr_rxn.add_metabolites(d)
                all_rxns.append(curr_rxn)
    
    model.add_reactions(all_rxns)
    if new_mids:
        print('Warning: following metabolites were not found and added as a stub:')
        print(', '.join(new_mids))
    return model           
        
# check how many (and which, if verb) metabolites of the given bio_rxn can be produced by the model
# build a biomass reaction for each of the components,
# optimize, check if component can be produced, then remove the reactions again
def dissect_bio(model, bio_rxn, verb=False):
    ori_obj = model.objective.expression
    if type(bio_rxn) == cobra.Reaction:
        bio_rxn = bio_rxn.metabolites
    producable = dict()
    num_consumed = 0
    num_produced = 0
    for m in bio_rxn:
        biomass_reaction = Reaction('BIOMASS_tmp')
        biomass_reaction.name = 'BIOMASS_tmp'
        if bio_rxn[m] < 0:
            num_consumed += 1
            biomass_reaction.add_metabolites({m: -1})
        else:
            num_produced += 1
            biomass_reaction.add_metabolites({m: 1})
        model.add_reactions([biomass_reaction])
        model.objective = 'BIOMASS_tmp'
        producable[m.id] = model.slim_optimize()
        model.remove_reactions(['BIOMASS_tmp'])

    # print compounds that can not be produced
    num_no_production = 0 
    num_no_consumption = 0
    for mid in producable:
        if producable[mid] <= 0.00001 or producable[mid] != producable[mid]:
            if bio_rxn[get_mid(model, mid)] < 0:
                num_no_production += 1
            else:
                num_no_consumption += 1
        
    if num_no_production:
        print(f'{num_no_production} of {num_consumed} compounds cannot be produced')
        if verb:
            print('Below is their quantity in the biomass reaction:')
            for mid in producable:
                if (producable[mid] <= 0.00001 or producable[mid] != producable[mid]) and bio_rxn[get_mid(model, mid)] < 0:
                    print(f'{mid}: {bio_rxn[model.metabolites.get_by_id(mid)]}')

    if num_no_consumption:
        print(f'{num_no_consumption} of {num_produced} compounds cannot be disposed of')
        if verb:
            print('Below is their quantity in the biomass reaction:')
            for mid in producable:
                if (producable[mid] <= 0.00001 or producable[mid] != producable[mid]) and bio_rxn[get_mid(model, mid)] > 0:
                    print(f'{mid}: {bio_rxn[model.metabolites.get_by_id(mid)]}')

    model.objective = ori_obj

# return_val can be one of:
#   'flux': only return number of maximum flux
#   'sol': return cobraPy solution object
#   'model': return solved model object
def check_production(model, mid, add_export = list(), add_import = list(), exclude_rxns = list(), 
                     consumption = False, return_val = 'flux', perform_pfba = False):
    if type(mid) == str:
        target_met = get_mid(model, mid)
    with model:
        for rid in exclude_rxns:
            get_rid(model, rid).knock_out()

        exp_rxns = list()
        for exp in add_export:
            rxn = Reaction('export' + exp)
            rxn.add_metabolites({get_mid(model, exp): 1 if consumption else -1})
            exp_rxns.append(rxn)
        model.add_reactions(exp_rxns)

        imp_rxns = list()
        for imp in add_import:
            rxn = Reaction('import' + imp)
            rxn.add_metabolites({get_mid(model, imp): -1 if consumption else 1})
            imp_rxns.append(rxn)
        model.add_reactions(imp_rxns)

        biomass_reaction = Reaction('BIOMASS_tmp')
        biomass_reaction.name = 'BIOMASS_tmp'
        biomass_reaction.add_metabolites({target_met: 1 if consumption else -1})
        model.add_reactions([biomass_reaction])
        model.objective = 'BIOMASS_tmp'
        if return_val == 'sol':
            result = cobra.flux_analysis.pfba(model) if perform_pfba else model.optimize()
        elif return_val == 'model':
            if perform_pfba:
                cobra.flux_analysis.pfba(model)
            else:
                model.optimize()
            result = model.copy()
        else:
            result = model.slim_optimize()
    
    return result

# check which imports are required for growth. If find_combinations, also tests all possible combinations of import reactions for essentiality.
# the import reactions are are to be tested need to have an id that starts with IM_, EX_ (or any lower case) and need to only add metabolites to the model, not consume any 
def analyze_uptake(model, eps = 1e-6, include_excretion = False, find_combinations = False):
    growth = model.slim_optimize()

    # find import reactions
    imports = set()
    for r in model.reactions:
        if ((r.id.lower().startswith('im_') or r.id.lower().startswith('ex_')) and
            (include_excretion or all([f > 0 for f in r.metabolites.values()]))):
            imports.add(r)
    
    # check which imports are vital
    lethal = list()
    detrimental = list()
    for r in imports:
        with model:
            r.knock_out()
            curr_growth = model.slim_optimize()
        if curr_growth < eps:
            lethal.append(r)
        elif (growth - curr_growth) > eps:
            detrimental.append(r)

    print('Uptake reactions that are necessary for growth:')
    for r in lethal:
        print_rxn(r)

    print('Uptake reactions that are necessary for the current amount of growth:')
    for r in detrimental:
        print_rxn(r)

    print('Uptake reactions without effect:')
    for r in imports.difference(lethal + detrimental):
            print_rxn(r)

    if find_combinations:
        
        def powerset(iterable):
            "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
            s = list(iterable)
            return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
        
        print('Finding lethal sets...')
        for to_test in powerset(imports.difference(lethal)):
            with model:
                for r in to_test:
                    r.knock_out()
                    if model.slim_optimize() < eps:
                        print('lethal set:', ', '.join(map(lambda r: r.id, to_test)))

# determine whether two reactions are equal (1), only different in reaction bounds (-1), or not equal (0)
# ignores metabolites defined by H_id
def reactions_equal(r1, r2, H_id = 'MNXM1[s]', custom_equivalences = dict()):
    # check whether both reactions involve the same set of metabolites
    # generate dict like r.metabolites but with m.id instead of m for comparability
    # also translate the m.id with custom_equivalences 
    mets1 = {(m.id if m.id not in custom_equivalences else custom_equivalences[m.id]): r1.metabolites[m] for m in r1.metabolites if m.id != H_id}
    mets2 = {(m.id if m.id not in custom_equivalences else custom_equivalences[m.id]): r2.metabolites[m] for m in r2.metabolites if m.id != H_id}
    
    if mets1.keys() != mets2.keys():
        return 0
    
    if mets1 == mets2:
        return 1
    
    # determine whether one reaction is simply a multiplicative of the other (includes backwards written rxns with factor -1)
    poss_factor = None
    for mid in mets1:
        if poss_factor is None:
            poss_factor = mets1[mid] / mets2[mid]
        elif mets2[mid] * poss_factor != mets1[mid]:   # not all metabolites have the same factor
            return 0
    
    # either one reaction is reversible and the other is not, or the reactions are swapped in direction but not their bounds
    if not (r1.lower_bound == poss_factor * r2.lower_bound and r1.upper_bound == poss_factor * r2.upper_bound):
        return -1
    
    return 1

# basically dfs/bfs with no guaranteed order. Finds cliques.
# excludes dead-end reactions
def find_connected(model, start_rxn_id, include_dead_end = False):
    dead_ends_mets = set()
    for m in model.metabolites:
        if len(m.reactions) <= 1:
            dead_ends_mets.add(m.id)

    to_search = set([start_rxn_id])
    found = set()

    while to_search:
        if len(found) % int(len(model.reactions) / 100) == 0:
            print(f'Found {len(found)} connected reactions', end='\r')
        curr_rid = to_search.pop()
        found.add(curr_rid)
        for curr_m in get_rid(model, curr_rid).metabolites:
            for r in curr_m.reactions:
                # exclude reactions containing metabolites that only contain one (this) reaction
                if r.id not in found and (include_dead_end or all([m not in dead_ends_mets for m in r.metabolites])):
                    to_search.add(r.id)
    
    print(f'Found {len(found)} connected reactions')
    return found

# determines whether a reaction is the only reaction for any of its metabolites
# the user can define a list of reactions that do not count
def immediate_dead_end(r, excluded_rxns = list()):
    for m in r.metabolites:
        if len([r for r in m.reactions if r not in excluded_rxns]) <= 1:
            return True
    return False

# asks whether the two reactions are the same, just defined backwards (stoichiometric coefficients of r1 are the negative of those of r2.
# completely ignores the metabolite given by H_id
def defined_backwards(r1, r2, H_id = 'MNXM1[s]'):
    mets1 = {m.id: r1.metabolites[m] for m in r1.metabolites if m.id != H_id}
    mets2 = {m.id: r2.metabolites[m] for m in r2.metabolites if m.id != H_id}
    return mets1 == {mid: -mets2[mid] for mid in mets2} and mets1 != mets2

# an attempt at a very naive gapfilling algorithm that is neither optimal nor efficient, but as opposed to the cobrapy function actually works. In my experience however, manual gapfilling using KEGG or Metacyc pathway map is usually the best option.  
def my_gapfill(input_model, universal, eps = 0.0001):
    model = input_model.copy()

    to_add_first = []
    to_bound_first = []
    i = 0
    while model.slim_optimize() < eps:
        i += 1
        next_it = False
        model_backup = model.copy()
        present_rxns = {r.id for r in model.reactions}
        added_counter = 0
        bounds_counter = 0

        # add / change bounds of reactions known to cause growth
        model.add_reactions(to_add_first)
        for r in to_bound_first:
            model.reactions.get_by_id(r.id).lower_bound = r.lower_bound
            model.reactions.get_by_id(r.id).upper_bound = r.upper_bound
        if model.slim_optimize() > eps:
            print(f'\nFound a solution with {len(to_add_first)} reactions to add and {len(to_bound_first)} reactions to change bounds of:')
            print(f'Reactions to add:')
            for r in to_add_first:
                print(r)
                print(r.name)
            print(f'\nReactions to change reversibility of:')
            for r in to_bound_first:
                print(r)
                print(r.name)
            break

        for r in universal.reactions:
            added, bounds = False, False

            if r.id not in present_rxns:
                tmp_r = Reaction(id=r.id, name=r.name, lower_bound=r.lower_bound, upper_bound=r.upper_bound)
                met_dict = dict()
                for met, v in r.metabolites.items():
                    if met.id in model.metabolites:
                        met_dict[model.metabolites.get_by_id(met.id)] = v
                    else:
                        new_met = cobra.Metabolite(id=met.id, name=met.name, formula=met.formula, charge=met.charge, compartment=met.compartment)
                        met_dict[new_met] = v
                tmp_r.add_metabolites(met_dict)
                model.add_reactions([tmp_r])
                added_counter += 1
                added = True

            # differing in reversibility
            if reactions_equal(model.reactions.get_by_id(r.id), r) == -1:
                model.reactions.get_by_id(r.id).lower_bound = r.lower_bound
                model.reactions.get_by_id(r.id).upper_bound = r.upper_bound
                bounds_counter += 1
                bounds = True
            
            if model.slim_optimize() > eps:
                break

        print(f'Iteration {i}:')
        print(f'Added {added_counter} reactions from universal model')
        print(f'Changed bounds of {bounds_counter} reactions to those of universal model')
        print(f'Now, the model does {"not " if model.slim_optimize() < eps else ""}grow')
        print('------------------------------------------------')
        if model.slim_optimize() > eps:
            # add reaction to actions to perform first
            if bounds:
                to_bound_first.append(r)
            if added:
                to_add_first.append(tmp_r)
            # reset model to original state
            model = model_backup
    return model

def namespace_stats(model):
    mc_rxns = 0
    mc_rxns_mets = {'mnx': 0, 'kegg': 0, 'metacyc': 0}
    kg_rxns = 0
    kg_rxns_mets = {'mnx': 0, 'kegg': 0, 'metacyc': 0}
    mx_rxns = 0
    mx_rxns_mets = {'mnx': 0, 'kegg': 0, 'metacyc': 0}
    for r in model.reactions:
        if 'RXN' in r.id:
            mc_rxns += 1
            for m in r.metabolites:
                if 'MNX' in m.id:
                    mc_rxns_mets['mnx'] += 1
                elif (m.id.startswith('C') or m.id.startswith('D')) and all([c.isnumeric() for c in m.id[1:6]]):
                    mc_rxns_mets['kegg'] += 1
                else:
                    mc_rxns_mets['metacyc'] += 1
        elif r.id.startswith('R') and all([c.isnumeric() for c in r.id[1:6]]):
            kg_rxns += 1
            for m in r.metabolites:
                if 'MNX' in m.id:
                    kg_rxns_mets['mnx'] += 1
                elif (m.id.startswith('C') or m.id.startswith('D')) and all([c.isnumeric() for c in m.id[1:6]]):
                    kg_rxns_mets['kegg'] += 1
                else:
                    kg_rxns_mets['metacyc'] += 1
        elif 'MNXR' in r.id:
            mx_rxns += 1
            for m in r.metabolites:
                if 'MNX' in m.id:
                    mx_rxns_mets['mnx'] += 1
                elif (m.id.startswith('C') or m.id.startswith('D')) and all([c.isnumeric() for c in m.id[1:6]]):
                    mx_rxns_mets['kegg'] += 1
                else:
                    mx_rxns_mets['metacyc'] += 1
        else:
            print('Could not attribute to any namespace:', r)
    if mc_rxns:
        print(f'MetaCyc reactions: {mc_rxns}, of these metabolites are from',
          ', '.join([db + ': ' + str(round(mc_rxns_mets[db]/sum(mc_rxns_mets.values()), 3)) for db in ['kegg', 'metacyc', 'mnx']]))
    if kg_rxns:
        print(f'KEGG reactions: {kg_rxns}, of these metabolites are from',
          ', '.join([db + ': ' + str(round(kg_rxns_mets[db]/sum(kg_rxns_mets.values()), 3)) for db in ['kegg', 'metacyc', 'mnx']]))
    if mx_rxns:
        print(f'MNX reactions: {mx_rxns}, of these metabolites are from',
          ', '.join([db + ': ' + str(round(mx_rxns_mets[db]/sum(mx_rxns_mets.values()), 3)) for db in ['kegg', 'metacyc', 'mnx']]))

# translate a model to MNX namespace. Returns a copy with metabolites and reaction ID translated
# split functions need to take a metabolite or rxn id and return a tuple of the rxn/met ID, the compartment str, and:
#   for split_rid a bool whether the rxn is an exchange. If true, returned rid should be only the exchanged metabolite (or the new name)
# if single_comp is not None, builds a model without any compartment information, and the given str is appended to each metabolite ID
# Split functions still need to work!
# collision_behavior can be one of 'choose', 'print', {rid:{mid:factor}}. In the latter case the correct stoichiometry is retrieved from given dict
#   print also chooses one of the definitions. Giving a dict is currently only supported with single_comp mode
# TODO THIS FUNCTION NEEDS A REWORK!
#   Why is transport reaction MNXR96495 translated to generate Cd out of nothing?
#       Because the underlying MetaCyc reactions 3.6.3.3-RXN and RXN-21035 are already imbalanced
#   translating to correct compartments even in case of collision
#   Single comp should be migrated to its own function
#       Re-think what to do with reactions across multiple compartments
#           Ignore and remove in later stage; add transporters with more sophisticated system
#   The new gene rules might contain some genes multiple times. This is not easily solvable if we also have AND relationships
#   Think of manual translation I did so far in model_improvements

def translate_model(model, split_mid, split_rid, rxns_to_mnx, mets_to_mnx, single_comp = None, comp_transform = lambda x: x, collision_behavior = 'print'):

    counter = 0
    model_mnx = cobra.Model(model.id + '_mnx')
    rxns = dict()  # {rid: cobra.Reaction}
    gprs = dict()  # {rid: str(gpr)}
    r_mapping = dict()
    m_mapping = dict()
    not_in_collision_behavior = list()
    stoichiometry_to_correct = list()
    if type(collision_behavior) == dict:
        if single_comp is None:
            print('Error in translate_model: Giving a dict as collision_behavior is currently only supported in single_comp mode')
            return
        all_translatable_mets = {mets_to_mnx.loc[split_mid(m.id)[0]].item(): m.id for m in model.metabolites if split_mid(m.id)[0] in mets_to_mnx.index}

    for r in model.reactions:

        # generate reaction with translated id
        rid, comp, is_exch = split_rid(r.id)
        
        # try to translate reaction id as well as possible
        if is_exch:
            if rid in mets_to_mnx.index:
                rid = mets_to_mnx.loc[rid].item()
            rid = 'EX_' + rid
        else:
            if rid in rxns_to_mnx.index:
                rid = rxns_to_mnx.loc[rid].item()
            
        # put it into a compartment
        if single_comp is None and comp:
            rid += '[' + comp_transform(comp) + ']'
        elif single_comp is not None:
            rid += '[' + single_comp + ']'
        
        # assemble reaction
        curr_r = Reaction(id=rid, name=r.name, lower_bound=r.lower_bound, upper_bound=r.upper_bound, subsystem=r.subsystem)
        r_mapping[r.id] = curr_r.id

        # add metabolites with translated ids
        curr_mets = dict()
        curr_mids = set()
        empty_transport = 0
        for m, v in r.metabolites.items():
            mid, comp = split_mid(m.id)
            if mid in mets_to_mnx.index:
                mid = mets_to_mnx.loc[mid].item()
            if single_comp is None and comp:
                mid += '[' + comp_transform(comp) + ']'
            elif single_comp is not None:
                mid += '[' + single_comp + ']'
            curr_m = cobra.Metabolite(id=mid, formula=m.formula, name=m.name, charge=m.charge,
                                      compartment=single_comp if single_comp is not None else comp_transform(m.compartment))

            # for transport reactions, the same metabolite could be already present if mapping to single compartment
            if mid in curr_mids:
                for m_it in curr_mets:
                    if mid == m_it.id:
                        curr_mets[m_it] += v
            else:
                curr_mets[curr_m] = v
            curr_mids.add(curr_m.id)
            m_mapping[m.id] = curr_m.id

        if all([factor == 0 for factor in curr_mets.values()]):
            empty_transport += 1
            continue

        curr_r.add_metabolites(curr_mets)
        if curr_r.id in rxns:

            counter += 1
            # join the two gpr with OR; in case one of the reactions had no gpr, it is omitted to not break the gpr string
            gprs[curr_r.id] = ' or '.join([gpr_str for gpr_str in [gprs[curr_r.id], str(r.gpr)] if gpr_str])
            
            if not reactions_equal(curr_r, rxns[curr_r.id]) in [-1, 1]:

                # if the reactions are different but have the same ID, have a closer look
                if collision_behavior == 'print':
                    print('Encountered reactions with the same name but different definitions:')
                    print(curr_r)
                    print(rxns[curr_r.id])
                    print()

                elif type(collision_behavior) == dict:
                    if curr_r.id.removesuffix(f'[{single_comp}]') in collision_behavior:
                        new_stoich = {}
                        for mid in collision_behavior[curr_r.id.removesuffix(f'[{single_comp}]')]:

                            # check whether mid can be translated from a metabolite in the old model
                            # this makes searching in metabolites already translated obsolete
                            if mid in all_translatable_mets:
                                # if found, use name, formula etc from the original metabolite
                                curr_m = get_mid(model, all_translatable_mets[mid]).copy()
                                curr_m.compartment = single_comp if single_comp is not None else comp_transform(m.compartment)
                                curr_m.id = mid

                            else:
                                curr_m = cobra.Metabolite(id=mid, compartment=single_comp)

                            new_stoich[curr_m] = collision_behavior[curr_r.id.removesuffix(f'[{single_comp}]')][mid]
                        
                        # change reaction.metabolites
                        curr_r.subtract_metabolites(curr_r.metabolites)
                        curr_r.add_metabolites(new_stoich)
                    else:
                        not_in_collision_behavior.append(curr_r.id.removesuffix(f'[{single_comp}]'))
        else:
            # only save GPR to be added as a string; if reactions are to be merged it will be more time efficient to build the whole gpr in the end
            gprs[curr_r.id] = str(r.gpr)

            rxns[curr_r.id] = curr_r

    # add gpr rules
    for rid in gprs:
        rxns[rid].gpr = cobra.core.gene.GPR.from_string(gprs[rid])

    model_mnx.add_reactions(list(rxns.values()))

    # carry over groups
    grps = list()
    for group in model.groups:
        curr_gr = cobra.core.Group(group.id, name = group.name, kind=group.kind)
        curr_gr.add_members([model_mnx.reactions.get_by_id(r_mapping[r.id]) for r in group.members])
        grps.append(curr_gr)
    model_mnx.add_groups(grps)

    print('Encountered', counter, 'duplicated reactions')
    print('Encountered', empty_transport, 'now empty reactions')
    if not_in_collision_behavior:
        print('Encountered', len(not_in_collision_behavior), 'reactions not in the provided collision_behavior dictionary:')
        for rid in not_in_collision_behavior:
            print(rid)
    return model_mnx, m_mapping, r_mapping

# Define functions to read in result from BLAST pipeline as model

# This might be subject to change, depending on whether the rxns.csv will still change
# maybe so that this can be equivalent to ast.literal_eval. Currently, the ' are missing
def parse_genes(s):
    return ast.literal_eval(s)

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

# Read in result from BLAST pipeline as model
# This function takes a long time, the majority is spent adding (any) gene_reaction_rule to reactions for some reason
# If return_python_gpr, GPR are not generated to save a lot (!) of time, but instead returned as python dict
def build_cobra(wd, acc, tool='diamond', return_python_gpr=False):

    model = cobra.Model(f"{acc}_BLAST_KEGG")
    df = pd.read_csv(f'{wd}{tool}_models/{acc}.rxns.csv')

    rxns = []
    rids = set()
    incomplete_rxn = []

    if return_python_gpr:
        rid2genes = dict()

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
                        if return_python_gpr:
                            rid2genes[rid] = rid2genes[rid].union(genes)
                        else:
                            genes_prev = rxn.gene_reaction_rule.split(' or ')
                            rxn.gene_reaction_rule = f'( {" or ".join(set(genes_prev + genes))} )'
                    break
            continue

        # Build reaction
        rxn = cobra.Reaction(rid, lower_bound=0, upper_bound=1000)  # KEGG assumes each reaction to be reversible;
        rxn.name = row['name']
        curr_mets = ast.literal_eval(row['stoichiometry'])
        rxn.add_metabolites({cobra.Metabolite(m, compartment='s'): curr_mets[m] for m in curr_mets})
        genes = parse_gpr(row['comment'].split('\n'))
        if return_python_gpr:
            rid2genes[rid] = set(genes)
        else:
            rxn.gene_reaction_rule = f'( {" or ".join(genes)} )'

        rxns.append(rxn)
        rids.add(rid)

    model.add_reactions(rxns)
    if return_python_gpr:
        return model, incomplete_rxn, rid2genes
    return model, incomplete_rxn
