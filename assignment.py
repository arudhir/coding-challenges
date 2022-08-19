#!/usr/bin python3

from itertools import product, combinations, chain
from collections import Counter, defaultdict
from pprint import pprint, pformat
import logging
logging.basicConfig(
    level=logging.INFO,
    format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)

translation_table = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                     'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                     'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                     'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                     'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                     'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                     'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                     'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                     'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                     'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                     'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                     'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                     'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                     'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                     'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                     'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
                     }


# nomenclature for degenerate codons
expanded_code = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
                 'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T'],
                 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
                 'N': ['A', 'C', 'G', 'T']
                 }

# {amino_acid: [codons]}
codon_table = defaultdict(list)
for codon, aa in translation_table.items():
    codon_table[aa].append(codon)

# {unambiguous_nucleotides: ambiguous_nucleotides}
degenerate_code = {
    frozenset(v): k
    for k, v in expanded_code.items()
}

# helpful for validating input
valid_nucleotides = 'ACGTWSMKRYBDHVN'
valid_aa = 'GAVLIMFWPSTCYNQDEKRH*'

def reduce_codons(codons):
    """Reduces a list of codons to the ambiguous representation"""
    pos1 = set()
    pos2 = set()
    pos3 = set()
    
    for codon in codons:
        logging.debug(f'The current codon is: {codon}')
        pos1.add(codon[0])
        pos2.add(codon[1])
        pos3.add(codon[2])
    expanded_ambiguous_codon = [frozenset(pos1), frozenset(pos2), frozenset(pos3)]
    ambiguous_codon = ''.join([degenerate_code[p] for p in expanded_ambiguous_codon])
    logging.debug(f'The reduced representation is: {ambiguous_codon}')
    return ambiguous_codon

def expand_codon(codon):
    """Expands degenerate codon to a list of their canonical form
    
    ::doctest::
    >>> expand_codon('RYA')
    ['ACA', 'ATA', 'GCA', 'GTA']
    """
    expansion = []
    for nt in codon:
        nts = expanded_code[nt]
        logging.debug(f'The expansion of {nt} is {nts}')
        expansion.append(nts)
    expansion_product = product(*expansion)
    codons = list(map(lambda x: ''.join(x), expansion_product))
    logging.debug(f'The codons are: {codons}')
    # codons = list(map(lambda x: ''.join(x), (product(*[expanded_code[nt] for nt in codon]))))
    return codons

def codons_to_amino_acids(codons):
    """Turn a list of codons to the amino acids they code for"""
    amino_acids = [translation_table[c] for c in codons]
    return amino_acids

def degenerate_codon_to_amino_acids(codon):
    """Turn a degenerate codon into the amino acids made"""
    return frozenset(codons_to_amino_acids(expand_codon(codon)))

def flatten(l):
    return [x for y in l for x in y]

def get_codon_for_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set, float
        returns two values the set of most efficient codons for the input set list, 
        e.g. {'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'} and the achieved efficiency e.g. 0.75
    """
    codons = list(product(*[codon_table[aa] for aa in amino_acids]))
    coding_efficiencies = {}

    # this can't be the right way to do this, but I didn't start this early enough...
    # and I couldn't figure out a better way to solve the situation where 
    # {M, F} gave WTB _and_ WTK + WTS, each of these operands solved for itx
    reduced_codons = set([reduce_codons(c) for c in codons] + [reduce_codons(flatten(codons))])
    for dcodon in reduced_codons:
        aas = degenerate_codon_to_amino_acids(dcodon)
        coding_efficiencies[dcodon] = len(amino_acids) / len(aas)
    max_efficiency = max(coding_efficiencies.values())
    efficient_codons = set([k for k, v in coding_efficiencies.items() if v >= max_efficiency])
    logging.info(f'For {amino_acids} the most efficient codons: {efficient_codons} with efficiency {max_efficiency}')
    return efficient_codons, max_efficiency

# Stolen from https://docs.python.org/3/library/itertools.html
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def truncate_list_of_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set
        the set of sets of amino acids that can be coded with 100% efficiency, i.e. {frozenset({'V', 'A'}), frozenset({'V', 'I'})}
    """
    ret = []

    # I mean, theoretically a single amino acid can be made with 100% efficency but the assert doesn't allow that
    for aas in [aas for aas in powerset(amino_acids) if len(aas) > 1]:
        codons, eff = get_codon_for_amino_acids(aas)
        if eff == 1:
            ret.append(frozenset(aas))
    logging.info(f'Truncated list of amino acids: {ret}')
    return set(ret)


if __name__ == "__main__":
    # using sets instead of lists throughout the code since the order doesn't matter and all items should be unique
    assert get_codon_for_amino_acids({'A', 'I', 'V'}) == ({'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'}, 0.75)
    assert get_codon_for_amino_acids({'M', 'F'}) == ({'WTS', 'WTK', "WTB"}, 0.5)

    # "frozenset" here since this seems to be the only way to get a set of sets - see https://stackoverflow.com/questions/5931291/how-can-i-create-a-set-of-sets-in-python
    assert truncate_list_of_amino_acids({'A', 'V', 'I'}) == {frozenset({'V', 'A'}), frozenset({'V', 'I'})}
