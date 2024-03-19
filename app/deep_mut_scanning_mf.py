#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Written for python 3, not tested under 2.
"""
Deep mutational scan primer design for Carlos Acevedo-Rocha based on this paper: http://nar.oxfordjournals.org/content/43/2/e12.long
"""
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = ""
__version__ = "1.0.8"

N = "\n"
T = "\t"
# N = "<br/>
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import math
from warnings import warn
import random, re

# MF copied the class over from direvo/mutagenesis.py to not having to import it 
class Mutation:
    """Accepts as arguments:
    * a mutation string
    * (opt) Seq object
    * (opt) forceDNA flag (def. False), if DNA is used but with protein notation
    * (opt) coding flag (def. True) to use the ref sequence as coding.
    It has the following groups of attributes:
    * from_nuc, to_nuc, num_nuc: nucleotide from, to and number.
    * from_aa, to_aa, num_aa: protein from, to and number
    * from_codon, to_codon: codon from and to
    * type: synonymous, non-synonymous and nonsense, and frameshift
    * is_substitution: true if a substitution
    It checks whether the mutation is legittimate if the position is not zero. If it is not it will raise an Error.
    It does not change the sequence passed as argument. To do that use:
    >>> dna = MutationDNASeq('ATG')
    >>> dna.mutate('1A>T') #alters the MutationDNASeq object itself
    >>> dna.variant('1A>T') #returns a copy of the MutationDNASeq object


    Has also the method apply which returns a string where the mutation is applied to the Seq object (unchanged).
    """
    codon_codex = {
        'ATG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'YTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'ACC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TAA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TAA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'CAG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'CTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'CGG', 'M': 'ATG'},
        'CGC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TGC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TGA', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TCT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GAT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'AAC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TAC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'CCA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'CTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'CGA', 'M': 'ATG'},
        'TGA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'MGA', 'M': 'ATG'},
        'CAT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TTG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'GGA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'MGA', 'M': 'ATG'},
        'TTA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'ACG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'GTG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'YTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'AGC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'AGR', 'M': 'ATG'},
        'GCG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'GTC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'GAC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TGG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'MGG', 'M': 'ATG'},
        'AGG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'AGY', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'GGT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'CCC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'AGA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'AGY', 'Y': 'TAY',
                'V': 'GTA', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'GTT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TCG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'CAA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TAA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'CTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'CGA', 'M': 'ATG'},
        'CCG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'CTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'CGG', 'M': 'ATG'},
        'ACT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GGG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'MGG', 'M': 'ATG'},
        'ATT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GGC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'GCC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'GTA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'YTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'ACA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'CGG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'CTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'CGG', 'M': 'ATG'},
        'CCT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GAA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TAA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'ATA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'YTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'AAG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'TTC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'TTR',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TCA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'CTC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TAG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'CGA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'CTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'CGA', 'M': 'ATG'},
        'TTT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'TTR',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TGT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TGA', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TAT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TCC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'AGT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'AGR', 'M': 'ATG'},
        'CGT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'AAA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TAA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'ATC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'CTG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'CTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'CGG', 'M': 'ATG'},
        'CTA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'CTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'CGA', 'M': 'ATG'},
        'GCA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'GAG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'CTT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'AAT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GCT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'CAC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'}}

    def __init__(self, mutation, seq=None, forceDNA=False, coding=True):
        # TODO frameshift.
        # regarding frameshifts and co. there are lots of notations (http://www.hgmd.cf.ac.uk/docs/mut_nom.html seems helpful).dels are marked with 76-78delACT or 76_78del 83^84insTG 76_77insT
        # I'll implement one first.
        # TODO check how unicode in code is handled when not on my machine... delta and omega would be cool.
        # TODO seq should be a weak reference.
        self.from_aa = None
        self.to_aa = None
        self.num_aa = None
        self.from_codon = None
        self.to_codon = None
        self.from_nuc = None
        self.to_nuc = None
        self.num_nuc = None
        self.is_substitution = False
        self.type = "ERROR"
        mutation = mutation.replace("_", "-")  # not implemented yet
        mutation = mutation.replace("del", "\u0394")  # \u0394 is uppercase delta
        rexprotsub = re.match("([A-Z])(\d+)([A-Z])", mutation)  # A23T
        rexnuclsub = re.match("(\d+)([A-Z])\>([A-Z])", mutation)  # 234A>T
        rexprotdel = re.match("([A-Z])(\d+)\u0394", mutation)  # A23del
        rexnucldel = re.match("(\d+)\u0394([A-Z]?)", mutation)  # 234delA
        rexprotmanydel = re.match("([A-Z])(\d+)\-([A-Z])(\d+)\u0394", mutation)  # A23-D24del
        rexnuclmanydel = re.match("(\d+)\-(\d+)\u0394([A-Z]+?)", mutation)  # 234-235delAT
        # deal with forceDNA flag
        if forceDNA:  # a hack...
            if rexprotsub:
                mutation = str(rexprotsub.group(2)) + str(rexprotsub.group(1)) + ">" + str(rexprotsub.group(3))
                rexnuclsub = re.match("(\d+)(\w)\>(\w)", mutation)
            elif rexprotdel:
                mutation = str(rexprotsub.group(2)) + "\u0394" + str(rexprotsub.group(1))
                rexnucldel = re.match("(\d+)\u0394(\w?)", mutation)  # 234delA
            elif mutation.find(">") != -1:  # 234A>T
                warn('forceDNA flag called even if DNA mutation given')
            else:
                MutationFormatError()
        # NUCLEOTIDE
        if rexnuclsub:
            self.is_substitution = True
            self.from_nuc = rexnuclsub.group(2)
            self.to_nuc = rexnuclsub.group(3)
            self.num_nuc = int(rexnuclsub.group(1))
            if seq:
                assert seq[self.num_nuc - 1] == self.from_nuc, str(self.num_nuc) + " is " + seq[
                    self.num_nuc - 1] + ", not " + self.from_nuc
            if seq and coding:
                translation = seq.translate()._data
                r = math.floor((self.num_nuc - 1) / 3)
                self.num_aa = r + 1
                self.from_codon = seq[r * 3:r * 3 + 3]._data
                self.to_codon = seq[r * 3:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc:r * 3 + 3]._data
                self.from_aa = translation[r]
                self.to_aa = Seq(self.to_codon).translate()._data
                if self.from_aa == self.to_aa:
                    self.type = "synonymous"
                elif self.to_aa == "*":
                    self.type = "nonsense"
                else:
                    self.type = "non-synonymous"
        elif rexnucldel:  # rexnucldel = re.match("(\d+)\u0394(\w?)", mutation)  # 234delA
            self.from_nuc = rexnucldel.group(2)
            self.to_nuc = ''
            self.num_nuc = int(rexnucldel.group(1))
            if seq:
                if self.from_nuc:
                    assert seq[self.num_nuc - 1] == self.from_nuc, str(self.num_nuc) + " is " + seq[
                        self.num_nuc - 1] + ", not " + self.from_nuc
                else:
                    self.from_nuc = seq[self.num_nuc - 1]
            if seq and coding:
                translation = seq.translate()._data
                r = math.floor((self.num_nuc - 1) / 3)
                self.num_aa = r + 1
                self.from_codon = seq[r * 3:r * 3 + 3]._data
                self.to_codon = seq[r * 3:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc:r * 3 + 3]._data
                self.from_aa = translation[r]
                self.to_aa = Seq(self.to_codon).translate()._data  # TODO check if it is a frameshift
                self.type = "deletion"
        # PROTEIN
        elif rexprotsub:
            self.is_substitution = True
            self.from_aa = rexprotsub.group(1)
            self.to_aa = rexprotsub.group(3).replace("X", "*")
            self.num_aa = int(rexprotsub.group(2))
            if self.to_aa == self.from_aa:
                self.type = "synonymous"  # no questions asked.
            elif self.to_aa == "*":
                self.type = "nonsense"
            else:
                self.type = "non-synonymous"
            if seq and coding:
                assert seq.translate()[self.num_aa - 1] == self.from_aa, str(self.num_aa) + " is " + seq.translate()[
                    self.num_aa - 1] + ", not " + self.from_aa
                self.from_codon = seq._data[(self.num_aa - 1) * 3: (self.num_aa - 1) * 3 + 3]
                self.to_codon = self.codon_codex[self.from_codon][self.to_aa]
                if self.from_aa == self.to_aa:  # avoid raising errors...
                    self.from_nuc = self.from_codon[0]
                    self.to_nuc = self.from_nuc
                    self.num_nuc = self.num_aa * 3
                # crap. what if there are two or three mutations to make an aa change?
                diff = [i for i in range(3) if self.to_codon[i] != self.from_codon[i]]
                self.from_nuc = self.from_codon[diff[0]:diff[-1] + 1]
                self.to_nuc = self.to_codon[diff[0]:diff[-1] + 1]
                self.num_nuc = self.num_aa * 3 - 2 + diff[0]
        else:
            raise MutationFormatError(str(mutation))
            # TODO handle other cases

    def apply(self, seq):
        return seq[0:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc + len(self.from_nuc) - 1:]._data

    def __str__(self):
        text = str(self.num_nuc) + self.from_nuc + ">" + self.to_nuc
        if self.num_aa:
            text += " (" + self.type + ": " + self.from_aa + str(self.num_aa) + self.to_aa + ")"
        return text

    def shortform(self):
        return self.from_nuc + ">" + self.to_nuc




def parse_AAmutation(mutation, sequence, offset=0,check=True):
    """
    Coverts a AA mutation to codon
    :param mutation: AA mutation, eg. A23K or A23[NNK]
    :param sequence: the sequence
    :param offset: the start codon is n (from zero)
    :param check: bol whether to check if the original codon matches the mutated-from AA.
    :return:
    """
    sequence=str(sequence) #for now..
    AA_choice='QWERTYIPASDFGHKLCVNM*X'
    rex=re.match('(\w)(\d+)(.*)',mutation)
    if not rex:
        raise ValueError('{0} is not a valid mutation, unlike say A23K'.format(mutation))
    (start_AA,num,target)=rex.groups()
    num=int(num)
    assert start_AA in AA_choice, ValueError('{0} is not an amino acid letter'.format(start_AA))
    pos = (num - 1) * 3 + offset
    if check:
        if isinstance(sequence,Seq):
            actual_AA = str(sequence[pos:pos + 3].translate())
        else:
            actual_AA=str(Seq(sequence[pos:pos+3]).translate())
        assert actual_AA == start_AA, ValueError('In {mutation}, Pos {pos}-{pos2} ({seq}) encodes a {actual} not a {start}. Neighbourhood: {neigh}'.format(pos=pos, pos2=pos+3, seq=str(Seq(sequence[pos:pos+3])), actual=actual_AA, start=start_AA, mutation=mutation,neigh=str(Seq(sequence[pos-9:pos+12]))))
    if target in AA_choice and target != 'X': #AA target
        ori_codon=str(Seq(sequence[pos:pos+3]))
        #temp solution. Not using the Mutation class itself
        return (pos,Mutation.codon_codex[ori_codon][target])
    elif target.find('[') != -1: # nucleotide target A23[NNK]
        return (pos,re.search('\[(\w{3})\]', target).groups()[0])
        # assert if real codon?
    else:
        raise ValueError('{0} is an unrecognised target'.format(target))

def deep_mutation_scan(region, section, target_temp=55, overlap_len=22, primer_range=(None, 60), mutation='NNK',
                       GC_bonus=1, Tm_bonus=2.8, staggered=True, count_from_one=True, task='DS'):
    """
    Designs primers for quikchange for deep mutation scanning.
    Based on the overlap principle of http://nar.oxfordjournals.org/content/43/2/e12.long
    In terms of calculating the melting temperature module is used.
    For now the calculations are based on everything after the NNK. The problem is the missing thermodynamics for the various mismatches.
    :param region: the sequence of a region including neighbouring parts and not only the seq of interest.
    :param section: the range of bases to make primers for as a list/tuple of two items or a slice
    :param target_temp: the temp threshold. Remember that Phusion HF buffer is secret, but they say the Tm is +3. For salt correction see below
    :param overlap_len: the length of the overlap of the two primers. ignored if stagged = False
    :param primer_range: the min and max len of the primer.
    :param mutation: str of the desired mutatated codon
    :param GC_bonus: 5' GC clamp. This number is added to the Tm when checking if is greater than the set threshold target_temp.
    :param Tm_bonus: Temp increase due to salt. Distilled water (+0). To match the IDT oligoanalyser with 50 mM Na (+2.8&deg;C). Taq buffer (+4.9&deg;C), Phusion buffer (+11.6&deg;C), Q5 buffer (+13.3&deg;C)
    :param staggered: Boolean. If true (default) staggered primers will be designed. If not, Agilent-QC primers with full overlap will be designed.
    :param count_from_one: Boolean or int. If true (default) the first amino acid mutated will be the first one, else whatever number it was in the sequence.
    :param mode: DS the normal way deepscan. MP make primers based on mutation list.
    :return: a list of dictionaries with the following keys: base codon primer len_homology len_anneal len_primer homology_start homology_stop homology_Tm anneal_Tm

    Regarding salts. check out mt.salt_correction at http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html#salt_correction
    if needed.
    """
    #
    # sanity
    if isinstance(region, Seq):
        pass
    elif isinstance(region, str):
        region = Seq(region)
    else:
        raise TypeError('Sequence is neither string or Seq.')
    # region = region.ungap(' ')  # spaces!
    if isinstance(section, slice):
        pass
    else:  # try if it is an iterable.
        section = slice(*section)
    gene = region[section]
    #if len(gene) % 3 != 0:
    #    warn('The length of the region is not a multiple of three: ' + str(len(region)))
    # max size
    if not staggered and primer_range[0]:  # in case someone gives a rather odd input.
        overlap_len = primer_range[0]
    if not primer_range[0] or primer_range[0] < overlap_len:
        primer_range = (overlap_len, primer_range[1])
    # iterate across codons.
    geneball = []
    if count_from_one is True:
        offset = int(section.start/3 -1)
    elif isinstance(count_from_one,int):
        offset=count_from_one
    else:
        offset=0
    mutagenplan=[] #list of tuple of two: nt position, codon
    if task =='DS':
        for x in range(section.start, section.stop, 3):
            mutagenplan.append((x,mutation)) #mutation is the same for all and is a codon
    elif task == 'MP':
        for AA_mut in mutation.replace("\n"," ").split(): # in MP mode mutation is a list A45K
            mutagenplan.append(parse_AAmutation(AA_mut, region,section.start))
    else:
        raise NotImplementedError
    for (x,mutcodon) in mutagenplan:
        start = x - int(overlap_len / 2)
        stop = x + overlap_len - int(overlap_len / 2)
        codon = {'codon': region[x:x + 3],
                 'AA': str(region[x:x + 3].translate()) + str(int(x/3)-offset),
                 'base': x,
                 'homology_start': start,
                 'homology_stop': stop,
                 'len_homology': overlap_len,
                 'homology_Tm': round(mt.Tm_NN(region[start:stop]), 1)
                 }
        # iterate to find best fw primer.
        for dir in ('fw', 'rv'):
            for i in range(primer_range[0] - int(overlap_len / 2) - 3, primer_range[1] - int(overlap_len / 2) - 3):
                # the length of the annealing part of the primer is i+3 (the end of the mutated codon)
                # so the region prior to the mutation does not count: -int(overlap_len/2)-3

                # This cannot be done:
                # mut=region[start:x]+Seq('NNK')+region[x+3:x+i]
                # ori=region[start:x+i]
                # t= mt.Tm_NN(mut, c_seq=ori.complement())
                # ValueError: no thermodynamic data for neighbors 'GG/TA' available

                # this seems to pick the weakest
                # mut = region[start:x] + Seq('NNK') + region[x + 3:x + i]
                # t = mt.Tm_NN(mut)

                # so ignoring forepart
                if dir == 'fw':
                    mut = region[x + 3:x + i]
                elif dir == 'rv':
                    mut = region[x - i:x].reverse_complement()
                else:
                    raise Exception
                try:
                  t = mt.Tm_NN(mut) + float(Tm_bonus)

                  # check if the tms are good...
                  if mut[-1].upper() in ['C', 'G'] and t > target_temp - GC_bonus:
                      break
                  elif t > target_temp:
                      break
                except:
                  #print(dir)
                  #print('mut', mut)
                  #print('region', region)
                  #print('x - i', x - i)
                  #print('x', x)
                  #print('i', i)
                  #print('region[x - i:x]', region[x - i:x])
                  #t = 0
                  #mut = None
                  #break
                  warn('Problem creating primer')
                  return 'NA'
            else:
                warn('Target temperature not met. {0}C > {1}C'.format(target_temp, t))
            if dir == 'fw':
                if staggered:
                    codon[dir + '_primer'] = region[start:x].upper() + mutcodon.lower() + mut.upper()
                else:
                    # placeholder
                    codon[dir + '_primer'] = mut.upper()
            else:  # dir == 'rv'
                if staggered:
                    codon[dir + '_primer'] = region[x + 3:stop].reverse_complement().upper() + Seq(
                        mutcodon).reverse_complement().lower() + mut.upper()
                else:
                    codon[dir + '_primer'] = mut.upper()
            codon[dir + '_len_primer'] = len(codon[dir + '_primer'])
            codon[dir + '_anneal_Tm'] = round(t, 1)
            codon[dir + '_len_anneal'] = i
        if not staggered:
            codon['fw_primer'] = codon['rv_primer'].reverse_complement().upper() + mutcodon.lower() + codon['fw_primer']
            codon['rv_primer'] = codon['fw_primer'].reverse_complement()
        geneball.append(codon)
    return geneball


def randomer(n):
    """
    Generate random DNA (not a randomer (NNNNN) which is a mix of random DNA)
    :param n: length
    :return: string
    """
    alphabet = ['A', 'T', 'G', 'C']
    return ''.join([random.choice(alphabet) for x in range(n)])


def test():
    """Diagnostic!"""
    print('Testing deep_mutation_scan...')
    n = 30
    m = 21
    query = randomer(n).lower() + randomer(m).upper() + randomer(n).lower()
    print('sequence:', query)
    import csv
    f=open('out.csv', 'w', newline='')
    w = csv.DictWriter(f,
                       fieldnames='base AA codon fw_primer rv_primer len_homology fw_len_anneal rv_len_anneal fw_len_primer rv_len_primer homology_start homology_stop homology_Tm fw_anneal_Tm rv_anneal_Tm'.split())
    w.writeheader()
    w.writerows(deep_mutation_scan(query, (n, n + m)))
    f.close()
    print(open('out.csv').read())
    ##
    print('Testing parse_AAmutation...')
    print('Y2K',parse_AAmutation('Y2K', 'ATGTATGGT', 0,True)[1])


def cmdline():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="the fasta file with the sequence")
    parser.add_argument("outfile", help="the csv file to save the data into")
    parser.add_argument("section_start", type=int, help="the start of the mutagenised region")
    parser.add_argument("section_end", type=int, help="the end of the mutagenised region")
    # parser.add_argument("-T", "--target_temp", type=int, nargs = 1, dest = "target_temp", help="Target temperature")
    args = parser.parse_args()
    seq = ''.join(open(args.infile).read().split('\n')[1:])  # crude fasta reading...
    deep_mutation_scan(seq, (args.section_start, args.section_stop))


if __name__ == "__main__":
    #cmdline()
    test()
