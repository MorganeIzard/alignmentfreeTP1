#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@teacher: yoann dufresne
@students: camillo and izard
"""
from alive_progress import alive_bar

def kmer2str(val, k):
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def stream_kmers(text, k):
    # Precompute the (k-1)-mer (and reverse)
    kmer = 0
    rkmer = 0
    for letter in text[:k-1]:
        # A = 00, C = 01, T = 10, G = 11
        # Forward kmer
        kmer <<= 2
        letter_value = (ord(letter) >> 1) & 0b11
        kmer += letter_value
        # Reverse kmer
        rkmer >>= 2
        rev_letter_value = (letter_value + 2) & 0b11
        rkmer += rev_letter_value << (2 * (k - 1))

    # Stream kmers
    mask = (1 << (2 * k)) - 1
    for letter in text[k-1:]:
        # Forward kmer
        kmer <<= 2
        letter_value = (ord(letter) >> 1) & 0b11
        kmer += letter_value
        kmer &= mask
        # Reverse kmer
        rkmer >>= 2
        rev_letter_value = (letter_value + 2) & 0b11
        rkmer += rev_letter_value << (2 * (k - 1))

        yield kmer, rkmer

def stream_kmers_ini(text : str, k : int, bar) -> list:
    '''
    Stream kmers from a text
    ------------
    parameters
    text : a string
    k : the length of the kmers
    ------------
    output
    a list of kmers
    '''
    list_kmer, kmer, rev_kmer = [], 0, 0

    for i in range(k):
        kmer = kmer << 2
        kmer += encode(text[i])

        rev_kmer += (3 - encode(text[i])) * 2**(2*i)
    bar()
    list_kmer.append(min(kmer, rev_kmer))

    mask = (1<<(k-1)*2)-1
    for nucl in text[k:]:
        kmer = kmer & mask
        kmer = kmer << 2
        kmer = kmer + encode(nucl)

        rev_kmer = rev_kmer >> 2
        rev_kmer += (3 - encode(nucl)) * (2**(2*k))
        list_kmer.append(min(kmer, rev_kmer))

        bar()

    return list_kmer

def encode(x : str) -> int:
    '''
    Encode a nucleotide into a number
    ------------
    parameter
    x : a nucleotide
    ------------
    output
    the number corresponding to the nucleotide
    '''
    dico = {'A':0, 'C':1, 'T':2, 'G':3}
    if x not in dico.keys() :
        return 0
    else :
        return dico[x]
