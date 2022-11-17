
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def stream_kmers(text, k):
    list_kmer = []
    kmer = 0
    rev_kmer = 0

    for i in range(k):
        kmer = kmer << 2
        kmer += encore(text[i])

        rev_kmer += (3 - encode(text[i])) * 2**(2*i)

    list_kmer.append(min(kmer, rev_kmer))

    mask = (1 << (k-1) * 2) - 1
    for nucl in text[k:]:
        kmer = kmer & mask
        kmer = kmer << 2
        kmer = kmer + encode(text[nucl])

        rev_kmer = reverse_kmer >> 2
        rev_kmer += (3 - encode(nucl)) * (2**(2*k))
        list_kmer.append(min(kmer, rev_kmer))

    return list_kmer

def encode(x) :
    dico = {'A':0, 'C':1, 'T':2, 'G':3}
    if x not in dico.keys() :
        return 0
