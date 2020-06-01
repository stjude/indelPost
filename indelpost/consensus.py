#!/usr/bin/env python3
import re
import numpy as np
from collections import OrderedDict
#from .pileup import read_updater
from .utilities import most_common, to_flat_list
#from .utilities import *

cigar_ptrn = re.compile(r"[0-9]+[MIDNSHPX=]")


def make_consensus(target, targetpileup):
    
    target_pos, target_type, target_len = (
        target.pos,
        target.variant_type,
        len(target.indel_seq),
    )
    
    lt_indexed = [
        index_bases(
            read["read_start"],
            target_pos,
            target_type,
            target_len,
            read["lt_cigar"],
            read["lt_flank"],
            read["lt_ref"],
            read["lt_qual"],
            left=True,
        )
        for read in targetpileup
    ]
    
    rt_indexed = [
        index_bases(
            read["read_start"],
            target_pos,
            target_type,
            target_len,
            read["rt_cigar"],
            read["rt_flank"],
            read["rt_ref"],
            read["rt_qual"],
            left=False,
        )
        for read in targetpileup
    ]

    lt_consensus = consensus_data(lt_indexed, left=True)
    rt_consensus = consensus_data(rt_indexed, left=False)
    
    return lt_consensus, rt_consensus


def index_bases(
    read_pos, target_pos, target_type, target_len, cigar, flank, ref, qual, left
):
    indexedbases = {}
    
    if left:
        current_pos = read_pos
    else:
        if "I" in cigar[0] or "D" in cigar[0]:
            cigar = cigar[1:]
        if target_type == "I":
            current_pos = target_pos + 1
        else:
            #ref = ref[target_len:]
            current_pos = target_pos + target_len + 1

    for c in cigar:
        
        event, event_len = c[-1], int(c[:-1])

        if event == "M" or event == "S":
            for i in range(event_len):
                if ref and event == "M":
                    indexedbases[current_pos] = (ref[0], flank[0], qual[0])
                    ref = ref[1:]
                else:
                    indexedbases[current_pos] = ("", flank[0], qual[0])
                
                flank = flank[1:]
                qual = qual[1:]
                current_pos += 1

        elif event == "I":
            padding_ref, padding_qual = (
                indexedbases[current_pos - 1][0],
                indexedbases[current_pos - 1][2],
            )
            ins_seq, flank, ins_qual, qual = (
                flank[:event_len],
                flank[event_len:],
                qual[:event_len],
                qual[event_len:],
            )
            indexedbases[current_pos - 1] = (
                padding_ref,
                padding_ref + ins_seq,
                np.median([padding_qual] + list(ins_qual)),
            )
        elif event == "D":
            padding_ref, padding_qual = (
                indexedbases[current_pos - 1][0],
                indexedbases[current_pos - 1][2],
            )
            del_seq, ref = ref[:event_len], ref[event_len:]
            indexedbases[current_pos - 1] = (
                padding_ref + del_seq,
                padding_ref,
                padding_qual,
            )

            current_pos += event_len

        elif event == "N":
            current_pos += event_len
    
    
    return indexedbases


def consensus_data(indexedbases_list, left):

    consensus_index = OrderedDict()
    skip_loci = []
    for locus in locus_list(indexedbases_list, left):
        ref, consensus_base, consensus_score = get_consensus_base(indexedbases_list, locus)
        
        if len(ref) > len(consensus_base):
            del_len = len(ref) - len(consensus_base)
            skip_loci += [locus + i for i in range(1, del_len + 1)]
            
        consensus_index[locus] = (ref, consensus_base, consensus_score)
    
    for locus in skip_loci:
        if locus in consensus_index:
            del consensus_index[locus]
    
    conseq, refseq = "", ""
    scores = []
    prev_ref = ""
    prev_locus = -1
    ref_end = -1
    for locus, data in consensus_index.items():
        
        ref, consensus_base, consensus_score = data[0], data[1], data[2]
        
        if left and len(ref) != len(consensus_base):
            ref = ref[::-1]
            consensus_base = consensus_base[::-1]

        refseq += ref
        conseq += consensus_base
        scores += [consensus_score] * len(consensus_base)
        
        if prev_ref and not ref:
            ref_end = prev_locus
        
        prev_locus = locus
        prev_ref = ref

    if left:
        conseq = conseq[::-1]
        refseq = refseq[::-1]
        scores = scores[::-1]
     
    return consensus_index, ref_end, refseq, conseq, scores


def locus_list(dict_list, left):
    loci = to_flat_list([[*d] for d in dict_list])
    loci = list(set(loci))
    loci.sort(reverse=left)
    return loci


def get_consensus_base(indexedbases_list, locus, qual_lim=23):
    refs = [
        indexedbases[locus][0].upper()
        for indexedbases in indexedbases_list
        if indexedbases.get(locus, False)
    ]
    bases = [
        indexedbases[locus][1]
        for indexedbases in indexedbases_list
        if indexedbases.get(locus, False) and indexedbases[locus][1] != "N"
    ]
    quals = [
        indexedbases[locus][2]
        for indexedbases in indexedbases_list
        if indexedbases.get(locus, False)
    ]
    
    
    hq_bases = [base for base, qual in zip(bases, quals) if qual > qual_lim]

    ref = most_common(refs) if refs else ""

    if ref and bases:
        consensus_base = most_common(bases)
        consensus_score = bases.count(consensus_base) / len(bases)
        if ref != consensus_base:
            if not consensus_base in hq_bases:
                consensus_base = "N"
                consensus_score = 0.0
    elif not bases:
         consensus_base = "N"
         consensus_score = 0.0
    else:
        if hq_bases:
            consensus_base = most_common(hq_bases)
            consensus_score = bases.count(consensus_base) / len(bases)
        else:
            consensus_base = "N"
            consensus_score = 0.0

    return ref, consensus_base, consensus_score


def consensus_refseq(refseq_lst, left=False):

    if left:
        refseq_lst = [seq[::-1].upper() for seq in refseq_lst]
    else:
        refseq_lst = [seq.upper() for seq in refseq_lst]
         
    consensus_seq = ""
    consensus_rates = []
    for i in range(len(max(refseq_lst, key=len))):
        ith_chars = [ith_char(seq, i) for seq in refseq_lst if ith_char(seq, i)]
        cosensus_base = most_common(ith_chars)
            
        if cosensus_base == "N":
            consensus_rate = 0.0
        else:
            consensus_rate = ith_chars.count(cosensus_base) / len(ith_chars)
        
        consensus_seq += cosensus_base
        consensus_rates.append(consensus_rate)
    
    if left:
        consensus_seq = consensus_seq[::-1]
        consensus_rates = consensus_rates[::-1]

    return consensus_seq, consensus_rates


def ith_char(seq, i):
    try:
        return seq[i]
    except:
        return None


def is_compatible(query, subject, indel_type, partial_match=True):
    """Check if query indel looks same as subject indel

    Args:
        query (dict): dictized read
        subject (dict): indel template 
        indel_type (str): "I" for ins "D" for del
        patial_match (bool): True to allow partial match for longer insertions
    Returns:
        True/False (bool): True if query indel is found the same as subject indel
    """
    query_lt_flank, query_indel, query_del, query_rt_flank = (
        query["lt_flank"],
        query["indel_seq"],
        query.get("del_seq", ""),
        query["rt_flank"],
    )
     
    query_indel_seq = query_indel if query_indel else query_del

    # left-align check (if not, an alternative alignment of something else)
    if query_indel_seq and query_lt_flank and query_lt_flank[-1] == query_indel_seq[-1]:
        return False

    subject_lt_flank, subject_lt_scores, subject_indel, subject_rt_flank, subject_rt_scores = (
        subject.lt_consensus_seq,
        subject.lt_consensus_scores,
        subject.indel_seq,
        subject.rt_consensus_seq,
        subject.rt_consensus_scores,
    )

    # check flanking sequences for similarity
    lt_len = min(len(query_lt_flank), len(subject_lt_flank))
    rt_len = min(len(query_rt_flank), len(subject_rt_flank))

    if lt_len > 0:
        lt_query, lt_subject, lt_scores = (
            query_lt_flank[-lt_len:],
            subject_lt_flank[-lt_len:],
            subject_lt_scores[-lt_len:],
        )
    else:
        lt_query, lt_subject, lt_scores = "", "", [0]

    rt_query, rt_subject, rt_scores = (
        query_rt_flank[:rt_len],
        subject_rt_flank[:rt_len],
        subject_rt_scores[:rt_len],
    )
    
    if lt_query and not is_almost_same(
        lt_query[::-1], lt_subject[::-1], lt_scores[::-1]
    ):
        return False

    if rt_query and not is_almost_same(rt_query, rt_subject, rt_scores):
        return False

    # repeat boundary check (assuming left-aligned bam)
    rt_check = contains_repeat_end(subject_indel, rt_query, subject_rt_flank)
    if rt_query and subject_rt_flank and not rt_check:
        return False

    # check inserted/deleted sequences for similarity
    if query_indel and indel_type == "I":
        subject_len = len(subject_indel)
        query_len = len(query_indel)

        if subject_len < query_len:
            return False
        elif subject_indel == query_indel:
            return True
        elif 4 <= subject_len <= 6 and partial_match:
            return identical_for_end_n_bases(query_indel, subject_indel, 3)
        elif 7 <= subject_len <= 8 and partial_match:
            return identical_for_end_n_bases(query_indel, subject_indel, 4)
        elif 9 <= subject_len <= 10 and partial_match:
            return identical_for_end_n_bases(query_indel, subject_indel, 5)
        elif 11 <= subject_len and partial_match:
            return identical_for_end_n_bases(query_indel, subject_indel, 6)
        else:
            return False
    elif not query_indel and indel_type == "D":
        return True
    else:
        return False


def contains_repeat_end(indel_seq, query_flank, subject_flank):
    """Check if repeat boundary is contained
    """
    tmp = subject_flank.replace(indel_seq, "")

    if tmp:
        repeat_end = tmp[0]
    else:
        return False

    if repeat_end == "N":
        return False

    tmp2 = query_flank.replace(indel_seq, "")
    if tmp2:
        repeat_end2 = tmp2[0]
        return repeat_end == repeat_end2
    else:
        return False


def identical_for_end_n_bases(query_str, subject_str, n):
    return (query_str[:n] == subject_str[:n]) or (query_str[-n:] == subject_str[-n:])


def is_almost_same(
    query_seq,
    subject_seq,
    consensus_score,
    consensus_lim=0.7,
    len_lim=2,
    mismatch_lim=2,
):
    """Check if query and subject sequences are same for high consensus bases
    
    Args:
        query_seq (str): query flanking seq
        subject_seq (str): indel template flanking seq constructed by consensus
        consensus_score (list): consensus score for template flanking seq
        consensus_lim (float): threshold to define "high consensus"
        len_lim (int): seq shorter than len_lim may not contain any mismatch
        mismatch_lim (int): allow (mismatch_lim)x mismatches for seq longer than len_lim 
    Returns:
        True/False (bool): True if query and subject look same
    """
    seq_len = len(query_seq)
    if seq_len > 0 and query_seq[0] != subject_seq[0]:
        return False

    # mismatches for high-consensus bases
    mismatches = [
        (query_seq[i] != subject_seq[i] and consensus_score[i] > consensus_lim)
        for i in range(seq_len)
    ]

    if seq_len < len_lim:
        return sum(mismatches) == 0
    else:
        near_mismatches = mismatches[:len_lim]
        mid_mismatches = mismatches[len_lim : (10 * len_lim)]
        far_mismatches = mismatches[(10 * len_lim) :]
        mismatch_score = (
            sum(near_mismatches) * 2 + sum(mid_mismatches) + sum(far_mismatches) * 0.5
        )
        return mismatch_score < mismatch_lim
