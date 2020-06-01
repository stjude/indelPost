#!/usr/bin/env python3

import numpy as np
from collections import OrderedDict, Counter

# from statistics import median

from difflib import SequenceMatcher

# from .localn import make_aligner, align, findall_indels
from .utilities import *

# from .consensus import is_compatible
from .localn import findall_mismatches
from .variant import Variant


#def hard_phase_nearby_variants(target, contig, snv_nearby, indel_nearby, mapq_thresh, dbsnp)
def suggest_from_alignment(
    target, contig, pileup, dbsnp=None, mapq_lim=20, max_common_str_len=7
):

    if contig.failed or contig.mapq < mapq_lim:
        return None
     
    pileup_mapq = np.percentile([read["mapq"] for read in pileup], 15)
    if pileup_mapq < mapq_lim:
        return None

    variants_list = contig.mismatches + contig.non_target_indels

    if not variants_list:
        return None

    indexed_contig = contig.genome_indexed_contig
    
    indexed_contig, variants_list = precleaning(
        indexed_contig, variants_list, target.pos
    )

    
    deletables = deletable_variants(pileup, target)
    
    lt_loci, rt_loci, tmp = [], [], variants_list
    for var in tmp:
        if is_deletable(var, indexed_contig[var.pos][2], deletables, dbsnp):
            if var.pos < target.pos:
                lt_loci.append(var.pos)
            elif var.pos > target.pos:
                rt_loci.append(var.pos)

            variants_list.remove(var)
    
    if not variants_list:
        return None

    lt_end = max(lt_loci) if lt_loci else -np.inf
    rt_end = min(rt_loci) if rt_loci else np.inf

    if are_scattered_snvs(indexed_contig, variants_list, target):
        return None
    
    remove_deletables(indexed_contig, lt_end, target.pos, rt_end)
    
    remove_common_substrings(indexed_contig, target.pos, max_common_str_len)

    cpos = 0
    cref = ""
    calt = ""
    for k, v in indexed_contig.items():

        if not cpos:
            cpos = k
        cref += v[0]
        calt += v[1]

    cvar = Variant(target.chrom, cpos, cref, calt, target.reference).normalize()
    
    if SequenceMatcher(None, cvar.ref, cvar.alt).ratio() > 0.75:
        return None

    if cvar != target:
        return cvar
    else:
        return None


def precleaning(genome_indexed_contig, variants_list, target_pos):
    lt_loci, rt_loci = [], []

    for k, v in genome_indexed_contig.items():
        ref, alt, score = v[0], v[1], v[2]
        if not ref or not alt:
            if k < target_pos:
                lt_loci.append(k)
            elif k > target_pos:
                rt_loci.append(k)

        elif "N" in ref or "N" in alt:
            if k < target_pos:
                lt_loci.append(k)
            elif k > target_pos:
                rt_loci.append(k)

        elif score < score_lim(ref, alt):
            if k < target_pos:
                lt_loci.append(k)
            elif k > target_pos:
                rt_loci.append(k)

    lt_lim = max(lt_loci) if lt_loci else -np.inf
    rt_lim = min(rt_loci) if rt_loci else np.inf

    tmp = genome_indexed_contig.copy()
    for k, v in genome_indexed_contig.items():
        if k <= lt_lim or rt_lim <= k:
            del tmp[k]

    variants_list = [var for var in variants_list if lt_lim < var.pos < rt_lim]

    return tmp, variants_list

def score_lim(ref, alt):
    if len(ref) == len(alt) == 1:
        return 0.7
    elif len(ref) > len(alt):
        return 0.65
    elif len(alt) > 6:
        return 0.5
    else:
        return 0.6



def are_scattered_snvs(indexed_contig, variants_list, target, neighborhood=10):

    mismatches = [var for var in variants_list if not var.is_indel]
    
    if not mismatches:
        return False

    lt_score, lt_cluster_center = cluster_center(indexed_contig, mismatches, target, left=True)
    rt_score, rt_cluster_center = cluster_center(indexed_contig, mismatches, target, left=False)
    
    if lt_score and target.pos - neighborhood <= lt_cluster_center < target.pos:
        if rt_cluster_center == np.inf:
            pass
        elif rt_score and target.pos < rt_cluster_center <= target.pos + neighborhood:
            pass
        else:
            return True
    elif rt_score and target.pos < rt_cluster_center <= target.pos + neighborhood:
        if lt_cluster_center == -np.inf:
            pass
        elif lt_score and target.pos - neighborhood <= lt_cluster_center < target.pos:
            pass
        else:
            return True
    else:
        return True
    
    if is_enriched_in_near_five(mismatches, target):
        return False
    else:
        return True
    

def cluster_center(indexed_contig, mismatches, target, left):
    target_pos = target.pos
    if left:
        loci = [k for k, v in indexed_contig.items() if k <= target_pos][::-1]
        snv_loci = [var.pos for var in mismatches if var.pos < target_pos]      
    else:
        loci = [k for k, v in indexed_contig.items() if k > target_pos]
        snv_loci = [var.pos for var in mismatches if var.pos > target_pos]
   
    score, gain = 0.0, 1.0
    max_locus = -np.inf if left else np.inf
    
    if not snv_loci:
        return score, max_locus
    
    indel_len = len(target.indel_seq)
    max_score = score
    for i, locus in enumerate(loci):
        
        if locus in snv_loci:
            score += gain
        else:
            score += loss(i, indel_len)

        if score > max_score:
            max_score = score
            max_locus = locus
    
    return max_score, max_locus


def loss(i, indel_len):
    if i < 4:
        return -0.1
    elif 4 <= i and indel_len < 10:
        return -0.4
    elif 4 <= i and indel_len >= 10:
        return -0.2


def is_enriched_in_near_five(mismatches, target):
    lt_near_snvs = [var for var in mismatches if target.pos - 5 <= var.pos < target.pos]
    lt_far_snvs = [var for var in mismatches if var.pos < target.pos - 5]   
    
    rt_margin = 0 if target.is_ins else len(target.indel_seq)
    rt_near_snvs = [var for var in mismatches if target.pos < var.pos <= target.pos + rt_margin +5]
    rt_far_snvs = [var for var in mismatches if target.pos + rt_margin + 5 < var.pos]
     
    if len(lt_near_snvs) >= 3:
        return True
    elif len(lt_near_snvs) > len(lt_far_snvs):
        return True

    if len(rt_near_snvs) >= 3:
        return True
    elif len(rt_near_snvs) >  len(rt_far_snvs):
        return True

    return False


def deletable_variants(pileup, target):
    nontarget_pileup = [
        findall_mismatches(read, end_trim=10)
        for read in pileup
        if not read["is_target"] and not read["is_dirty"]
    ]

    if not nontarget_pileup:
        return []

    margin = max(10, min(20, len(target.indel_seq) * 2))
    indels = [
        v[-1]
        for read in nontarget_pileup
        for v in read["I"] + read["D"]
        if not "S" in read["cigar_string"]
        and read["covering_subread"]
        and read["covering_subread"][0] + margin
        < target.pos
        < read["covering_subread"][1] - margin
    ]

    mismatches = [
        Variant(target.chrom, v[0], v[1], v[2], target.reference)
        for read in nontarget_pileup
        for v in read["mismatches"] if v[3] > 24
    ]
    
    non_target_depth = len(nontarget_pileup)
    mismatches = [var for var, cnt in Counter(mismatches).items() if cnt / non_target_depth > 0.05]
    
    return set(indels + mismatches)


def is_deletable(variant, consensus_score, deletable_variants, dbsnp):

    if variant in deletable_variants:
        return True

    if variant.is_indel:
        rep = repeats(variant)
        if len(variant.indel_seq) < 4:
            if rep > 1:
                return True

            multiallelic = [
                True for var in deletable_variants if var.pos == variant.pos
            ]
            if multiallelic:
                return True

    if dbsnp:
        hits = variant.query_vcf(dbsnp)
        for hit in hits:

            info = hit["INFO"]
            if info.get("COMMON", False):
                return True

            if (
                info.get("G5A", False)
                or info.get("G5A", False)
                or info.get("HD", False)
            ):
                return True

            if info.get("KGPhase1", False) or info.get("KGPhase1", False):
                return True

            topmedfreq = get_freq(info.get("TOPMED", 0.0))
            caffreq = get_freq(info.get("CAP", 0.0))
            if min(topmedfreq, caffreq) < 0.999:
                return True

    return False


def repeats(indel):
    unit = to_minimal_repeat_unit(indel.indel_seq)
    return repeat_counter(unit, indel.right_flank())  # leftaligned


def get_freq(freqinfo):
    try:
        wildtype_freq = float(freqinfo.split(",")[0])
    except:
        wildtype_freq = 1.0

    return wildtype_freq


def remove_deletables(indexed_contig, lt_end, target_pos, rt_end):
    tmp = indexed_contig.copy()
    
    for k, v in tmp.items():
        if k <= lt_end < target_pos:
            del indexed_contig[k]
        elif lt_end < k < target_pos:
            if v[0] == v[1]:
                del indexed_contig[k]
            else:
                break
    
    tmp = OrderedDict(reversed(list(tmp.items())))
    for k, v in tmp.items():
        if target_pos < rt_end <= k:
            del indexed_contig[k]
        elif target_pos < k < rt_end:
            if v[0] == v[1]:
                del indexed_contig[k]
            else:
                break

    return indexed_contig


def remove_common_substrings(indexed_contig, target_pos, max_common_str_len):

    common_sub_strs = profile_common_substrings(indexed_contig)

    lt_commons = [sub_str for sub_str in common_sub_strs if sub_str[1] < target_pos]
    rt_commons = [sub_str for sub_str in common_sub_strs if target_pos < sub_str[0]]

    trim_common(indexed_contig, lt_commons, max_common_str_len, left=True)
    trim_common(indexed_contig, rt_commons, max_common_str_len, left=False)

    return indexed_contig


def trim_common(indexed_contig, commons, max_common_str_len, left):
    if not left:
        commons[::-1]
    
    deletable_commons = []
    for sub_str in commons:

        if sub_str[0] == sub_str[-1]:
            start = sub_str[0]
        else:
            start = search_nearest_lt_locus(indexed_contig, sub_str[0])

        end = sub_str[-1]
        start_event = indexed_contig[start]
        end_event = indexed_contig[end]

        event = start_event if left else end_event

        if len(event[0]) != len(event[1]):
            break
        else:
            sub_str_len = end - start
            if sub_str_len >= max_common_str_len:
                if left:
                    deletable_commons.append(end) 
                else:
                    deletable_commons.append(start)
     
   
    if deletable_commons:
        loci = [item[0] for item in list(indexed_contig.items())]
        if left:
            lim = max(deletable_commons)
            for locus in loci:
                if locus < lim:
                    del indexed_contig[locus]
        else:
            lim = min(deletable_commons)
            for locus in loci:
                if locus > lim:
                    del indexed_contig[locus]        
                   
          
def search_nearest_lt_locus(indexed_contig, pos):
    not_found = True
    while not_found:
        pos -= 1

        if indexed_contig.get(pos, False):
            not_found = False

    return pos


def profile_common_substrings(indexed_contig):

    commons = []

    items = list(indexed_contig.items())

    contig_pos = items[0][0]
    contig_end = items[-1][0]
    
    while contig_pos < contig_end:
        common_sub_str = extend_sub_str(contig_pos, indexed_contig)
        end = common_sub_str[-1]
        commons.append(common_sub_str)
        contig_pos = find_next_rt_locus(indexed_contig, end, contig_end)
        
    return commons


def find_next_rt_locus(indexed_contig, pos, contig_end):
    found = False

    while not found and pos < contig_end:
        pos += 1
        found = indexed_contig.get(pos, False)
    
    return pos

def extend_sub_str(start, indexed_contig):
    common_start, common_end = start, start

    common_sub_str = []
    for k, v in indexed_contig.items():
        if k > start and v[0] == v[1]:
            common_start = k
            common_sub_str.append(k)
        elif k > common_start > start and v[0] != v[1]:
            common_end = k
            common_sub_str.append(k)
            break

    if not common_sub_str:
        common_sub_str = [common_start, common_end]

    return common_sub_str
