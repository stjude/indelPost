#!/usr/bin/env python3

import numpy as np
from difflib import SequenceMatcher
from collections import OrderedDict, Counter

from .utilities import *
from .variant import Variant
from .localn import findall_mismatches


def phase_nearby_variants(
    target,
    contig,
    pileup,
    mapq_thresh,
    low_qual_frac_thresh,
    basequalthresh,
    snv_neighborhood,
    indel_neighborhood,
    indel_repeat_thresh,
    sequence_complexity_thresh,
    dbsnp,
    hard=False,
):

    if hard:
        if contig.failed:
            return target
        else:
            variants_list = []
            cleaned, variant_list = precleaning(contig.contig_dict, variants_list, target.pos, pileup, limit_to_target_exon=True)
            return greedy_phasing(target, cleaned)
    
    if contig.failed or contig.mapq < mapq_thresh or contig.is_target_right_aligned:
        return None
    
    if (
        seq_complexity(contig, snv_neighborhood, indel_neighborhood)
        < sequence_complexity_thresh
    ):
        return None

    if low_qual_fraction(pileup) > low_qual_frac_thresh:
        return None
    
    pileup_mapq_low_20 = np.percentile([read["mapq"] for read in pileup], 20)
    if pileup_mapq_low_20 < mapq_thresh:
        return None

    variants_to_phase = contig.mismatches + contig.non_target_indels
    if not variants_to_phase:
        return None

    indexed_contig = contig.contig_dict
    indexed_contig, variants_to_phase = precleaning(
        indexed_contig, variants_to_phase, target.pos, pileup
    )
    
    if not variants_to_phase:
        return None
    else:
        variants_in_non_targets, mut_frac = variants_in_non_target_pileup(
            pileup, target, basequalthresh
        )
        if mut_frac > 0.05:
            return None
    
    lt_loci, rt_loci, tmp = [], [], variants_to_phase.copy()
    for var in tmp:
        if is_deletable(var, variants_in_non_targets, indel_repeat_thresh, dbsnp):
            if var.pos < target.pos:
                lt_loci.append(var.pos)
            elif var.pos > target.pos:
                rt_loci.append(var.pos)
            
            variants_to_phase.remove(var)
    
    if not variants_to_phase:
        return None

    lt_end = max(lt_loci) if lt_loci else -np.inf
    rt_end = min(rt_loci) if rt_loci else np.inf

    remove_deletables(indexed_contig, lt_end, target.pos, rt_end)
    
    mismatches_to_phase = [var for var in variants_to_phase if not var.is_indel]
    non_target_indels_to_phase = [var for var in variants_to_phase if var.is_indel]
    
    if variants_to_phase:
        if not non_target_indels_to_phase:
            peak_locs = locate_mismatch_cluster_peaks(
                indexed_contig, mismatches_to_phase, target, snv_neighborhood
            )

            if peak_locs:
                remove_deletables(
                    indexed_contig, peak_locs[0], target.pos, peak_locs[1]
                )
            else:
                return None
        else:
            target_len = len(target.indel_seq)
            non_target_max_len = max(
                [len(var.indel_seq) for var in non_target_indels_to_phase]
            )

            if max(target_len, non_target_max_len) < 4:
                indel_neighborhood = int(indel_neighborhood / 2)

            remove_common_substrings(indexed_contig, target.pos, indel_neighborhood)
     
    cvar = greedy_phasing(target, indexed_contig)
    
    if (
        len(cvar.ref) > 14
        and len(cvar.alt) > 14
        and SequenceMatcher(None, cvar.ref, cvar.alt).ratio() > 0.75
    ):
        return None
    
    if cvar != target:
        return cvar
    else:
        return None


def greedy_phasing(target, indexed_contig):

    cpos = 0
    cref = ""
    calt = ""
    
    for k, v in indexed_contig.items():
        if not cpos:
            cpos = k
         
        cref += v[0]
        calt += v[1]

    return Variant(target.chrom, cpos, cref, calt, target.reference).normalize()


def seq_complexity(contig, snv_neighborhood, indel_neighorhood):

    splits = contig.get_reference_seq(split=True)
    lt_flank, rt_flank = splits[0], splits[2]
    neighorbood = min(snv_neighborhood, indel_neighorhood, len(lt_flank), len(rt_flank))

    return min(
        linguistic_complexity(lt_flank[-neighorbood:]),
        linguistic_complexity(rt_flank[:neighorbood]),
    )


def precleaning(genome_indexed_contig, variants_list, target_pos, pileup, limit_to_target_exon=True):
    lt_loci, rt_loci = [], []

    # filter low qual loci
    for k, v in genome_indexed_contig.items():
        ref, alt, score, cov = v[0], v[1], v[2], v[3]
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

        elif score < score_thresh(ref, alt, cov):
            if k < target_pos:
                lt_loci.append(k)
            elif k > target_pos:
                rt_loci.append(k)

    lt_lim = max(lt_loci) if lt_loci else -np.inf
    rt_lim = min(rt_loci) if rt_loci else np.inf

    if limit_to_target_exon:
        # within the same exon
        spliced_subreads = [
            read["covering_subread"]
            for read in pileup
            if read["is_target"] and read["covering_subread"]
        ]

        if spliced_subreads:
            lt_exon_end = min([subread[0] for subread in spliced_subreads])
            rt_exon_end = max([subread[1] for subread in spliced_subreads])
            lt_lim = max(lt_lim, lt_exon_end)
            rt_lim = min(rt_lim, rt_exon_end)

        tmp = genome_indexed_contig.copy()
        for k, v in genome_indexed_contig.items():
            if k <= lt_lim or rt_lim <= k:
                del tmp[k]

    variants_list = [var for var in variants_list if lt_lim < var.pos < rt_lim]

    return tmp, variants_list


def score_thresh(ref, alt, cov):
    if len(ref) == len(alt) == 1:
        if ref == alt:
            return 0.0
        else:
            if cov > 4:
                return 0.7 if ref == alt else 0.79
            elif 2 < cov <= 4:
                return 0.65 
            else:
                return 1.0
    elif len(ref) > 6 or len(alt) > 6:
        return 0.6
    else:
        return 0.67


def locate_mismatch_cluster_peaks(
    indexed_contig, mismatches_to_phase, target, snv_neighborhood
):

    lt_peak, lt_peak_pos = calc_peak(
        indexed_contig, mismatches_to_phase, target, snv_neighborhood, left=True
    )
    rt_peak, rt_peak_pos = calc_peak(
        indexed_contig, mismatches_to_phase, target, snv_neighborhood, left=False
    )
    
    if lt_peak > 0:
        if rt_peak > 0 or rt_peak_pos == np.inf:
            pass
        else:
            return None
    elif rt_peak > 0:
        if lt_peak > 0 or lt_peak_pos == -np.inf:
            pass
        else:
            return None
    else:
        return None

    if is_tight_cluster(mismatches_to_phase, target, snv_neighborhood):
        return (lt_peak_pos - 1, rt_peak_pos + 1)
    else:
        return None


def calc_peak(indexed_contig, mismatches, target, snv_neighborhood, left):
    target_pos = target.pos

    if left:
        loci = [k for k, v in indexed_contig.items() if k <= target_pos][::-1]
        snv_loci = [var.pos for var in mismatches if var.pos < target_pos]
    else:
        loci = [k for k, v in indexed_contig.items() if k > target_pos]
        snv_loci = [var.pos for var in mismatches if var.pos > target_pos]
    
    score, gain = 0.0, 1.0
    peak_locus = -np.inf if left else np.inf

    if not snv_loci:
        return score, peak_locus

    indel_len = len(target.indel_seq)
    scores = []
    for i, locus in enumerate(loci):

        if locus in snv_loci:
            score += gain
        else:
            score += loss(i, indel_len, snv_neighborhood)
        
        scores.append(score)
    
    peak_score = max(scores)
    if peak_score > 0.0:
        peak_idx = [i for i, j in enumerate(scores) if j == peak_score][-1]
        peak_locus = loci[peak_idx]
        score = peak_score

    return score, peak_locus


def loss(i, indel_len, snv_neighborhood):
    if indel_len < 10:
        return -1 * min(i * 1 / snv_neighborhood, 1.0)
    else:
        return -1 * min(i * 1 * 0.6 / snv_neighborhood, 1.0)


def is_tight_cluster(mismatches, target, snv_neighborhood):
    neigborhood = snv_neighborhood / 2

    lt_near_snvs = [
        var for var in mismatches if target.pos - neigborhood <= var.pos < target.pos
    ]
    lt_far_snvs = [var for var in mismatches if var.pos < target.pos - neigborhood]

    rt_margin = 0 if target.is_ins else len(target.indel_seq)
    rt_near_snvs = [
        var
        for var in mismatches
        if target.pos < var.pos <= target.pos + rt_margin + neigborhood
    ]
    rt_far_snvs = [
        var for var in mismatches if target.pos + rt_margin + neigborhood < var.pos
    ]

    if len(lt_near_snvs) < len(lt_far_snvs):
        return False

    if len(rt_near_snvs) < len(rt_far_snvs):
        return False

    return True


def variants_in_non_target_pileup(pileup, target, basequalthresh):
    nontarget_pileup = [
        findall_mismatches(read, end_trim=10)
        for read in pileup
        if not read["is_target" ] and read["is_covering"]
    ]
            
    if not nontarget_pileup:
        return [], 0.0

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
        for v in read["mismatches"] if v[3] > basequalthresh
    ]

    nontarget_pileup_vol = sum(
        max(0, len(read["ref_seq"]) - 20) for read in nontarget_pileup
    ) + 1

    mutation_frac = (len(mismatches) + len(indels)) / nontarget_pileup_vol
    
    mismatches = [
        var
        for var, cnt in Counter(mismatches).items()
        if (cnt > 2 and cnt / len(nontarget_pileup) > 0.15)
        or cnt > 5
    ]
    
    return set(indels + mismatches), mutation_frac


def is_deletable(variant, deletable_variants, indel_repeat_thresh, dbsnp):
    
    if variant in deletable_variants:
        return True

    if variant.is_indel:
        if repeats(variant) >= indel_repeat_thresh:
            return True

    if dbsnp:
        hits = variant.query_vcf(dbsnp, indel_only=False)
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
    
    #if lt_end == -np.inf:
    #    lt_end = target_pos - 1
    #if rt_end == np.inf:
    #    rt_end = target_pos + 1 
    
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
            start = search_nearest_lt_locus(indexed_contig, sub_str[0], left)

        end = sub_str[-1]
        start_event = indexed_contig[start]
        end_event = indexed_contig[end]

        # event = start_event if left else end_event

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


def search_nearest_lt_locus(indexed_contig, pos, left=True):
    if left:
        not_found = True
    else:
        not_found = False if indexed_contig.get(pos, None) else True

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
