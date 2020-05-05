#!/usr/bin/env python3

import re
import numpy as np
from collections import namedtuple
from pysam.libcbcf cimport VariantRecord, VariantFile

cigar_ptrn = re.compile(r"[0-9]+[MIDNSHPX=]")


def most_common(lst):
    alst = list(set(lst))
    alst.sort()
    return max(alst, key=lst.count)


cpdef list to_flat_list(list lst_of_lst):
    cdef list lst
    return [i for lst in lst_of_lst for i in lst]


cpdef list to_flat_vcf_records(VariantRecord record):
    
    cdef str alt  
    
    VcfRec = namedtuple(
        "VcfRec", "chrom pos id ref alt qual filter info format samples orig"
    )

    flat_record = [
        VcfRec(
            chrom=record.chrom,
            pos=record.pos,
            id=record.id,
            ref=record.ref,
            alt=alt,
            qual=record.qual,
            filter=record.filter,
            info=record.info,
            format=record.format,
            samples=record.samples,
            orig=record,
        )
        for alt in record.alts
    ]
    
    return flat_record


cpdef dict to_dict(pysam_vcf_rec):
    cdef str k
    
    d = {} 
    for k, v in pysam_vcf_rec.items():
        if isinstance(v, tuple):
            d[k] = ",".join([str(i) for i in v])
        else:
            d[k] = v

    if d:
        return d


def match_indels(query, subject, matchby):
    if matchby != "equivalence" and not query.is_indel:
        return False

    if matchby == "equivalence":
        return query == subject

    elif matchby == "locus":
        if query.chrom != subject.chrom:
            return False

        query.normalize(inplace=True)
        subject.normalize(inplace=True)
        
        return query.pos == subject.pos

    elif matchby == "exact":
        return (query.chrom == subject.chrom) and (query.pos == subject.pos) and (query.ref == subject.ref) and (query.alt == subject.alt)


def to_minimal_repeat_unit(seq):
    """Find repeat unit in indel sequence
    """
    mid = int(len(seq) / 2)
    min_unit = seq
    
    j = 1
    found = False
    
    while j <= mid and not found:
        tandems = [seq[i : i + j] for i in range(0, len(seq), j)]
        if len(set(tandems)) == 1:
            found = True
            min_unit = list(set(tandems))[0]
        j += 1
    
    return min_unit


def repeat_counter(query_seq, flank_seq):
    """
    """
    qlen, flen = len(query_seq), len(flank_seq)
    count = 0 

    if flen < qlen:
        return count
    
    for i in range(0, flen, qlen):
        if flank_seq[i : i + qlen] == query_seq:
            count += 1
        else:
            break
    
    return count


cpdef list get_mapped_subreads(str cigarstring, int aln_start_pos, int aln_end_pos):
    
    cdef int event_len, current_pos
    cdef str cigar, event
    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef list res = []
    
    current_pos = aln_start_pos
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        if event == "M":
            res.append((current_pos, (current_pos + event_len - 1)))
            current_pos += event_len
        elif event == "I" or event == "S":
            pass
        else:
            current_pos += event_len 
    
    return res


cpdef list get_spliced_subreads(str cigarstring, int read_start_pos, int read_end_pos):
    
    cdef int i = 0
    cdef int event_len
    cdef str cigar, event, prev_event
    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef list pos_lst = [read_start_pos]
    cdef list res = []
     
    if not "N" in cigarstring:
        return [(read_start_pos, read_end_pos)]
    
    prev_event = "A"
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])
        if event == "I":
            pass
        else:
            if event == "N":
                pos_lst.append(read_start_pos - 1)
            elif prev_event == "N":
                pos_lst.append(read_start_pos)
          
            read_start_pos += event_len
        
        prev_event = event
    
    pos_lst.append(read_end_pos) 
    
    while i < len(pos_lst):
        res.append(pos_lst[i : i+2])
        i += 2
        
    return res


cpdef int get_end_pos(int read_start_pos, str lt_flank, str cigarstring):
    read_start_pos -= 1

    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef str cigar, event
    cdef int event_len, i = 0, flank_len = len(lt_flank)
    
    while flank_len > 0:
        cigar = cigar_lst[i]
        event, event_len = cigar[-1], int(cigar[:-1])
        
        if event == "D" or event == "N":
            read_start_pos += event_len
        elif event == "I":
            flank_len -= event_len
        else: 
            flank_len -= event_len
            read_start_pos += event_len

        i += 1

    return read_start_pos + flank_len


cpdef tuple locate_indels(str cigarstring, int aln_start_pos):
    aln_start_pos -= 1 

    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef str cigar, event
    cdef int event_len
    
    cdef list ins = []
    cdef list dels = []
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])
        if event == "I":
            ins.append((aln_start_pos, event_len))
        elif event == "D":
            dels.append((aln_start_pos, event_len))
            aln_start_pos += event_len
        else:
            aln_start_pos += event_len

    return ins, dels


cpdef tuple split_cigar(str cigarstring, int target_pos, int start):
     
    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef str cigar, event
    cdef int event_len, move, diff
    cdef list lt_lst = []
    cdef list rt_lst = cigar_lst

    start -= 1
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])
        
        move = 0 if event == "I" else event_len
        start += move
        rt_lst = rt_lst[1 :]

        if target_pos <= start:
            diff = start - target_pos
            lt_cigar = str(event_len - diff) + event
            lt_lst.append(lt_cigar)
            
            if diff:
                    rt_lst = [str(diff) + event] + rt_lst
            
            return lt_lst, rt_lst
        else:
            lt_lst.append(cigar)


cpdef tuple split(data, str cigarstring, int target_pos, int string_pos, bint is_for_ref, bint reverse):
    
    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef int _size = len(cigar_lst)

    cdef str cigar, event
    cdef int event_len, d_move, g_move
    
    cdef double [:] data_moves = np.zeros((_size,))
    cdef double [:] genome_moves = np.zeros((_size,))
    
    cdef int i = 0
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])
        
        if event == "N":
            d_move = 0
            g_move = event_len
        elif event == "I":
            g_move = 0
            d_move = 0 if is_for_ref else event_len
        elif event == "D":
            g_move = event_len
            d_move = event_len if is_for_ref else 0
        else:
            g_move, d_move = event_len, event_len
        
        data_moves[i] = d_move
        genome_moves[i] = g_move
        i += 1

    if reverse:
        string_pos += 1
        data = data[::-1]
        data_moves = data_moves[::-1]
        genome_moves = genome_moves[::-1]
    else:
        string_pos -= 1

    cdef int j = 0
    for d_move, g_move in zip(data_moves, genome_moves):
        if reverse:
            if target_pos < string_pos:
                string_pos -= g_move
            else:
                break
        else:
            if string_pos < target_pos:
                string_pos += g_move
            else:
                break
        j += d_move
     
    diff = string_pos - (target_pos + 1)if reverse else target_pos - string_pos
    if reverse:
        lt = data[j + diff :]
        lt = lt[::-1]
        rt = data[: j + diff]
        rt = rt[::-1]
    else:
        lt = data[: j + diff]
        rt = data[j + diff :]
     
    return lt, rt


def get_local_reference(target, pileup):

    chrom, pos, reference = target.chrom, target.pos, target.reference
    splice_patterns = [read["splice_pattern"] for read in pileup if read["splice_pattern"] != ("", "")]
   
    ref_len = reference.get_reference_length(chrom)
    
    if splice_patterns:
        lt_patterns = [ptrn[0] for ptrn in splice_patterns if ptrn[0]]
        if lt_patterns:
            lt_pattern = most_common(lt_patterns)
            lt_spl_pos = []
            for span in lt_pattern.split(":"):
                lt_spl_pos += [int(i) for i in span.split("-")] 
        else:
            lt_spl_pos = []

        rt_patterns = [ptrn[1] for ptrn in splice_patterns if ptrn[1]]
        if rt_patterns:
            rt_pattern = most_common(rt_patterns)
            rt_spl_pos = []
            for span in rt_pattern.split(":"):
                rt_spl_pos += [int(i) for i in span.split("-")]
        else:
            rt_spl_pos = []

        spl_pos = lt_spl_pos + rt_spl_pos
        last_idx = len(spl_pos) - 1
        
        left_len = 0
        first_pass = False
        local_reference = ""
        
        for i, x in enumerate(spl_pos):
            if i == 0:
                lt_end = max(0, x - 100)
                local_reference += reference.fetch(chrom, lt_end, x -1)
                rt_end = x - 1
            elif i % 2 == 1 and i != last_idx:
                local_reference += reference.fetch(chrom, x, spl_pos[i+1] - 1)
                rt_end = spl_pos[i+1] - 1
            elif i % 2 == 0:
                pass
            elif i == last_idx:
                rt_end = min(x + 100, ref_len)
                local_reference += reference.fetch(chrom, x, rt_end)
            
            if pos <= rt_end and not first_pass:
                left_len = len(local_reference) - (rt_end - pos)
                first_pass = True

    else:   
        local_reference = reference.fetch(chrom, max(0, pos - 150), min(pos + 150, ref_len))
        left_len = pos - max(0, pos - 150)
    
    return local_reference, left_len    
