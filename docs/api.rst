API
===

indel processing on the reference
---------------------------------
.. autoclass:: indelpost.Variant()
   :members:
   :exclude-members: chrom, pos, ref, alt, reference

.. autoclass:: indelpost.NullVariant()

indel processing on alignment files (BAM)
-----------------------------------------
.. autoclass:: indelpost.VariantAlignment()
   :members:

.. autoclass:: indelpost.Contig()
   :members:
   :exclude-members: contig_dict, failed, gaps, indel_seq, is_target_right_aligned, 
                     low_qual_mapping_rate, lt_consensus_scores, lt_consensus_seq,
                     lt_target_block_consensus_scores, lt_target_block_consensus_seq,
                     mapq, mismatches, non_target_indels, qc_passed, qc_stats, rt_consensus_scores,
                     rt_consensus_scores, rt_consensus_seq, rt_target_block_consensus_seq, rt_target_block_consensus_scores,
                     splice_pattern, lt_end_pos

.. autoclass:: indelpost.FailedContig()
   :members:
   :exclude-members: target_not_found, is_low_quality, failed_anyway
