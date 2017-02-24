# ----------------------------------------------------------------------------
# Copyright (c) 2016--,  Calour development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import numpy as np
from logging import getLogger

logger = getLogger(__name__)


def is_sample_v4(exp, region_seq='TACG', frac_have=0.4, min_reads=10):
    '''Test which samples in the experiment are not from the region.
    Based on the consensus sequence at the beginning of the region.

    Parameters
    ----------
    region : str (optional)
        The nucelotide sequence which is the consensus
    frac_have : float (optional)
        The fraction (per sample) of sequences containing the consensus in order to be from the region
    min_reads : float
        test only sequences with at least total min_reads

    Returns
    -------
    good_samples : list of str
        List of samples which have at least frac_have of sequences matching region_seq
    bad_samples : list of str
        List of samples which don't have at least frac_have of sequences matching region_seq
    '''

    newexp = exp.filter_min_abundance(min_reads)
    seqs_ok = newexp.feature_metadata.index.str.startswith(region_seq)
    num_seqs_ok = np.sum(newexp.data[:, seqs_ok] > 0, axis=1)
    num_seqs = np.sum(newexp.data > 0, axis=1)
    frac_ok = num_seqs_ok / num_seqs
    ok_samples = np.where(frac_ok >= frac_have)[0]
    bad_samples = np.where(frac_ok < frac_have)[0]
    return list(newexp.sample_metadata.index[ok_samples]), list(newexp.sample_metadata.index[bad_samples])
