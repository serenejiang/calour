# calour functions for doing differential abundance

from logging import getLogger
from collections import defaultdict

import numpy as np
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests


import calour.pbfdr


logger = getLogger(__name__)


def getpbfdr(exp, field, val1=None, val2=None, method='meandiff', transform='rankdata', numperm=1000, alpha=0.1, fdrmethod='pbfdr'):
    """
    test the differential expression between 2 groups (val1 and val2 in field field)
    using permutation based fdr (pbfdr)
    for bacteria that have a high difference.
    input:
    exp : Experiment
    field - the field for the 2 categories
    val1 - values for the first group
    val2 - value for the second group or false to compare to all other or None for all except val1
    method - the test to compare the 2 groups:
        mean - absolute difference in mean frequency
        binary - abs diff in binary presence/absence
        ranksum - abs diff in rank order (to ignore outliers)
    transform : how to transform the data
    numperm - number of random permutations to run
    maxfval - the maximal f-value (FDR) for a bacteria to keep
    usepfdr : bool
        True to use the new permutation fdr, False to use bh-fdr

    output:
    newexp - the experiment with only significant (FDR<=maxfval) difference, sorted according to difference
    """

    # if val2 is not none, need to get rid of all other samples (not val1/val2)
    if not isinstance(val1, (list, tuple)):
        val1 = [val1]
    if val2 is not None:
        if not isinstance(val2, (list, tuple)):
            val2 = [val2]
        cexp = exp.filter_samples(field, val1+val2, negate=False)
        logger.warn('%d samples left with both values' % cexp.get_num_samples())
    else:
        cexp = exp

    # prepare the labels. If correlation method, get the values, otherwise the group
    labels = np.zeros(len(cexp.sample_metadata))
    if method in ['spearman', 'pearson', 'nonzerospearman', 'nonzeropearson']:
        labels = cexp.sample_metadata[field].values
    else:
        logger.warn(np.sum(cexp.sample_metadata[field].isin(val1).values))
        labels[cexp.sample_metadata[field].isin(val1).values] = 1
        logger.warn('%d samples for value 1' % np.sum(labels))
    data = cexp.get_data(getcopy=True, sparse=False).transpose()
    keep, odif = calour.pbfdr.pbfdr(data, labels, method=method, transform=transform, alpha=alpha, numperm=numperm, fdrmethod=fdrmethod)

    keep = np.where(keep)

    if len(keep[0]) == 0:
        logger.warn('no significant features found')
        return None

    newexp = exp.reorder(keep[0], axis=1)
    odif = odif[keep[0]]
    si = np.argsort(odif, kind='mergesort')
    odif = odif[si]
    logger.warn('method %s. number of higher in %s : %d. number of higher in %s : %d. total %d' % (method, val1, np.sum(odif > 0), val2, np.sum(odif < 0), len(odif)))
    newexp = newexp.reorder(si, axis=1)
    newexp.feature_metadata['odif'] = odif
    return newexp


def get_group_terms(seqs, sequence_annotations):
    '''Get dict of number of each ontology term in seqs

    Parameters
    ----------
    seqs : list of str sequences (ACGT)
    sequence_annotations : dict of (sequence, list of ontology terms)

    Returns
    -------
    seq_annotations : dict of (ontology_term : count)
        number of times each ontology term appears in the annotations for all sequences in seqs combined
    '''
    seq_annotations = defaultdict(int)
    num_seqs_no_annotations = 0
    for cseq in seqs:
        if cseq not in sequence_annotations:
            num_seqs_no_annotations += 1
            continue
        for cterm in sequence_annotations[cseq]:
            seq_annotations[cterm] += 1
    if num_seqs_no_annotations > 0:
        logger.warn('found %d sequences with no annotations out of %d' % (num_seqs_no_annotations, len(seqs)))
    return seq_annotations


def get_term_seqs(seqs, sequence_annotations):
    '''Get dict of number of appearances in each sequence keyed by term

    Parameters
    ----------
    seqs : list of str sequences (ACGT)
    sequence_annotations : dict of (sequence, list of ontology terms)

    Returns
    -------
    seq_annotations : dict of (ontology_term : num per sequence)
        number of times each ontology term appears in each sequence in seqs
    '''
    # get all terms
    terms = set()
    for ctermlist in sequence_annotations.values():
        for cterm in ctermlist:
            terms.add(cterm)

    seq_annotations = {}
    for cterm in terms:
        seq_annotations[cterm] = np.zeros([len(seqs)])

    num_seqs_no_annotations = 0
    for idx, cseq in enumerate(seqs):
        if cseq not in sequence_annotations:
            num_seqs_no_annotations += 1
            continue
        for cterm in sequence_annotations[cseq]:
            seq_annotations[cterm][idx] += 1
    if num_seqs_no_annotations > 0:
        logger.warn('found %d sequences with no annotations out of %d' % (num_seqs_no_annotations, len(seqs)))
    return seq_annotations


def enrichment(seqs1, seqs2, sequence_annotations):
    '''Test annotation enrichment for 2 groups of sequences

    Parameters
    ----------
    exp : Experiment
    seqs1 : list of str sequences (ACGT)
    seqs2 : list of str sequences (ACGT)
    sequence_annotations : dict of (sequence, list of ontology terms)

    Returns
    -------
    '''
    logger.debug('enrichment. number of sequences in group1, 2 is %d, %d' % (len(seqs1), len(seqs2)))
    group1_terms = get_group_terms(seqs1, sequence_annotations)
    group2_terms = get_group_terms(seqs2, sequence_annotations)

    all_terms = set(group1_terms.keys()).union(set(group2_terms.keys()))
    len1 = len(seqs1)
    len2 = len(seqs2)
    tot_len = len1 + len2
    # the list of p-values for fdr
    allp = []
    # list of info per p-value
    pv = []
    for cterm in all_terms:
        num1 = group1_terms.get(cterm, 0)
        num2 = group2_terms.get(cterm, 0)
        pval = float(num1 + num2) / tot_len
        pval1 = 1 - stats.binom.cdf(num1, len1, pval)
        pval2 = 1 - stats.binom.cdf(num2, len2, pval)
        p = np.min([pval1, pval2])
        # store the result
        allp.append(p)
        cpv = {}
        cpv['pval'] = p
        cpv['observed'] = num1
        cpv['expected'] = pval * num1
        cpv['description'] = cterm
        pv.append(cpv)

    reject, _, _, _ = multipletests(allp, method='fdr_bh')
    keep = np.where(reject)[0]
    plist = []
    rat = []
    for cidx in keep:
        plist.append(pv[cidx])
        rat.append(np.abs(float(pv[cidx]['observed'] - pv[cidx]['expected'])) / np.mean([pv[cidx]['observed'], pv[cidx]['expected']]))
    si = np.argsort(rat)
    si = si[::-1]
    newplist = []
    for idx, crat in enumerate(rat):
        newplist.append(plist[si[idx]])

    return(newplist)


def enrichment2(seqs1, seqs2, sequence_annotations):
    '''Test annotation enrichment for 2 groups of sequences

    Parameters
    ----------
    exp : Experiment
    seqs1 : list of str sequences (ACGT)
    seqs2 : list of str sequences (ACGT)
    sequence_annotations : dict of (sequence, list of ontology terms)

    Returns
    -------
    '''
    logger.debug('enrichment. number of sequences in group1, 2 is %d, %d' % (len(seqs1), len(seqs2)))
    group1_terms = get_term_seqs(seqs1, sequence_annotations)
    group2_terms = get_term_seqs(seqs2, sequence_annotations)

    len1 = len(seqs1)
    len2 = len(seqs2)

    all_terms = set(group1_terms.keys()).union(set(group2_terms.keys()))
    # the list of p-values for fdr
    allp = []
    # list of info per p-value
    pv = []
    for cterm in all_terms:
        # t, p = stats.ranksums(group1_terms[cterm], group2_terms[cterm])
        t, p = stats.mannwhitneyu(group1_terms[cterm]+np.random.normal(size=len(group1_terms[cterm]))*0.001, group2_terms[cterm]+np.random.normal(size=len(group2_terms[cterm]))*0.001, alternative='two-sided')
        # store the result
        allp.append(p)
        cpv = {}
        cpv['tstat'] = t
        cpv['pval'] = p
        cpv['group1'] = np.sum(group1_terms[cterm])
        cpv['group2'] = np.sum(group2_terms[cterm])
        cpv['description'] = cterm
        pv.append(cpv)

    reject, _, _, _ = multipletests(allp, method='fdr_bh')
    keep = np.where(reject)[0]
    plist = []
    rat = []
    newpvals = []
    for cidx in keep:
        plist.append(pv[cidx])
        crat = np.abs(pv[cidx]['group1'] / len1 - pv[cidx]['group2'] / len2) / np.mean([pv[cidx]['group1'] / len1, pv[cidx]['group2'] / len2])
        rat.append(crat)
        pv[cidx]['ratio'] = crat
        newpvals.append(pv[cidx]['pval'])
    rat = np.array(rat)
    # si = np.argsort(rat)
    si = np.argsort(newpvals)[::-1]
    rat = rat[si]
    si = si[::-1]
    newplist = []
    for idx, crat in enumerate(rat):
        newplist.append(plist[si[idx]])

    return(newplist)
