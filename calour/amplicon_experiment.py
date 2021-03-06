'''
amplicon experiment (:mod:`calour.amplicon_experiment`)
=======================================================

.. currentmodule:: calour.amplicon_experiment

Classes
^^^^^^^^
.. autosummary::
   :toctree: generated

   AmpliconExperiment
'''

# ----------------------------------------------------------------------------
# Copyright (c) 2016--,  Calour development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

import numpy as np
import skbio

from .experiment import Experiment
from .util import _get_taxonomy_string, _to_list


logger = getLogger(__name__)


class AmpliconExperiment(Experiment):
    '''This class contains the data for a experiment or a meta experiment.

    Parameters
    ----------
    data : :class:`numpy.ndarray` or :class:`scipy.sparse.csr_matrix`
        The abundance table for OTUs, metabolites, genes, etc. Samples
        are in row and features in column
    sample_metadata : :class:`pandas.DataFrame`
        The metadata on the samples
    feature_metadata : :class:`pandas.DataFrame`
        The metadata on the features
    description : str
        name of experiment
    sparse : bool
        store the data array in :class:`scipy.sparse.csr_matrix`
        or :class:`numpy.ndarray`

    Attributes
    ----------
    data : :class:`numpy.ndarray` or :class:`scipy.sparse.csr_matrix`
        The abundance table for OTUs, metabolites, genes, etc. Samples
        are in row and features in column
    sample_metadata : :class:`pandas.DataFrame`
        The metadata on the samples
    feature_metadata : :class:`pandas.DataFrame`
        The metadata on the features
    exp_metadata : dict
        metadata about the experiment (data md5, filenames, etc.)
    shape : tuple of (int, int)
        the dimension of data
    sparse : bool
        store the data as sparse matrix (scipy.sparse.csr_matrix) or numpy array.
    description : str
        name of the experiment

    See Also
    --------
    Experiment
    '''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.heatmap_feature_field = 'taxonomy'
        self.heatmap_databases = ('dbbact',)

    def __repr__(self):
        '''Return a string representation of this object.'''
        return 'AmpliconExperiment %s with %d samples, %d features' % (
            self.description, self.data.shape[0], self.data.shape[1])

    def filter_taxonomy(exp, values, negate=False, inplace=False, substring=True):
        '''filter keeping only observations with taxonomy string matching taxonomy

        if substring=True, look for partial match instead of identity.
        Matching is case insensitive

        Parameters
        ----------
        values : str or list of str
            the taxonomy string/strings to filter (can be partial if substring is True)
        negate : bool (optional)
            False (default) to keep matching taxonomies, True to remove matching taxonomies
        inplace : bool (optional)
            do the filtering on the original ``Experiment`` object or a copied one.
        substring : bool (optional)
            True (default) to do partial (substring) matching for the taxonomy string,
            False to do exact matching

        Returns
        -------
        ``AmpliconExperiment``
            With only features with matching taxonomy
        '''
        if 'taxonomy' not in exp.feature_metadata.columns:
            logger.warn('No taxonomy field in experiment')
            return None

        values = _to_list(values)

        taxstr = exp.feature_metadata['taxonomy'].str.lower()

        select = np.zeros(len(taxstr), dtype=bool)
        for cval in values:
            if substring:
                select += [cval.lower() in ctax for ctax in taxstr]
            else:
                select += [cval.lower() == ctax for ctax in taxstr]

        if negate is True:
            select = ~ select

        logger.warn('%s remaining' % np.sum(select))
        return exp.reorder(select, axis=1, inplace=inplace)

    def filter_fasta(exp, filename, negate=False, inplace=False):
        '''Filter features from experiment based on fasta file

        Parameters
        ----------
        filename : str
            the fasta filename containing the sequences to use for filtering
        negate : bool (optional)
            False (default) to keep only sequences matching the fasta file, True to remove sequences in the fasta file.
        inplace : bool (optional)
            False (default) to create a copy of the experiment, True to filter inplace

        Returns
        -------
        newexp : Experiment
            filtered so contains only sequence present in exp and in the fasta file
        '''
        logger.debug('filter_fasta using file %s' % filename)
        okpos = []
        tot_seqs = 0
        for cseq in skbio.read(filename, format='fasta'):
            tot_seqs += 1
            cseq = str(cseq).upper()
            if cseq in exp.feature_metadata.index:
                pos = exp.feature_metadata.index.get_loc(cseq)
                okpos.append(pos)
        logger.debug('loaded %d sequences. found %d sequences in experiment' % (tot_seqs, len(okpos)))
        if negate:
            okpos = np.setdiff1d(np.arange(len(exp.feature_metadata.index)), okpos, assume_unique=True)

        newexp = exp.reorder(okpos, axis=1, inplace=inplace)
        return newexp

    @Experiment._record_sig
    def sort_taxonomy(exp, inplace=False):
        '''Sort the features based on the taxonomy

        Sort features based on the taxonomy (alphabetical)

        Parameters
        ----------
        inplace : bool (optional)
            False (default) to create a copy
            True to Replace data in exp
        Returns
        -------
        exp : Experiment
            sorted by taxonomy
        '''
        logger.debug('sorting by taxonomies')
        taxonomy = _get_taxonomy_string(exp, remove_underscore=True)
        sort_pos = np.argsort(taxonomy, kind='mergesort')
        exp = exp.reorder(sort_pos, axis=1, inplace=inplace)
        return exp

    @Experiment._record_sig
    def filter_orig_reads(exp, minreads, **kwargs):
        '''Filter keeping only samples with >= minreads in the original reads column
        Note this function uses the _calour_original_abundance field rather than the current sum of sequences per sample.
        So if you start with a sample with 100 reads, normalizing and filtering with other functions with not change the original reads column
        (which will remain 100).
        If you want to filter based on current total reads, use ``filter_by_data()`` instead

        Parameters
        ----------
        minreads : numeric
            Keep only samples with >= minreads reads (when loaded - not affected by normalization)

        Returns
        -------
        ``AmpliconExperiment`` - with only samples with enough original reads
        '''
        origread_field = '_calour_original_abundance'
        if origread_field not in exp.sample_metadata.columns:
            raise ValueError('%s field not initialzed. Did you load the data with calour.read_amplicon() ?' % origread_field)

        good_pos = (exp.sample_metadata[origread_field] >= minreads).values
        newexp = exp.reorder(good_pos, axis=0, **kwargs)
        return newexp

    def collapse_taxonomy(exp, level='genus', inplace=False):
        '''Collapse all features sharing the same taxonomy up to level into a single feature

        Sums abundances of all features sharing the same taxonomy up to level.

        Parameters
        ----------
        level: str or int (optional)
            the level to bin the taxonmies. can be int (0=kingdom, 1=phylum,...6=species)
            or a string ('kingdom' or 'k' etc.)
        inplace : bool (optional)
            False (default) to create a copy
            True to Replace data in exp
        '''
        level_dict = {'kingdom': 0, 'k': 0, 'phylum': 1, 'p': 1, 'class': 2, 'c': 2, 'order': 3, 'o': 3,
                      'family': 4, 'f': 4, 'genus': 5, 'g': 5, 'species': 6, 's': 6}
        if not isinstance(level, int):
            if level not in level_dict:
                raise ValueError('Unsupported taxonomy level %s. Please use out of %s' % (level, list(level_dict.keys())))
            level = level_dict[level]
        if inplace:
            newexp = exp
        else:
            newexp = exp.copy()

        def _tax_level(tax_str, level):
            # local function to get taxonomy up to given level
            ctax = tax_str.split(';')
            level += 1
            if len(ctax) < level:
                ctax.extend(['other'] * (level - len(ctax)))
            return ';'.join(ctax[:level])

        newexp.feature_metadata['_calour_tax_group'] = newexp.feature_metadata['taxonomy'].apply(_tax_level, level=level)
        newexp.merge_identical('_calour_tax_group', method='sum', axis=1, inplace=True)
        newexp.feature_metadata['taxonomy'] = newexp.feature_metadata['_calour_tax_group']
        return newexp
