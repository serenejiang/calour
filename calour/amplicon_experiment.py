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


logger = getLogger(__name__)


class AmpliconExperiment(Experiment):
    def plot(self, databases=('dbbact',), feature_field='taxonomy', **kwargs):
        # plot the experiment using taxonmy field and dbbact database
        super().plot(feature_field=feature_field, databases=databases, **kwargs)

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

        if not isinstance(values, (list, tuple)):
            values = [values]

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
