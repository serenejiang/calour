import requests
from logging import getLogger
import operator
import scipy.stats

logger = getLogger(__name__)


class DBSponge:
    def __init__(self):
        # Web address of the bact server
        # self.dburl = 'http://127.0.0.1:5000'
        self.dburl = 'http://amnonim.webfactional.com/spongeworld'

    def get_name(self):
        '''Get the name of the database.
        Used for displaying when no annotations are found

        Returns
        -------
        dbname : str
            nice name of the database
        '''
        return 'SpongeWorld'

    def _post(self, api, rdata):
        '''POST a request to dbBact using authorization parameters

        Parameters
        ----------
        api : str
            the location in the dbBact REST API to post the request to
        rdata : dict
            parameters to pass to the dbBact REST API

        Returns
        -------
        res : request
            the result of the request
        '''
        res = requests.post(self.dburl + '/' + api, json=rdata)
        if res.status_code != 200:
            logger.warn('REST error %s enountered when accessing spongeworld %s' % (res.reason, api))
        return res

    def _get(self, api, rdata):
        '''GET a request to dbBact using authorization parameters

        Parameters
        ----------
        api : str
            the location in the dbBact REST API to post the request to
        rdata : dict
            parameters to pass to the dbBact REST API

        Returns
        -------
        res : request
            the result of the request
        '''
        res = requests.get(self.dburl + '/' + api, json=rdata)
        if res.status_code != 200:
            logger.warn('REST error %s enountered when accessing spongeworld %s' % (res.reason, api))
        return res

    def get_seq_annotations(self, sequence):
        '''Get the annotations for a sequence

        Parameters
        ----------
        sequence : str (ACGT)

        Returns
        -------
        annotations : dict (see spongeworld sequences/info documentation)
        '''
        rdata = {}
        rdata['sequence'] = sequence
        res = self._get('sequence/info', rdata)
        if res.status_code != 200:
            logger.warn('error getting annotations for sequence %s' % sequence)
            return None
        annotations = res.json()
        return annotations

    def get_annotation_string(self, annotations, pval=0.1):
        '''Get nice string summaries of annotations

        Parameters
        ----------
        annotations : dict (see get_sequence_annotations)
            'total_samples' : int
                the total amount of samples in the database
            'total_observed' : int
                the total number of samples where the sequence is present
            'info' : dict of {field(str): information(dict)}
                the frequency of the sequence in each field.
                information is a dict of {value(str): distribution(dict)}
                distribution contains the following key/values:
                    'total_samples': int
                        the total number of samples having this value
                    'observed_samples': int
                        the number of samples with this value which have the sequence present in them

        Returns
        -------
        desc : list of str
            a short summary of each annotation, sorted by importance
        '''
        keep = []
        total_observed = annotations['total_observed']
        if total_observed is None:
            logger.debug('sequence %s not found in database')
            return []
        total_samples = annotations['total_samples']
        null_pv = 1 - (total_observed / total_samples)
        for cfield in annotations['info'].keys():
            for cval, cdist in annotations['info'][cfield].items():
                observed_val_samples = cdist['observed_samples']
                total_val_samples = cdist['total_samples']
                cfrac = observed_val_samples / total_val_samples
                cpval = scipy.stats.binom.cdf(total_val_samples - observed_val_samples, total_val_samples, null_pv)
                if cpval <= pval:
                    cdesc = '%s:%s (%d/%d)' % (cfield, cval, observed_val_samples, total_val_samples)
                    keep.append([cdesc, cfrac, cpval])
        logger.debug('found %d significant annotations' % len(keep))
        # sort first by p-value and then by fraction (so fraction is more important)
        keep = sorted(keep, key=operator.itemgetter(2), reverse=False)
        keep = sorted(keep, key=operator.itemgetter(1), reverse=True)
        desc = [[{'annotationtype': 'dbSponge'}, ckeep[0]] for ckeep in keep]
        desc = [[{'annotationtype': 'dbSponge'}, 'Found in %f samples (%d / %d)' % (total_observed / total_samples, total_observed, total_samples)]] + desc
        return desc

    def get_seq_annotation_strings(self, sequence):
        '''Get nice string summaries of annotations for a given sequence

        Parameters
        ----------
        sequence : str (ACGT)
            the sequence to query the annotation strings about

        Returns
        -------
        shortdesc : list of (dict,str) (annotationdetails,annotationsummary)
            a list of:
                annotationdetails : dict
                    'annotationid' : int, the annotation id in the database
                    'annotationtype : str
                    ...
                annotationsummary : str
                    a short summary of the annotation
        '''
        annotations = self.get_seq_annotations(sequence)
        if annotations is None:
            return []
        shortdesc = self.get_annotation_string(annotations)
        return shortdesc
