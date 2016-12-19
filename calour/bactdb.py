import requests
from logging import getLogger


logger = getLogger(__name__)


class BactDB:
    def __init__(self):
        # Web address of the bact server
        self.dburl = 'http://amnonim.webfactional.com/scdb_main'

    def get_seq_annotations(self, sequence):
        '''Get the annotations for a sequence

        Parameters:
        -----------
        sequence : str (ACGT)

        Returns:
        --------
        curs : list of list of (curation dict,list of [Type,Value] of curation details)
        '''
        rdata = {}
        rdata['sequence'] = sequence
        res = requests.get(self.dburl + '/sequences/get_annotations', json=rdata)
        if res.status_code != 200:
            logger.warn('error getting annotations for sequence %s' % sequence)
            return []
        annotations = res.json()['annotations']
        logger.debug('Found %d annotations for sequence %s' % (len(annotations), sequence))
        return annotations

    def get_annotation_strings(self, sequence):
        '''Get nice string summaries of annotations for a given sequence

        Parameters:
        -----------
        sequence : str (ACGT)
            the sequence to query the annotation strings about

        Returns:
        --------
        shortdesc : list of (dict,str) (annotationdetails,annotationsummary)
            a list of:
                annotationdetails : dict
                    'annotationid' : int, the annotation id in the database
                    'annotationtype : str
                    ...
                annotationsummary : str
                    a short summary of the annotation
        '''
        shortdesc = []
        annotations = self.get_seq_annotations(sequence)
        for cann in annotations:
            annotationdetails = cann
            cdesc = ''
            if cann['description']:
                cdesc += cann['description']+' ('
            if cann['annotationtype'] == 'diffexp':
                chigh = []
                clow = []
                call = []
                for cdet in cann['details']:
                    if cdet[0] == 'all':
                        call.append(cdet[1])
                        continue
                    if cdet[0] == 'low':
                        clow.append(cdet[1])
                        continue
                    if cdet[0] == 'high':
                        chigh.append(cdet[1])
                        continue
                cdesc += ' high in '
                for cval in chigh:
                    cdesc += cval+' '
                cdesc += ' compared to '
                for cval in clow:
                    cdesc += cval+' '
                cdesc += ' in '
                for cval in call:
                    cdesc += cval+' '
            elif cann['annotationtype'] == 'isa':
                cdesc += ' is a '
                for cdet in cann['details']:
                    cdesc += 'cdet,'
            elif cann['annotationtype'] == 'contamination':
                cdesc += 'contamination'
            else:
                cdesc += cann['annotationtype']+' '
                for cdet in cann['details']:
                    cdesc = cdesc + ' ' + cdet[1] + ','
            shortdesc.append((annotationdetails, cdesc))
        return shortdesc

    def find_study_id(self,datamd5='',mapmd5='',getall=False):
        """
        find the data id for the data/map md5 (which are calculated on load)
        note the md5s don't change following filtering/normalization/etc... - only the original data
        input:
        datamd5 : str
            from Experiment.datamd5
        mapmd5 : str
            from Experiment.mapmd5
        getall : bool (optional)
            False (default) to get only 1st id, True to get a list of all

        output:
        expids: int (if getall=False - default) or list of int (if getall=True)
            an id or a list of ids of matching dataID indices (or None if no match)
        """
        logger.debug('findexpid for datamd5 %s mapmd5 %s' % (datamd5,mapmd5))
        details = []
        if datamd5:
            details.append(['DataMD5', datamd5])
        if mapmd5:
            details.append(['MapMD5', mapmd5])
        if len(details) == 0:
            logger.warn('Error. MapMD5 and DataMD5 both missing from finddataid')
            return None

        rdata={}
        rdata['details'] = details
        res = requests.get(self.dburl+'/experiments/get_id', json=rdata)
        if res.status_code == 200:
            expids=res.json()['expId']
            if not getall:
                if len(expids) > 1:
                    logger.warn('Problem. Found %d matches for data' % len(expids))
                logger.debug('Found study id %d' % expids[0])
                return expids[0]
            logger.debug("Found %d matches to data" % len(expids))
            return expids
        logger.error('Error getting expid from details')
        return None

    def get_study_info(self, expid):
        """
        get the information about a given study dataid
        input:
        dataid : int
            The dataid on the study (DataID field)

        output:
        info : list of (str,str)
            list of tuples for each entry in the study:
            (type,value) about dataid (i.e. ('PubMedID','100234'))
            empty if dataid not found
        """
        logger.debug('get experiment details for expid %d' % expid)
        rdata = {}
        rdata['expId'] = expid
        res = requests.get(self.dburl+'/experiments/get_details', json=rdata)
        if res.status_code == 200:
            details = res.json()['details']
            logger.debug('Found %d details for experiment %d' % (len(details), expid))
            return details
        return []


    def get_study_annotations(self, expid):
        """
        get the list of annotations for study studyid

        input:
        expid : int
            The dataid of the study

        output:
        info: list of str
            the list of curations for this study (1 item per curation)
        """
        logger.debug('get experiment annotations for expid %d' % expid)
        rdata = {}
        rdata['expId'] = expid
        res = requests.get(self.dburl+'/experiments/get_annotations', json=rdata)
        if res.status_code != 200:
            logger.warn('error getting annotations for experiment %d' % expid)
            return []
        annotations = res.json()['annotations']
        logger.debug('Found %d annotations for experiment %d' % (len(annotations),expid))
        # make it into a nice list of str
        info = []
        for cann in annotations:
            cstr = 'date:%s description:%s user:%s private:%s' % (cann['date'],cann['description'],cann['userid'],cann['private'])
            info.append(cstr)
        return info

    def add_study_data(self,data,studyid=None):
        """
        add new data entries (for a new study)
        input:
        data : list of tuples (Type:Value)
            a list of tuples of (Type,Value) to add to Data table (i.e. ("PUBMEDID","322455") etc)
        studyid : list of int
            the ids in which this study appears (from finddataid)

        output:
        suid : int
            the value of DataID for the new study (from Data table)
        """
        # we need to get a new identifier for all entries in the study
        # there should be a more elegant way to do it
        logger.debug("add_study_data for %d enteries" % len(data))
        if studyid is None:
            # add new study
            logger.debug("add_study_data for a new study")
        else:
            logger.debug('add_study_data for existing study %d' % studyid)
        rdata = {}
        rdata['expId'] = studyid
        rdata['details'] = data
        res = requests.post(self.dburl+'/experiments/add_details',json=rdata)
        if res.status_code == 200:
            newid = res.json()['expId']
            logger.debug('experiment added. id is %d' % newid)
            return newid
        else:
            logger.debug('error adding experiment. msg: %s' % res.content)
            return None
