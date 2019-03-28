#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys
import bz2
import gzip
import logging

'''
Files parse of Fastx, GTF, NCBI GFF, contain following utility function

    Fasta_Reader
        parse fasta

    Fastq_Reader
        parse single fastq

    GTF_Reader
        parse universal GTF/GFF file, return Transcript object, could convert
        annotation infor as GTF, BED, GenePred format, and extract genome,
        transcript, CDS and UTR sequence with reference genome file

    RefSeq_GFF_Reader
        specific parse GFF file from NCBI, ref above
'''

auther = 'chengchaoze@gmail.com'
update_data = '2019-03'

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s',
                    stream=sys.stderr)
logger = logging.getLogger(__name__)
# stdlog = logging.StreamHandler(sys.stdout)
# stdlog.setLevel(logging.INFO)
# formatter = logging.Formatter('%(asctime)-12s: %(levelname)-8s\n    %(message)s')
# stdlog.setFormatter(formatter)
# logger.addHandler(stdlog)


class Files(object):
    '''
        open file as file handle for iterator
    '''

    def __init__(self, File):
        self._fos = File

    def __iter__(self):
        if self._fos.lower().endswith(('.gz', '.gzip')):
            lines = gzip.open(self._fos)
        elif self._fos.lower().endswith('.bz2'):
            lines = bz2.BZ2File(self._fos)
        else:
            lines = open(self._fos)
        for index, line in enumerate(lines):
            if index % 50000 == 0 and index != 0:
                logger.info('working on line NO. of file: {:,}'.format(index))
            try:
                # if isinstance(line, bytes):
                line = line.decode('utf-8')
            except:
                pass
            yield line
        lines.close()


class Sequence(object):
    '''
        nucleic acid sequence object, contain seq, name and seq description
    '''

    def __init__(self, name, seq, descr=''):
        self._name = name
        self._seq = seq
        self._descr = descr

    @property
    def name(self):
        return self._name

    @property
    def seq(self):
        return self._seq

    @property
    def descr(self):
        return self._descr

    @property
    def length(self):
        return len(self._seq)

    def reverse_complement(self):
        paired_rc = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C',
            'N': 'N',
            'a': 't',
            't': 'a',
            'c': 'g',
            'g': 'c',
            'n': 'n',
        }
        seq = ''.join([paired_rc[i] for i in self._seq[::-1]])
        return Sequence(self._name, seq, self._descr)

    def write_to_fasta_file(self, fp, chara=80):
        if self._descr:
            fp.write('>{} {}\n'.format(self._name, self._descr))
        else:
            fp.write('>{}\n'.format(self._name))
        fp.write('{}'.format(self.__class__._seq_formater(self._seq, chara)))

    def write_to_tab_file(self, fp):
        fp.write('{}\t{}\n'.format(self._name, self._seq))

    @classmethod
    def _seq_formater(cls, seq, length):
        tmp = ''
        for i in range(0, len(seq), length):
            tmp += seq[i : (i + length)] + '\n'
        return tmp


class SequenceWithQual(Sequence):
    '''
        nucleic acid sequence object, contain seq, name and qual
    '''

    def __init__(self, name, seq, qual):
        if len(seq) != len(qual):
            raise Exception('length is Inconsistent seq and qual string')
        self._name = name
        self._seq = seq
        # super().__init__(name, seq)
        # Sequence.__init__(self, name, seq)
        self._qual = qual

    @property
    def name(self):
        return self._name

    @property
    def seq(self):
        return self._seq

    @property
    def qual(self):
        return [ord(i) for i in self._qual]

    @property
    def qualstr(self):
        return self._qual

    @property
    def length(self):
        return len(self._seq)

    def write_to_fastq_file(self, fp):
        fp.write('@{}\n'.format(self._name))
        fp.write('{}\n'.format(self._seq))
        fp.write('+\n{}\n'.format(self._qual))

    def write_to_fasta_file(self, fp):
        fp.write('>{}\n'.format(self._name))
        fp.write('{}\n'.format(self._seq))


class Fasta_Reader(Files):
    '''
        File parser for fasta file of nucleic acid
    '''

    def __init__(self, fasta):
        Files.__init__(self, fasta)

    def __iter__(self):
        seq = None
        for line in Files.__iter__(self):
            if line.startswith('>'):
                if seq:
                    yield Sequence(seqid, seq, descr)
                seqid, _, descr = line.strip('>\t \n').partition(' ')
                seq = ''
            else:
                assert seq is not None, 'FASTA file does not start with \'>\'.'
                seq += line.strip()
        yield Sequence(seqid, seq, descr)


class Fastq_Reader(Files):
    '''
        File parser for fastq file of nucleic acid
    '''

    def __init__(self, fastq):
        self._fastq = fastq
        Files.__init__(self, fastq)

    def __iter__(self):
        iter = Files.__iter__(self)
        try:
            while True:
                seqid = next(iter)
                try:
                    seq, _, qual = next(iter), next(iter), next(iter)
                except StopIteration:
                    logger.error('incompelete Fastq file: {}\n'.format(self._fastq))
                    sys.exit(1)
                seqid, seq, qual = [i.strip('\t \n') for i in [seqid[1:], seq, qual]]
                yield SequenceWithQual(seqid, seq, qual)
        except StopIteration:
            pass

    def phred_judge(self):
        iter = Files.__iter__(self)
        _, _, _, qual = next(iter), next(iter), next(iter), next(iter)
        phred = min([ord(i) for i in qual.strip()])
        # phred = fromstring(qual, dtype=byte)
        return 'phred+64' if phred > 64 else 'phred+33'


class Transcript(object):
    '''
        Transcript structure annotation information

        Parameter
        ---------
        Tid:   sring
            Transcript IDs
        chro:   sting
            chromesome ID
        start, end:   int
            0-based location
        strand:   string
            choices from {-, +}
        exon:  int list
            [[exonstart1, exonend1], [exonstart2, exonend2], ]

        cds:   int list
            [[CDSstart1, CDSend1], [CDSstart2, CDSend2], ]
        infor:    dict
            transcript_type:   string
            gene_id:   string
            gene_type:   string
    '''

    def __init__(self, Tid, chro, start, end, strand, exon, cds=None, infor=None, suffix=None):
        self._name = self.__del_suffix(Tid, suffix)
        self._chro = chro
        self._start = start
        self._end = end
        self._strand = strand
        # msg = 'missing "exon" and "CDS" feature of transcript "{}".'.format(Tid)
        # assert any([exon, cds]), msg
        self._exon = sorted(exon, key=lambda x: x[0], reverse=False)
        self._cds = sorted(cds, key=lambda x: x[0], reverse=False) if cds else None
        infor = infor if infor else {}
        geneid = infor.get('gene_id', Tid)
        self._gene_id = self.__del_suffix(geneid, suffix)
        self._type = infor.get('transcript_type', None)  # 'protein_coding'
        self._gene_type = infor.get('gene_type', None)
        self._transcript_name = infor.get('transcript_name', self._name)
        self._gene_name = infor.get('gene_name', None)
        self._supply = infor
        self.__valid_data()

    def __del_suffix(self, string, suffix):
        if suffix:
            suffix_len = len(suffix)
            assert len(string) > suffix_len, 'Gene suffix is illegal'
            if suffix == string[-suffix_len:]:
                return string[:-suffix_len]
        return string

    def __valid_data(self):
        assert self._strand in ['-', '+'], 'unsupported strand Symbol, {-, +}'
        if self._exon:
            start, end = self._exon[0][0], self._exon[-1][1]
            msg = ('{}: Exon region Over range of the transcriptional region.\n'
                   '    transcript region:  {}:{}-{}\n    Exon region:  {}-{}'
                  ).format(self._name, self._chro, self._start, self._end, start, end)
            assert self._start <= start < end <= self._end, msg

            for each in self._exon:
                assert each[1] > each[0], 'End loc must more than Start of Exon.'
        if self._cds:
            start, end = self._cds[0][0], self._cds[-1][1]
            msg = ('{}: CDS region Over range of the transcriptional region.\n'
                   '    transcript region:  {}:{}-{}\n    coding sequence region:  {}-{}'
                  ).format(self._name, self._chro, self._start, self._end, start, end)
            assert self._start <= start < end <= self._end, msg
            for each in self._cds:
                assert each[1] > each[0], 'End loc must more than Start of CDS.'

    @property
    def id(self):
        return self._name

    @property
    def chro(self):
        return self._chro

    @property
    def start(self):
        '''
            1-based
        '''
        return self._start + 1

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def basic_info(self):
        '''
            return list object, 1-based postion
                [name, chro, start, end, strand]
        '''
        return (self._name, self._chro, self._start + 1, self._end, self._strand)

    @property
    def type(self):
        return self._type

    @property
    def name(self):
        return self._transcript_name

    @property
    def exon(self):
        return self._exon

    @property
    def CDS(self):
        return self._cds

    @property
    def gene_id(self):
        return self._gene_id

    @property
    def gene_type(self):
        return self._gene_type

    @property
    def gene_name(self):
        return self._gene_name

    @property
    def length_transcript(self):
        EXON = self._exon if self._exon else self._cds
        return sum([i[1] - i[0] for i in EXON])

    @property
    def length(self):
        return self._end - self._start

    @property
    def exon_num(self):
        EXON = self._exon if self._exon else self._cds
        return len(EXON)

    def __supply_annot_for_gtf(self, hash):
        defaultkey = {
            'gene_id',
            'gene_name',
            'gene_type',
            'transcript_name',
            'transcript_type',
        }
        keys = set(hash.keys()) - defaultkey
        supply_info = {i: hash[i] for i in keys}
        tmp = None
        if keys:
            tmp = ' '.join(
                ['{} "{}";'.format(i[0], i[1]) for i in supply_info.items() if i[1]]
            )
        return tmp

    def __UTR_identify(self):
        utr = []
        if self._cds and self._exon and (self._cds != self._exon):
            coding_start, coding_end = self._cds[0][0], self._cds[-1][1]
            pos_s = '5UTR' if self._strand == '+' else '3UTR'
            pos_l = '3UTR' if self._strand == '+' else '5UTR'
            for each in self._exon:
                start, end = each
                if start < end < coding_start:
                    utr.append([pos_s, start, end])
                elif start < coding_start < end:
                    utr.append([pos_s, start, coding_start])
                if start < coding_end < end:
                    utr.append([pos_l, coding_end, end])
                elif coding_end < start < end:
                    utr.append([pos_l, start, end])
        return utr

    @classmethod
    def __list2str(cls, lst):
        return '{}\n'.format('\t'.join([str(i) for i in lst]))

    def to_gtf(self, fp):
        '''
            parameter
            ---------
            fp:    file handle for output standrand GTF file
        '''
        type_t = 'transcript_type "{}"; '.format(self._type) if self._type else ''
        type_g = 'gene_type "{}"; '.format(self._gene_type) if self._gene_type else ''
        name_t = 'transcript_name "{}"; '.format(self._transcript_name) if self._transcript_name else ''
        name_g = 'gene_name "{}"; '.format(self._gene_name) if self._gene_name else ''
        identifier = 'transcript_id "{}"; gene_id "{}"; '.format(self._name, self._gene_id)
        supply_info = self.__supply_annot_for_gtf(self._supply)
        supply_info = supply_info if supply_info else ''
        attr = ''.join([identifier, type_t, type_g, name_t, name_g, supply_info])

        transcript = [
            self._chro,
            '.',
            'transcript',
            self._start + 1,
            self._end,
            '.',
            self._strand,
            '.',
            attr,
        ]
        fp.write(self.__class__.__list2str(transcript))

        Feature = []
        order = False if self._strand == '+' else True
        if self._exon:
            for i, each in enumerate(self._exon):
                start, end = each
                Feature.append(
                    [self._chro, '.', 'exon', start + 1, end, '.', self._strand, '.', attr]
                )
        if self._cds:
            for i, each in enumerate(self._cds):
                start, end = each
                Feature.append(
                    [self._chro, '.', 'CDS', start + 1, end, '.', self._strand, '.', attr]
                )
        UTR = self.__UTR_identify()
        for utr in UTR:
            labels, start, end = utr
            Feature.append(
                [self._chro, '.', labels, start + 1, end, '.', self._strand, '.', attr]
            )
        Feature = sorted(Feature, key=lambda x: x[3], reverse=order)
        for each in Feature:
            fp.write(self.__class__.__list2str(each))

    def to_bed(self, fp):
        '''
            parameter
            ---------
            fp:    file handle for output 12 columns bed file
        '''
        if self._cds:
            cstart = self._cds[0][0]
            cend = self._cds[-1][1]
            # coding transcript
        else:
            cstart, cend = self._end, self._end
            # non coding transcript
        exon_num = len(self._exon)
        exon_len = ''.join(['{},'.format(x[1] - x[0]) for x in self._exon])
        exon_start = ''.join(['{},'.format(x[0] - self._start) for x in self._exon])

        transcript = [
            self._chro,
            self._start,
            self._end,
            self._name,
            0,
            self._strand,
            cstart,
            cend,
            0,
            exon_num,
            exon_len,
            exon_start,
        ]
        fp.write(self.__class__.__list2str(transcript))

    def to_genePred(self, fp, CIRCexplorer2=False):
        '''
            parameter
            ---------
            fp:     file handle for output GenePred file,
            CIRCexplorer2:   bool
                whether create CIRCexplorer2 style genePred, {True, False}

            reference
            ---------
            url: https://genome.ucsc.edu/FAQ/FAQformat#format9
            CIRCexplorer2_format:
                http://circexplorer2.readthedocs.io/en/latest/tutorial/setup/
        '''
        if self._cds:
            cstart = self._cds[0][0]
            cend = self._cds[-1][1]
            # coding transcript
        else:
            cstart, cend = self._end, self._end
            # non coding transcript
        exon_num = len(self._exon)
        estart = ''.join(['{},'.format(i[0]) for i in self._exon])
        eend = ''.join(['{},'.format(i[1]) for i in self._exon])

        transcript = [
            self._name,
            self._chro,
            self._strand,
            self._start,
            self._end,
            cstart,
            cend,
            exon_num,
            estart,
            eend,
        ]
        if CIRCexplorer2:
            transcript.insert(0, self._gene_id)
        fp.write(self.__class__.__list2str(transcript))

    def extract_genomic_seq(self, seq_dict, fp):
        '''
            parameter
            ---------
            seq_dict: dict
                {chro: seq, }
            fp:         file handle for output sequence file
        '''
        seq = seq_dict[self._chro][self._start : self._end]
        wanted = Sequence(
            '{}_genomic'.format(self._name),
            seq,
            'gene_id:{}'.format(self._gene_id)
        )
        if self._strand == '-':
            wanted = wanted.reverse_complement()
        wanted.write_to_fasta_file(fp)

    def extract_isoform_seq(self, seq_dict, fp):
        '''
            parameter
            ---------
            seq_dict: dict
                {chro: seq, }
            fp:         file handle for output sequence file
        '''
        if self._exon is None:
            logger.info('{}: -- not annot exon --'.format(self._name))
            return None
        seqtmp = ''
        for each in self._exon:
            start, end = each
            seqtmp += seq_dict[self._chro][start:end]
        wanted = Sequence(self._name, seqtmp, 'gene_id:{}'.format(self._gene_id))
        if self._strand == '-':
            wanted = wanted.reverse_complement()
        wanted.write_to_fasta_file(fp)

    def extract_CDS_seq(self, seq_dict, fp):
        '''
            parameter
            ---------
            seq_dict: dict
                {chro: seq, }
            fp:         file handle for output sequence file
        '''
        if self._cds is None:
            logger.info('{}: -- not annot CDS --'.format(self._name))
            return None
        seqtmp = ''
        for each in self._cds:
            start, end = each
            seqtmp += seq_dict[self._chro][start:end]
        wanted = Sequence(self._name, seqtmp, 'gene_id:{}'.format(self._gene_id))
        if self._strand == '-':
            wanted = wanted.reverse_complement()
        wanted.write_to_fasta_file(fp)

    def extract_exon_seq(self, seq_dict, fp):
        '''
            parameter
            ---------
            seq_dict: dict
                {chro: seq, }
            fp:         file handle for output sequence file
        '''
        if self._exon is None:
            logger.info('{}: -- not annot exon --'.format(self._name))
            return None
        exon_num = len(self._exon)
        for index, each in enumerate(self._exon):
            start, end = each
            order = index + 1 if self._strand == '+' else exon_num - index
            wanted = Sequence(
                '{}_exon{}'.format(self._name, order),
                seq_dict[self._chro][start:end],
                'gene_id:{}'.format(self._gene_id),
            )
            wanted = wanted.reverse_complement() if self._strand == '-' else wanted
            wanted.write_to_fasta_file(fp)

    def extract_intron_seq(self, seq_dict, fp):
        '''
            parameter
            ---------
            seq_dict: dict
                {chro: seq, }
            fp:         file handle for output sequence file
        '''
        if self._exon is None:
            logger.info('{}: -- not annot exon --'.format(self._name))
            return None
        intron_num = len(self._exon) - 1
        if intron_num >= 1:
            for index in range(intron_num):
                _, start = self._exon[index]
                end, _ = self._exon[index + 1]
                order = index + 1 if self._strand == '+' else intron_num - index
                wanted = Sequence(
                    '{}_intron{}'.format(self._name, order),
                    seq_dict[self._chro][start:end],
                    'gene_id:{}'.format(self._gene_id),
                )
                wanted = wanted.reverse_complement() if self._strand == '-' else wanted
                wanted.write_to_fasta_file(fp)

    def extract_UTR_seq(self, seq_dict, fp):
        '''
            parameter
            ---------
            seq_dict: dict
                {chro: seq, }
            fp:         file handle for output sequence file
        '''
        UTR = self.__UTR_identify()
        if UTR:
            for index, utr in enumerate(UTR):
                labels, start, end = utr
                wanted = Sequence(
                    '{}_{}_{}'.format(self._name, labels, index),
                    seq_dict[self._chro][start:end],
                    'gene_id:{}'.format(self._gene_id),
                )
                wanted = wanted.reverse_complement() if self._strand == '-' else wanted
                wanted.write_to_fasta_file(fp)
        else:
            logger.info('{}: -- not annot exon and CDS --'.format(self._name))
            return None

    # def structure_dat(self):
    #     '''
    #         specific function for plot gene structure
    #
    #         return
    #         ------
    #         list, 0-based postion
    #             [start, end, [[es,ee],], [[cs,ce],], strand]
    #     '''
    #     return [self._start, self._end, self._exon, self._cds, self._strand]
    #
    #
    # def gene_structure_for_itol(self):
    #     '''
    #     '''
    #     start, end, strand = self._start, self._end, self._strand
    #     if self._cds:
    #         cds = [['CDS', i[0], i[1]] for i in self._cds]
    #         utr = self.__UTR_identify()
    #         feature = sorted(cds+utr, key=lambda x: x[1], reverse=False)
    #     else:
    #         feature = sorted([['EXON', i[0], i[1]] for i in self._exon],
    #                          key=lambda x: x[1], reverse=False)
    #     length = end-start
    #     if strand == '+':
    #         feature = [[i[0], i[1]-start, i[2]-start] for i in feature]
    #     else:
    #         feature = feature[::-1]
    #         feature = [[i[0], end-i[2], end-i[1]] for i in feature]
    #     # SHAPE|START|END|COLOR|LABEL
    #     element = [self._name, str(length)]
    #     for i in feature:
    #         label, start, end = i
    #         color = '#ff9800' if label=='CDS' or label=='EXON' else '#0067ff'
    #         element.append('|'.join(['RE', str(start), str(end), color, label]))
    #     return '\t'.join(element)


class GTF_Reader(Files):
    '''
        File parser for GTF/GFF file of gene annotation
    '''

    def __init__(self, gtf, suffix=None):
        self.suffix = suffix
        self._transcript_feature = {'mRNA', 'transcript'}
        self._keep_feature = {
            'CDS',
            'exon',
            'transcript',
            'mRNA'
        }
        self._skip_feature = {
            'Selenocysteine',
            'start_codon',
            'stop_codon',
            'UTR',
            'gene',
            'five_prime_UTR',
            'three_prime_UTR',
            'five_prime_utr',
            'three_prime_utr',
            '3UTR',
            '5UTR',
        }
        Files.__init__(self, gtf)

    def __iter__(self):
        logger.info('Skip known annotation feature: ({})'.format(', '.join(self._skip_feature)))

        t_id, t_chro, t_gene = None, None, None
        t_exon, t_cds, t_info = [], [], {}
        for line in Files.__iter__(self):
            if line.startswith('#'):
                continue
            chro, _, feature, start, end, _, strand, _, attr = line.strip().split('\t')
            if feature not in self._keep_feature:
                if feature not in self._skip_feature:
                    logger.warning('skip novel annotation feature: {}'.format(feature))
                continue
            start, end = int(start) - 1, int(end)
            if '=' in attr:
                attr = [i.strip().partition('=') for i in attr.split(';') if i.strip()]
            else:
                attr = [i.strip().partition(' ') for i in attr.split(';') if i.strip()]
            attr = {i[0]: i[2].strip('\'\t\n\r"') for i in attr}
            try:
                line_id = attr['transcript_id']
            except KeyError:
                line_id = attr['ID'] if feature in self._transcript_feature else attr['Parent']

            logger.debug(line)
            if t_id and ((t_id != line_id) or ((t_id == line_id) and (t_chro != chro))):
                yield Transcript(t_id, t_chro, t_start, t_end, t_strand, t_exon, t_cds, t_info, self.suffix)
                t_exon, t_cds, t_info = [], [], {}

            if feature in self._transcript_feature:
                t_id = line_id
                t_chro, t_start, t_end, t_strand = chro, start, end, strand
                try:
                    t_info['gene_id'] = attr['gene_id']
                except KeyError:
                    t_info['gene_id'] = attr['Parent']
                # t_gene = t_info['gene_id']
                try:
                    t_info['transcript_type'] = attr['transcript_type']
                    t_info['gene_type'] = attr['gene_type']
                    t_info['transcript_name'] = attr['transcript_name']
                    t_info['gene_name'] = attr['gene_name']
                except KeyError:
                    pass
            elif feature == 'exon':
                t_exon.append([start, end])
            elif feature == 'CDS':
                t_cds.append([start, end])
        yield Transcript(t_id, t_chro, t_start, t_end, t_strand, t_exon, t_cds, t_info, self.suffix)

    def Saft_Model(self):
        '''
            parse GTF/GFF file and sava transcript structure data as python
            object to memory, Saft_Model methods is Applicable to unsorted files,
            the same transcript annotation line is no longer in the same block
        '''
        logger.info('Skip known annotation feature: ({})'.format(', '.join(self._skip_feature)))
        logger.info('Start Read gff file to python object...')

        transcriptlst, annot = [], {}
        for line in Files.__iter__(self):
            if line.startswith('#'):
                continue
            chro, _, feature, start, end, _, strand, _, attr = line.strip().split('\t')
            if feature not in self._keep_feature:
                if feature not in self._skip_feature:
                    logger.warning('skip novel annotation feature: {}'.format(feature))
                continue

            start, end = int(start) - 1, int(end)
            if '=' in attr:
                attr = [i.strip().partition('=') for i in attr.split(';') if i.strip()]
            else:
                attr = [i.strip().partition(' ') for i in attr.split(';') if i.strip()]
            attr = {i[0]: i[2].strip('\'\t\n\r"') for i in attr}
            try:
                t_id = attr['transcript_id']
            except KeyError:
                t_id = attr['ID'] if feature in self._transcript_feature else attr['Parent']
            transcriptlst.append(t_id)

            if feature == 'transcript' or feature == 'mRNA':
                annot.setdefault(t_id, {})['pos'] = [chro, start, end, strand]
                try:
                    annot[t_id]['gene_id'] = attr['gene_id']
                except KeyError:
                    annot[t_id]['gene_id'] = attr['Parent']
                try:
                    annot[t_id]['transcript_type'] = attr['transcript_type']
                    annot[t_id]['gene_type'] = attr['gene_type']
                    annot[t_id]['transcript_name'] = attr['transcript_name']
                    annot[t_id]['gene_name'] = attr['gene_name']
                except KeyError:
                    pass
            elif feature == 'exon':
                annot.setdefault(t_id, {}).setdefault('exon', []).append([start, end])
            elif feature == 'CDS':
                annot.setdefault(t_id, {}).setdefault('cds', []).append([start, end])

        logger.info('Done of parse gff, sort transcript list...')
        transcriptlst = sorted(
            set(transcriptlst),
            key=lambda x: transcriptlst.index(x),
            reverse=False
        )
        logger.info('Done of sort transcript list.')
        for iso in transcriptlst:
            chro, start, end, strand = annot[iso]['pos']
            exon = annot[iso].get('exon', None)
            cds = annot[iso].get('cds', None)
            if not any([exon, cds]):
                exon = [[start, end]]
            elif not exon:
                exon = cds.copy()
            t_info = {
                'gene_id': annot[iso].get('gene_id', None),
                'gene_name': annot[iso].get('gene_name', None),
                'gene_type': annot[iso].get('gene_type', None),
                'transcript_name': annot[iso].get('transcript_name', None),
                'transcript_type': annot[iso].get('transcript_type', None),
            }
            t_info = annot[iso]
            yield Transcript(iso, chro, start, end, strand, exon, cds, t_info, self.suffix)


class RefSeq_GFF_Reader(Files):
    '''
        File parser for GFF file of gene annotation from NCBI RefSeq or Genome database
    '''

    def __init__(self, gtf, chrom=None):
        Files.__init__(self, gtf)
        self._chrom = self.__convert_chrom_id(chrom) if chrom else {}

    def __iter__(self):
        genelst, annot = self.__parse_gff()
        for gene in genelst:
            for iso in annot[gene]:
                chro, start, end, strand = annot[gene][iso]['infor']
                chro = self._chrom.get(chro, chro)
                gene_name = annot[gene][iso]['gene_name']
                gene_type = annot[gene][iso]['gene_type']
                transcript_name = annot[gene][iso]['transcript_name']
                transcript_type = annot[gene][iso]['transcript_type']
                exon = annot[gene][iso].get('exon', None)
                cds = annot[gene][iso].get('CDS', None)
                if not any([exon, cds]):
                    exon = [[start, end]]
                elif not exon:
                    exon = cds.copy()
                t_info = {
                    'gene_id': gene,
                    'gene_name': gene_name,
                    'gene_type': gene_type,
                    'transcript_name': transcript_name,
                    'transcript_type': transcript_type,
                    'EntrezID': annot[gene][iso]['EntrezID']
                }
                yield Transcript(iso, chro, start, end, strand, exon, cds, t_info)

    def __parse_gff(self):
        SKIP = {
            'cDNA_match',
            'repeat_region',
            'D_loop',
            'enhancer',
            'match',
            'promoter',
            'region',
            'origin_of_replication',
            'centromere',
            'sequence_feature',
            'biological_region',    # Igl
            'sequence_alteration',
        }
        logger.info('Skip annotation feature: ({})'.format(', '.join(SKIP)))

        TRANSCRIPT = {
            'mRNA',
            'lnc_RNA',
            'ncRNA',
            'antisense_RNA',
            'transcript',  # non coding
            'RNase_MRP_RNA',
            'RNase_P_RNA',
            'Y_RNA',
            'tRNA',
            'rRNA',
            'snoRNA',
            'snRNA',
            'primary_transcript',  # miRNA precursor seq
            'miRNA',
            'C_gene_segment',
            'D_gene_segment',
            'V_gene_segment',
            'J_gene_segment',
            'SRP_RNA',
            'telomerase_RNA',
            'vault_RNA',
            'guide_RNA',
            'scRNA',
        }
        logger.info('the following annotation feature as transcript: ({})'.format(', '.join(TRANSCRIPT)))

        logger.info('Start Read gff file to python object...')
        pattern_ENTREZID = re.compile(r'GeneID:(\d+)')
        pattern_HGNC = re.compile(r'HGNC:(\d+)')
        geneidlst, gene_annot = [], {}
        transcript_annot, transcript_structure = {}, {}
        for line in Files.__iter__(self):
            if line.startswith('#'):
                continue
            chro, _, feature, start, end, _, strand, _, attr = line.strip().split('\t')
            if feature in SKIP:
                continue

            start, end = int(start) - 1, int(end)
            attr = [i.strip().partition('=') for i in attr.split(';') if i.strip()]
            attr = {i[0]: i[2].strip('\'\t\n\r"') for i in attr}

            try:
                ENTREZID = pattern_ENTREZID.search(attr['Dbxref']).group(1)
            except:
                ENTREZID = None
            try:
                HGNCID = pattern_HGNC.search(attr['Dbxref']).group(1)
            except:
                HGNCID = None
            if feature == 'gene' or feature == 'pseudogene':
                line_id = attr['ID']
                geneidlst.append(line_id)
                try:
                    gchro, _, _, _ = gene_annot[line_id]['infor']
                    if gchro == chro:
                        logger.error('Find duplicate gene({}) annotation line: {}'.format(line_id, line))
                        sys.exit(1)
                    else:
                        logger.warning('Duplicate gene_id({}) on different chromosome'.format(line_id))
                except KeyError:
                    pass
                gene_annot.setdefault(line_id, {})['infor'] = [chro, start, end, strand]
                gene_annot[line_id]['gene_name'] = attr['Name']
                gene_annot[line_id]['gene_type'] = attr['gene_biotype']
                gene_annot[line_id]['EntrezID'] = ENTREZID
                gene_annot[line_id]['HGNC'] = HGNCID

            elif feature in TRANSCRIPT:
                g_id = attr['Parent']
                t_id = attr['ID']
                geneidlst.append(g_id)
                transcript_annot.setdefault(g_id, {}).setdefault(t_id, {})['infor'] = [chro, start, end, strand]
                if feature == 'miRNA':
                    transcript_id = attr['product']
                    transcript_name = attr['product']
                    transcript_type = 'miRNA'
                elif feature == 'tRNA':
                    transcript_id = attr['ID']
                    transcript_name = attr['ID']
                    transcript_type = 'tRNA'
                else:
                    transcript_id = attr.get('transcript_id', t_id)
                    transcript_name = attr.get('Name', transcript_id)
                    tmp = attr['gbkey']
                    transcript_type = 'protein_coding' if tmp == 'mRNA' else tmp
                transcript_annot[g_id][t_id]['transcript_id'] = transcript_id
                transcript_annot[g_id][t_id]['transcript_name'] = transcript_name
                transcript_annot[g_id][t_id]['transcript_type'] = transcript_type
                transcript_annot[g_id][t_id]['EntrezID'] = ENTREZID
                transcript_annot[g_id][t_id]['HGNC'] = HGNCID
            elif feature == 'exon':
                t_id = attr['Parent']
                transcript_structure.setdefault(t_id, {}).setdefault('exon', []).append([start, end])
            elif feature == 'CDS':
                t_id = attr['Parent']
                transcript_structure.setdefault(t_id, {}).setdefault('CDS', []).append([start, end])
            else:
                logger.warning('Skip A novel annotation feature: {}'.format(feature))

        logger.info('Done of Read gff, sort transcript list...')
        geneidlst = sorted(set(geneidlst), key=lambda x: geneidlst.index(x), reverse=False)
        logger.info('Done of sort transcript list, convert raw annot infor to appropriate transcript structure annotation...')
        annot = {}
        for index, gene in enumerate(geneidlst):
            if index % 5000 == 0 and index != 0:
                logger.info('Already processed gene NO. : {}'.format(index))

            infor_gene = gene_annot.get(gene, None)
            infor_transcript = transcript_annot.get(gene, None)
            flag_gene = True
            if infor_gene is None:
                flag_gene = False

            if infor_transcript is None:
                logger.debug('Missing transcript feature line for gene: {}'.format(gene))
                chro, start, end, strand = infor_gene['infor']
                gene_name = infor_gene['gene_name']
                gene_type = infor_gene['gene_type']
                annot.setdefault(gene, {}).setdefault(gene, {})['infor'] = [chro, start, end, strand]
                annot[gene][gene]['gene_name'] = gene_name
                annot[gene][gene]['gene_type'] = gene_type
                annot[gene][gene]['EntrezID'] = infor_gene['EntrezID']
                annot[gene][gene]['HGNC'] = infor_gene['HGNC']
                annot[gene][gene]['transcript_name'] = gene_name
                annot[gene][gene]['transcript_type'] = gene_type
                # exon, CDS
                try:
                    exon = transcript_structure[gene]['exon']
                except KeyError:
                    exon = None
                try:
                    cds = transcript_structure[gene]['CDS']
                except KeyError:
                    cds = None
                if not any([exon, cds]):
                    logger.debug('Missing exon, CDS feature line for gene: {}'.format(gene))
                    annot[gene][gene]['exon'] = [[start, end], ]
                else:
                    annot[gene][gene]['exon'] = exon
                    annot[gene][gene]['CDS'] = cds
            else:
                for iso in infor_transcript:
                    chro, start, end, strand = infor_transcript[iso]['infor']
                    transcript_id = infor_transcript[iso]['transcript_id']
                    if flag_gene:
                        EntrezID = infor_gene['EntrezID']
                        HGNC = infor_gene['HGNC']
                        gene_name = infor_gene['gene_name']
                        gene_type = infor_gene['gene_type']
                        G_chro, G_start, G_end, G_strand = infor_gene['infor']
                    else:
                        EntrezID = infor_transcript[iso]['EntrezID']
                        HGNC = infor_transcript[iso]['HGNC']
                        msg = ('Missing gene feature annotation for transcript({}), '
                               'use transcript infor as gene infor').format(gene)
                        logger.warning(msg)
                        gene_name = infor_transcript[iso]['transcript_name']
                        gene_type = infor_transcript[iso]['transcript_type']
                        G_chro, G_start, G_end, G_strand = infor_transcript[iso]['infor']
                    msg = 'transcript({}) region over range of gene({}) region'.format(transcript_id, gene)
                    assert G_start <= start < end <= G_end, msg
                    msg = 'chro is inconsistent of transcript({}) and gene({})'.format(transcript_id, gene)
                    assert chro == G_chro, msg
                    msg = ('Transcription direction(strand) is inconsistent of '
                           'transcript({}) and gene({})').format(transcript_id, gene)
                    assert strand == G_strand, msg

                    transcript_type = infor_transcript[iso]['transcript_type']
                    if gene_type == 'other':
                        gene_type = transcript_type
                    annot.setdefault(gene, {}).setdefault(transcript_id, {})['infor'] = [chro, start, end, strand]
                    annot[gene][transcript_id]['gene_name'] = gene_name
                    annot[gene][transcript_id]['gene_type'] = gene_type
                    annot[gene][transcript_id]['transcript_name'] = infor_transcript[iso]['transcript_name']
                    annot[gene][transcript_id]['transcript_type'] = transcript_type
                    annot[gene][transcript_id]['EntrezID'] = EntrezID
                    annot[gene][transcript_id]['HGNC'] = HGNC
                    # CDS, exon
                    try:
                        exon = transcript_structure[iso]['exon']
                    except KeyError:
                        exon = None
                    try:
                        cds = transcript_structure[iso]['CDS']
                    except KeyError:
                        cds = None
                    if not any([exon, cds]):
                        logger.info('Missing exon, CDS feature line for transcript: {}'.format(transcript_id))
                        annot[gene][transcript_id]['exon'] = [[start, end], ]
                    else:
                        annot[gene][transcript_id]['exon'] = exon
                        annot[gene][transcript_id]['CDS'] = cds
        logger.info('Done of parse gff.')
        return geneidlst, annot

    def __convert_chrom_id(self, tab):
        with open(tab) as f:
            tab = [i.strip().split()[:2] for i in f if not i.startswith('#')]
        return {i[0]: i[1] for i in tab}


class Bed_Reader(Files):
    '''
    '''
    def __init__(self, bed):
        Files.__init__(self, bed)

    def __iter__(self):
        for line in Files.__iter__(self):
            chro, start, end, name, _, strand, cstart, cend, _, _, elen, estart = line.strip().split('\t')
            start, end, cstart, cend = map(int, [start, end, cstart, cend])
            elen = [int(i) for i in elen.split(',') if i]
            estart = [int(i) + start for i in estart.split(',') if i]
            eend = [estart[i] + elen[i] for i, j in enumerate(elen)]
            EXON = [[estart[i], eend[i]] for i, j in enumerate(estart)]
            if cstart == cend:
                CDS = []
            else:
                CDS = [i.copy() for i in EXON if i[1] > cstart]
                CDS = [i for i in CDS if i[0] < cend]
                CDS[0][0] = cstart
                CDS[-1][1] = cend
            yield Transcript(name, chro, start, end, strand, EXON, CDS)


# if __name__ == '__main__':
#     _, gff, output = sys.argv
#     with open(output, 'w') as f:
#         for i in GTF_Reader(gff):
#             i.to_gtf(f)
