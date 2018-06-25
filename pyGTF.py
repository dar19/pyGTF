#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import bz2
import gzip
from collections import defaultdict

'''
High Frequency Function of Files parse
    Sequence File:
        Fastq, Fasta
      output:
        Fastq, Fasta

    Gene Annotation File:
        GTF, GFF
      output:
        GTF, BED, GenePred

    update data: 20180427
    Bug report: chengchaoze@163.com
'''

class Files(object):
    '''
    '''
    def __init__(self, File):
        self.fos = File

    def __iter__(self):
        if self.fos.lower().endswith((".gz", ".gzip")):
            lines = gzip.open(self.fos)
        elif self.fos.lower().endswith(".bz2"):
            lines = bz2.BZ2File(self.fos)
        else:
            lines = open(self.fos)
        for line in lines:
            try:
            #if isinstance(line, bytes):
                line = line.decode('utf-8')
            except:
                pass
            yield line
        lines.close()


class Sequence(object):
    '''
    '''
    def __init__(self, name, seq, descr=''):
        self.name = name
        self.seq = seq
        self.descr = descr

    def reverse_complement(self):
        paired_rc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                     'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
        tmp = Sequence(self.name,
                       ''.join([paired_rc[i] for i in self.seq[::-1]]),
                       self.descr)
        return tmp

    def write_to_fasta_file(self, fp, chara=80):
        if self.descr:
            fp.write('>{} {}\n'.format(self.name, self.descr))
        else:
            fp.write('>{}\n'.format(self.name))
        fp.write('{}'.format(self.__class__.__seq_formater(self.seq, chara)))

    def write_to_tab_file(self, fp,):
        fp.write('{}\t{}\n'.format(self.name, self.seq))

    @classmethod
    def __seq_formater(cls, seq, length):
        tmp = ''
        for i in range(0, len(seq), length):
            tmp += (seq[i:(i+length)]+'\n')
        return tmp


class SequenceWithQual(Sequence):
    '''
    '''
    def __init__(self, name, seq, qual):
        if len(seq) != len(qual):
            raise Exception('length is Inconsistent seq and qual string')
        self.name = name
        self.seq = seq
        # super().__init__(name, seq)
        # Sequence.__init__(self, name, seq)
        self.qual = qual

    def write_to_fastq_file(self, fp):
        fp.write('@{}\n'.format(self.name))
        fp.write('{}\n'.format(self.seq))
        fp.write('+\n{}\n'.format(self.qual))

    def write_to_fasta_file(self, fp):
        fp.write('>{}\n'.format(self.name))
        fp.write('{}\n'.format(self.seq))


class Fasta_Reader(Files):
    '''
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
                assert seq is not None, "FASTA file does not start with '>'"
                seq += line.strip()
        yield Sequence(seqid, seq, descr)


class Fastq_Reader(Files):
    '''
    '''
    def __init__(self, fastq):
        self.fastq = fastq
        Files.__init__(self, fastq)

    def __iter__(self):
        iter = Files.__iter__(self)
        try:
            while True:
                seqid = next(iter)
                try:
                    seq, _, qual = next(iter), next(iter), next(iter)
                except StopIteration:
                    print('incompelete Fastq file: {}'.format(self.fastq))
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
    transcript annotation information
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

    [option argument]: dict
        transcript_type:   string
        gene_id:   string
        gene_type:   string
    '''
    def __init__(self, Tid, chro, start, end, strand, exon, cds, infor, suffix=''):
        self.name = self.__class__.__del_suffix(Tid, suffix) if suffix else Tid
        self.chro = chro
        self.start = start
        self.end = end
        self.strand = strand
        assert any([exon, cds]), 'missing "exon" and "CDS" feature of transcript "{}"'.format(Tid)
        self.__exon = sorted(exon, key=lambda x: x[0], reverse=False) if exon else None
        self.__cds = sorted(cds, key=lambda x: x[0], reverse=False) if cds else None
        Gid = infor.get('gene_id', Tid)
        self.gene_id = self.__class__.__del_suffix(Gid, suffix) if suffix else Gid
        self.__type = infor.get('transcript_type', None)    # 'protein_coding'
        self.__Gtype = infor.get('gene_type', None)
        self.__Tname = infor.get('transcript_name', None)
        self.__Gname = infor.get('gene_name', None)
        self.__valid_data()

    @classmethod
    def __del_suffix(cls, string, suffix):
        suffix_len = len(suffix)
        assert len(string) > suffix_len, 'Gene suffix is illegal'
        if suffix:
            if suffix == string[-suffix_len:]:
                return string[:-suffix_len]
        return string
        #loc = string.find(suffix)
        #if loc == -1:
        #    return string
        #return string[:loc]

    def __valid_data(self):
        assert self.strand in ['-', '+'], 'unsupported strand Symbol, {-, +}'
        if self.__exon:
            start, end = self.__exon[0][0], self.__exon[-1][1]
            info = ('Error: Over range of the transcriptional region.\n',
                    '    location:  {}:{}-{}\n'.format(self.chro, self.start, self.end),
                    '    transcriptional region:  {}-{}\n'.format(start, end))
            assert self.start == start < end == self.end, info
            for each in self.__exon:
                assert each[1]>each[0], 'End loc must more than Start of Exon.'
        if self.__cds:
            start, end = self.__cds[0][0], self.__cds[-1][1]
            info = ('Error: Over range of the transcriptional region.\n',
                    '    location:  {}:{}-{}\n'.format(self.chro, self.start, self.end),
                    '    coding sequence region:  {}-{}\n'.format(start, end))
            assert self.start <= start < end <= self.end, info
            for each in self.__cds:
                assert each[1]>each[0], 'End loc must more than Start of CDS.'

    def basic_info(self):
        '''
        return list object, 1-based postion
            [name, chro, start, end, strand]
            value attribute: [str, str, int, int, str]
        '''
        return [self.name, self.chro, self.start+1, self.end, self.strand]
        # return self.__class__.__list2str([self.name, self.chro, self.start, self.end, self.strand])

    def to_gtf(self, fp):
        '''
        fp:    file handle for output standrand GTF file
        '''
        type_t = 'transcript_type "{}"; '.format(self.__type) if self.__type else ''
        type_g = 'gene_type "{}"; '.format(self.__Gtype) if self.__Gtype else ''
        name_t = 'transcript_name "{}"; '.format(self.__Tname) if self.__Tname else ''
        name_g = 'gene_name "{}"; '.format(self.__Gname) if self.__Gname else ''
        identifier = 'transcript_id "{}"; gene_id "{}"; '.format(self.name, self.gene_id)
        attr = ''.join([identifier, type_t, type_g, name_t, name_g])

        transcript = [self.chro, '.', 'transcript', self.start+1, self.end,
                      '.', self.strand, '.', attr]
        fp.write(self.__class__.__list2str(transcript))

        Feature = []
        if self.__exon:
            for i, each in enumerate(self.__exon):
                start, end = each
                Feature.append([self.chro, '.', 'exon', start+1, end,
                                '.', self.strand, '.', attr])
        if self.__cds:
            for i, each in enumerate(self.__cds):
                start, end = each
                Feature.append([self.chro, '.', 'CDS', start+1, end,
                                '.', self.strand, '.', attr])
        UTR = self.__UTR_identify()
        for utr in UTR:
            labels, start, end = utr
            Feature.append([self.chro, '.', labels, start+1, end,
                            '.', self.strand, '.', attr])
        order = False if self.strand == '+' else True
        Feature = sorted(Feature, key=lambda x: x[3], reverse=order)
        for each in Feature:
            fp.write(self.__class__.__list2str(each))

    def to_bed(self, fp):
        '''
        fp:    file handle for output 12 columns bed file
        '''
        if self.__cds:
            cstart = self.__cds[0][0]
            cend = self.__cds[-1][1]
        else:
            cstart, cend = self.start, self.end
        EXON = self.__exon if self.__exon else self.__cds
        exon_num = len(EXON)
        exon_len = ''.join(['{},'.format(x[1]-x[0]) for x in EXON])
        exon_start = ''.join(['{},'.format(x[0]-self.start) for x in EXON])

        transcript = [self.chro, self.start, self.end, self.name, 1000, self.strand,
                      cstart, cend, 0, exon_num, exon_len, exon_start]
        fp.write(self.__class__.__list2str(transcript))

    def to_genePred(self, fp, CIRCexplorer2=False):
        '''
        fp:     file handle for output GenePred file,
                ref: https://genome.ucsc.edu/FAQ/FAQformat#format9
        CIRCexplorer2:   bool
                whether create CIRCexplorer2 style genePred, {True, False}
        '''
        if self.__cds:
            cstart = self.__cds[0][0]
            cend = self.__cds[-1][1]
        else:
            cstart, cend = self.start, self.end
        EXON = self.__exon if self.__exon else self.__cds
        exon_num = len(EXON)
        estart = ''.join(['{},'.format(i[0]) for i in EXON])
        eend = ''.join(['{},'.format(i[1]) for i in EXON])

        transcript = [self.gene_id, ] if CIRCexplorer2 else []
        transcript += [self.name, self.chro, self.strand, self.start, self.end,
                       cstart, cend, exon_num, estart, eend]
        fp.write(self.__class__.__list2str(transcript))

    @classmethod
    def __list2str(cls, lst):
        return '{}\n'.format('\t'.join([str(i) for i in lst]))

    def length_transcript(self):
        EXON = self.__exon if self.__exon else self.__cds
        return sum([i[1]-i[0] for i in EXON])

    def length(self):
        return self.end-self.start

    def exon_num(self):
        EXON = self.__exon if self.__exon else self.__cds
        return len(EXON)

    def extract_genomic_seq(self, seq_dict, fp):
        '''
        seq_dict:    ref dict object
        fp:         file handle for output sequence file
        '''
        seq = seq_dict[self.chro][self.start:self.end]
        wanted = Sequence('{}_genomic'.format(self.name), seq, 'gene_id:{}'.format(self.gene_id))
        if self.strand == '-':
            wanted = wanted.reverse_complement()
        wanted.write_to_fasta_file(fp)
        # return wanted

    def extract_isoform_seq(self, seq_dict, fp):
        '''
        seq_dict:    ref dict object
        fp:         file handle for output sequence file
        '''
        if self.__exon is None:
            sys.stderr.write('{}: -- not annot exon --\n'.format(self.name))
            return None
            # return Sequence('{}_transcript'.format(self.name), '-- not annot exon --')
        seqtmp = ''
        for each in self.__exon:
            start, end = each
            seqtmp += seq_dict[self.chro][start:end]
        wanted = Sequence(self.name, seqtmp, 'gene_id:{}'.format(self.gene_id))
        if self.strand == '-':
            wanted = wanted.reverse_complement()
        wanted.write_to_fasta_file(fp)
        # return wanted

    def extract_CDS_seq(self, seq_dict, fp):
        '''
        seq_dict:    ref dict object
        fp:         file handle for output sequence file
        '''
        if self.__cds is None:
            sys.stderr.write('{}: -- not annot CDS --\n'.format(self.name))
            return None
            # return Sequence('{}_coding_sequence'.format(self.name), '-- not annot CDS --')
        seqtmp = ''
        for each in self.__cds:
            start, end = each
            seqtmp += seq_dict[self.chro][start:end]
        wanted = Sequence(self.name, seqtmp, 'gene_id:{}'.format(self.gene_id))
        if self.strand == '-':
            wanted = wanted.reverse_complement()
        wanted.write_to_fasta_file(fp)
        # return wanted

    def extract_exon_seq(self, seq_dict, fp):
        '''
        seq_dict:    ref dict object
        fp:         file handle for output sequence file
        '''
        if self.__exon is None:
            sys.stderr.write('{}: -- not annot exon --\n'.format(self.name))
            return None
            # return [Sequence('{}_exon'.format(self.name), '-- not annot exon --'),]
        # seqlist = []
        exon_num = len(self.__exon)
        for index, each in enumerate(self.__exon):
            start, end = each
            order = index+1 if self.strand == '+' else exon_num-index
            wanted = Sequence('{}_exon{}'.format(self.name, order),
                              seq_dict[self.chro][start:end],
                              'gene_id:{}'.format(self.gene_id))
            # seqlist.append(wanted)
            wanted = wanted.reverse_complement() if self.strand == '-' else wanted
            wanted.write_to_fasta_file(fp)
        # if self.strand == '-':
        #     seqlist = [i.reverse_complement() for i in seqlist]
        # return seqlist

    def extract_intron_seq(self, seq_dict, fp):
        '''
        seq_dict:    ref dict object
        fp:         file handle for output sequence file
        '''
        if self.__exon is None:
            sys.stderr.write('{}: -- not annot exon --\n'.format(self.name))
            return None
            # return [Sequence('{}_intron'.format(self.name), '-- not annot exon --'), ]
        # seqlist = []
        intron_num = len(self.__exon)-1
        if intron_num >= 1:
            for index in range(intron_num):
            #for index, each in enumerate(self.__exon):
                _, start = self.__exon[index]
                end, _ = self.__exon[index+1]
                order = index+1 if self.strand == '+' else intron_num-index
                wanted = Sequence('{}_intron{}'.format(self.name, order),
                                  seq_dict[self.chro][start:end],
                                 'gene_id:{}'.format(self.gene_id))
                # seqlist.append(wanted)
                wanted = wanted.reverse_complement() if self.strand == '-' else wanted
                wanted.write_to_fasta_file(fp)
            # if self.strand == '-':
            #     seqlist = [i.reverse_complement() for i in seqlist]
        # return seqlist

    def extract_UTR_seq(self, seq_dict, fp):
        '''
        seq_dict:    ref dict object
        fp:         file handle for output sequence file
        '''
        UTR = self.__UTR_identify()
        if UTR:
            # seqlist = []
            for index,utr in enumerate(UTR):
                labels, start, end = utr
                # seqlist.append(Sequence('{}_{}_{}'.format(self.name, labels, index),
                #                seq_dict[self.chro][start:end]))
                wanted = Sequence('{}_{}_{}'.format(self.name, labels, index),
                                  seq_dict[self.chro][start:end],
                                  'gene_id:{}'.format(self.gene_id))
                wanted = wanted.reverse_complement() if self.strand == '-' else wanted
                wanted.write_to_fasta_file(fp)
            # if self.strand == '-':
            #     seqlist = [i.reverse_complement() for i in seqlist]
            # return seqlist
        else:
            sys.stderr.write('{}: -- not annot exon and CDS --\n'.format(self.name))
            return None
            # return [Sequence('{}_UTR'.format(self.name), '-- not annot exon and CDS --'), ]

    def __UTR_identify(self):
        utr = []
        if self.__cds and self.__exon and (self.__cds != self.__exon):
            coding_start, coding_end = self.__cds[0][0], self.__cds[-1][1]
            pos_s = '5UTR' if self.strand == '+' else '3UTR'
            pos_l = '3UTR' if self.strand == '+' else '5UTR'
            for each in self.__exon:
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

    def structure_dat(self):
        '''
        specific function for plot gene structure
            return list object, 0-based postion
                [start, end, [[es,ee],], [[cs,ce],], strand]
                value attribute: [int, int, list, list, str]
        '''
        return [self.start, self.end, self.__exon, self.__cds, self.strand]


class GTF_Reader(Files):
    '''
    '''
    def __init__(self, gtf, suffix=''):
        self.suffix = suffix
        Files.__init__(self, gtf)

    def __iter__(self):
        keep = {'CDS', 'exon', 'transcript', 'mRNA'}
        skip = set(['Selenocysteine', 'start_codon', 'stop_codon', 'UTR',
                    'gene', 'five_prime_UTR', 'three_prime_UTR'])
        t_id, t_exon, t_cds, t_info = None, [], [], {}
        for line in Files.__iter__(self):
            if line.startswith( "#" ):
                continue
            chro, _, feature, start, end, _, strand, _, attr = line.strip().split('\t')
            if feature not in keep:
                continue
            start, end = int(start)-1, int(end)
            if '=' in attr:
                attr = [i.strip().partition('=') for i in attr.split(';') if i.strip()]
            else:
                attr = [i.strip().partition(' ') for i in attr.split(';') if i.strip()]
            attr = {i[0]:i[2].strip('"\t\n\r\'') for i in attr}
            try:
                line_id = attr['transcript_id']
            except KeyError:
                line_id = attr['ID'] if feature in ['transcript', 'mRNA'] else attr['Parent']

            if t_id and t_id != line_id:
                yield Transcript(t_id, t_chro, t_start, t_end, t_strand, t_exon, t_cds, t_info, self.suffix)
                t_exon, t_cds, t_info = [], [], {}

            if feature == 'transcript' or feature == 'mRNA':
                t_id = line_id
                t_chro, t_start, t_end, t_strand = chro, start, end, strand
                try:
                    t_info['gene_id'] = attr['gene_id']
                except KeyError:
                    t_info['gene_id'] = attr['Parent']
                try:
                    t_info['transcript_type'] = attr['transcript_type']
                    t_info['gene_type'] = attr['gene_type']
                    t_info['transcript_name'] = attr['transcript_name']
                    t_info['gene_name'] = attr['gene_name']
                except:
                    pass
            elif feature == 'exon':
                t_exon.append([start, end])
            elif feature == 'CDS':
                t_cds.append([start, end])
        yield Transcript(t_id, t_chro, t_start, t_end, t_strand, t_exon, t_cds, t_info, self.suffix)


# class GTF_Reader_Saft(Files):
#     '''
#     safe read GFF file
#     '''
#     def __init__(self, gtf, suffix=''):
#         self.suffix = suffix
#         Files.__init__(self, gtf)
#
#         keep = {'CDS', 'exon', 'transcript', 'mRNA'}
#         isoform = defaultdict(dict)
#         for line in Files.__iter__(self):
#             if line.startswith('#'):
#                 continue
#             chro, _, feature, start, end, _, strand, _, attr = line.strip().split('\t')
#             if feature not in keep:
#                 continue
#             start, end = int(start)-1, int(end)
#             if '=' in attr:
#                 attr = [i.strip().partition('=') for i in attr.split(';') if i.strip()]
#             else:
#                 attr = [i.strip().partition(' ') for i in attr.split(';') if i.strip()]
#             attr = {i[0]:i[2].strip('"\t\n\r\'') for i in attr}
#             try:
#                 t_id = attr['transcript_id']
#             except KeyError:
#                 t_id = attr['ID'] if feature in ['transcript', 'mRNA'] else attr['Parent']
#
#             if feature == 'transcript' or feature == 'mRNA':
#                 isoform[t_id]['pos'] = [chro, start, end, strand]
#                 try:
#                     isoform[t_id]['gene_id'] = attr['gene_id']
#                 except KeyError:
#                     isoform[t_id]['gene_id'] = attr['Parent']
#                 try:
#                     isoform[t_id]['transcript_type'] = attr['transcript_type']
#                     isoform[t_id]['gene_type'] = attr['gene_type']
#                     isoform[t_id]['transcript_name'] = attr['transcript_name']
#                     isoform[t_id]['gene_name'] = attr['gene_name']
#                 except:
#                     pass
#             elif feature == 'exon':
#                 isoform[t_id]['exon'] = isoform[t_id].get(t_id, []) + [[start, end], ]
#             elif feature == 'CDS':
#                 isoform[t_id]['cds'] = isoform[t_id].get(t_id, []) + [[start, end], ]
#         self.annot = isoform
#
#     def __iter__(self):
#         for t_id in self.annot:
#             chro, start, end, strand = self.annot[t_id]['pos']
#             exon = self.annot[t_id].get('exon', None)
#             cds = self.annot[t_id].get('cds', None)
#             t_type = self.annot[t_id].get('transcript_type', None)
#             g_id = self.annot[t_id].get('gene_id', None)
#             g_type = self.annot[t_id].get('gene_type', None)
#             t_info = self.annot[t_id]
#             yield Transcript(t_id, chro, start, end, strand, exon, cds, t_info, self.suffix)


class RefSeq_GFF_Reader(Files):
    '''
    '''
    def __init__(self, gtf):
        Files.__init__(self, gtf)

    def __iter__(self):
        SKIP = {'region', 'cDNA_match', 'repeat_region', 'D_loop',
                'promoter', 'enhancer', 'miRNA',
                'match', 'origin_of_replication', 'centromere', 'sequence_feature'}

        TRANSCRIPT = {'mRNA',
                      'lnc_RNA', 'ncRNA', 'antisense_RNA', 'transcript',    # non coding
                      'RNase_MRP_RNA', 'RNase_P_RNA', 'Y_RNA',
                      'tRNA', 'rRNA', 'snoRNA', 'snRNA',
                      'primary_transcript',    # miRNA precursor seq
                      'C_gene_segment', 'D_gene_segment', 'V_gene_segment', 'J_gene_segment',
                      'SRP_RNA', 'telomerase_RNA', 'vault_RNA'}

        t_id, t_exon, t_cds, t_info = None, [], [], {}
        for line in Files.__iter__(self):
            if line.startswith( "#" ):
                continue
            chro, _, feature, start, end, _, strand, _, attr = line.strip().split('\t')
            for feature in SKIP:
                continue

            start, end = int(start)-1, int(end)
            attr = [i.strip().partition('=') for i in attr.split(';') if i.strip()]
            attr = {i[0]: i[2].strip('"\t\n\r\'') for i in attr}
            if feature in SKIP:
                continue

            line_id = attr['ID'] if feature in TRANSCRIPT else attr['Parent']

            if t_id and t_id != line_id:
                yield Transcript(t_id, t_chro, t_start, t_end, t_strand, t_exon, t_cds, t_info, self.suffix)
                t_exon, t_cds, t_info = [], [], {}


            if feature == 'gene':
                pass
            elif feature in TRANSCRIPT:
                pass
            elif feature == 'exon':
                pass
            elif feature == 'CDS':
                pass
            else:
                print('Feature:({}) is ignore, please check line:\n    {}'.format(feature, line))
        yield Transcript(t_id, t_chro, t_start, t_end, t_strand, t_exon, t_cds, t_info, self.suffix)
