'''
Created on Jul 22, 2014

@author: Kunal Bhutani
'''
import pandas as pd

class VariantFile(object):

    index = ['CHROM', 'POS', 'REF', 'ALT']
    features = {'allele': ['HQ', 'EHQ', 'CGA_CEHQ', 'AD'],
                'binary': ['CGA_XR', 'CGA_RPT', 'multiallele'],
                'categorical': ['FT', 'vartype1', 'vartype2', 'phase', 'zygosity'],
                'numeric': ['FT', 'vartype1', 'vartype2', 'phase', 'zygosity']
                }

    all_features = [feature
                    for subfeatures in features.itervalues()
                    for feature in subfeatures]

    def __init__(self):
        self._df = None

    @property
    def df(self):
        return self._df

    @property
    def features(self):

        return self._df[features].values

    @staticmethod
    def _get_phase(series, gt_column):
        '''
        
        :param series:
        :type series:
        :param gt_column:
        :type gt_column:
        '''
        '''
        Returns phase status from a GT field
        '''
        if "|" in series[gt_column]:
            return "|"
        if "/" in series[gt_column]:
            return "/"
        else:
            return '-'

    @staticmethod
    def _get_allele(series, sample_column, phase_column,
                    chromosome_num=0):
        '''
        Using the field info field in a sample column and the phasing
        field, parses out the allele from the first or second chromosome.
        
        :param series:
        :type series:
        :param sample_column:
        :type sample_column:
        :param phase_column:
        :type phase_column:
        :param chromosome_num:
        :type chromosome_num:
        '''

        # Assumes first info field is GT
        gt_split = series[phase_column]
        allele_calls = series[sample_column].split(":")[0].split(gt_split)
        alleles = [series['REF']] + series['ALT'].split(",")
        try:
            allele_number = int(allele_calls[chromosome_num])
            allele = alleles[allele_number]
        except:
            allele = "."
        return allele

    @staticmethod
    def _get_zygosity(series):
        '''
        Based on the alleles present in sample, determines the zygosity
        of the individual in comparison to the reference.
        
        :param series:
        :type series:
        '''

        allele_set = set([series['a1'], series['a2'], series['REF']])
        num_of_alleles = len(allele_set)

        if num_of_alleles == 1:
            return 'hom-ref'

        elif (num_of_alleles == 2  and
            series['a1'] == series['a2']):
            return 'hom-alt'

        elif (num_of_alleles == 2 and
              series['a1'] != series['a2']):
            return 'het-ref'

        elif num_of_alleles == 3 :
            return 'het-alt'

        else:
            assert False

    @staticmethod
    def _get_vartype(series, allele_base_col):
        '''
        This function assigns the following vartypes to the 
        allele specified by allele_base_col: snp, mnp, ins, del, indel or SV
        
        '''
        ref = series['REF']
        alt = series[allele_base_col]
        if is_snp(ref, alt): return 'snp'
        elif is_sv(ref, alt): return 'sv'
        elif is_mnp(ref, alt): return 'mnp'
        elif is_insertion(ref, alt): return 'ins'
        elif is_del(ref, alt): return 'del'
        else: return 'indel or SV'

        def is_snp(ref, alt):
            """ Return whether or not the variant is a SNP """
            if len(ref) != len(alt): return False
            base_diff = [nt for i, nt in enumerate(alt) if ref[i] != alt[i]]
            if len(base_diff) <= 1:
                    return True
            return False

        def is_sv(ref, alt):
            if "[" in alt or "]" in alt:
                return True
            return False

        def is_mnp(ref, alt):
            """ Return whether or not the variant is an MNP"""
            if len(ref) > 1 and len(ref) == len(alt):
                return True
            return False

        def is_insertion(ref, alt):
            if len(ref) < len(alt):
                return True
            return False

        def is_del(ref, alt):
            if len(ref) > len(alt):
                return True
            return False


class VcfTsv(VariantFile):
    def __init__(self, filename):
        self._df = self._load_vcf_tsv(filename)
        self._process()
        return

    def _process(self):
        '''
        Main processing function. Initial processing steps that might include feature generation
        or alternate representation of the dataframe.
        
        '''

        df = self._df

        # Overall processing functions
        self._df['variant_type'] = self._df.apply(self._variant_type, axis=1)

        # Ignore CNV and Repeats
        discard_vars = ['<INS:ME:ALU>', '<INS:ME:L1>', '<CGA_CNVWIN>', '<INS:ME:SVA>']
        df = df[~df.ALT.isin(discard_vars)]

        # Phase genotype and remove edge cases, likely at sex chromosome sites
        df['phase'] = df.apply(self._get_phase, args=['GT'], axis=1)
        df = df[df.phase != "-"]

        # Get allele base sequences
        # TODO: Combine phasing and getting alleles into one step
        # i.e. df['a1'], df['a2'] = df['GT'].map(lambda x: x.replace('|','/').split('/'))
        df['a1'] = df.apply(self._get_allele, args=['GT', 'phase', 0], axis=1)
        df['a2'] = df.apply(self._get_allele, args=['GT', 'phase', 1], axis=1)

        # Multiallelic site if more than one alternative allele
        df['multiallele'] = df['ALT'].map(lambda x: 1 if "," in x else 0)

        # Zygosity: hom-ref, het-ref, hom-alt, het-alt
        df['zygosity'] = df.apply(self._get_zygosity, axis=1)

        # Set the variant type based on number of bases in reference versus allele
        df['vartype1'] = df.apply(self._get_vartype, args=['a1'], axis=1)
        df['vartype2'] = df.apply(self._get_vartype, args=['a2'], axis=1)

        # Ignore structural variants
        df = df[(df.vartype1 != 'sv') & (df.vartype2 != 'sv')]

        # Fill missing allelic depth with 0,0 reads
        df['AD'].fillna(value='0,0', inplace=True)

        # Get number of reads supporting the reference allele
        df['ref_read_depth'] = df.apply(ref_read_depth, axis=1)

        # Generate features
        self._generate_features()

        return

    def _generate_features(self):
        # Generate features
        self._df['DP'] = 0


    def _load_vcf(self, filename):
        '''
        Loads in a vcf file, aware of gzipped files.
        '''
        compression = None
        if filename.endswith('.gz'):
            compression = 'gzip'

        df = pd.read_table(filename, sep='\t', compression=compression)
        df.set_index(['#CHROM', 'POS', 'REF', 'ALT'], drop=False, inplace=True)
        return df


class MasterVar(VariantFile):

    def __init__(self, filename):
        self._df = self._load_master(filename)
        self._process()
        return

    def _load_master(self, filename):
        return None

    def _process(self):
        self._df = None




