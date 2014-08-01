'''
Created on Jul 22, 2014

@author: Kunal Bhutani
'''
import pandas as pd
import numpy as np


class VariantFile(object):

    index = ['CHROM', 'POS', 'REF', 'ALT']
    features_columns = {'allele': ['HQ', 'EHQ', 'CGA_CEHQ', 'AD'],
                'binary': ['CGA_XR', 'CGA_RPT', 'multiallele'],
                'categorical': ['FT', 'vartype1', 'vartype2',
                                'phase', 'zygosity'],
                'genotype': ['GL', 'CGA_CEGL'],
                'numeric': ['CGA_SDO', 'GQ', 'DP', 'CGA_RDP']
                }

    def __init__(self, variant_df=None, features_df=None):
        self._variant_df = variant_df
        self._features_df = features_df

    @property
    def variants(self):
        return self._variant_df

    @property
    def features(self):
        return self._features_df.astype(float)

    @property
    def _autosome_variants(self):
        chroms = map(str, range(1, 23))
        in_autosome = self._variant_df['CHROM'].astype(str).isin(chroms)
        return in_autosome

    @property
    def _not_half_calls(self):
        half_calls = ['0/.', './0', '0|.', '.|0', '1/.', './1', '1|.', '.|1']
        not_half_calls = ~self._variant_df['GT'].isin(half_calls)
        return not_half_calls

    def filter(self, autosomes=True, not_half_calls=True, inplace=True):
        '''

        Filters variant df to only include autosomes and genotypes
        that do not contain half calls. It is possible to do this
        in place or create a new VariantFile object

        :param autosomes:
        :type autosomes:
        :param not_half_calls:
        :type not_half_calls:
        :param inplace:
        :type inplace:
        '''
        assert autosomes or not_half_calls, "Nothing to filter on."
        indices = np.ones(len(self._variant_df), dtype=bool)
        if autosomes:
            indices = indices & self._autosome_variants
        if not_half_calls:
            indices = indices & self._not_half_calls

        variant_df = self._variant_df[indices]
        features_df = self._features_df[indices]

        if inplace:
            self._variant_df = variant_df
            self._features_df = features_df

        else:
            return VariantFile(variant_df, features_df)

    @staticmethod
    def _get_phase(genotype):
        '''
        Returns phase from genotype
        '''

        if "|" in genotype:
            return "|"
        if "/" in genotype:
            return "/"
        else:
            return '-'

    @staticmethod
    def _get_allele(sample, phase, ref, alt,
                    chromosome_num=0):
        '''
        Using the field info field in a sample column and the phasing
        field, parses out the allele from the first or second chromosome.

        '''

        # Assumes first info field is GT
        allele_calls = sample.split(":")[0].split(phase)
        alleles = [ref] + alt.split(",")
        try:
            allele_number = int(allele_calls[chromosome_num])
            allele = alleles[allele_number]
        except:
            allele = "."
        return allele

    @staticmethod
    def _get_zygosity(a1, a2, ref):
        '''
        Based on the alleles present in sample, determines the zygosity
        of the individual in comparison to the reference.

        :param series:
        :type series:
        '''

        allele_set = set([a1, a2, ref])
        num_of_alleles = len(allele_set)

        if num_of_alleles == 1:
            return 'hom-ref'

        elif (num_of_alleles == 2  and
              a1 == a2):
            return 'hom-alt'

        elif (num_of_alleles == 2 and
              a1 != a2):
            return 'het-ref'

        elif num_of_alleles == 3:
            return 'het-alt'

        else:
            assert False

    @staticmethod
    def _get_ref_read_depth(series):
        '''
        This function returns the number of reads
        supporting the reference allele

        Required columns: [DP, AD, zygosity, a1]
        '''

        dp = series['DP']  # total read depth
        ad = series['AD'].replace(".", '0')  # replace missing with 0
        ad1, ad2 = ad.split(",")  # allelic depth for allele1 and allele2

        if series['zygosity'] == 'het-ref':
            if series['a1'] != series['REF']:  # allele1 is non-reference
                return  dp - int(ad1)  # total read depth - allele1 read depth
            return dp - int(ad2)  # total read depth - allele2 read depth

        # allele1 and allele2 read depths are identical
        if series['zygosity'] == 'hom-alt':
            return  dp - int(ad2)

        if series['zygosity'] == 'het-alt':
            return  dp - (int(ad1) + int(ad2))

        return dp  # hom-ref read depth == total read depth

    @staticmethod
    def _get_vartype(ref, alt):
        '''
        This function assigns the following vartypes to the
        allele specified by allele_base_col: snp, mnp, ins, del, indel or SV

        '''

        def is_snp(ref, alt):
            """ Return whether or not the variant is a SNP """
            if len(ref) != len(alt):
                return False
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

        if is_snp(ref, alt):
            return 'snp'
        elif is_sv(ref, alt):
            return 'sv'
        elif is_mnp(ref, alt):
            return 'mnp'
        elif is_insertion(ref, alt):
            return 'ins'
        elif is_del(ref, alt):
            return 'del'
        else:
            return 'indel or SV'

    def _generate_features(self):

        def multiallele_GL_formatting(df, gl_col):
            data = [i.split(",") for i in df[gl_col]]
            columns = ['AA_' + gl_col,
                        'AB_' + gl_col,
                        'BB_' + gl_col,
                        'AC_' + gl_col,
                        'BC_' + gl_col,
                        'CC_' + gl_col]
            GL_df = pd.DataFrame.from_records(data,
                                              columns=columns)
            GL_df.index = df.index
            GL_df.fillna(value=-999, inplace=True)
            return GL_df

        def allele_col_formatting(df, al_col):
            col1 = al_col + "_1"
            col2 = al_col + "_2"
            # TODO beautify this using partial functions
            df[col1] = [int(i.split(",")[0]) if i.split(",")[0] != "." else -999
                        for i in df[al_col]]
            df[col2] = [int(i.split(",")[1]) if i.split(",")[1] != "." else -999
                        for i in df[al_col]]
            return df[[col1, col2]]

        def format_category(df, category):
            return pd.get_dummies(df[category], prefix=category)

        df = self._variant_df
        features = VariantFile.features_columns

        # Formatting Allelic Features
        df_allele = df[features['allele']]
        df_allele = pd.concat([allele_col_formatting(df_allele, col)
                               for col in df_allele.columns],
                              axis=1)

        # Formatting binary features
        binarize = lambda x: 1 if x != 0 else 0
        df_binary = df[features['binary']].fillna(value=0).applymap(binarize)

        # Formatting categorical features
        df_categ = pd.concat([format_category(df, category)
                              for category in features['categorical']],
                              axis=1)

        # Formatting genotype features
        df_geno = df[features['genotype']]
        df_geno = pd.concat([multiallele_GL_formatting(df_geno, col)
                             for col in df_geno.columns],
                            axis=1)

        # Formatting numeric features
        df_numeric = df[features['numeric']].astype(float)

        # Create feature dataframe
        self._features_df = df_numeric.join([df_categ, df_geno,
                                            df_allele, df_binary])

        # Get rid of NaNs
        self._features_df.CGA_SDO.fillna(value=0, inplace=True)
        self._features_df.GQ.fillna(value=-999, inplace=True)


class VcfTsv(VariantFile):
    def __init__(self, filename, autosome_only=False, no_halfs=False):

        # Load VCF
        self._variant_df = self._load_vcf_tsv(filename)

        # Filter
        if autosome_only:
            self._autosome_variants()

        if no_halfs:
            self._no_half_variants()

        # Process variants
        self._process()

        # Generate features
        self._generate_features()

    def _process(self):
        df = self._variant_df

        # Ignore CNV and Repeats
        discard_vars = ['<INS:ME:ALU>', '<INS:ME:L1>',
                        '<CGA_CNVWIN>', '<INS:ME:SVA>']
        df = df[~df.ALT.isin(discard_vars)]

        # Parse many fields
        parse_columns = ['phase', 'multiallele', 'zygosity',
                         'a1', 'a2', 'vartype1', 'vartype2',
                          'ref_read_depth']
        for column in parse_columns:
            df[column] = ''

        df['phase'] = df['GT'].map(self._get_phase)
        df = df[df['phase'] != '-']

        df = df.apply(self._parse, axis=1)

        # Ignore sex chromosomes with - as phaser
        df = df[df.phase != '-']

        # Ignore structural variants
        df = df[(df.vartype1 != 'sv') & (df.vartype2 != 'sv')]

        # Fill missing allelic depth with 0,0 reads
        df['AD'].fillna(value='0,0', inplace=True)

        self._variant_df = df

    def _parse(self, series):
        '''
        Main processing function. Initial processing steps that might include
        feature generation or alternate representation of the dataframe.

        '''

        ref, alt = series['REF'], series['ALT']
        sample = series['GT']

        # Phase genotype and remove edge cases, likely at sex chromosome sites
        # series['phase'] = self._get_phase(series['GT'])
        # if series['phase'] == '-':
        #    return series

        # Get allele base sequences
        series['a1'] = self._get_allele(sample, series['phase'], ref, alt,
                                    chromosome_num=0)
        series['a2'] = self._get_allele(sample, series['phase'], ref, alt,
                                    chromosome_num=1)

        # Multiallelic site if more than one alternative allele
        series['multiallele'] = 1 if ',' in alt else 0

        # Zygosity: hom-ref, het-ref, hom-alt, het-alt
        series['zygosity'] = self._get_zygosity(series['a1'],
                                                series['a2'],
                                                ref)

        # Set the variant type based on number of bases in reference vs allele
        series['vartype1'] = self._get_vartype(ref, series['a1'])
        series['vartype2'] = self._get_vartype(ref, series['a2'])

        # Get number of reads supporting the reference allele
        series['ref_read_depth'] = self._get_ref_read_depth(series)

        return series

    def _load_vcf_tsv(self, filename):
        '''
        Loads in a vcf file, aware of gzipped files.
        '''
        compression = None
        if filename.endswith('.gz'):
            compression = 'gzip'

        df = pd.read_table(filename, sep='\t', compression=compression)
        df.set_index(['CHROM', 'POS', 'REF', 'ALT'], drop=False, inplace=True)
        return df


class MasterVar(VariantFile):

    def __init__(self, filename):
        self._variant_df = self._load_master(filename)
        self._process()
        return

    def _load_master(self, filename):
        return None

    def _process(self):
        self._variant_df = None
