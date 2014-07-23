'''
Created on Jul 22, 2014

@author: Kunal Bhutani
'''
import pandas as pd

class VariantFile(object):

    index = ['CHROM', 'POS', 'REF', 'ALT']
    features_columns = {'allele': ['HQ', 'EHQ', 'CGA_CEHQ', 'AD'],
                'binary': ['CGA_XR', 'CGA_RPT', 'multiallele'],
                'categorical': ['FT', 'vartype1', 'vartype2', 'phase', 'zygosity'],
                'genotype': ['GL', 'CGA_CEGL'],
                'numeric': ['CGA_SDO', 'GQ', 'DP', 'CGA_RDP']
                }


    def __init__(self):
        self._variant_df = None
        self._features_df = None

    @property
    def variants(self):
        return self._variant_df

    @property
    def features(self):
        return self._features_df.astype(float)

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
            return  dp - sum([ int(ad1), int(ad2) ])

        return dp  # hom-ref read depth == total read depth

    @staticmethod
    def _get_vartype(series, allele_base_col):
        '''
        This function assigns the following vartypes to the 
        allele specified by allele_base_col: snp, mnp, ins, del, indel or SV
        
        '''

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

        ref = series['REF']
        alt = series[allele_base_col]
        if is_snp(ref, alt): return 'snp'
        elif is_sv(ref, alt): return 'sv'
        elif is_mnp(ref, alt): return 'mnp'
        elif is_insertion(ref, alt): return 'ins'
        elif is_del(ref, alt): return 'del'
        else: return 'indel or SV'




    def _generate_features(self):

        def multiallele_GL_formatting(df, gl_col):
            GL_df = pd.DataFrame.from_records([i.split(",") for i in df[gl_col]],
                                               columns=['AA_' + gl_col,
                                                        'AB_' + gl_col,
                                                        'BB_' + gl_col,
                                                        'AC_' + gl_col,
                                                        'BC_' + gl_col,
                                                        'CC_' + gl_col])
            GL_df.index = df.index
            GL_df.fillna(value=-999, inplace=True)
            df = df.join(GL_df)
            return df


        def allele_col_formatting(df, al_col):
            df[al_col + "_1"] = [int(i.split(",")[0]) if i.split(",")[0] != "." else -999
                                 for i in df[al_col]]
            df[al_col + "_2"] = [int(i.split(",")[1]) if i.split(",")[1] != "." else -999
                                 for i in df[al_col]]
            return df


        def categ_formatting(df, cat_col):
            temp_dummies = pd.get_dummies(df[cat_col], prefix=cat_col)
            df = df.join(temp_dummies, how='left')
            return df

        df = self._variant_df
        features = VariantFile.features_columns

        # Formatting Allelic Features
        df_allele = df[features['allele']]
        for a_cols in df_allele.columns:
            df_allele = allele_col_formatting(df_allele, a_cols)
            del df_allele[a_cols]


        # Formatting binary features
        df_binary = df[features['binary']]
        for b in df_binary.columns:
            df_binary[b].fillna(value=0, inplace=True)
            df_binary[b] = df_binary[b].map(lambda x: 1 if x != 0 else 0)


        # Formatting categorical features
        df_categ = df[features['categorical']]
        for c in df_categ.columns:
            df_categ = categ_formatting(df_categ, c)
            del df_categ[c]


        # Formatting genotype features
        df_geno = df[features['genotype']]
        for g_cols in df_geno.columns:
            df_geno = multiallele_GL_formatting(df_geno, g_cols)
            del df_geno[g_cols]


        # Formatting numeric features
        df_numeric = df[features['numeric']]
        for n in df_numeric.columns:
            df_numeric[n] = df_numeric[n].astype(float)


        # Create feature dataframe
        self._features_df = df_numeric.join([df_categ, df_geno, df_allele, df_binary])


class VcfTsv(VariantFile):
    def __init__(self, filename):

        # Load VCF
        self._variant_df = self._load_vcf_tsv(filename)

        # Process variants
        self._process()

        # Generate features
        self._generate_features()


    def _process(self):
        '''
        Main processing function. Initial processing steps that might include feature generation
        or alternate representation of the dataframe.
        
        '''

        df = self._variant_df

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
        df['ref_read_depth'] = df.apply(self._get_ref_read_depth, axis=1)

        self._variant_df = df


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




