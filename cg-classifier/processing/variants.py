'''
Created on Jul 22, 2014

@author: Kunal Bhutani
'''
import pandas as pd
import pybedtools
import numpy as np
import os



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
        chroms = range(1, 23)
        chroms = map(str, chroms)
        chroms = set(chroms)
        self.variants['CHROM'] = self.variants['CHROM'].astype(str)
        in_autosome = self.variants['CHROM'].map(lambda x: x in chroms)
        return in_autosome

    @property
    def _not_half_calls(self):
        not_half_calls = ~self.variants[self.sample_col].map(lambda x: '.' in x)
        return not_half_calls
    
    @property
    def _not_haploid(self):
        sv = ~self.variants[self.sample_col].map(lambda x: len(x.split(":")[0]) == 1 )
        return sv
    


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

    def create_bed(self):
        df = self._variant_df
        df['Start'] = df['POS'] - 1
        bed_index = ['CHROM', 'Start', 'POS', 'REF', 'ALT']
        bed_repr = df.to_string(header=False, index_names=False,
                                index=None,
                                sparsify=False,
                                columns=bed_index)
        bedtool = pybedtools.BedTool(bed_repr, from_string=True)
        return bedtool

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
        if ref == alt:
            return 'ref'
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

            GL_df = pd.DataFrame.from_records(data)
            num_of_genotypes = len(GL_df.columns)
            GL_df.columns = columns[:num_of_genotypes]
            for unused_col in columns[num_of_genotypes:]:
                GL_df.unused_col = -999
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
    def __init__(self, filename, generate_features=False,
                 autosome_only=True, no_halfs=True):

        # Load VCF
        self._variant_df = self._load_vcf_tsv(filename)

        # Filter
        if autosome_only:
            self._variant_df = self._variant_df[self._autosome_variants]

        if no_halfs:
            self._variant_df = self._variant_df[self._not_half_calls]

        # Process variants
        self._process()

        # Generate features
        if generate_features:
            self._generate_features()
        else:
            self._features_df = pd.DataFrame([])

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
        # series['ref_read_depth'] = self._get_ref_read_depth(series)

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



class Vcf(VariantFile):
    def __init__(self, filename, sample_col, generate_features=False,
                 autosome_only=True, no_halfs=True):

        # Load VCF
        self._variant_df = self._load_vcf(filename)
        
        # Sample Column Name
        self.sample_col = sample_col

        
        # Filter
        if autosome_only:
            self._variant_df = self._variant_df[self._autosome_variants]
            
        if no_halfs:
            self._variant_df = self._variant_df[self._not_half_calls]
            
        #remove haploid called variants, occurs with sex chromosomes, mitochondria, or some structural variants    
        self._variant_df = self._variant_df[self._not_haploid]
 
        # Process variants
        self._process()
 
        # Generate features
        if generate_features:
            self._generate_features()
        else:
            self._features_df = pd.DataFrame([])
            


    def _process(self):
        df = self._variant_df

        # Ignore CNV and Repeats
        discard_vars = ['<INS:ME:ALU>', '<INS:ME:L1>',
                        '<CGA_CNVWIN>', '<INS:ME:SVA>']
        df = df[~df.ALT.isin(discard_vars)]

        # Parse many fields
        df = self._parse(df)
        
        df = df[df['phase'] != '-']

        # Ignore sex chromosomes with - as phaser
        df = df[df.phase != '-']

        # Ignore structural variants
        df = df[(df.vartype1 != 'sv') & (df.vartype2 != 'sv')]

        # Fill missing allelic depth with 0,0 reads
        #df['AD'].fillna(value='0,0', inplace=True)

        self._variant_df = df


    def _get_allele(self, line, gt_col):
        '''
        Returns allelic base, handles multi-allelic variants
        '''
        alleles = [line['REF']]
        alleles.extend(line['ALT'].split(","))
        a1 = "."
        try:
            a1 = alleles[int(line[gt_col])]  #returns missing if gt_int_call is "."
        except:
            a1 = "."
        return a1
        
        
    def _get_GT_multisample_vcf(self, line, sample_col, gt_index):
        '''
        Slow parser for multisample vcf
        '''
        return int( line[sample_col].split(line['phase'])[int(gt_index)])



    def _get_allele_bases(self, df, sample_col, single_sample_vcf=True):
        '''
        Adds phase, GT1, GT2, a1, a2 to self._variants_df dataframe
        
        phase : { /, | }, unphased or phased call
        GT1: int, first allele call in the numeric sample genotype column
        GT2: int, second allele call in the numeric sample genotype column
        a1: {A, T, G, C, AA, etc}, nucleotide base representation for GT1
        a2: {A, T, G, C, AA, etc}, nucleotide base representation for GT2
        
        '''
        
        
        if single_sample_vcf:
            df['phase'] = df[sample_col].str[1]
            df['GT1'] = df[sample_col].str[0]
            df['GT1'] = df['GT1'].astype(int)
            df['GT2'] = df[sample_col].str[2]
            df['GT2'] = df['GT2'].astype(int)
        
        
        if not single_sample_vcf:
            df['phase'] = df.apply(get_phase, args=['GT'], axis=1)  #get phase
            df = df[df.phase != "-"]  #likley occurs at sex chromosome sites
            df['GT1'] = df.apply(self._get_GT_multisample, args=[sample_col, 0], axis=1)
            df['GT2'] = df.apply(self._get_GT_multisample, args=[sample_col, 1], axis=1)
            
            
        
        #SLOW PROCESS MULTIPLE ALLELE GENOTYPES
        df_multi = df[(df.GT1.astype(int)>1) | (df.GT2.astype(int)>1)] #select all multi-alleleic variants
        df_multi['a1'] = df_multi.apply(self._get_allele, args=['GT1'], axis=1)  #
        df_multi['a2'] = df_multi.apply(self._get_allele, args=['GT2'], axis=1)
        
        
        #FAST PROCESS SIMPLE ALLELE GENOTYPES
        df_simple = df[~df.index.isin(df_multi.index)][['REF', 'ALT', 'GT1', 'GT2']]  #dropping multiallele variants, minimize memory usage
        
        df_gt1_ref = df_simple[df_simple.GT1==0][['REF']]  #get a1 ref alleles
        df_gt1_ref.columns = ['a1']
        df_gt2_ref = df_simple[df_simple.GT2==0][['REF']]  #get a2 ref alleles
        df_gt2_ref.columns = ['a2']
        
        
        df_gt1_alt = df_simple[df_simple.GT1==1][['ALT']]  #get a1 alt alleles
        df_gt1_alt.columns = ['a1']
        df_gt2_alt = df_simple[df_simple.GT2==1][['ALT']]  #get a2 alt alleles
        df_gt2_alt.columns = ['a2']
        
        
        gt1_alleles = pd.concat([df_gt1_ref,df_gt1_alt])  #merging GT1 allele bases into a single df
 
        gt2_alleles = pd.concat([df_gt2_ref,df_gt2_alt])  #merging GT2 allele bases into a single df

        gt1_2_allele_df = gt1_alleles.join(gt2_alleles, how='outer')  #Joining the GT1 and GT2 simple allele bases 
        
        
        df = df.join(gt1_2_allele_df, how='inner')  #Adding simle allele a1 and a2 columns to original df
        df = df.append(df_multi)  #Adding multi-alleleic bases to original df
        
        return df



    def _zygosity_fast(self, df):
        '''
        This function quickly assigns zygosity states to 
        each variant using set logic.
        
        zygosity: {hom-ref, hom-miss, het-miss, het-alt, hom-alt, het-ref}
        examples:  REF     ALT
        hom-ref    A        A/A
        hom-miss   A        ./.
        het-miss   A        A/.
        het-alt    A        T/G
        hom-alt    A        T/T
        het-ref    A        A/T
        '''
        
        #A/A
        df_hom_ref = df[ (df['a1'] == df['REF']) & (df['a2'] == df['REF'])]
        df_hom_ref['zygosity'] = 'hom-ref'
        
        #./.
        df_hom_miss = df[ (df['a1'] == '.') & (df['a2'] == '.')]
        df_hom_miss['zygosity'] = 'hom-miss'
        
        #A/.
        df_het_miss = df[ (df['a1'] == '.') | (df['a2'] == '.') ]
        df_het_miss['zygosity'] = 'het-miss'
        
        # NO MISSING CALLS
        df_not_miss = df[~df.index.isin( set(df_hom_miss.index) | set(df_het_miss.index))  ]
        
        #T/G
        df_het_alt = df_not_miss[ ((df_not_miss['a1'] != df_not_miss['REF']) & (df_not_miss['a2'] != df_not_miss['REF'])) & (df_not_miss['a1'] != df_not_miss['a2']) ]
        df_het_alt['zygosity'] = 'het-alt'
        
        #T/T
        df_hom_alt = df_not_miss[ (((df_not_miss['a1'] != df_not_miss['REF']) & (df_not_miss['a2'] != df_not_miss['REF']))) & (df_not_miss['a1'] == df_not_miss['a2']) ]
        df_hom_alt['zygosity'] = 'hom-alt'
        
        #A/T
        df_het_ref = df_not_miss[ ((df_not_miss['a1'] == df_not_miss['REF']) & (df_not_miss['a2'] != df_not_miss['REF'])) | ((df_not_miss['a1'] != df_not_miss['REF']) & (df_not_miss['a2'] == df_not_miss['REF']))  ]
        df_het_ref['zygosity'] = 'het-ref'
        
        #Appending dataframes
        df_zygosity = pd.concat([df_hom_ref, df_hom_miss, df_het_miss, df_het_ref, df_het_alt, df_hom_alt])
        
        
        
        #print set(df_zygosity.index) - set(df.index)
        assert len(df_zygosity) == len(df) #ensures no variants were missed or double counted
        return df_zygosity



    def _vartype_map(self, ref_alt_bases):
        '''
        This function assigns the following vartypes to the 
        allele specified by allele_base_col: snp, mnp, ins, del, indel or SV
        
        '''
        ref, alt = ref_alt_bases
        len_diff = len(ref) - len(alt)
        
        if ref == alt: return 'ref' #Orderd by frequency of the variant to reduce complexity
        
        if len_diff == 0:
            base_diff = [nt for i,nt in enumerate(alt) if ref[i] != alt[i]]
            if len(base_diff) == 1: return 'snp' 
            else: return 'mnp'
        
        if len_diff > 0:
            base_diff = [nt for i,nt in enumerate(alt) if ref[i] != alt[i]]
            if len(base_diff) > 0: return 'indel'
            else: return 'del'
            
        if len_diff < 0:
            return 'ins'
        
        elif is_sv(ref,alt): return 'sv'
        
        else: return 'indel or SV'


    def _parse(self, df):
        '''
        Main processing function. Initial processing steps that might include
        feature generation or alternate representation of the dataframe.

        '''

        # Phase genotype and remove edge cases, likely at sex chromosome sites
        # series['phase'] = self._get_phase(series['GT'])
        # if series['phase'] == '-':
        #    return series

        # Get allele base sequences
        df = self._get_allele_bases(df, self.sample_col)  #adding 'phase', 'GT1', 'GT2', 'a1', 'a2' to df
        
        # Zygosity: hom-ref, het-ref, hom-alt, het-alt
        df = self._zygosity_fast(df)
        
        # Set the variant type based on number of bases in reference vs allele
        df['vartype1'] = map(self._vartype_map, df[['REF','a1']].values)
        df['vartype2'] = map(self._vartype_map, df[['REF','a2']].values)
        
        # Multiallelic site if more than one alternative allele
        df['multiallele'] = df.ALT.map(lambda x: 1 if "," in x else 0)
        
        return df


    def _get_header(self, filename):
        '''
        Creates Header object (list)
        '''
        if filename.endswith('.gz'):
            header_lines = os.popen('tabix -H ' + filename).readlines()
            return [l.replace('#CHROM','CHROM') for l in header_lines if l.startswith('#')]
        header_lines = os.popen('head -5000 ' + filename).readlines()
        return [l.replace('#CHROM','CHROM') for l in header_lines if l.startswith('#')]


    def _load_vcf(self, filename):
        '''
        Loads in a vcf file, aware of gzipped files.
        '''
        
        header = self._get_header(filename)
        #print header
        df = pd.read_table(filename, sep="\t", compression='gzip', skiprows=(len(header)-1))
        df.columns = header[-1].rstrip('\n').split("\t")
        df.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True, drop=False)
    
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
