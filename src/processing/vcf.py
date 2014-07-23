'''
Created on Jul 22, 2014

@author: Kunal Bhutani
'''

class VCF(object):
    '''
    pdVCF class. 
    
    Note: It's frowned upon / hacky to inherit dataframe class, so we will have to work with a 
    VCF class and its representation.
    '''
    features = ['DP']

    def __init__(self, filename, filetype='vcf'):
        self._df = load_vcf(filename)
        self._process()
        return

    @property
    def df(self):
        return self._df

    @property
    def features(self):
        return self._df[features].values

    def _process(self):
        '''
        Main processing function. Initial processing steps that might include feature generation
        or alternate representation of the dataframe.
        
        '''

        # Overall processing functions
        self._df['variant_type'] = self._df.apply(_variant_type, axis=1)

        # Add other processing functions

        # Generate features
        self._generate_features()

        return

    def _generate_features(self):
        # Generate features
        self._df['DP'] = 0


    def load_vcf(self, filename):
        '''
        Loads in a vcf file, aware of gzipped files.
        '''
        compression = None
        if filename.endswith('.gz'):
            IN = gzip.open(filename)
            compression = 'gzip'
        else:
            IN = open(filename)

        num_of_header = 0
        for line in IN:
            if line.startswith('#'):
                num_of_header += 1
            else:
                break
        IN.close()

        df = pd.read_table(filename, sep='\t', compression=compression,
                            skiprows=num_of_header - 1)
        df.set_index(['#CHROM', 'POS', 'REF', 'ALT'], drop=False, inplace=True)
        return df

    def _variant_type(series):
        def is_snp(ref, alt):
            """ Return whether or not the variant is a SNP """
            if len(ref) != len(alt): return False
            base_diff = [nt for i, nt in enumerate(alt) if ref[i] != alt[i]]
            if len(base_diff) <= 1:
                    return True
            return False


        def is_mnp(ref, alt):
            """ Return whether or not the variant is an MNP"""
            if len(ref) > 1 and len(ref) == len(alt):
                return True

        def is_insertion(ref, alt):
            if len(ref) < len(alt):
                return True
            return False

        def is_del(ref, alt):
            if len(ref) > len(alt):
                return True
            return False


        ref = series['REF']
        alt = series['ALT']
        if is_snp(ref, alt):
            return 'snp'
        elif is_mnp(ref, alt):
             return 'mnp'
        elif is_insertion(ref, alt):
            return 'ins'
        elif is_del(ref, alt):
            return 'del'
        else:
            return 'indel or SV'


