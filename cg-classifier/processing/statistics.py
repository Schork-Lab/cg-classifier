from collections import defaultdict
import cPickle

import pandas as pd

import pybedtools

from variants import VariantFile


class Statistics():

    def __init__(self, variants, exon_file=None, from_df=False):

        if from_df:
            self.df = variants

        else:
            self.variants = variants
            self.df = self.variants.variants.copy()

        self.df = self.df[~self.df.GT.map(lambda x: '.' in x)]
        '''
        if exon_file:
            self.exon_df = self._intersect_exon_bed(exon_file)
        else:
            self.exon_df = None
        '''
        self._generate_statistics(exon_file)
        return

    def from_df(self, df, exon_file=None):
        self.__init__(df, exon_file, True)

    def _generate_statistics(self, exon_file=None):

        df = self.df
        self.counts = self._overall_counts(df)
        '''
        self.base_changes = self._base_changes(df)
        self.allele_depths = self._read_depths_2(df)
        self.ts_tv = self._ts_tv(self.base_changes)
        self.indels = self._indel_lengths(df)

        if exon_file:
            exon_df = self.exon_df
            self.exon_counts = self._overall_counts(exon_df)
            self.exon_base_changes = self._base_changes(exon_df)
            self.exon_allele_depths = self._read_depths_2(exon_df)
            self.exon_ts_tv = self._ts_tv(self.exon_base_changes)
            self.exon_indels = self._indel_lengths(exon_df)
        '''
        return

    @staticmethod
    def _overall_counts(df):
        df_not_hom = df[df.zygosity != 'hom-alt']
        df_hom = df[df.zygosity == 'hom-alt']
        count_df = pd.concat([df_hom.vartype1.value_counts(),
                               df_not_hom.vartype2.value_counts(),
                               df_not_hom.vartype1.value_counts()],
                             axis=1)
        return count_df.sum(axis=1)

    def _indel_lengths(self, df):

        def get_lengths(var_df, allele):
            ref_length = var_df['REF'].map(len)
            allele_length = var_df[allele].map(len)
            return allele_length - ref_length

        lengths = defaultdict(lambda: defaultdict(list))

        hom_alt = df.zygosity == 'hom-alt'
        het_ref = df.zygosity == 'het-ref'
        het_alt = df.zygosity == 'het-alt'

        for vartype in ['ins', 'del']:

            var1 = df.vartype1 == vartype
            var2 = df.vartype2 == vartype

            # If 1/1, add just one read depth
            var_df = df[var1 & var2 & hom_alt]
            lengths[vartype]['hom-alt'].extend(get_lengths(var_df, 'a1'))

            # If 1/0, add first read depth
            var_df = df[var1 & het_ref]
            lengths[vartype]['het-ref'].extend(get_lengths(var_df, 'a1'))

            # If 1/2, add first read depth
            var_df = df[var1 & het_alt]
            lengths[vartype]['het-alt'].extend(get_lengths(var_df, 'a1'))

            # If 0/1, add second read depth
            var_df = df[var2 & het_ref]
            lengths[vartype]['het-ref'].extend(get_lengths(var_df, 'a2'))

            # If 2/1, add second read depth
            var_df = df[var2 & het_alt]
            lengths[vartype]['het-alt'].extend(get_lengths(var_df, 'a2'))

        return pd.DataFrame(lengths)

    def _read_depths_2(self, df):
        read_depths = defaultdict(lambda: defaultdict(list))
        df['AD_1'] = df['AD'].map(lambda x: int(x.split(',')[0]))
        df['AD_2'] = df['AD'].map(lambda x: int(x.split(',')[1]))

        hom_alt = df.zygosity == 'hom-alt'
        het_ref = df.zygosity == 'het-ref'
        het_alt = df.zygosity == 'het-alt'

        vartypes = df.vartype1.unique()

        for vartype in vartypes:

            var1 = df.vartype1 == vartype
            var2 = df.vartype2 == vartype

            # If 1/1, add just one read depth
            var_df = df[var1 & var2 & hom_alt]
            read_depths[vartype]['hom-alt'].extend(var_df['AD_1'].values)

            # If 1/0, add first read depth
            var_df = df[var1 & het_ref]
            read_depths[vartype]['het-ref'].extend(var_df['AD_1'].values)

            # If 1/2, add first read depth
            var_df = df[var1 & het_alt]
            read_depths[vartype]['het-alt'].extend(var_df['AD_1'].values)

            # If 0/1, add second read depth
            var_df = df[var2 & het_ref]
            read_depths[vartype]['het-ref'].extend(var_df['AD_2'].values)

            # If 2/1, add second read depth
            var_df = df[var2 & het_alt]
            read_depths[vartype]['het-alt'].extend(var_df['AD_2'].values)

        return pd.DataFrame(read_depths)

    def _read_depths(self, df):
        read_depths = defaultdict(lambda: defaultdict(list))
        df['AD_1'] = df['AD'].map(lambda x: x.split(',')[0])
        df['AD_2'] = df['AD'].map(lambda x: x.split(',')[1])
        for row in df.iterrows():
            # If 1/1, add just one read depth
            if row[1]['zygosity'] == 'hom-alt':
                read_depths[row[1]['vartype1']]['hom-alt'].append(row[1]['AD_1'])

            # If 1/2, add read depth for each type of variant
            if row[1]['zygosity'] == 'het-alt':
                read_depths[row[1]['vartype1']]['het-alt'].append(row[1]['AD_1'])
                read_depths[row[1]['vartype2']]['het-alt'].append(row[1]['AD_2'])

            # If 0/1, add read depth for the variant
            if row[1]['zygosity'] == 'het-ref':
                if row[1]['vartype1'] != 'ref':
                    read_depths[row[1]['vartype1']]['het-ref'].append(row[1]['AD_1'])
                else:
                    read_depths[row[1]['vartype2']]['het-ref'].append(row[1]['AD_2'])
        return pd.DataFrame(read_depths)

    def _intersect_exon_bed(self, exon_file):
        def create_bed(df):
            df['Start'] = df['POS'] - 1
            bed_index = ['CHROM', 'Start', 'POS', 'REF', 'ALT']
            bed_repr = df.to_string(header=False, index_names=False,
                                    index=None,
                                    sparsify=False,
                                    columns=bed_index)
            bedtool = pybedtools.BedTool(bed_repr, from_string=True)
            return bedtool

        exon_bed = pybedtools.BedTool(exon_file)
        variant_bed = create_bed(self.df)
        overlapping_variants = variant_bed.intersect(exon_bed)
        if self.df.CHROM.dtype == 'int64':
            indices = [(int(x[0]), int(x[2]), x[3], x[4])
                       for x in overlapping_variants]
        else:
            indices = [(x[0], int(x[2]), x[3], x[4])
                       for x in overlapping_variants]

        return self.df.ix[indices]

    def _base_changes(self, df):
        '''

         Ugly ugly ugly. Complicated by current representation;
         would be easier with allelic primitives or no indel/snp variants
         at same position.

        :param df:
        :type df:
        '''

        def get_snp(series, alt='a1'):
            """ Return snp """
            ref = series['REF']
            alt = series[alt]
            if len(ref) != len(alt):
                return False
            for ref_nt, alt_nt in zip(ref, alt):
                if ref_nt != alt_nt:
                    return ref_nt + '>' + alt_nt

        pattern = dict()
        pattern['A'] = 'T'
        pattern['C'] = 'G'
        pattern['T'] = 'A'
        pattern['G'] = 'C'
        pattern['.'] = '.'
        change = set(['A', 'G'])

        def change_orientation(idx):
            ref = idx[0]
            alt = idx[2]
            if ref in change:
                ref = pattern[ref]
                alt = pattern[alt]
            return ref + '>' + alt


        base_changes = []

        hom_snps = df[(df.vartype1 == 'snp') &
                      (df.zygosity == 'hom-alt')]
        base_change_hom = hom_snps['REF'] + '>' + hom_snps['a1']
        base_changes.append(base_change_hom.value_counts())

        var1_snps = df[(df.vartype1 == 'snp') &
                       (df.zygosity != 'hom-alt')]
        if len(var1_snps) > 0:
            base_change_var1 = var1_snps.apply(get_snp, alt='a1', axis=1)
            base_changes.append(base_change_var1.value_counts())

        var2_snps = df[(df.vartype2 == 'snp') &
                       (df.zygosity != 'hom-alt')]
        if len(var2_snps) > 0:
            base_change_var2 = var2_snps.apply(get_snp, alt='a2', axis=1)
            base_changes.append(base_change_var2.value_counts())

        base_change = pd.concat(base_changes, axis=1)

        base_change.index = base_change.index.map(change_orientation)
        base_change['base_change'] = base_change.index
        base_change = base_change.groupby('base_change').agg(sum).sum(axis=1)
        return base_change

    def _ts_tv(self, base_change):
        transition = base_change.ix[['C>T', 'T>C']].sum(axis=1)
        transversion = base_change.ix[['C>A', 'C>G', 'T>A', 'T>G']].sum(axis=1)
        return transition, transversion

    def save(self, filename):
        stats = {}
        stats['counts'] = self.counts
        '''
        stats['base_changes'] = self.base_changes
        stats['ts_tv'] = self.ts_tv
        stats['allele_depths'] = self.allele_depths
        stats['indels'] = self.indels

        if self.exon_df is not None:
            stats['exon_counts'] = self.exon_counts
            stats['exon_base_changes'] = self.exon_base_changes
            stats['exon_allele_depths'] = self.exon_allele_depths
            stats['exon_ts_tv'] = self.exon_ts_tv
            stats['exon_indels'] = self.exon_indels
        '''
        cPickle.dump(stats, open(filename, 'wb'))
