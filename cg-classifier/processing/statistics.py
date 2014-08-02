from collections import defaultdict
import cPickle

import pandas as pd

import pybedtools

from variants import VariantFile


class Statistics():

    def __init__(self, variant_file, exon_file=None):
        assert isinstance(variant_file, VariantFile), "VariantFile needed"
        self.variant_file = variant_file
        self.df = self.variant_file.variants.copy()
        self.df = self.df[~self.df.GT.map(lambda x: '.' in x)]
        if exon_file:
            self.exon_df = self._intersect_exon_bed(exon_file)
        else:
            self.exon_df = None
        self._generate_statistics(exon_file)
        return

    def _generate_statistics(self, exon_file=None):

        df = self.df
        self.counts = self._overall_counts(df)
        self.base_changes = self._base_changes(df)
        self.allele_depths = self._read_depths(df)
        self.ts_tv = self._ts_tv(self.base_changes)

        if exon_file:
            exon_df = self.exon_df
            self.exon_counts = self._overall_counts(exon_df)
            self.exon_base_changes = self._base_changes(exon_df)
            self.exon_allele_depths = self._read_depths(exon_df)
            self.exon_ts_tv = self._ts_tv(self.exon_base_changes)

        return

    def _overall_counts(self, df):
        df_not_hom = df[df.zygosity != 'hom-alt']
        df_hom = df[df.zygosity == 'hom-alt']
        count_df = pd.concat([df_not_hom.vartype1.value_counts(),
                               df_hom.vartype2.value_counts(),
                               df_hom.vartype1.value_counts()],
                             axis=1)
        return count_df.sum(axis=1)

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
        if self.variant_file.variants.CHROM.dtype == 'int64':
            indices = [(int(x[0]), int(x[2]), x[3], x[4])
                       for x in overlapping_variants]
        else:
            indices = [(x[0], int(x[2]), x[3], x[4])
                       for x in overlapping_variants]

        return self.variant_file.variants.ix[indices]

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
        stats['base_changes'] = self.base_changes
        stats['ts_tv'] = self.ts_tv
        stats['allele_depths'] = self.allele_depths

        if self.exon_df is not None:
            stats['exon_counts'] = self.exon_counts
            stats['exon_base_changes'] = self.exon_base_changes
            stats['exon_allele_depths'] = self.exon_allele_depths
            stats['exon_ts_tv'] = self.exon_ts_tv

        cPickle.dump(stats, open(filename, 'wb'))
