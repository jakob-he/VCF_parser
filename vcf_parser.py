#!/usr/bin/env python3
'''
Class used to parse and process a VCF file.
'''
import gzip
import numpy as np
import pandas


class VCF_parser:

    def __init__(self, vcf_file: str):
        self._df = self._read_vcf(vcf_file)

    def __repr__(self):
        return self._df

    def _read_vcf(self, file: str):
        """
        Reads a gziped VCF file and returns a dataframe containing
        the variant information
        """
        if file.endswith("gz"):
            with gzip.open(file, "rt") as vcf:
                return self._process_vcf(vcf)
        else:
            with open(file, "r") as vcf:
                return self._process_vcf(vcf)

    def _process_vcf(self, file_object) -> pandas.DataFrame:
        line = file_object.readline().rstrip('\n')
        while line.startswith('##'):
            line = file_object.readline()
        column_names = line.split('\t')
        df = pandas.read_csv(file_object, sep='\t', names=column_names)
        # In case of VCFs that only containt numerical chromsome values
        # all values in the CHROM column have to converted to strings
        df['#CHROM'] = df['#CHROM'].apply(str)
        return df

    def _get_mean_phred(self):
        return np.mean(self._df['QUAL'])

    def _get_statistics(self):
        """
        Parses the VCF dataframe for necessary statistics:

        - Number of transitions
        - Number of transversion
        - Number of heterozygous calls (for substiutions of one base)
        - Number of homozygous calls (for substitiotions of one base)
        - mean,sd and median of the Genotype quality
        """
        ti = 0
        tv = 0
        het = 0
        hom = 0
        gq = []

        for index, variant in self._df.iterrows():
            purine = ['A','G']
            pyrimidine = ['C','T']

            # check for transitions and transversions
            if variant['REF'] in purine and variant['ALT'] in purine:
                ti += 1
            elif variant['REF'] in purine and variant['ALT'] in purine:
                ti += 1
            elif variant['ALT'] in purine and variant['REF'] in pyrimidine:
                tv += 1
            elif variant['ALT'] in purine and variant['REF'] in pyrimidine:
                tv += 1

            # check for heterozygous/homozygous calls
            # TODO generalize for multisample vcfs
            genotype = variant[-1].split(":")

            if genotype[0] in ['0/1', '1/2', '0|1', '1|2']:
                het += 1
            elif genotype[0] in ['0/0', '1/1', '0|0', '1|1']:
                hom += 1
            gq.append(int(genotype[-2]))

        parameters = ['Variant call Phred score', 'Total variants called', 'Transitions',
                      'Transversions', 'Ti/Tv', 'Het calls', 'Hom calls', 'Het/Hom', 'mean GQ', 'std GQ', 'median GQ']
        values = [self._get_mean_phred(), self._df.shape[0], ti, tv, ti /
                  tv, het, hom, het / hom, np.mean(gq), np.std(gq), np.median(gq)]
        range = ['-', '-', '-', '-', '[2-3]',
                 '-', '-', '[1.3-2.0]', '-', '-', '-']
        statistics = pandas.DataFrame(
            {'Parameters': parameters, 'Values': values, 'Range': range})
        return statistics

    def show_statistics(self):
        self._statistics = self._get_statistics()
        print(self._statistics.to_string())


# test
def main():
    testvcf = "daughter.vcf"
    parsed = VCF_parser(testvcf)
    parsed.show_statistics()


if __name__ == "__main__":
    main()
