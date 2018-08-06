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
        self._ti = 0
        self._tv = 0
        self._het = 0
        self._hom = 0
        self._gq = []
        self._mu_types = {}

        for index, variant in self._df.iterrows():
            self._check_ti_tv(variant)
            self._check_genotype(variant)
            self._check_mutation_types(variant)

        parameters = ['Variant call Phred score', 'Total variants called', 'Transitions',
                      'Transversions', 'Ti/Tv', 'Het calls', 'Hom calls', 'Het/Hom', 'mean GQ', 'std GQ', 'median GQ']
        values = [self._get_mean_phred(), self._df.shape[0], self._ti, self._tv, self._ti /
                  self._tv, self._het, self._hom, self._het / self._hom, np.mean(self._gq), np.std(self._gq), np.median(self._gq)]
        range = ['-', '-', '-', '-', '[2-3]',
                 '-', '-', '[1.3-2.0]', '-', '-', '-']
        statistics = pandas.DataFrame(
            {'Parameters': parameters, 'Values': values, 'Range': range})

        # add mutation type information
        #!! only available after annotation with e.g. jannovar
        nonsynomous = self._mu_types['missense_variant'] + self._mu_types['stop_gained']
        statistics = statistics.append(pandas.DataFrame({'Parameters': ['Synonymous variants', 'Nonsynomous variants', 'ns/ss'], 'Values': [self._mu_types['synonymous_variant'], nonsynomous, nonsynomous / self._mu_types['synonymous_variant']], 'Range': ['-', '-', '[0.8-1.0]']}))
        statistics.reset_index(inplace=True)
        return statistics

    def _check_ti_tv(self, variant: pandas.Series):
        purine=['A', 'G']
        pyrimidine=['C', 'T']

        # check for transitions and transversions
        if variant['REF'] in purine and variant['ALT'] in purine:
            self._ti += 1
        elif variant['REF'] in purine and variant['ALT'] in purine:
            self._ti += 1
        elif variant['ALT'] in purine and variant['REF'] in pyrimidine:
            self._tv += 1
        elif variant['ALT'] in purine and variant['REF'] in pyrimidine:
            self._tv += 1

    def _check_genotype(self, variant: pandas.Series):
        genotype=variant[-1].split(":")

        if genotype[0] in ['0/1', '1/2', '0|1', '1|2']:
            self._het += 1
        elif genotype[0] in ['0/0', '1/1', '0|0', '1|1']:
            self._hom += 1
        self._gq.append(int(genotype[-2]))

    def _check_mutation_types(self, variant: pandas.Series):
        # hardcoded! TODO: dynamic approach
        annotation=variant['INFO'].split(";")[3].split("|")[1]
        if annotation in self._mu_types:
            self._mu_types[annotation] += 1
        else:
            self._mu_types[annotation]=1

    def statistic_report(self):
        print("=========== Statistic Report ============\n")
        self._statistics=self._get_statistics()
        print(self._statistics.to_string())
        print(self._get_variant_dist().to_string())

    def _get_variant_dist(self):
        print("\n=========== Variant Distribution ============")
        dist = pandas.DataFrame.from_dict(self._mu_types,columns = ["Count"], orient = 'index')
        return

    def filter_vcf(self,threshold:int = 0.1):
        print("\n========== Filter variants with MAF <= {} ===========".format(threshold))
        current_length = self._df.shape[0]
        for index, variant in self._df.iterrows():
            maf = self.get_maf(variant)
            if float(maf) <= threshold:
                self._df.drop(index, inplace = True)
        filtered_length = self._df.shape[0]
        print("\n filtered {} variants".format(current_length-filtered_length))

    def get_maf(self,variant: pandas.Series):
        if variant['ID'] != '.':
            afs = variant['INFO'].split('|')[-1].split(';')[2]
            if 'CAF' in afs:
                return(min(afs.split(',')))
            else:
                return(1)
        else:
            return(1)

# test
def main():
    testvcf="../data/daughter-annotated-all-dbsnp.vcf"
    parsed=VCF_parser(testvcf)
    parsed.filter_vcf()


if __name__ == "__main__":
    main()
