#!/usr/bin/env python3
'''
Class used to parse and process a VCF file.
'''
import gzip

import pandas


class VCF_parser:

    def __init__(self, vcf_file: str):
        self._df = self._read_vcf(vcf_file)

    def __repr__(self):
        return self._df

    def _read_vcf(self,file: str):
        """
        Reads a gziped VCF file and returns a dataframe containing
        the variant information
        """
        if file.endswith("gz"):
            with gzip.open(file, "rt") as vcf:
                return self._process_vcf(vcf)
        else:
            with open(file,"r") as vcf:
                return self._process_vcf(vcf)

    def _process_vcf(self,file_object) -> pandas.DataFrame:
        line = file_object.readline().rstrip('\n')
        while line.startswith('##'):
            line = file_object.readline()
        column_names = line.split('\t')
        df = pandas.read_csv(file_object, sep='\t', names=column_names)
        # In case of VCFs that only containt numerical chromsome values
        # all values in the CHROM column have to converted to strings
        df['#CHROM'] = df['#CHROM'].apply(str)
        return df


#test

def main():
    testvcf = "sample.vcf"
    parsed = VCF_parser(testvcf)


if __name__ == "__main__":
    main()
