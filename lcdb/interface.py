import re

import pandas as pd
from snakemake.io import expand, regex


class SampleHandler(object):
    """ """
    def __init__(self, config):
        self.config = config

        # Load sampleTable
        self._load_sample_table()

        # Build all regexs
        self._compile_regex()

        # Order Dict
        self.order = {
            'runLevel': 'rawLevel',
            'sampleLevel': 'runLevel',
            'aggLevel': 'sampleLevel'
        }

    def _compile(self, level):
        return list(set([self.config[level].format_map(x) for x in self.samples]))

    def _compile_regex(self):
        """ Build all possible prefixes for each level """
        self.raw = self._compile('rawLevel')
        self.run = self._compile('runLevel')
        self.sample = self._compile('sampleLevel')
        self.agg = self.config['aggLevel']

    def _load_sample_table(self):
        """ Import the sample table and make a dictionary version too. 

        This function imports the sample table and saves its information in two
        formats. The first is a pandas.DataFrame and the second is a list of
        dictionaries.

        """
        self.sampleTable = pd.read_table(self.config['sampletable'], sep='\t')
        self.sampleTable.set_index('sampleID', inplace=True)
        self.samples = self.sampleTable.reset_index().to_dict('records')

    def find_level(self, prefix):
        """ Figure out which regex the prefix matches. """
        if prefix in self.raw:
            return 'rawLevel'
        elif prefix in self.run:
            return 'runLevel'
        elif prefix in self.sample:
            return 'sampleLevel'
        elif prefix == self.config['aggLevel']:
            return 'aggLevel'
        else:
            raise ValueError

    def find_sample(self, prefix):
        """ Find which sample(s) to use """
        pass

    def make_input(self, suffix, agg=False):
        """ Generates Input Function based on wildcards """
        def _input(wildcards):
            if agg:
                level = self.find_level(wildcards['prefix'])
                _samples = self.find_sample(wildcards['prefix'])
                return expand('{prefix}{ext}', prefix=wildcards['prefix'], ext=suffix)
            else:
                return expand('{prefix}{ext}', prefix=wildcards['prefix'], ext=suffix)
        return _input
