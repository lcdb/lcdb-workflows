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
        self.levelMap = {
            'runLevel': 'rawLevel',
            'sampleLevel': 'runLevel',
            'aggLevel': 'sampleLevel'
        }

    def _compile(self, level):
        #return list(set([self.config[level].format_map(x) for x in self.samples]))
        return regex(self.config[level])

    def _compile_regex(self):
        """ Build all possible prefixes for each level """
        self.raw = self._compile('rawLevel')
        self.run = self._compile('runLevel')
        self.sample = self._compile('sampleLevel')
        self.agg = self._compile('aggLevel')

    def _load_sample_table(self):
        """ Import the sample table and make a dictionary version too. 

        This function imports the sample table and saves its information in two
        formats. The first is a pandas.DataFrame and the second is a list of
        dictionaries.

        """
        self.sampleTable = pd.read_table(self.config['sampletable'], sep='\t', dtype=str)
        self.sampleTable.set_index('sampleID', inplace=True)
        self.samples = self.sampleTable.reset_index().to_dict('records')

    def find_level(self, prefix):
        """ Figure out which regex the prefix matches. """
        if re.match(self.raw, prefix):
            return 'rawLevel', self.find_sample(self.raw, prefix)
        elif re.match(self.run, prefix):
            return 'runLevel', self.find_sample(self.run, prefix)
        elif re.match(self.sample, prefix):
            return 'sampleLevel', self.find_sample(self.sample, prefix)
        elif re.match(self.agg, prefix):
            return 'aggLevel'
        else:
            raise ValueError

    def find_sample(self, pattern, prefix):
        """ Find which sample(s) to use """
        m = re.match(pattern, prefix).groupdict()
        _samples = []
        for s in self.samples:
            if m.items() <= s.items():
                _samples.append(s)
        return _samples

    def make_input(self, suffix, lookup=False, agg=False):
        """ Generates Input Function based on wildcards """
        def _input(wildcards):
            if lookup:
                _suffix = wildcards[suffix]
            else:
                _suffix = suffix

            if agg:
                level, _samples = self.find_level(wildcards['prefix'])
                _sampleList = pd.DataFrame(_samples).to_dict('list')
                return expand(self.config[self.levelMap[level]] + '{ext}', ext=_suffix, **_sampleList)
            else:
                return expand('{prefix}{ext}', prefix=wildcards['prefix'], ext=_suffix)

        return _input

    def build_targets(self, patterns):
        """ Build target file names based on pattern namming scheme. """
        _targets = []
        for p in patterns:
            p = p.format_map(self.config)
            for s in self.samples:
                e = dict(s, **self.config)
                _targets.append(p.format_map(e))
        return list(set(_targets))
