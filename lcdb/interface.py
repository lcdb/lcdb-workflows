import re
import os

import pandas as pd
import yaml

from snakemake.io import expand, regex


class SampleHandler(object):
    """ Basic interface to help handle filenames in snakemake """
    def __init__(self, config):
        self.config = config

        # Load sampleTable
        self._load_sample_table()

        # Build all regexs
        self._compile_regex()

        # Order Dict to specify hierarchy of levels
        self.levelMap = {
            'runLevel': 'rawLevel',
            'sampleLevel': 'runLevel',
            'aggLevel': 'runLevel'
        }

    def _compile(self, level):
        """ Use snakemake regex to compile regex from format string.

        Snakemake provides a nice regex function that converts format strings
        to regex. To help I also add sampleTable values to the regex to limit
        the results.

        Parameters
        ----------
        level: str
            A string consiting of rawLevel, runLevel, sampleLevel, aggLevel

        Returns
        -------
        str
            regex pattern generated by snakemake.io.regex

        Example
        -------

        >>> SH = SampleHandler(test_config)
        >>> level = 'runLevel'
        >>> pattern = SH.config[level]

        >>> pattern
        'pasilla_sample/{sampleID}/{sampleID}_{treatment}_{replicate}_R1'

        >>> assert SH._compile(level) == (
        ... 'pasilla_sample\\\\/(?P<sampleID>treated1|treated2|untreated1'
        ... '|untreated2)\\\\/(?P=sampleID)_(?P<treatment>treated|untreated)'
        ... '_(?P<replicate>1|2)_R1'
        ... )
        """
        pattern = self.config[level]
        for name, values in self.sampleTable.reset_index().to_dict('list').items():
            # Subsitute the first instance of each sampleTable column name and
            # add the unique list of column values. This will help narrow down
            # regex. NOTE: This may not be needed, but thought it might be useful.
            pattern = re.sub(
                '{{{name}}}'.format(name=name),
                '{' + '{name}, {res}'.format(
                    name=name, res='|'.join(sorted(set(values)))
                ) + '}', pattern, count=1)
        # Retrun regex removing the '$' off of the end to allow partial matches
        return regex(pattern)[:-1]

    def _compile_regex(self):
        """ Compile regex for each level

        Attributes
        ----------
        raw: str
            compiled regex for raw
        run: str
            compiled regex for run
        sample: str
            compiled regex for sample
        agg: str
            compiled regex for agg

        """
        self.raw = self._compile('rawLevel')
        self.run = self._compile('runLevel')
        self.sample = self._compile('sampleLevel')
        self.agg = self._compile('aggLevel')

    def _load_sample_table(self):
        """ Import the sample table and make a dictionary version too.

        This function imports the sample table and saves its information in two
        formats. The first is a pandas.DataFrame and the second is a list of
        dictionaries.

        Attributes
        ----------
        sampleTable: pandas.DataFrame
        samples: list of dict

        """
        self.sampleTable = pd.read_table(self.config['sampletable'], sep='\t', dtype=str)
        self.sampleTable.set_index('sampleID', inplace=True)
        self.samples = self.sampleTable.reset_index().to_dict('records')

    def find_level(self, prefix):
        """ Figure out which regex the prefix matches.

        Scans each level and tries to identify which level the prefix matches.

        Parameters
        ----------
        prefix: str
            string that you want to match to the prefix, this string would have
            sample information filled in.

        Returns
        -------
        tuple:
            [0] is a string indicating which level
            [1] is a dict with sample information

        Raises
        ------
        ValueError
            Too lazy to make a custom exception, `raises` ValueError if the
            prefix does not match any of the patterns.

        Example
        -------
        >>> SH = SampleHandler(test_config)
        >>> prefix = 'pasilla/treated1/treated1_treated_1_0001_R1'
        >>> level, attrs = SH.find_level(prefix)

        >>> level
        'rawLevel'

        >>> assert attrs == [
        ... {'replicate': '1', 'sampleID': 'treated1', 'treatment': 'treated'}
        ... ]

        This prefix doesn't resolve to an actual sample defined in the
        sampletable (note mismatch between "1" in everything but the last
        "_2_". Should this raise ValueError?

        >>> prefix = 'pasilla_sample/treated1/treated1_treated_2_R1'
        >>> level, attrs = SH.find_level(prefix)
        >>> level
        'runLevel'
        >>> assert attrs == []

        Run-level prefix:

        >>> prefix = 'pasilla_sample/treated2/treated2_treated_2_R1'
        >>> level, attrs = SH.find_level(prefix)
        >>> level
        'runLevel'
        >>> assert attrs == [
        ... {'replicate': '2', 'sampleID': 'treated2', 'treatment': 'treated'}
        ... ]

        Agg-level prefix (note the multiple sets of attributes returned):

        >>> prefix = 'pasilla_agg/treated'
        >>> level, attrs = SH.find_level(prefix)
        >>> level
        'aggLevel'
        >>> assert attrs == [
        ... {'treatment': 'treated', 'replicate': '1', 'sampleID': 'treated1'},
        ... {'treatment': 'treated', 'replicate': '2', 'sampleID': 'treated2'}
        ... ]

        """
        if re.match(self.raw, prefix):
            return 'rawLevel', self.find_sample(self.raw, prefix)
        elif re.match(self.run, prefix):
            return 'runLevel', self.find_sample(self.run, prefix)
        elif re.match(self.sample, prefix):
            return 'sampleLevel', self.find_sample(self.sample, prefix)
        elif re.match(self.agg, prefix):
            return 'aggLevel', self.find_sample(self.agg, prefix)
        else:
            raise ValueError("Can't find a match for %s" % prefix)

    def find_sample(self, pattern, prefix):
        """ Find which sample(s) to use.

        Parameters
        ----------
        pattern: str
            A regex expression for the level

        prefix: str
            The prefix that is being matched to the pattern

        Returns:
        --------
        dict:
            Sample attributes corresponding to the sample(s) that the prefix
            contains.

        Example
        -------
        >>> SH = SampleHandler(test_config)

        >>> pattern = (
        ... 'pasilla_sample\\/(?P<sampleID>treated1|treated2|untreated1|untreated2)'
        ... '\\/(?P=sampleID)_(?P<treatment>treated|untreated)_(?P<replicate>1|2)_R1'
        ... )

        >>> prefix = 'pasilla_sample/treated2/treated2_treated_2_R1'

        >>> assert SH.find_sample(pattern, prefix) == [
        ... {'treatment': 'treated', 'sampleID': 'treated2', 'replicate': '2'}
        ... ]
        """
        m = re.match(pattern, prefix).groupdict()
        _samples = []
        for s in self.samples:
            if m.items() <= s.items():
                _samples.append(s)
        return _samples

    def make_input(self, prefix='prefix', midfix='', suffix='', agg=False):
        """ Generates Input Function based on wildcards.

        Notes
        -----

        example: '{first_element}{second_element}{third_element}'

        Parameters
        ----------
        prefix: str
            This can either be the entire prefix, or the name of the string
            formatting group that contains the prefix. This would be the
            'first_element' in the example.

        midfix: str
            This can either be the entire midfix, or the name of the string
            formating group that contains the midfix. This would be the
            'second_element' in the example.

        suffix: str
            This can either be the entire suffix, or the name of the string
            formating group that contains the suffix. This would be the
            'third_element' in the example.

        agg: bool
            True if you want file names from the level above the current
            prefix.

        Returns
        -------
        function:
            Retruns a snakemake input function that generates a list of files.

        Examples
        --------

        >>> SH = SampleHandler(test_config)
        >>> func = SH.make_input('prefix', '', '.fastq')
        >>> func({'prefix': 'pasilla'})
        ['pasilla.fastq']

        Instead of directly specifying the midfix, we pull it from the
        wildcards dict entry for "asdf":

        >>> func = SH.make_input('prefix', 'asdf', '.fastq')
        >>> func({'prefix': 'pasilla', 'asdf': '.cutadapt'})
        ['pasilla.cutadapt.fastq']

        We want to aggregate into "pasilla_agg/treated.fastq", but the
        following doesn't work because we need the directory slash:

        >>> func = SH.make_input('pasilla_agg', 'treated', '.fastq', agg=True)
        >>> assertRaises(ValueError, func, {})

        Also doesn't work:

        >>> func = SH.make_input('pasilla_agg/', midfix='treated', suffix='.fastq', agg=True)
        >>> assertRaises(ValueError, func, {})

        However if we provide most of the path in the prefix it works:

        >>> func = SH.make_input('pasilla_agg/treated', suffix='.fastq', agg=True)
        >>> assert sorted(func({})) == [
        ... 'pasilla_sample/treated1/treated1_treated_1_R1.fastq',
        ... 'pasilla_sample/treated1/treated1_treated_2_R1.fastq',
        ... 'pasilla_sample/treated2/treated2_treated_1_R1.fastq',
        ... 'pasilla_sample/treated2/treated2_treated_2_R1.fastq'
        ... ]
        """
        def _input(wildcards):

            try:
                _prefix = wildcards[prefix]
            except:
                _prefix = prefix

            try:
                _midfix = wildcards[midfix]
            except:
                _midfix = midfix

            try:
                _suffix = wildcards[suffix]
            except:
                _suffix = suffix

            if agg:
                level, _samples = self.find_level(_prefix)
                _sampleList = pd.DataFrame(_samples).to_dict('list')

                # Combine midfix and suffix and expand out any format strings
                _suffix = expand(_midfix + _suffix, **_sampleList, **self.config)

                # Retrun a list of sample ids by comining format string for higher level along with the full suffix.
                return list(set(expand(self.config[self.levelMap[level]] + '{suffix}', suffix=_suffix, **_sampleList, **self.config)))
            else:
                return expand('{prefix}{midfix}{suffix}', prefix=_prefix, midfix=_midfix, suffix=_suffix)

        return _input

    def build_targets(self, patterns):
        """ Build target file names based on pattern naming scheme.

        Given a list of string formatted patterns will use config information
        to fill in the patterns and generate a list of file targets.

        Parameters
        ----------
        patterns: list
            List of files with string formating marks that can be filled in
            from the config or the sampleTable.

        Returns
        -------
        list:
            Filled in list of file names.

        """
        _targets = []
        for p in patterns:
            p = p.format_map(self.config)
            for s in self.samples:
                e = dict(s, **self.config)
                _targets.append(p.format_map(e))
        return list(set(_targets))

if __name__ == "__main__":
    import doctest
    from unittest import TestCase

    def assertRaises(*args, **kwargs):
        return TestCase.assertRaises(None, *args, **kwargs)

    from textwrap import dedent
    test_sampletable = dedent("""\
    sampleID	treatment	replicate
    treated1	treated	1
    treated2	treated	2
    untreated1	untreated	1
    untreated2	untreated	2
    """)
    sampletable = os.path.join(
        os.path.dirname(__file__), 'test', 'test_sampletable.tsv')
    with open(sampletable, 'w') as fout:
        fout.write(test_sampletable)

    test_config = yaml.load("""
    sampletable: test/test_sampletable.tsv
    sample_dir: pasilla
    assembly: dm6
    data_dir: references

    rawLevel: 'pasilla/{sampleID}/{sampleID}_{treatment}_{replicate}_0001_R1'
    runLevel: 'pasilla_sample/{sampleID}/{sampleID}_{treatment}_{replicate}_R1'
    sampleLevel: 'pasilla_sample/{sampleID}/{sampleID}_{treatment}'
    aggLevel: 'pasilla_agg/{treatment}'
    aggLevel2:
        - '{treatment}'
    """)

    doctest.testmod()
    os.unlink('test/test_sampletable.tsv')
