import unittest
import sys
import os
import yaml

# This may not be needed if tests are always ran from project folder, but
# guarantees that paths are correct.
HERE = os.path.dirname(os.path.realpath(__file__))
os.chdir(HERE)

configYaml = """
sampletable: ../../test/pasilla_sampletable.tsv
sample_dir: ../../pasilla
assembly: dm6
data_dir: /data/LCDB/references

rawLevel: '../../test/passilla/{sampleID}/{sampleID}_R1'
runLevel: '../../test/passilla_sample/{sampleID}/{sampleID}_R1'
sampleLevel: '../../test/passilla_sample/{sampleID}/{sampleID}'
aggLevel: '../../test/passilla_agg'
aggLevel2:
    - '../../test/{treatment}'
"""

samples = [
    {'sampleID': 'treated1', 'replicate': 1, 'treatment': 'treated'}, 
    {'sampleID': 'treated2', 'replicate': 2, 'treatment': 'treated'}, 
    {'sampleID': 'untreated1', 'replicate': 1, 'treatment': 'untreated'}, 
    {'sampleID': 'untreated2', 'replicate': 2, 'treatment': 'untreated'} 
]


class TestSampleHandler(unittest.TestCase):
    def setUp(self):
        from interface import SampleHandler

        # Load config
        self.config = yaml.load(configYaml)

        # Initialize the class
        self.SH = SampleHandler(self.config)

    def tearDown(self):
        pass

    def test_init(self):
        self.assertEqual(self.config, self.SH.config)
        self.assertEqual(samples, self.SH.samples)

    def test_find_level_raw(self):
        prefix = self.config['rawLevel'].format_map(samples[0])
        self.assertEqual('rawLevel', self.SH.find_level(prefix))

    def test_find_level_run(self):
        prefix = self.config['runLevel'].format_map(samples[0])
        self.assertEqual('runLevel', self.SH.find_level(prefix))

    def test_find_level_sample(self):
        prefix = self.config['sampleLevel'].format_map(samples[0])
        self.assertEqual('sampleLevel', self.SH.find_level(prefix))

    def test_find_level_bad_pattern(self):
        prefix = '../../test/{sampleID}/{sampleID}/{sampleID}'.format_map(samples[0])
        with self.assertRaises(ValueError):
            self.SH.find_level(prefix)

    def test_find_sample(self):
        prefix = self.config['runLevel'].format_map(samples[0])
        print(prefix)
        _samples = self.SH.find_sample(prefix)
        print(_samples)
        self.assertEqual(_samples, samples[0])

    def test_make_input(self):
        wildcards = {'prefix': self.config['runLevel'].format_map(self.SH.samples[0])}
        _input = self.SH.make_input(suffix='.fastq', agg=True)
        self.assertEqual(_input(wildcards), [wildcards['prefix'] + '.fastq'])

