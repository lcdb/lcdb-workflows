import re

rawLevel = '/data/LCDB/users/fearjm/lcdb-workflow/test/passilla'
runLevel = '/data/LCDB/users/fearjm/lcdb-workflow/test/passilla_sample/{sample}/{}'
sampleLevel = '/data/LCDB/users/fearjm/lcdb-workflow/test/passilla_sample/{sample}/{}'
aggLevel = '/data/LCDB/users/fearjm/lcdb-workflow/test/passilla_sample/{}'

class SampleHandler(object):
    """ """
    def __init__(self, rawLevel, runLevel, sampleLevel, aggLevel):
        self.raw = re.compile(rawLevel)
        self.run = re.compile(runLevel)
        self.sample = re.compile(sampleLevel)
        self.agg = re.compile(aggLevel)

    def generatePatter(self):
        pass
