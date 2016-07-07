import re
from blaze import resource, DataFrame
import pandas as pd

@resource.register('.+\.bwt2.log')
def resource_bowtie2_log(uri, **kwargs):
    with open(uri) as fh:
        data = "".join(fh)

    index = ['Number of Reads', 'Number Unpaired',
             'Number Unaligned', 'Number Uniquely Aligned',
             'Number Ambiguously Aligned']

    values = []
    for row in data.strip().split('\n'):
        if not row.startswith('Warn'):
            values.append(re.sub(r' *(\d+)[ %]+.+', r'\1', row))

    df = DataFrame([int(x) for x in values[:-1]], index=index)
    df.index.name = 'statistic'
    df.columns = ['counts']
    return df
