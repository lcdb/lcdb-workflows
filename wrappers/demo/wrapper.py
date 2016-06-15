__author__ = "Your Name"
__copyright__ = "Copyright 2016, Your Name"
__email__ = "your@email.edu"
__license__ = "MIT"

# include any imports that will be needed, example below
from snakemake.shell import shell

# Use this block to support arbitrary arguments to be passed in to the shell
# call without raising an error if no params.extra were provided. If there
# aren't any additional params possible then you can delete this.
try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""

# Use this block to redirect stdout and stderr to a log if it was provided.
# Adjust the stdout/stderr redirection as needed.
if snakemake.log:
    log = "> {} 2>&1".format(snakemake.log)
else:
    log = ""

# Formatting: put logical blocks of arguments together on one line, but use
# multiple lines for clarity. Note the space at the end of each line for when
# the lines are concatenated. Also note the use of {extra} and {log} as
# specified above.
shell(
    "program-name "
    "-i {snakemake.input} "
    "-o {snakemake.output} "
    "{extra} "
    "{log} ")
