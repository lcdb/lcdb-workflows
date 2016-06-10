from jsonschema import Draft4Validator, validators
import helpers
import yaml
from collections import OrderedDict
from textwrap import wrap as _wrap


def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    """
    Load YAML into an ordered dictionary to maintain key sorting.
    """
    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)


def represent_odict(dump, tag, mapping, flow_style=None):
    """
    Dump an ordered dictionary to YAML, maintaining the key order but making it
    look like a normal dictionary (without the !!python/object extra stuff).

    From https://gist.github.com/miracle2k/3184458
    """
    value = []
    node = yaml.MappingNode(tag, value, flow_style=flow_style)
    if dump.alias_key is not None:
        dump.represented_objects[dump.alias_key] = node
    best_style = True
    if hasattr(mapping, 'items'):
        mapping = mapping.items()
    for item_key, item_value in mapping:
        node_key = dump.represent_data(item_key)
        node_value = dump.represent_data(item_value)
        if not (isinstance(node_key, yaml.ScalarNode) and not node_key.style):
            best_style = False
        if (
            not (isinstance(node_value, yaml.ScalarNode) and not
                 node_value.style
                 )
        ):
            best_style = False
        value.append((node_key, node_value))
    if flow_style is None:
        if dump.default_flow_style is not None:
            node.flow_style = dump.default_flow_style
        else:
            node.flow_style = best_style
    return node

yaml.SafeDumper.add_representer(
    OrderedDict, lambda dumper, value: represent_odict(
        dumper, u'tag:yaml.org,2002:map', value))


def access(dct, keys):
    """
    Access a value from an arbitrarily-nested dictionary, given a set of keys.

    If any key doesn't exist, returns None.

    >>> access({'a':{'z': {'y': {'b': 1, 'c': 2}}}}, keys=['a', 'z', 'y', 'b'])
    1
    """
    o = dct
    for k in keys:
        o = o.get(k)
        if o is None:
            return None
    return o


def follow_ref(ref, dct):
    """
    Follow a "$ref" JSON Schema reference
    """
    ref = ref.lstrip('#/')
    keys = ref.split('/')
    return access(dct, keys)


def _indent(s, amount):
    """
    Indents string `s` by `amount` doublespaces.

    Assumes, e.g., vim sw=2 ts=2 indentation.
    """
    pad = '  ' * amount
    return '\n'.join([pad + i for i in s.splitlines(False)])


def create_config(schema, genome, site_file, context=None, fout=None):
    """
    Generates an example config file based on the schema alone and writes it to
    `fout`, which can then be validated.

    The goal is to have the schema be the sole source of config files.

    "description" fields in the schema will be printed out as YAML comments,
    directly above the entry, thereby serving as documentation.  "default"
    fields in the schema will be used to fill in the generated YAML. "use_site"
    means that the site.yaml file will be accessed for that chain of keys (see
    the `access()` function above), which allows us to fill in info like
    indexes and GTF files in a genome-specific manner.

    Parameters
    ----------
    schema: str
        Filename of config schema

    genome : str
        Used as a top-level key to access the site info

    site_file : str
        Path to "site" info, which is a YAML file describing default locations
        of files in a genome-specific manner.

        Currently assumes a format like this, where values are typically
        filenames::

            dm6:
                fastq_screen:
                    config: fastq_screen.conf
                fasta: /data/dm6.fa
                indexes:
                    bowtie2: /data/indexes/bowtie2
                    hisat2: /data/indexes/hisat2
                annotations:
                    gene: /data/annotations/dm6.genes.gtf
                    rRNA_fasta: /data/rDNA.fa
                    refflat: /data/annotations/dm6.refflat
            hg38:
                fastq_screen:
                    config: fastq_screen.conf
                fasta: /data/hg38.fa
                indexes:
                    bowtie2: /data/indexes/bowtie2
                    hisat2: /data/indexes/hisat2
                annotations:
                    gene: /data/annotations/hg38.genes.gtf
                    rRNA_fasta: /data/rDNA.fa
                    refflat: /data/annotations/hg38.refflat


    context : dict
        Context dictionary. This is used as a way of getting custom information
        into the schema. If a schema key has the property "use_context: True",
        then that key will be used to get a value from `context`.

    fout : file-like
        Output file to write YAML to.
    """
    if context is None:
        context = {}

    context['genome'] = genome

    d = ordered_load(open(schema), yaml.SafeLoader)

    site = yaml.load(open(site_file))

    def props(path, v, fout=None, print_key=True):
        """
        Recursively print out a filled-in config file based on the schema

        Parameters
        ----------
        path : list
            List of keys to access from `v`

        v : dict
            Dictionary from which to access `path`

        fout : file handle
            Output will be written here

        print_key : boolean


        """

        # Set the full original and output file objects, but only on the first
        # time.
        try:
            props.orig
        except AttributeError:
            # so now we can always access the full original YAML dict
            props.orig = v

        try:
            props.out
        except AttributeError:
            props.out = fout

        # Key into schema. If a path is provided, we start with the last one.
        if path:
            k = path[-1]
        else:
            k = None

        # First thing to do is write the comments for whatever value we're on.
        # This will include the description field (if provided) as well as any
        # enum options. Much of the complexity here is in making the comments
        # be indented at the proper level, which is tracked by the `level`
        # attribute of the function

        def wrap(s):
            "wrap a comment string at the current indent level"
            return ("\n%s# " % (indent)).join(_wrap(s))

        indent = '  ' * props.level

        # Print out the description as a comment.
        if 'description' in v:
            props.out.write('\n%s# %s\n' % (indent, wrap(v['description'])))

        # Describe the possible values of any enums.
        if 'enum' in v:
            if len(v['enum']) > 1:
                enum = '\n'.join(['%s# - "%s"' % (indent, i) for i in v['enum']])
                props.out.write(
                    '{indent}# options for "{k}" are:\n{enum}\n'
                    .format(**locals()))

        default = v.get('default')

        # If "use_site" key exists, override any defaults with the
        # corresponding entry in site.yaml.  "Corresponding" means having the
        # same "path" of keys, with 'genome' prepended. So if the schema has an
        # entry like::
        #
        #   dm6:
        #       indexes:
        #           bowtie:
        #               type: string
        #               description: bowtie index
        #               use_site: True
        #
        #
        # Then we will make the default value in the generated config be the
        # value of::
        #
        #   site['dm6']['indexes']['bowtie']
        #
        #
        # If "use_context" is True, then we assume that this key was passed in
        # as part of the parent function and is available in the current
        # context dictionary.
        if v.get('use_site'):
            default = access(site, [genome] + path)
        if v.get('use_context'):
            default = context[k]

        # Create an appropriate prefix to be inserted right after the ":". This
        # avoids trailing spaces for keys with no value; this also sets strings
        # with no default to the empty string  "" instead of a blank space
        # (since a blank space would fail string validation).
        prefix = " "
        if default is None:
            if v['type'] == "string":
                default = '""'
            else:
                default = ""
                prefix = ""

        else:
            # It's tricky to get arrays to be generated using $ref references;
            # this allows a default to be written into the schema and printed
            # nicely in the generated output.
            if v['type'] in ['object', 'array']:
                default = '\n' + _indent(
                    yaml.safe_dump(default, default_flow_style=False, indent=2,
                                   line_break=False), props.level + 2)

        if not print_key:
            k = ""
            colon = ""
            indent = ""
        else:
            colon = ":"

        # Write out what we have so far.
        props.out.write(
            '{indent}{k}{colon}{prefix}{default}\n'.format(**locals()))

        # Recursively call this function for everything in properties.
        # These are the only things that are printed to the output YAML;
        # everything else is only used for schema validation.
        if 'properties' in v:
            for k, v in v['properties'].items():

                # indent the next output by one level
                props.level += 1

                # stick the key onto the end of the path
                path.append(k)

                # recursively call
                props(path, v)

                # now we back out, decreasing the indent level and removing the
                # key we just provided.
                props.level -= 1
                path.pop()

        # Follow references if needed for an array
        if v['type'] == 'array':
            if 'default' not in v:

                # TODO: is this really needed?
                props.level += 1
                props.out.write('%s-' % ('  ' * props.level))
                props.level -= 1

                # TODO: should follow $ref for non-array items as well.
                if '$ref' in v['items']:
                    props.level += 1
                    props(path,
                          follow_ref(v['items']['$ref'],
                                     props.orig),
                          print_key=False)
                    props.level -= 1

        return

    # Set the default level and then recursively call props to generate the
    # YAML file.
    #
    props.level = -1
    props([], d, fout, print_key=False)
    return

if __name__ == "__main__":
    import argparse
    import sys
    from io import StringIO

    ap = argparse.ArgumentParser()
    ap.add_argument('--site')
    ap.add_argument('--schema')
    ap.add_argument('--genome')
    ap.add_argument('--out')
    args = ap.parse_args()

    with open(args.out, 'w') as fout:
        create_config(
            schema=args.schema,
            genome=args.genome,
            site_file=args.site,
            fout=fout)

    helpers.validate_config(fout.name, args.schema)
