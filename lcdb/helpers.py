import yaml
from jsonschema import validate, ValidationError


def validate_config(config, schema):
    schema = yaml.load(open(schema))
    cfg = yaml.load(open(config))
    try:
        validate(cfg, schema)
    except ValidationError as e:
        msg = '\nPlease fix %s: %s\n' % (config, e.message)
        raise ValidationError(msg)
