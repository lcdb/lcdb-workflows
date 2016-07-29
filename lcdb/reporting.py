#!/usr/bin/env python
""" Reporting System """

import os
import sys
import base64
from textwrap import dedent

import markdown
import jinja2

TEMPLATES_DIR = os.path.join(os.path.dirname(__file__), 'templates')
loader = jinja2.FileSystemLoader(TEMPLATES_DIR)
ENV = jinja2.Environment(loader=loader)

def report(string, **kwargs):
    """ Build a report """
    md = markdown.Markdown(output_format='html5')
    tp = ENV.get_template('base.html')
    return tp.render(content=md.convert(dedent(string).format(**kwargs)))


class Report(object):
    """ Basic report class """
    pass


class Section(object):
    """ Basic section class """
    pass


class Image(object):
    """ Image class """

    def __init__(self, png, name='', caption=''):
        with open(png, 'rb') as fh:
            uri = base64.b64encode(fh.read()).decode('utf-8').replace('\n', '')
            self.image = 'data:image/png;base64,{uri}'.format(uri=uri)
            self.image_html = '<a href={0}><img src="{0}"/></a>'.format(self.image)

            self.name = name
            self.caption = caption

    def __repr__(self):
        return self.image

    def __str__(self):
        return self.image


class JinjaPanel(object):
    """ Image panel class """

    def __init__(self):
        self.tp = ENV.get_template('panel.html')
        self.panel = ''

    def add_row(self, images, caption=''):
        self.panel += self.tp.render(caption=caption, images=images)

    def __repr__(self):
        return self.panel

    def __str__(self):
        return self.panel


if __name__ == '__main__':
    import doctest
    from textwrap import dedent
    from unittest import TestCase

    def assertRaises(*args, **kwargs):
        return TestCase.assertRaises(None, *args, **kwargs)
