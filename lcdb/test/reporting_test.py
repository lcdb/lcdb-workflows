import os
import sys
import unittest
from urllib.request import urlretrieve
from textwrap import dedent

imgURL = 'https://avatars1.githubusercontent.com/u/7584817?v=3&s=200'

imgB64 = ('data:image/png;base64,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')

mdownText = dedent("""\
        # Heading Level 1

        This is just a test text document that will be translated into HTML
        upon rendering. *Some example italics* and maybe some **bold**. 

        * simple bullets 1
        * simple bullets 2
        * simple bullets 3

        ## Heading Level 2

        1. now some numbers
        2. now some numbers
        3. now some numbers

        ### Heading level 3

        Some more text.
        """)


def renderHtml(string):
        from lcdb.reporting import ENV
        with open('tmp_render.html', 'w') as OUT:
            tp = ENV.get_template('base.html')
            OUT.write(tp.render(content=string))
            

class TestReportImage(unittest.TestCase):
    def setUp(self):
        self.imageFile, _ = urlretrieve(imgURL)

    def tearDown(self):
        os.unlink(self.imageFile)

    def test_Image(self):
        from lcdb.reporting import Image
        img = Image(self.imageFile)
        self.assertEqual(str(img), imgB64)


class TestReportPanel(unittest.TestCase):
    def setUp(self):
        self.imageFile, _ = urlretrieve(imgURL)

    def tearDown(self):
        os.unlink(self.imageFile)

    def test_jinja_panel_with_caption(self):
        from lcdb.reporting import JinjaPanel, Image, ENV
        img1 = Image(self.imageFile, name='Sample1', caption='male')
        img2 = Image(self.imageFile, name='Sample2', caption='female')
        img3 = Image(self.imageFile, name='Sample3', caption='male')

        examplePanel = dedent("""\
        <div class="row">
            <div class="legend">Row caption</div>
 
                <div class="col-sm-3 thumbnail figure">
                    <div class="legend">Sample1</div>
                    <a href="{0}"><img src="{0}" /></a>
                    <div class="legend">male</div>
                </div>
 
                <div class="col-sm-3 thumbnail figure">
                    <div class="legend">Sample2</div>
                    <a href="{0}"><img src="{0}" /></a>
                    <div class="legend">female</div>
                </div>
 
                <div class="col-sm-3 thumbnail figure">
                    <div class="legend">Sample3</div>
                    <a href="{0}"><img src="{0}" /></a>
                    <div class="legend">male</div>
                </div>
 
        </div><div class="row">
            <div class="legend">Row caption2</div>
    
                <div class="col-sm-3 thumbnail figure">
                    <div class="legend">Sample1</div>
                    <a href="{0}"><img src="{0}" /></a>
                    <div class="legend">male</div>
                </div>
 
                <div class="col-sm-3 thumbnail figure">
                    <div class="legend">Sample2</div>
                    <a href="{0}"><img src="{0}" /></a>
                    <div class="legend">female</div>
                </div>
 
                <div class="col-sm-3 thumbnail figure">
                    <div class="legend">Sample3</div>
                    <a href="{0}"><img src="{0}" /></a>
                    <div class="legend">male</div>
                </div>

        </div>""".format(str(img1)))

        table = JinjaPanel()
        table.add_row([img1, img2, img3], caption='Row caption')
        table.add_row([img1, img2, img3], caption='Row caption2')
        self.assertEqual(dedent(table.panel), examplePanel)


class TestReportText(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_initialize(self):
        from lcdb.reporting import Text
        txt = Text(mdownText)

        exampleHtml = dedent("""\
                <h1>Heading Level 1</h1>
                <p>This is just a test text document that will be translated into HTML
                upon rendering. <em>Some example italics</em> and maybe some <strong>bold</strong>. </p>
                <ul>
                <li>simple bullets 1</li>
                <li>simple bullets 2</li>
                <li>simple bullets 3</li>
                </ul>
                <h2>Heading Level 2</h2>
                <ol>
                <li>now some numbers</li>
                <li>now some numbers</li>
                <li>now some numbers</li>
                </ol>
                <h3>Heading level 3</h3>
                <p>Some more text.</p>""")

        self.assertEqual(dedent(txt.html), exampleHtml)


if __name__ == '__main__':
    unittest.main()
