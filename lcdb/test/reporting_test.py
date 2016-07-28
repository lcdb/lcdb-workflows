import os
import sys
import unittest
from urllib.request import urlretrieve
from textwrap import dedent

imgURL = 'https://avatars1.githubusercontent.com/u/7584817?v=3&s=200'
imgB64 = ('data:image/png;base64,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')

def writeHtml(string):
    with open('tmp.html', 'w') as OUT:
        OUT.write(dedent("""<html>
        <body>
        {0}
        </body>
        </html>""".format(string)))

def renderHtml(string):
        from lcdb.reporting import ENV
        with open('tmp_render.html', 'w') as OUT:
            tp = ENV.get_template('base.html')
            OUT.write(tp.render(content=string))
            

class TestReport(unittest.TestCase):

    def setUp(self):
        self.imageFile, _ = urlretrieve(imgURL)

    def tearDown(self):
        os.unlink(self.imageFile)

    def test_Image(self):
        from lcdb.reporting import Image
        img = Image(self.imageFile)
        self.assertEqual(str(img), imgB64)

    def test_panel_add_row(self):
        from lcdb.reporting import Panel
        table = Panel()
        table.add_row([imgB64, imgB64, imgB64])

        exampleTable = dedent("""<table>
        <tr>
        <td>{img}</td> <td>{img}</td> <td>{img}</td>
        </tr>
        </table>""".format(img=imgB64))

        self.assertEqual(str(table), exampleTable)

    def test_panel_with_caption(self):
        from lcdb.reporting import Panel
        table = Panel()
        table.add_row([imgB64, imgB64, imgB64], caption='Row caption')
        table.add_row([imgB64, imgB64, imgB64], caption='Row caption2')

        exampleTable = dedent("""<table>
        <tr>
        <td>{img}</td> <td>{img}</td> <td>{img}</td>
        </tr>
        </table>""".format(img=imgB64))

        writeHtml(str(table))

        #self.assertEqual(str(table), exampleTable)

    def test_jinja_panel_with_caption(self):
        from lcdb.reporting import JinjaPanel, Image, ENV
        img1 = Image(self.imageFile, name='Sample1', caption='male')
        img2 = Image(self.imageFile, name='Sample2', caption='female')
        img3 = Image(self.imageFile, name='Sample3', caption='male')
        table = JinjaPanel()
        table.add_row([img1, img2, img3], caption='Row caption')
        table.add_row([img1, img2, img3], caption='Row caption2')
        renderHtml(table.panel)



if __name__ == '__main__':
    unittest.main()
