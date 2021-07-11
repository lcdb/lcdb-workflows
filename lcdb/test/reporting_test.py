import os
import sys
import unittest
from urllib.request import urlretrieve
from textwrap import dedent

import pandas as pd

test_imgURL = 'https://avatars1.githubusercontent.com/u/7584817?v=3&s=200'

test_imgB64 = ('data:image/png;base64,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')

test_imgB64_html = '<a href={0}><img src="{0}"/></a>'.format(test_imgB64)

test_imagePanel = dedent("""\
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

    </div>""".format(test_imgB64))

test_report = dedent("""\
    # The title of the report

    This is a test report

    Lorem ipsum dolor sit amet, alii noluisse vim no, mazim mediocritatem
    quo et, vis modo forensibus interpretaris ea. Et sea probatus atomorum
    mediocrem, ex mei saperet percipit. Nec vitae aeterno maiestatis ea,
    eripuit instructior ex vel. Est ad legere everti deseruisse, usu
    feugiat pertinax ut, pri dicam omnesque epicurei ut.

    Exerci laudem option pri te, vix modus inani signiferumque an, magna
    conclusionemque has ex. Vitae graeco quo ex, labore qualisque vel in, at mei
    eius nihil postea. Eum iudico audire eligendi ea, at eam malis aliquip. Vix at
    feugiat postulant ullamcorper, ius no prodesset repudiandae. Hendrerit
    complectitur mea ex.

    ## Section 1

    Exerci laudem option pri te, vix modus inani signiferumque an, magna
    conclusionemque has ex. Vitae graeco quo ex, labore qualisque vel in, at mei
    eius nihil postea. Eum iudico audire eligendi ea, at eam malis aliquip. Vix at
    feugiat postulant ullamcorper, ius no prodesset repudiandae. Hendrerit
    complectitur mea ex.

    {Panel1}

    ## Section 2

    ### Sub section 2.1

    {Image1}

    Exerci laudem option pri te, vix modus inani signiferumque an, magna
    conclusionemque has ex. Vitae graeco quo ex, labore qualisque vel in, at mei
    eius nihil postea. Eum iudico audire eligendi ea, at eam malis aliquip. Vix at
    feugiat postulant ullamcorper, ius no prodesset repudiandae. Hendrerit
    complectitur mea ex.

    {Table1}

    ### Sub section 2.2

    Bulleted list

    * one
    * two
    * thee

    Sed sumo posse liberavisse ex. Eu ponderum neglegentur cum, ius eu tractatos
    partiendo. Albucius explicari vim et. Id maluisset argumentum vim. Ex mei summo
    conceptam, ut adipisci dignissim disputando per. At detracto intellegebat vim,
    mea nihil everti at. Wisi error vis no, pri choro ceteros ne.

    ### Sub section 2.3

    Exerci laudem option pri te, vix modus inani signiferumque an, magna
    conclusionemque has ex. Vitae graeco quo ex, labore qualisque vel in, at mei
    eius nihil postea. Eum iudico audire eligendi ea, at eam malis aliquip. Vix at
    feugiat postulant ullamcorper, ius no prodesset repudiandae. Hendrerit
    complectitur mea ex.""")

test_table = pd.DataFrame({'sampleID': ['treated1', 'treated2'], 
                           'treatment': ['trt1', 'trt2'], 
                           'replicate': ['1', '2']
                           })[['sampleID', 'treatment', 'replicate']]

def renderHtml(string):
        from lcdb.reporting import ENV
        with open('tmp_render.html', 'w') as OUT:
            tp = ENV.get_template('base.html')
            OUT.write(tp.render(content=string))
            

class TestReportImage(unittest.TestCase):
    def setUp(self):
        self.imageFile, _ = urlretrieve(test_imgURL)

    def tearDown(self):
        os.unlink(self.imageFile)

    def test_Image(self):
        from lcdb.reporting import Image
        img = Image(self.imageFile)
        self.assertEqual(str(img), test_imgB64)


class TestReportPanel(unittest.TestCase):
    def setUp(self):
        self.imageFile, _ = urlretrieve(test_imgURL)

    def tearDown(self):
        os.unlink(self.imageFile)

    def test_jinja_panel_with_caption(self):
        from lcdb.reporting import JinjaPanel, Image, ENV
        img1 = Image(self.imageFile, name='Sample1', caption='male')
        img2 = Image(self.imageFile, name='Sample2', caption='female')
        img3 = Image(self.imageFile, name='Sample3', caption='male')

        table = JinjaPanel()
        table.add_row([img1, img2, img3], caption='Row caption')
        table.add_row([img1, img2, img3], caption='Row caption2')
        self.assertEqual(table.panel, test_imagePanel)


class Test_report(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_report(self):
        from lcdb.reporting import report


        values = {'Panel1': dedent(test_imagePanel),
                  'Image1': dedent(test_imgB64_html),
                  'Table1': test_table.to_html(index=False, classes='table table-striped')
                }

        with open('tmp_render.html', 'w') as OUT:
            OUT.write(report(test_report, **values))

if __name__ == '__main__':
    unittest.main()
