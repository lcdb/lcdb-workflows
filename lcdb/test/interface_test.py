import unittest

class TestSampleHandler(unittest.TestCase):
    def setUp(self):
        from interface import SampleHander
        SH = SampleHander()
        pass

    def tearDown(self):
        pass

    def test_bob(self):
        self.assertFalse(True)
