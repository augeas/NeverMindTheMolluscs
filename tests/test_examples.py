
import unittest

from mollusc.examples import MolluscExample, sp_pars


class ExamplesTestCase(unittest.TestCase):
    
    def setUp(self):
        self.all_equations = set(('eg_{}'.format(pars['ke'])
            for pars in sp_pars.values()))
        self.mollusc_dir = dir(MolluscExample)


    def test_all_equations_implemented(self):
        for eq in self.all_equations:
            self.assertTrue(eq in self.mollusc_dir)
            print('{} is implemented'.format(eq))
            
