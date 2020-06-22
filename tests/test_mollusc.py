
from itertools import count
import unittest

import numpy as np

from mollusc.mollusc import Mollusc


class MolluscTestCase(unittest.TestCase):
    def setUp(self):
        self.mollusc = Mollusc()
        
        
    def test_getset_pars(self):
        self.mollusc['ga'] = 9
        self.assertEqual(self.mollusc['ga'], 9)
        
    
    def test_argfilter(self):
        gen_conc = dict(zip(['g'+c for c in  'abcdef'], count(1)))      
        
        # Empty Mollusc:
        filtered_args = self.mollusc.argfilter('g', **gen_conc)
        expected_args = np.array([[]]).T
        
        self.assertEqual(filtered_args.shape, expected_args.shape)
        
        # Set number of substances:
        self.mollusc['kn'] = 6
        
        filtered_args = self.mollusc.argfilter('g', **gen_conc)
        expected_args = np.array([[1],[2],[3],[4],[5],[6]])
        
        self.assertEqual(filtered_args.shape, expected_args.shape)
        self.assertTrue((filtered_args==expected_args).all())
        

        
            
if __name__ == '__main__':
    unittest.main()
