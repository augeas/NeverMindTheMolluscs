
from io import BytesIO
from itertools import islice

from imageio import imwrite
import numpy as np

from .equations import MolluscMixin


class Mollusc(MolluscMixin, object):
    """Implementation of the code from "The Algorithmic Beauty of Sea Shells"
    (https://doi.org/10.1007/978-3-662-03617-4) by the late Hans Meinhardt
    (https://www.eb.tuebingen.mpg.de/emeriti/hans-meinhardt).
    
    Attributes
    ----------
    num_vars: int
        The number of substances represented by the simulation.
    
    concentrations: numpy.ndarray
        A numpy array of the concentrations of the substances in the
        simulation, of dimension width * n_vars * height, if populated by calls
        to trace or render. This can be accessed by indexing the Mollusc object
        directly.
    
    Methods
    ------
    trace(height=None)    
        Populates the "concentrations" array. 
    
    render(height=None, redraw=False, substances=[])
        Returns a bytes object containing a .png rendering of the
        concentrations.
    """
    
    def argfilter(self, initial, **kwargs):
        linitial = initial.lower()
        args = filter(lambda k: len(k)==2 and k.lower().startswith(linitial), kwargs.keys())
        return np.array([[kwargs[k] for k in sorted(args)[0:self.num_vars]]]).T
    
    
    def __init__(self,width=640, **kwargs):
        """
        """
        self.width = width
        
        self.num_vars = kwargs['kn']
        self.title = kwargs.get('title')
        
        self.diff_coeffs = self.argfilter('d', **kwargs)
        self.decay_rates = self.argfilter('r', **kwargs)
        self.basic_prod = self.argfilter('b', **kwargs)
        self.saturation = self.argfilter('s', **kwargs)
        self.coupling = self.argfilter('c', **kwargs)
        self.init_conc = self.argfilter('a', **kwargs)
        self.gen_conc = self.argfilter('g', **kwargs)

        self.activator_saturation = self.saturation[0]

        self.substrate = np.ones(width)

        self.decay_diff_consts = 1.0 - self.decay_rates - 2 * self.diff_coeffs

        eq_methods = {method:getattr(self, method)
            for method in dir(self) if method.startswith('eq_')}

        self.eq = eq_methods['eq_'+str(kwargs['ke'])]
        
        self.skip = kwargs['kp']
        
        self.perturb = kwargs.get('kr', 2) / 100.0
        self.fluctuations = self.fluctuate()   
        self.s = self.decay_rates[0] * self.fluctuate()
    
        self.old = 0
        self.new = 1
        
        self.condition_pars = [kwargs.get('k'+str(i),0) for i in range(1,5)]
        
        self.cells = np.zeros((2,self.num_vars,self.width))
        self.cells[0,:,:] = self.gen_conc

        self.concentrations = None
        
        initial_conditions =  kwargs.get('ki', 3)
        
        if initial_conditions == 1:
            self.cells[self.old,:,0] = self.init_conc.T
        elif initial_conditions == 2:
            self.cells[self.old,:,0] = self.init_conc.T
            self.cells[self.old,:,-1] = self.init_conc.T
        elif initial_conditions == 6:
            selected = np.random.random((self.width,)) > 1 / 20
            self.cells[self.old,:,selected] = self.init_conc.T
        else:
            self.cells[self.old,:,self.width//2] = self.init_conc.T
    

    def __getitem__(self, key):
        if self.concentrations is None:
            raise IndexError()
        else:
            return self.concentrations.__getitem__(key)


    def __str__(self):
        return self.title.strip()


    def fluctuate(self):
        return 1.0 + self.perturb * (np.random.random(self.width) - 0.5)


    def stabilize(self,*args):
        substances = np.array(sorted(args))
        to_zero = self.cells[self.new,substances,:] < 0.0
        self.cells[self.new,substances,:][to_zero] = 0.0

        
    def step(self):
        old_cells = self.cells[self.old,:,:]
        new_cells = self.cells[self.new,:,:]


        l_bound = np.copy(self.cells[self.old,:,0])
        r_bound = np.copy(self.cells[self.old,:,-1])

        left_cells = np.roll(old_cells, -1, axis=1)
        right_cells = np.roll(old_cells, 1, axis=1)

        left_cells[:,0] = l_bound
        right_cells[:,-1] = r_bound

        decay = old_cells * self.decay_diff_consts + self.diff_coeffs * (left_cells + right_cells)
        
        a_sq = np.square(old_cells[0,:])
        
        self.eq(old_cells, new_cells, decay, a_sq)
        
        self.old, self.new = self.new, self.old
        

    def __iter__(self):
        return self
    
    
    def __next__(self):
        for i in range(self.skip):
            self.step()
        return np.copy(self.cells[self.new,:,:])


    def trace(self, height=None):
        """Populates the "concentrations" array with dimensions
        width * height * n_vars, defaulting to a square array.
        
        Parameters
        ----------
            height: int, optional
                The height of the concentrations array, defaulting to thw width
                of simulation.
        """
        if not height:
            rows = self.width
        else:
            rows = height
        self.concentrations = np.stack(list(islice(self, 0, rows)))
    
    
    def render(self, height=None, redraw=False, substances=[]):
        if redraw or self.concentrations is None:
            self.trace(height)
        
        width, _, height = self.concentrations.shape
        
        image = np.zeros(shape=(height, width, 3))

        if not substances:
            substances = tuple(range(min(self.num_vars, 3)))
        
        for i, col in enumerate(substances[0:3]):
            if col:
                image[:,:,i] = np.tanh(self.concentrations[:, col, :])
                image[:,:,i] = 255 * image[:,:,i] / image[:,:,i].max()
        
        im_bytes = BytesIO()
        imwrite(im_bytes, image.astype(np.uint8), format='png')
        
        return im_bytes.getvalue()

