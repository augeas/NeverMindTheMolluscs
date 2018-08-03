import numpy as np

class Mollusc(object):
    
    def argfilter(self,initial,**kwargs):
        linitial = initial.lower()
        args = filter(lambda k: len(k)==2 and k.lower().startswith(linitial), kwargs.keys())
        return np.array([[kwargs[k] for k in sorted(args)[0:self.num_vars]]]).T
    
    
    def __init__(self,width=640,**kwargs):
        self.width = width
        
        self.num_vars = kwargs['kn']
        
        self.diff_coeffs = self.argfilter('d', **kwargs)
        self.decay_rates = self.argfilter('r', **kwargs)
        self.basic_prod = self.argfilter('b', **kwargs)
        self.saturation = self.argfilter('s', **kwargs)
        self.coupling = self.argfilter('c', **kwargs)
        self.init_conc = self.argfilter('a', **kwargs)
        self.gen_conc = self.argfilter('g', **kwargs)

        self.activator_saturation = 0.0

        self.decay_diff_consts = 1.0 - self.decay_rates - 2 * self.diff_coeffs

        eq_methods = {method:getattr(self, method)
            for method in dir(self) if method.startswith('eq_')}

        self.eq = eq_methods['eq_'+str(kwargs['ke'])]
        
        perturb = kwargs.get('kr',2) / 100.0
        fluctuations = 1.0 + perturb * (np.random.random(self.width) - 0.5)
    
        self.s = self.decay_rates[0] * fluctuations
    
        self.tick = 0
        
        self.cells = np.zeros((2,self.num_vars,self.width))
        self.cells[0,:,:] = self.gen_conc
        
        self.cells[0,:,self.width//2] = self.init_conc[:,0]
        
    def decay_diff(self, old_cells):
        left_cells = np.roll(old_cells,-1,axis=1)
        right_cells = np.roll(old_cells,1,axis=1)
        return old_cells * self.decay_diff_consts + self.diff_coeffs * (left_cells + right_cells)   
        
    def step(self):
        
            if self.tick:
                old = 1
                new = 0
                self.tick = 0
            else:
                old = 0
                new = 1
                self.tick = 1

            old_cells = self.cells[old]
            new_cells = self.cells[new]
            decay = self.decay_diff(old_cells)
            self.eq(old_cells,new_cells,decay)

    def __iter__(self):
        while True:
            yield self.cells[self.tick,:,:]
            self.step()
    
    
    def eq_21(self,old_cells,new_cells,decay):
        """- activator - inhibitor mechanism: B is inhibitor"""
        a_sq = np.square(old_cells[0])
        new_cells[0] = decay[0] + self.s * (a_sq / old_cells[1] + self.basic_prod[0])
        new_cells[1] = decay[1] + self.s * a_sq + self.basic_prod[1]

#CASE 21 '- activator - inhibitor mechanism: B is inhibitor --------------
#   FOR i = ja TO js: GOSUB olddecay:
#     axt(1, i) = olddecaydiffA + s * (a * a / b + ba)
#     axt(2, i) = olddecaydiffB + s * a * a + bb
#   NEXT i

        
    def eq_61(self,old_cells,new_cells,decay):
        """Hormone (c) changes lifetime of the inhibitor"""
        a_sq = np.square(old_cells[0])
        aq = self.s * a_sq / (1 + self.activator_saturation * a_sq) + self.basic_prod[0]
        new_cells[0] = decay[0] + aq / (self.saturation[1] + old_cells[1])
        new_cells[1] = decay[1] + aq + self.basic_prod[1]
        #mean_hormone = self.decay_rates[2] * old_cells[0]
        
        #self.decay_diff_consts[1] = 1.0 - 2 * 
        
#CASE 61 '-- Branches controlled by a hormone  : Olivia Porphyria ----------
#     '  Hormone (c) changes lifetime of the inhibitor
#   FOR i = ja TO js: GOSUB olddecay:
#     aq = s * a * a / (1! + sA * a * a) + ba
#     axt(1, i) = olddecaydiffA + aq / (sb + b)
#     axt(2, i) = olddecaydiffB + aq + bb
#     ahorm = ahorm + rc * a 'hormone production by a
#     IF i = js THEN 'averaging
#     CALL hormone(3, ahorm, ja, js)
#      rbb = rB / C '---- effective inhibitor decay rate
#      drb = 1! - 2! * db - rbb
#     END IF
#   NEXT i

#SUB hormone (ila, ahorm, ja, js)
#rx = flv(ila,2)
#	  axt(ila, 1) = axt(ila, 1) * (1! - rx) + ahorm / (js - ja + 1)'C() represents hormone
#	  FOR iic = ja TO js 
#	    axt(ila, iic) = axt(ila, 1)
#	  NEXT iic
#	  ahorm = 0
#  EXIT SUB
#END SUB