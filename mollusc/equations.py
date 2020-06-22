

import numpy as np


class MolluscMixin(object):

    def eq_21(self, old_cells, new_cells, decay, a_sq):
        """- activator - inhibitor mechanism: B is inhibitor"""
        new_cells[0] = decay[0] + self.s * (a_sq / old_cells[1] + self.basic_prod[0])
        new_cells[1] = decay[1] + self.s * a_sq + self.basic_prod[1]


    def eq_61(self, old_cells, new_cells, decay, a_sq):
        """Hormone (c) changes lifetime of the inhibitor"""
        aq = self.s * a_sq / (1 + self.activator_saturation * a_sq) + self.basic_prod[0]
        new_cells[0] = decay[0] + aq / (self.saturation[1] + old_cells[1])
        new_cells[1] = decay[1] + aq + self.basic_prod[1]


    '''
    Here is a short "fair-use" extract from the original FreeBASIC code for
    equation 86, global state, GOSUBs and all. This is what we're up against...

    CASE 86 'three inhibitors (b,e,f) and oscillation - L. hieroglyphica
    FOR i = ja TO js
        C = axt(3, i): D = axt(4, i)
        rbb = rB * (1 + cB * D + cd * C) / (1 + cc * C)'---- effective deccay
        drb = 1! - 2! * db - rbb           '     rate of the inhibitor
    GOSUB olddecay:
        aq = s * (a * a / (1! + sA * a * a) + ba)
        axt(1, i) = olddecaydiffA + aq / (sb + sf * f + se * e + b)
        axt(2, i) = b * (1 - rbb - 2 * db) + db * (bl + axt(2, i + 1)) + aq + bb
        cq = rc * arandomxt(i) * (C * C + bc)         'c and d - > oscillation
        axt(3, i) = olddecaydiffC + cq / (sd + D)  '= A-I system
        axt(4, i) = olddecaydiffD + cq
        axt(5, i) = olddecaydiffE + re * a
        axt(6, i) = olddecaydiffF + rf * a
    NEXT i
    '''


    def eq_86(self, old_cells, new_cells, decay, a_sq):
        '''three inhibitors (b,e,f) and oscillation - L. hieroglyphica'''
        c = old_cells[2]
        d = old_cells[3]
        
        rbb = self.decay_rates[1] * (
            1 + self.coupling[1] * d + self.coupling[3] * c) / (1
            + self.coupling[2] * c) 

        drb = 1 - 2 * self.diff_coeffs[1] - rbb

        aq = self.s * (a_sq / (1 + self.saturation[0] * a_sq)
            + self.basic_prod[0])
        
        divisor = (self.saturation[1] + self.saturation[5] * old_cells[5]
            + self.saturation[4] * old_cells[4] + old_cells[1])
        
        valid = divisor > 0.0
        
        '''The rather ill-defined behaviour of FreeBASIC for division by zero
        copes with this quite quietly...'''
        new_cells[0][valid] = decay[0][valid] + aq[valid] / divisor[valid]
        
        new_cells[1] = (old_cells[1] * (1 - rbb - 2 * self.diff_coeffs[1])
            + self.diff_coeffs[1] * (np.roll(old_cells[1], -1)
            + np.roll(old_cells[1], 1)) + aq + self.basic_prod[1]) 

        cq = self.decay_rates[2] * self.fluctuate() * (np.square(c)
            + self.basic_prod[2])
        
        new_cells[2] = decay[2] + cq / (self.saturation[3] + d)
        new_cells[3] = decay[3] + cq
        new_cells[4] = decay[4] + self.decay_rates[4] * old_cells[0]
        new_cells[5] = decay[5] + self.decay_rates[5] * old_cells[0]


    def eq_711(self, old_cells, new_cells, decay, a_sq):
        '''activation (ab) and extinguishing  (e,f) reaction as 71 without
        enhancing reaction'''
        aq = self.s * old_cells[1] * (a_sq /
            (1 + self.activator_saturation * a_sq) + self.basic_prod[0])
        new_cells[0] = (decay[0] + aq - self.coupling[4] * old_cells[0] *
            old_cells[4])
        new_cells[1] = decay[1] - aq + self.basic_prod[1]
        e_sq = np.square(old_cells[4])
        eq = self.decay_rates[4] * old_cells[5] * (e_sq /
            (1.0 + self.saturation[4] * e_sq) + self.basic_prod[4])
        new_cells[4] = decay[4] + eq + self.saturation[5] * old_cells[0]
        new_cells[5] = decay[5] - eq + self.basic_prod[5] * old_cells[0] + self.coupling[5] * self.substrate
        self.stabilize(1,5)
        
