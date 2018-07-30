import numpy as np


class Oliva(object):
    def __init__(self,width=640, skip=12, act_diff=0.015, act_decay=0.1,
        act_prod=0.1, sat=0.25, in_diff=0.0, in_decay=0.014, in_mm=0.1,
        h_decay=0.1, hormone=0.5):
         
        self.width = width
        self.cells = np.zeros((2,2,self.width))
        
        self.skip = skip
        
        self.act_diff = act_diff
        self.act_decay = act_decay
        self.act_prod = act_prod
        self.sat = sat
         
        self.in_diff = in_diff
        self.in_decay = in_decay
        self.in_mm = in_mm
        self.h_decay = h_decay
        self.h_fac = 1-self.h_decay
        self.hormone = hormone
        
        self.tick = False
         
        self.cells[0,1,:] = 0.1
         
        self.fluct = self.act_decay * (0.96 +
            0.08 *np.random.random(self.width))
            
        seeds = np.random.choice(np.arange(self.width),30,replace=False)
        self.cells[0,0,seeds] = 1.0
        
        self.act_diff_const = 1.0 - self.act_decay -2*self.act_diff
        self.in_diff_const = 1.0 - self.in_decay -2*self.in_diff


    def step(self):
        if self.tick:
            old = 1
            new = 0
        else:
            old = 0
            new = 1
        
        l_bound = np.copy(self.cells[old,:,0])
        r_bound = np.copy(self.cells[old,:,-1])
                
        act_sq = np.square(self.cells[old,0,:])
        auto_cat = self.fluct * act_sq / (1 + self.sat * act_sq)
        
        left_cells = np.roll(self.cells[old,:,:],-1,axis=1)
        right_cells = np.roll(self.cells[old,:,:],1,axis=1)

        left_cells[:,0] = l_bound
        right_cells[:,-1] = r_bound
                
        self.cells[new,0,:] = self.cells[old,0,:] * self.act_diff_const + self.act_diff * (left_cells[0,:] + right_cells[0,:]) + auto_cat / (self.in_mm + self.cells[old,1,:])
            
        self.cells[new,1,:] = self.cells[old,1,:] * self.in_diff_const + self.in_diff * (left_cells[1,:] + right_cells[1,:]) + auto_cat
            
        hormone_prod = (self.cells[old,0,:] * self.h_decay).sum()
        
        self.hormone = self.hormone * self.h_fac + hormone_prod / self.width
            
        self.in_diff_const = 1.0 - 2 * self.in_diff - self.in_decay / self.hormone
        
        self.tick = not self.tick
        
        
    def __iter__(self):
        return self
    
    
    def __next__(self):
        self.step()
        if self.tick:
            out = np.copy(self.cells[0,:,:])
        else:
            out = np.copy(self.cells[1,:,:])
        for i in range(self.skip):
            self.step()
        return out
        
        