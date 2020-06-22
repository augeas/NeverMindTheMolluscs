#! /usr/bin/python

import json
import os
import re
import string
import sys

line1_pars = ['kt','kp','kx','ky','kd','ki','ke','kr','kn','kg']
line2_pars = ['k1','k2','k3','k4','dx','dy','dz','dw']

substance_pars = 'drbscag'


def parse_line(line,pars,par_type):
    return dict(zip(pars,map(par_type,re.split(',\s?',line))))


def parse_file(fname):
    with open(fname) as f:
        lines = f.readlines()
    ref = int(''.join(filter(lambda c: c in string.digits, fname)))
    pars = {'title':lines[0], 'ref':ref}
    pars.update(parse_line(lines[1],line1_pars,int))
    pars.update(parse_line(lines[2],line2_pars,float))

    for key in line2_pars[0:4]:
        pars[key] = int(pars[key])

    for line, sub in zip(range(3,3+pars['kn']),string.ascii_lowercase):
        keys = [''.join([char,sub]) for char in substance_pars]
        pars.update(parse_line(lines[line],keys,float))
        
    return pars

    
def read_sp_pars(path):
    all_files = sorted(os.listdir(path))
    prm_files = filter(lambda n:n.endswith('prm'),all_files)
    prm_paths = ('/'.join([path,f]) for f in prm_files)
    dump = json.dumps([parse_file(f) for f in prm_paths])
    with open('sp_pars.json','w') as dumpfile:
        dumpfile.write(dump)
        
    
if __name__ == "__main__":
    read_sp_pars(sys.argv[1])
