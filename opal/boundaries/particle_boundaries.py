__author__ = 'swebb'

import numpy as np

class periodic_boundaries:


    def __init__(self, pd):

        self.boundary = []


    def add_boundary(self, bc):

        self.boundary.append(bc)


    def apply_boundary(self, ptcls):

        for bc in self.boundary:
            bc.apply_boundary(ptcls)