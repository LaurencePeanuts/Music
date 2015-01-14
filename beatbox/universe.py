import numpy as np

class Universe(object):
    """
    A simple model universe in a box. Coordinates are comoving cartesians!
    """
    def __init__(self):
        self.size = 100
        self.phi = np.zeros([3,self.size])
        return

    def __str__(self):
        return "a model universe, containing potential map phi"
