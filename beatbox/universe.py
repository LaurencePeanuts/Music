import numpy

class Universe(object):

    def __init__(self):
        self.size = 100
        self.phi = np.zeros(3,100)
        return

    def __str__(self):
        return "a model universe, containing potential map phi"
