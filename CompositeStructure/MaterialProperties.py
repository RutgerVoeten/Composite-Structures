import numpy as np

class Props:
    def __init__(self):
        """ Input the material properties of a composite lamina """
        self.E1 = 210e3
        self.E2 = 210e3
        self.G12 = 210e3/(2*(1+0.3))
        self.v12 = 0.3
        self.v21 = self.v12/self.E1 * self.E2
        self.t = 0.125

    def InputFile(self, file):
        #TODO Add the possibility to use an input file


        return file
