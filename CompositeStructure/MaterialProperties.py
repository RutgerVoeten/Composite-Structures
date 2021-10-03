import numpy as np

class MaterialProperties:
    def __init__(self):
        """ Input the material properties of a composite lamina """
        self.E1 = 180e3
        self.E2 = 70e3
        self.G12 = 45e3
        self.v12 = 0.3
        self.v21 = self.v12/self.E1 * self.E2
        self.t = 0.125

    def InputFile(self, file):
        print(file)
