
import CompositeStructure.MaterialProperties as MP
import CompositeStructure.BasicFunctions

""" Example for calculating the stiffness matrices for each lamina in x-y direction and
the laminate's ABD matrix """

stack = [45,90,0,0,90,45]

stack = CompositeStructure.BasicFunctions.Createdictfromlist(stack)

MaterialProperties = MP.Props()
StiffnessMatrices = CompositeStructure.BasicFunctions.ABDMatrix(stack,MaterialProperties)

