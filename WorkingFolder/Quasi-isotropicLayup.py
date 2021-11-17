import CompositeStructure.MaterialProperties as MP
import CompositeStructure.BasicFunctions

MaterialProperties = MP.Props()

stack = [45,-45,90,0,0,90,-45,45]

stack = CompositeStructure.BasicFunctions.Createdictfromlist(stack)

ABD = CompositeStructure.BasicFunctions.ABDMatrix(stack, MaterialProperties)

A11 = ABD.ABD[0,0]
A12 = ABD.ABD[0,1]
A22 = ABD.ABD[1,1]

print('A12/A11 = ' + str(A12/A11))
print('second = {}'.format(A11/(A11 + 3 * A22)))



