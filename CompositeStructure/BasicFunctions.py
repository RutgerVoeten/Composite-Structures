""" ABD Functions"""
import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi

class ABDMatrix:
    """
    A class used to evaluate the properties of the ABD matrix of a composite laminate
    ...

    """
    def __init__(self, stack , MaterialProps, fails = {}, degrade = {}, minreal = True):
        """

        Parameters
        ----------
        stack : dict
            a dictionary where each index has the corresponding angle in degrees
        E1 : float
            Stiffness of the ply on the first principle axis
        E2 : float
            Stiffness of the ply on the second principle axis
        v12 : float
            poisson ratio 12
        v21 : float
            poisson ratio 21 can be evaluated using v21 = E2 * v12/E1
        G12 : float
            Shear modulus of the ply in the 1-2 direction
        t : float
            Thickness of each ply
        fails : dict (optional)
            Dictionary of the plies which are completely failed
        degrade : dict (optional)
            Dictionary of the plies which are degraded in properties

        Attributes
        ----------
        ABD : arr
            ABD matrix

        Qmats : dict
            dictionary containing Q matrices for each ply

        stack : dict
            Returned stack excluding the failed plies

        z : arr
            z coordinates of original stacking sequence"""

        E1 = MaterialProps.E1
        E2 = MaterialProps.E2
        v12 = MaterialProps.v12
        v21 = MaterialProps.v21
        G12 = MaterialProps.G12
        t = MaterialProps.t

        self.CalcABD(stack, E1, E2, v12, v21, G12, t, fails = fails, degrade = degrade, minreal = minreal)


    def CalcABD(self, stack, E1, E2, v12, v21, G12, t, fails = {}, degrade = {}, minreal = True):

        z = np.arange(-len(stack)/2 * t, len(stack)/2*t+t/2,t)

        Axx = 0
        Axy = 0
        Ayy = 0
        Axs = 0
        Ays = 0
        Ass = 0
        Axs = 0


        Dxx = 0
        Dxy = 0
        Dxs = 0
        Dyy = 0
        Dys = 0
        Dxs = 0
        Dys = 0
        Dxx = 0
        Dss = 0

        Bxx = 0
        Bxy = 0
        Bxs = 0
        Byy = 0
        Bys = 0
        Bss = 0

        def mf(phi):
            return np.cos(phi)

        def nf(phi):
            return np.sin(phi)

        self.Qmats = {}

        for m in fails:
            del stack[m]


        for k in stack:

            if k in degrade:
                E2 = E2*degrade[k]
                G12 = G12*degrade[k]

            #Defining constants
            #For plane stress conditions
            Q = 1 - v12*v21
            Q11 = E1 * Q**-1
            Q22 = E2*Q**-1
            Q12 = v12*E2*Q**-1
            Q66 = G12
            Q21 = Q12

            phi = stack[k]*np.pi/180

            m = mf(phi)
            n = nf(phi)

            #Calculate QMatrix properties

            Qxx = Q11*m**4 + 2*(Q12+ 2*Q66)*m**2*n**2 + Q22*n**4
            Qxy = (Q11 + Q22 - 4*Q66)*m**2*n**2 + Q12*(m**4+n**4)
            Qyy = Q11*n**4 + 2*(Q12 + 2*Q66)*m**2*n**2 + Q22*m**4
            Qxs = (Q11 - Q12 - 2*Q66)*n*m**3 + (Q12-Q22+2*Q66)*n**3*m
            Qys = (Q11-Q12-2*Q66)*m*n**3+(Q12-Q22+2*Q66)*m**3*n
            Qss = (Q11 + Q22 -2*Q12 - 2*Q66)*n**2*m**2 + Q66*(n**4 + m**4)

            Qmat = np.array([[Qxx,Qxy,Qxs],[Qxy,Qyy,Qys],[Qxs,Qys,Qss]])


            self.Qmats[k] = Qmat

            d1 = (z[k] - z[k-1])
            d2 = (z[k]**2 - z[k-1]**2)
            d3 = (z[k]**3 - z[k-1]**3)

            Axx = Axx + Qxx * d1
            Axy = Axy + Qxy * d1
            Ayy = Ayy + Qyy * d1
            Axs = Axs + Qxs * d1
            Ays = Ays + Qys * d1
            Ass = Ass + Qss * d1

            Dxx = Dxx + 1/3*Qxx * d3
            Dxy = Dxy + 1/3*Qxy * d3
            Dxs = Dxs + 1/3*Qxs * d3
            Dyy = Dyy + 1/3*Qyy * d3
            Dys = Dys + 1/3*Qys * d3
            Dss = Dss + 1/3*Qss * d3

            Bxx = Bxx + 1/2*Qxx * d2
            Bxy = Bxy + 1/2*Qxy * d2
            Bxs = Bxs + 1/2*Qxs * d2
            Byy = Byy + 1/2*Qyy * d2
            Bys = Bys + 1/2*Qys * d2
            Bss = Bss + 1/2*Qss * d2


        self.ABD = np.array([[Axx,Axy,Axs,Bxx,Bxy,Bxs], [Axy,Ayy,Ays,Bxy,Byy,Bys],[Axs,Ays,Ass, Bxs, Bys, Bss],
                        [Bxx, Bxy, Bxs, Dxx, Dxy, Dxs] ,[Bxy, Byy, Bys, Dxy, Dyy, Dys],[Bxs, Bys, Bss, Dxs, Dys, Dss]])

        #Remove rounding errors:
        if minreal:
            self.ABD[np.where(np.abs(self.ABD) < 1e-5)] = 0

        # Stack excluding failed plies
        self.stack = stack

        # Z coordinates including failed plies
        self.z = z

def Createdictfromlist(stack):
    stack = dict(zip(np.arange(1,len(stack)+1),stack))
    return stack

class StressCalculations:
    """
    Calculate the Stresses of a composite laminate and Q matrices for each lamina
    """

    def __init__(self, StiffnessMatrices, MaterialProperties, force_vector):
        self.z_interfaces = StiffnessMatrices.z
        self.stack = StiffnessMatrices.stack
        self.force_vector = force_vector
        self.ABD = StiffnessMatrices.ABD
        self.Q_mats = StiffnessMatrices.Qmats
        self.t = MaterialProperties.t

        #
        self.StressCalc(force_vector)

    def StressCalc(self, force_vector):
        """
        Calculate the Stresses of a composite laminate and Q matrices for each lamina

        **Rutger Voeten 4973682 25-03-2021**

        Parameters
        ----------
        Stackor: Array with all the original plies without failure
        Stack: Sequence of angles starting from -Z to + Z compiled in a dictionary starting from 1
        force_vector: Vector containing the forces -> np.array([[Nx],[Ny],[Ns],[Mx],[My],[Ms]])
        ABD_matrix: The ABD matrix of the total laminate
        t: Thickness of each lamina, constant thickness for each lamina is assumed

        Returns
        -------
        Z_array for plotting, sigma_dict, Stresses per layer in the form: array([[Sigmax],[Sigmay],[Sigmas]]), global strains vector
        """

        z_laminate = self.z_interfaces

        # Calculate global strains
        glob_strains = np.linalg.inv(self.ABD) @ force_vector

        # z increment for stress calculations
        dz = 0.0001  # mm

        # total arrays
        Z_array = np.array([])

        sigmax = np.array([])
        sigmay = np.array([])
        sigmas = np.array([])

        epsx = np.array([])
        epsy = np.array([])
        epss = np.array([])

        Layer_stress = {}

        for i in self.stack:

            # Ply coordinates in Z
            z_layer = np.arange(z_laminate[i - 1], z_laminate[i] + dz, dz)

            # Combine in total Z Coordinates
            Z_array = np.append(Z_array, z_layer)

            # Total strains on ply level
            strain = glob_strains[:3] + z_layer * glob_strains[3:]

            # Stresses on ply level
            tot_stress = self.Q_mats[i] @ strain

            # Save stresses in dictionary in order to recall this
            Layer_stress[i] = tot_stress

            # Sigmay, Sigmax, and Sigmas
            sigmax = np.append(sigmax, tot_stress[0])
            sigmay = np.append(sigmay, tot_stress[1])
            sigmas = np.append(sigmas, tot_stress[2])

            epsx = np.append(epsx, strain[0])
            epsy = np.append(epsy, strain[1])
            epss = np.append(epss, strain[2])

        sigma = {"Sigma_x": sigmax, "Sigma_y": sigmay, "Sigma_xy": sigmas}
        eps = {"Eps_1": epsx, "Eps_2": epsy, "Eps_12": epss}

        self.Z_array = Z_array
        self.sigma = sigma
        self.eps = eps
        self.Layer_stress = Layer_stress

    def plot_stresses(self):

        plt.figure(figsize=(15,10))

        sigma_names = ['$\sigma_{x}$', '$\sigma_{y}$', '$\\tau_{xy}$']

        n = 0
        for i in self.sigma:
            plt.plot(self.sigma[i], self.Z_array, label = '{}'.format(sigma_names[n]))
            n = n+1

        minval = np.min([np.min(self.sigma['Sigma_x']), np.min(self.sigma['Sigma_y']), np.min(self.sigma['Sigma_xy'])])
        maxval = np.max([np.max(self.sigma['Sigma_x']), np.max(self.sigma['Sigma_y']), np.max(self.sigma['Sigma_xy'])])
        for m in self.z_interfaces:
            plt.hlines(m, minval, maxval, linestyles='dashed', color = 'k')

        plt.xlabel('Stresses')
        plt.ylabel('Z location')
        plt.legend()
        plt.grid()
        plt.show()




    def CalcFPF(self, dNx, stackor, stackdict, props, ABD, Qmats):
        r''' Calculate first ply failure due to axial stress '''
        Nx = 0
        failscount = 0

        while failscount == 0:
            force_vector = np.array([[Nx],[0],[0],[0],[0],[0]])
            #Calculated stresses in principle direction
            Z_array, sigma, eps, Layer_stress, strains = self.StressCalc(stackor, stackdict, force_vector, ABD, Qmats, props['tply'])
            #Calculate if failure happened
            failscount, ply_stress = self.failure(stackdict, Layer_stress, Qmats, props['Xt'], props['Yt'], props['Sxz'], props['Xc'], props['Yc'], props['vyx'], props['Gxy'], props['tply'])
            Nx += dNx

        sigmaXFPF = Nx/(len(stackdict)* props['tply'])

        return sigmaXFPF


    def Check_hashin(self, ply_stress, Xt,Yt,S12, Xc, Yc):

        #Tension fibres
        if len(ply_stress[0,ply_stress[0,:] >=0]) > 0:
            sigma1 = ply_stress[0,ply_stress[0,:] >=0]
            tau12 = ply_stress[2,ply_stress[0,:] >= 0]
            C1 = (sigma1/Xt)**2 + tau12**2/S12**2

        else:
            C1 = 0

        #Compressive fibre
        if len(ply_stress[0,ply_stress[0,:] < 0]) > 0:
            sigma1 = ply_stress[0, ply_stress[0,:] < 0]
            C2 = (sigma1/Xc)**2

        else:
            C2 = 0

        #Tension matrix
        if len(ply_stress[1,ply_stress[1,:] >= 0]) > 0:
            sigma2 = ply_stress[1,ply_stress[1,:] >= 0]
            tau12 = ply_stress[2,ply_stress[1,:] >= 0]
            C31 = (sigma2)**2/Yt**2 + tau12**2/S12**2

        else:
            C31 = 0

        #Compression Matrix
        if len(ply_stress[1,ply_stress[1,:] < 0]) > 0:
            sigma2 = ply_stress[1,ply_stress[1,:] < 0]
            tau12 = ply_stress[2,ply_stress[1,:] < 0]
            C32 = 1/Yc * ((Yc/(2*S12)**2)-1)*sigma2 + 1/(4*S12**2)*sigma2**2 + tau12**2/S12**2

        else:
            C32 = 0


        C1 = np.max(np.abs(C1))
        C2 = np.max(np.abs(C2))
        C31 = np.max(np.abs(C31))
        C32 = np.max(np.abs(C32))

        C3 = max([C31,C32])

        if C1 >= 1 or C2 >= 1:
            print('Fibre failure:')

        if C31 >= 1:
            print('Matrix failure tension:')
        elif C32 >= 1:
            print('Matrix failure compression:')

        return  C1, C2, C3

    def check_TsaiHill(self, ply_stress, Xt,Yt,S12, Xc, Yc):

        #Tension sigmax
        if len(ply_stress[0,ply_stress[0,:] >=0]) > 0:
            sigma1 = ply_stress[0,ply_stress[0,:] >=0]
            sigma2 = ply_stress[1,ply_stress[0,:] >=0]
            tau12 = ply_stress[2,ply_stress[0,:] >= 0]


            X = Xt
            Y = Yt

            signs = np.sign(sigma2)
            sigma1A = sigma1[signs == 1]
            sigma2A = sigma2[signs ==1]
            tau12A = tau12[signs==1]

            C1A = sigma1A**2/X**2 - sigma1A*sigma2A/X**2 + sigma2A**2/Y + tau12A**2/S12**2

            Y = Yc

            sigma1B = sigma1[signs == -1]
            sigma2B = sigma2[signs == -1]
            tau12B = tau12[signs== -1]

            C1B = sigma1B**2/X**2 - sigma1B*sigma2B/X**2 + sigma2B**2/Y + tau12B**2/S12**2

            if len(C1A) > 0:
                C1A = np.max(np.abs(C1A))
                print('C1A =', C1A)
            else:
                C1A = 0
            if len(C1B) > 0:
                C1B = np.max(np.abs(C1B))
                print('C1B =', C1B)
            else:
                C1B = 0
        else:
            C1A, C1B = 0,0


        #Compression sigmax
        if len(ply_stress[0,ply_stress[0,:] < 0]) > 0:
            sigma1 = ply_stress[0,ply_stress[0,:] < 0]
            sigma2 = ply_stress[1,ply_stress[0,:] < 0]
            tau12 = ply_stress[2,ply_stress[0,:] < 0]


            X = Xc
            Y = Yt

            signs = np.sign(sigma2)
            sigma1A = sigma1[signs == 1]
            sigma2A = sigma2[signs ==1]
            tau12A = tau12[signs==1]

            C2A = sigma1A**2/X**2 - sigma1A*sigma2A/X**2 + sigma2A**2/Y + tau12A**2/S12**2

            Y = Yc

            sigma1B = sigma1[signs == -1]
            sigma2B = sigma2[signs == -1]
            tau12B = tau12[signs== -1]

            C2B = sigma1B**2/X**2 - sigma1B*sigma2B/X**2 + sigma2B**2/Y + tau12B**2/S12**2

            if len(C2A) > 0:
                C2A = np.max(np.abs(C2A))
                print('C2A =', C2A)
            else:
                C2A = 0
            if len(C2B) > 0:
                C2B = np.max(np.abs(C2B))
                print('C2B =', C2B)
            else:
                C2B = 0

        else:
            C2A, C2B = 0,0

        Crits = [np.abs(C1A), np.abs(C1B), np.abs(C2A), np.abs(C2B)]


        return np.max(Crits)

    def check_maxStress(self, ply_stress, Xt,Yt,S12, Xc, Yc):

        r"""
        Calculate the Stresses of a composite laminate and Q matrices for each lamina

        **Rutger Voeten 4973682 25-03-2021**

        Parameters
        ----------
        ply_stress: The stress in the corresponding ply
        Xt
        Yt
        S12
        Xc

        Returns
        -------
        Failure criteria --> C1-tensionx, C1-compressionx, C2-tensiony, C2-compressiony, C3-shear
        """

        #Tension sigmax
        if len(ply_stress[0,ply_stress[0,:] >=0]) > 0:

            sigma1 = np.max(np.abs(ply_stress[0,ply_stress[0,:] >=0]))

            C1A = sigma1/Xt

        else:
            C1A = 0

        #Compression sigmax
        if len(ply_stress[0,ply_stress[0,:] < 0]) > 0:
            sigma1 = np.max(np.abs(ply_stress[0,ply_stress[0,:] < 0]))

            C1B = sigma1/(Xc)

        else:
            C1B = 0

        #Tension sigmay
        if len(ply_stress[1,ply_stress[1,:] >= 0]) > 0:
            sigma2 = np.max(np.abs(ply_stress[1,ply_stress[1,:] >= 0]))

            C2A = sigma2/Yt
        else:
            C2A = 0

        if len(ply_stress[1,ply_stress[1,:] < 0]) > 0:
            sigma2 = np.max(np.abs(ply_stress[1,ply_stress[1,:] < 0]))

            C2B = sigma2/Yc
        else:
            C2B = 0

        tau12 = np.max(np.abs(ply_stress[2,:]))

        C3 = np.abs(tau12/S12)

        return C1A, C1B, C2A, C2B, C3

    def check_HashinRottem(self, ply_stress, Xt,Yt,S12, Xc, Yc):
        #Tension fibres
        if len(ply_stress[0,ply_stress[0,:] >=0]) > 0:
            sigma1 = ply_stress[0,ply_stress[0,:] >=0]
            tau12 = ply_stress[2,ply_stress[0,:] >= 0]
            C1 = (sigma1/Xt)**2

        else:
            C1 = 0

        #Compressive fibre
        if len(ply_stress[0,ply_stress[0,:] < 0]) > 0:
            sigma1 = ply_stress[0, ply_stress[0,:] < 0]
            C2 = (sigma1/Xc)**2

        else:
            C2 = 0

        #Tension matrix
        if len(ply_stress[1,ply_stress[1,:] >= 0]) > 0:
            sigma2 = ply_stress[1,ply_stress[1,:] >= 0]
            tau12 = ply_stress[2,ply_stress[1,:] >= 0]
            C31 = (sigma2)**2/Yt**2 + tau12**2/S12**2

        else:
            C31 = 0

        #Compression Matrix
        if len(ply_stress[1,ply_stress[1,:] < 0]) > 0:
            sigma2 = ply_stress[1,ply_stress[1,:] < 0]
            tau12 = ply_stress[2,ply_stress[1,:] < 0]
            C32 = sigma2**2/Yc**2 + tau12**2/S12**2

        else:
            C32 = 0


        C1 = np.max(np.abs(C1))
        C2 = np.max(np.abs(C2))
        C31 = np.max(np.abs(C31))
        C32 = np.max(np.abs(C32))

        C3 = max([C31,C32])

        # if C1 >= 1 or C2 >= 1:
        #     print('Fibre failure:')

        # if C31 >= 1:
        #     print('Matrix failure tension:', C31)

        # elif C32 >= 1:
        #     print('Matrix failure compression:', C32)

        return  C1, C2, C3

    def failure(self, stack ,layer_stress,Qmats,Xt,Yt,S12,Xc,Yc, MaxStress = True, \
                Hashin = False, HashinRottem=False, Tsai_Hill = False):

        r"""
        Calculate the Stresses of a composite laminate and Q matrices for each lamina

        **Rutger Voeten 4973682 25-03-2021**

        Parameters
        ----------
        Stack: Stacking sequence dict
        Layer stresses: Stresses in each layer of a laminate - dict
        Qmats: Qmatrices for each ply on principle axis
        Xt
        Yt
        S12
        Xc
        Yc
        Failure criteria, standard = MaxStress failure criterion

        Returns
        -------
        Failure criteria --> C1-tensionx, C1-compressionx, C2-tensiony, C2-compressiony, C3-shear
        """
        failscount = 0
        for check in Qmats:
            theta = stack[check]*np.pi/180 #Calculate corresponding angle
            ply_stress = layer_stress[check] #Stress in each ply R[x]

            m  = cos(theta)
            n = sin(theta)

            # rotate stresses to ply orientation for each ply in ply direction
            rot = np.array([[m**2, n**2, 2*m*n],
                              [n**2, m**2, -2*m*n],
                              [-m*n, m*n, m**2-n**2]])

            #Take ply stresses in the orientation of the fibres
            ply_stress = rot @ ply_stress

            #Implement failure criteria
            if Hashin:
                Crit1, Crit2, Crit3 = self.Check_hashin(ply_stress, Xt, Yt, S12, Xc, Yc)

                if abs(Crit1) >= 1 or abs(Crit2) >= 1 or abs(Crit3) >= 1:
                    failscount = failscount + 1

            if Tsai_Hill:
                Crit, sigma1A, sigma1B, stressaa = self.check_TsaiHill(ply_stress, Xt, Yt, S12, Xc, Yc)

                if abs(Crit) >= 1:
                    failscount = failscount + 1


            if HashinRottem:
                Crit1, Crit2, Crit3 = self.check_HashinRottem(ply_stress, Xt, Yt, S12, Xc, Yc)

                if abs(Crit1) >= 1 or abs(Crit2) >= 1 or abs(Crit3) >= 1:
                    failscount = failscount + 1

            if MaxStress:
                C1, C2, C3, C4, C5 = self.check_maxStress(ply_stress, Xt, Yt, S12, Xc, Yc)

                if C1 >= 1 or C2 >= 1 or C3 >= 1 or C4 >=1 or C5 >= 1:
                    failscount = failscount+1

            if failscount >= 1:
                break

        return failscount, ply_stress
