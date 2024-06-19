# from curses.ascii import isdigit #uncomment if using linux
import os

import matplotlib.pyplot as plt
import numpy as np


"""
Disclaimer: Make sure to mave file to appropriate location in order to read it.
This file also assume that the WADAM file is on .txt format and that the format
of the text file is equal to standard WADAM files

Exciting Force have to be added manually in HydroDReader

"""


def str2num(num_str, Beta = False):
    """
    Takes in string on specific format and transforms it to float.


    @param num_str    Format 0.0000E-00 (string)

    @return float     The string num_str is mapped to float. (float)
    """
    decimal = 0
    num = 0.0
    sign = 1

    if Beta:
        return float(num_str)

    # Save sign of number
    if num_str[0] == "-":
        sign = -1.0
        num_str = num_str[1:]

    # Add together numbers left of E
    for i in range(len(num_str)-4):
        if i == 1:
            continue
        num += float(num_str[i])*10**decimal
        decimal -= 1

    # Take the sign of the input into consideration
    if num_str[-3] == "+":
        return sign*num*10**int(num_str[-2:])
    return sign*num*10**(-int(num_str[-2:]))


def addDimentions_A(A, rho, V, L):
    """
    Adding dimensions to added mass matrix or mass matrix.
    They are non-dimensionalized by the same factors.

    @param A       Added mass matrix dimensionless.   (6x6 matrix)
    @param rho     Density of fluid.    (float)
    @param V       Volume of structure. (float)
    @param L       Length of structure. (float)

    @return A      Added mass matrix.   (6x6 matrix)
    """
    A[0:3, 0:3] = A[0:3, 0:3]*rho*V
    A[3:6, 3:6] = A[3:6, 3:6]*rho*V*L*L
    A[0:3, 3:6] = A[0:3, 3:6]*rho*V*L
    A[3:6, 0:3] = A[3:6, 0:3]*rho*V*L
    return A


def addDimentions_D(D, rho, V, L, G):
    """
    Adding dimensions to damping matrix.

    @param D       Damping matrix dimensionless.   (6x6 matrix)
    @param rho     Density of fluid.    (float)
    @param V       Volume of structure. (float)
    @param L       Length of structure. (float)
    @param G       Gravity constant.    (float)

    @return D      Damping matrix.      (6x6 matrix)
    """
    D[0:3, 0:3] = D[0:3, 0:3]*rho*V*np.sqrt(G/L)

    D[3:6, 3:6] = D[3:6, 3:6]*rho*V*L*np.sqrt(G*L)

    D[0:3, 3:6] = D[0:3, 3:6]*rho*V*np.sqrt(G*L)
    D[3:6, 0:3] = D[3:6, 0:3]*rho*V*np.sqrt(G*L)
    return D

def addDimentions_C(C, rho, V, L, G):
    """
    Adding dimentions to restoring matrix.

    @param C       Restoring matrix dimensionless.   (6x6 matrix)
    @param rho     Density of fluid.    (float)
    @param V       Volume of structure. (float)
    @param L       Length of structure. (float)
    @param G       Gravity constant.    (float)

    @return C      Restoring matrix.      (6x6 matrix)
    """
    C[0:3, 0:3] = C[0:3, 0:3]*rho*V

    C[3:6, 3:6] = C[3:6, 3:6]*rho*V*L*G

    C[0:3, 3:6] = C[0:3, 3:6]*rho*V*G
    C[3:6, 0:3] = C[3:6, 0:3]*rho*V*G
    return C


def addDimentions_Fexc(Fexc, rho, V,G,WA,L):
    """
    Adding dimentions to Exciting forces

    @param Fexc    Exciting forces dimensionless.   (6x6 matrix)
    @param rho     Density of fluid.    (float)
    @param V       Volume of structure. (float)
    @param L       Length of structure. (float)
    @param G       Gravity constant.    (float)
    @param WA      Wave amplitude of the incoming waves (float)

    @return Fexc      Exciting force matrix.      (9x(6x4) matrix)
    """

    for i in range(len(Fexc)):
        Fexc[i][0:3, 0:3] = Fexc[i][0:3,0:3]*rho*V*G*WA/L
        Fexc[i][3:6,0:3] = Fexc[i][3:6,0:3]*rho*V*G*WA

    return Fexc


def addDimentions_motions(Motion, WA, L):
    """
    Adding dimentions to Exciting forces

    @param Motion    Exciting forces dimensionless.   (6x6 matrix)
    @param L       Length of structure. (float)
    @param WA      Wave amplitude of the incoming waves (float)

    @return Fexc      Exciting force matrix.      (9x(6x4) matrix)
    """

    for i in range(len(Motion)):
        Motion[i][0:3, 0:3] = Motion[i][0:3, 0:3] * WA
        Motion[i][3:6, 0:3] = Motion[i][3:6, 0:3] * WA / L

    return Motion


def readWadam(filename):
    """
    Collects constants (rho, g, V, L, WA),mass inertia matrix (M) ,added mass matrix (A), Hydostatic restoring matrix (C) and Damping matrix (D) from wadam file.


    Parameter:
    @param filename    Name of Wadam file (string)

    @return constants          dictionary that maps key[constant(str)] to value[constant(float)]
    @return periodToAddedMass    dictionary that maps key[period(float)] to value[added mass(6x6 matrix)]
    @return periodToDamping    dictionary that maps key[period(float)] to value[damping(6x6 matrix)]
    """


    #NB: Endret til os.path.join() for å bedre ta hensyn til fremtidige endringer. Sjekker også om main kjører, ellers
    # så antas det at multimodule.py kjører


    if (os.getcwd() == 'C:\\Python Prosjekter\\multimoduleWithBranch'):
        f = open(os.path.abspath(os.path.join(os.getcwd(), filename)), "r")
    else:
        f = open(os.path.abspath(os.path.join(os.getcwd(), os.pardir, filename)), "r")



    # Set decimals when printing matrices
    np.set_printoptions(precision=9)

    # Temporary variables
    T_i = 0
    HeadingFexc_i = []
    HeadingMotion_i = []

    M_i = np.zeros((6, 6))
    A_i = np.zeros((6, 6))
    C_i = np.zeros((6, 6))
    D_i = np.zeros((6, 6))
    Fexc_i = np.zeros((6,4))
    Motion_i = np.zeros((6,4))


    num_str = ""
    row_str = ""

    # Dictionarys
    constants = {"RO": 0.0,
                 "G": 0.0,
                 "VOL": 0.0,
                 "L": 0.0,
                 "WA": 1.0}

    periodToMass_noDim = 0
    periodToAddedMass = {}
    periodToAddedMass_noDim = {}
    HeadingAngles = {}
    periodToRestoring_noDim = 0
    periodToDamping = {}
    periodToDamping_noDIm = {}
    periodToHeadingToFexc = {}
    periodToHeadingToMotion = {}
    periodToHeadingToFexcDim = {}
    periodToHeadingToMotionDim = {}
    # Iterate through file
    for line in f:
        # Extract wave period
        if "::  W A V E  P E R I O D" in line:
            # find number str
            for c_i in range(len(line)):
                if line[c_i] == "=":
                    num_str = str(line[c_i + 1:]).replace(":", "").strip()
            # convert to number
            T_i = str2num(num_str)
            HeadingFexc_i = []
            HeadingMotion_i = []


        if ":  H E A D I N G   A N G L E" in line:
            #find number str
            for c_i in range(len(line)):
                if line[c_i] == "=":
                    num_str = str(line[c_i +1:]).replace(":", "").strip()

            #convert to number
            Beta_i = str2num(num_str,True)

        # extract mass matrix
        if "    MASS INERTIA COEFFICIENT MATRIX        " in line:
            f.readline()
            f.readline()
            f.readline()

            # fill in matrix row by row
            for i in range(6):
                line = f.readline()

                row_str = line.split()[1:]
                for j in range(len(row_str)):
                    M_i[i, j] = str2num(row_str[j])
                    try:
                        M_i[i, j] = str2num(row_str[j])
                    except:
                        print("ERROR: ", row_str[j],
                              " is not valid for C[", i, ",", j, "]!!")

            # Add mass matrix with wave period as key
            periodToMass_noDim = np.copy(M_i)

        # extract added mass matrix
        if "    ADDED MASS MATRIX        " in line:
            f.readline()
            f.readline()
            f.readline()

            # fill in matrix row by row
            for i in range(6):
                line = f.readline()

                row_str = line.split()[1:]
                for j in range(len(row_str)):
                    A_i[i, j] = str2num(row_str[j])
                    try:
                        A_i[i, j] = str2num(row_str[j])
                    except:
                        print("ERROR: ", row_str[j],
                              " is not valid for A[", i, ",", j, "]!!")

            # Add added mass matrix with wave period as key
            periodToAddedMass_noDim[T_i] = np.copy(A_i)

        if "    HYDROSTATIC RESTORING COEFFICIENT MATRIX        " in line:
            f.readline()
            f.readline()
            f.readline()

            # fill in matrix row by row
            for i in range(6):
                line = f.readline()

                row_str = line.split()[1:]
                for j in range(len(row_str)):
                    C_i[i, j] = str2num(row_str[j])
                    try:
                        C_i[i, j] = str2num(row_str[j])
                    except:
                        print("ERROR: ", row_str[j],
                              " is not valid for C[", i, ",", j, "]!!")

            # Add restoring matrix with wave period as key
            periodToRestoring_noDim = np.copy(C_i)

        # extract damping mass matrix
        if "    TOTAL DAMPING MATRIX         " in line:
            f.readline()
            f.readline()
            f.readline()

            # fill in matrix row by row
            for i in range(6):
                line = f.readline()

                row_str = line.split()[1:]
                for j in range(len(row_str)):
                    D_i[i, j] = str2num(row_str[j])
                    try:
                        D_i[i, j] = str2num(row_str[j])
                    except:
                        print("ERROR: ", row_str[j],
                              " is not valid for D[", i, ",", j, "].")

            # Add damping matrix with wave period as key
            periodToDamping_noDIm[T_i] = np.copy(D_i)

        # Extract Exciting forces and moments from integration pressure
        if "    EXCITING FORCES AND MOMENTS FROM INTEGRATION OF PRESSURE" in line:
            f.readline()
            f.readline()
            f.readline()
            f.readline()

            for i in range(6):
                line = f.readline()

                row_str = line.split()[1:]
                for j in range(len(row_str)):
                    Fexc_i[i,j] = str2num(row_str[j],True)
                    try:
                        Fexc_i[i,j] = str2num(row_str[j],True)
                    except:
                        print("Error: ", row_str[j],
                              " is not valid for Fexc[",i,",",j,"].")
                f.readline()
            HeadingFexc_i.append(np.copy(Fexc_i))
            periodToHeadingToFexc[T_i] = HeadingFexc_i

        if "    MOTION " in line:
            f.readline()
            f.readline()
            f.readline()
            f.readline()

            for i in range(6):
                line = f.readline()

                row_str = line.split()[1:]
                for j in range(len(row_str)):
                    Motion_i[i,j] = str2num((row_str[j]),True)
                    try:
                        Motion_i[i,j] = str2num((row_str[j]),True)
                    except:
                        print("Error: ", row_str[j],
                              " is not valid for Motion[",i,",",j,"].")
                f.readline()
            HeadingMotion_i.append(np.copy(Motion_i))
            periodToHeadingToMotion[T_i] = HeadingMotion_i

        # Extract constants(rho, g, wa, L, V) to constants map
        if "    NON-DIMENSIONALIZING FACTORS:" in line:
            # REMOVE 10 LINES
            c = 10
            for line in f:
                c -= 1
                if c == 0:
                    break

            # extract 4 parameters
            c = 4
            for line in f:
                l = str(line).split()
                constants[l[0]] = str2num(l[2])

                c -= 1
                if c == 0:
                    break

        if "IN DEGREES   IN RADIANS" in line:
            f.readline()
            for i in range(9):
                line = f.readline().strip().split("  ")
                HeadingAngles[float(line[1])] = int(line[0])-1  #Index in wadam starts on 1


    # Adding dimentions to added mass damping matrices
    periodToMass = addDimentions_A(periodToMass_noDim,constants["RO"],constants["VOL"],constants["L"])
    periodToRestoring = addDimentions_C(periodToRestoring_noDim, constants["RO"], constants["VOL"],
                                           constants["L"], constants["G"])

    for T in periodToAddedMass_noDim:
        periodToAddedMass[T] = addDimentions_A(
            periodToAddedMass_noDim[T], constants["RO"], constants["VOL"], constants["L"])
        periodToDamping[T] = addDimentions_D(
            periodToDamping_noDIm[T], constants["RO"], constants["VOL"], constants["L"], constants["G"])
        periodToHeadingToFexcDim[T] = addDimentions_Fexc(periodToHeadingToFexc[T],constants["RO"],constants["VOL"],constants["G"],constants["WA"],constants["L"])
        periodToHeadingToMotionDim[T] = addDimentions_motions(periodToHeadingToMotion[T],constants["WA"],constants["L"])
    """
    # Printing result
    print("___________________________________")
    print("Mass: ")
    print(periodToMass_noDim)
    # print(periodToAddedMass[T])

    print("___________________________________")
    print("Restoring: ")
    print(periodToRestoring_noDim)
    # print(periodToDamping)
    for T in periodToAddedMass:


        print("___________________________________")
        print("Period:", T)
        print("Added Mass: ")
        print(periodToAddedMass_noDim[T])
        # print(periodToAddedMass[T])

        print("___________________________________")
        print("Damped: ")
        print(periodToDamping_noDIm[T])
        # print(periodToDamping)



    print(constants)"""

    return constants, periodToDamping, periodToAddedMass, periodToMass, periodToRestoring,periodToHeadingToMotionDim,periodToHeadingToFexcDim, HeadingAngles

"""
print("constants :", constants)
print("periodToDamping :", type(periodToDamping))
print("periodToAddedMass :",periodToAddedMass)
print("periodToMass: ", type(periodToMass))
print("periodToRestoring: ", type(periodToRestoring))
print("periodToHeadingToMotionDim: ",type(periodToHeadingToMotionDim))
print("periodToHeadingToFexcDim: ", periodToHeadingToFexcDim[2.0])
print("HeadingAngles: ", type(HeadingAngles))
"""


"""
constants, periodToDamping, periodToAddedMass,periodToMass, periodToRestoring,motions, Fexc = readWadam("WADAMFILES/WADAM1LOTOFRUNS.txt")

A11 = []
A22 = []
A33 = []
Period = []
for T,A in periodToAddedMass.items():
    A11.append(A[0][0])
    A22.append(A[1][1])
    A33.append(A[2][2])
    Period.append(T)

plt.figure("Period")

plt.plot(Period,np.array(A22)/(1025*np.pi*(1.6/2)**2), label = "A22")
plt.plot(Period,np.array(A33)/(1025*np.pi*(1.6/2)**2),label = "A33")
plt.plot(Period,np.array(A11)/(1025*np.pi*(1.6/2)**2),linestyle = ":", label = "A11")
plt.hlines(4*1025*np.pi *(1.6/2)**3/(1025*np.pi*(1.6/2)**2),xmin=0,xmax=20,colors="orange",label = "A33 -theo")
plt.hlines((1025 * np.pi*(1.6/2)**2 *4 * 1.13/(1025*np.pi*(1.6/2)**2)),0,20,colors="blue",label= "A22 - theo")
plt.xlabel("Period T [s]")
plt.legend()


Period = 2*np.pi/np.array(Period)
plt.figure("Frequence")
plt.plot(Period,np.array(A22)/(1025*np.pi*(1.6/2)**2), label = "A22")
plt.plot(Period,np.array(A33)/(1025*np.pi*(1.6/2)**2),label = "A33")
plt.plot(Period,np.array(A11)/(1025*np.pi*(1.6/2)**2),linestyle = ":" ,label = "A11")
plt.xlabel("Frequncy")
plt.hlines((1025 * np.pi*(1.6/2)**2 *4 * 1.13/(1025*np.pi*(1.6/2)**2)),0,20,colors="blue", label = "A22 - theo")
plt.hlines(4*2 / 3 *1025*np.pi *(1.6/2)**3/(1025*np.pi*(1.6/2)**2),xmin=0,xmax=20,colors="orange" ,label = "A33 - theo")
plt.legend()
plt.show()
"""
