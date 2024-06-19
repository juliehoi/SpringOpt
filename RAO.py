import numpy as np
import matplotlib.pyplot as plt
import general_2D_model
from tqdm import tqdm


def RAO(steepn,WavePerVect,combinedSpringStiff,n):

    springStiffn = combinedSpringStiff*np.ones(n-1) #spring stiffness.
    RAO = np.zeros((n,len(WavePerVect))) #n x len(WavePerVect) which holds the RAO values for each module per wave period for a set spring stiffness

    for i in tqdm(range(len(WavePerVect))):

        x,spring_Force_Vector,global_positions,w_time,WA,m_inclAM,Mk,A,B,C = general_2D_model.main(WavePerVect[i],steepn,springStiffn,n)
        print("WA: ",WA)
        xposVectSteadyState = x[:n, len(x[0])//2:] #n x num_steps/2 vector. Eacho row is the position of each module in a timeseries after the transient period is over.
        for module in range(n):
            RAO[module][i] = np.max(xposVectSteadyState[module])/WA

    return(RAO)

#RAOtest = RAO(1/60,np.array([4,6,8,10]),10000,2)

"""
plt.plot(wavePeriodVect,RAOtest[0],marker='o', linestyle='-',label="module 1")
plt.plot(wavePeriodVect,RAOtest[1],marker='o', linestyle='-',label="module 2")
plt.xlabel("Wave period T [s]")
plt.ylabel("eta_surge/WA")
plt.legend()
plt.title("RAO")
plt.show()
"""

"""

TestingPeriods = np.array([2,3,4,5,6,7,8])

#Testing two different stiffnesses for module 1 and module
#S1
RAOtestS12mods = RAO(1/80,TestingPeriods,2*29332,2) #case 2
RAOtestS13mods = RAO(1/80,TestingPeriods,2*29332,3) #case 3
RAOtestS15mods = RAO(1/80,TestingPeriods,2*29332,5) #case 4
#S3
RAOtestS32mods = RAO(1/80,TestingPeriods,2*14899,2) #case 2
RAOtestS33mods = RAO(1/80,TestingPeriods,2*14899,3) #case 3
RAOtestS35mods = RAO(1/80,TestingPeriods,2*14899,5) #case 4

#plot 1
plt.plot(TestingPeriods,RAOtestS12mods[0],marker='o', linestyle='-',label="Case 2")
plt.plot(TestingPeriods,RAOtestS13mods[0],marker='o', linestyle='-',label="Case 3")
plt.plot(TestingPeriods,RAOtestS15mods[0],marker='o', linestyle='-',label="Case 4")
plt.title("Module 1 - Surge response - Spring S1")
plt.xlabel("Wave period T [s]")
plt.ylabel(r"$\eta_{1a}$/ $\zeta_{a}$ [-]")
plt.legend()
plt.show()

#plot 2
plt.plot(TestingPeriods,RAOtestS32mods[0],marker='o', linestyle='-',label="Case 2")
plt.plot(TestingPeriods,RAOtestS33mods[0],marker='o', linestyle='-',label="Case 3")
plt.plot(TestingPeriods,RAOtestS35mods[0],marker='o', linestyle='-',label="Case 4")
plt.title("Module 1 - Surge response - Spring S3")
plt.xlabel("Wave period T [s]")
plt.ylabel(r"$\eta_{1a}$/ $\zeta_{a}$ [-]")
plt.legend()
plt.show()

#plot 3
plt.plot(TestingPeriods,RAOtestS13mods[1],marker='o', linestyle='-',label="Case 3 - Module 2")
plt.plot(TestingPeriods,RAOtestS15mods[1],marker='o', linestyle='-',label="Case 4 - Module 2")
plt.plot(TestingPeriods,RAOtestS15mods[2],marker='o', linestyle='-',label="Case 4 - Module 3")
plt.plot(TestingPeriods,RAOtestS15mods[3],marker='o', linestyle='-',label="Case 4 - Module 4")
plt.title("Middle modules - Surge response - Spring S1")
plt.xlabel("Wave period T [s]")
plt.ylabel(r"$\eta_{1a}$/ $\zeta_{a}$ [-]")
plt.legend()
plt.show()

#plot 4
plt.plot(TestingPeriods,RAOtestS33mods[1],marker='o', linestyle='-',label="Case 3 - Module 2")
plt.plot(TestingPeriods,RAOtestS35mods[1],marker='o', linestyle='-',label="Case 4 - Module 2")
plt.plot(TestingPeriods,RAOtestS35mods[2],marker='o', linestyle='-',label="Case 4 - Module 3")
plt.plot(TestingPeriods,RAOtestS35mods[3],marker='o', linestyle='-',label="Case 4 - Module 4")
plt.title("Middle modules - Surge response - Spring S3")
plt.xlabel("Wave period T [s]")
plt.ylabel(r"$\eta_{1a}$/ $\zeta_{a}$ [-]")
plt.legend()
plt.show()
"""
