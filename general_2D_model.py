import numpy as np
import WadamReader
import waveForce
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.io import savemat

#simulation setup
dt = 0.005  #Time step
t_final = 200  #Final time
num_steps = int(t_final / dt) + 1
t = np.linspace(0, t_final, num_steps) #Time array

modules = 8
states = 2*modules
dx = 13.2 #[m], distance between the mass centers of the modules (Onsrud 2019)
m = 9088 #Per module
km = 4280 #Mooring stiffness page 42 (Sølvberg,2023)

def systemMatrices(n,m_inclAm,d):
    #Setting the sizes of the matrices
    C = np.zeros((n,2*n))
    B = np.zeros((2*n,n))
    A0 = np.zeros((2*n,2*n))
    Mk = [] #Holds the n-1 Mk matrices that will be multiplied by the stiffnesses between the modules
    A0[n][0] = -km/m_inclAm
    A0[(2*n)-1][n-1] = A0[n][0]

    for i in range(n):
        #C[i][i] = 1 #If z = motion of each modules

        if i == n-1: #If z = relative motion between each module
            break
        else:
            C[i][i] = -1
            C[i][i+1] = 1


    for i in range(2*n):
        if i > n-1:
            B[i][i-n] = 1/m_inclAm

    for i in range(2*n):
        for j in range(2*n):
            if i < n:
                if j == i:
                    A0[i][j+n] = 1

            else:
                if j > n-1:
                    A0[i][i] = -d/m_inclAm

    for i in range(n-1):
        Mk.append(np.zeros((2*n,2*n)))
        for j in range(2):
            if j == 0:
                Mk[i][n+j+i][i] = -1/m_inclAm
                Mk[i][n+j+i][i+1] = 1/m_inclAm
            else:
                Mk[i][n+j+i][i] = 1/m_inclAm
                Mk[i][n+j+i][i+1] = -1/m_inclAm
    systemMatricesMatlab = {"A0":A0,"B":B,"C":C,"Mk":Mk}
    savemat("C:/Users/Julie/OneDrive - NTNU/Vår2024/Masteroppgave/Matlab/SystMatr.mat",systemMatricesMatlab) #Save in matlab directory
    print("A0, B, C and Mk for a system with ",modules," modules has been saved in the matlab folder as a workspace")
    return(A0,B,C,Mk)

constants, periodToDamping, periodToAddedMass, periodToMass, periodToRestoring,periodToHeadingToMotionDim,periodToHeadingToFexcDim, HeadingAngles = WadamReader.readWadam("WADAM2.txt")

RO   =  constants['RO']
G    =  constants['G']
VOL  =  constants['VOL']
L    =  constants['L']
WA   =  constants['WA'] #(IF NOT GIVEN AS INPUT)
#EXCITING_FORCES_FACTOR = (RO*VOL*G*WA)/L #To re-dimensionalize
#print("EXCITING_FORCES_FACTOR: ",EXCITING_FORCES_FACTOR)

#Wave periods to choose from: 2s-20s
def main(WavePeriod,WaveSteepness,springVector,n): #send in wave period on the format 2.0, 3.0, 4.0, etc. wave stiffness = H/lambda
    print("Wave period: ",WavePeriod, ". Wave steepness: ",WaveSteepness)

    WA = (1/(4*np.pi))*WaveSteepness*G*(WavePeriod*WavePeriod) #Formula obtained from (Faltinsen,1990)
    print("Wave amplitude found from the wave steepness: ",WA)

    mA_11 = periodToAddedMass[WavePeriod][0][0] #Added mass in direction 11
    print("added_mass: ",mA_11)
    m_inclAM = m + mA_11

    #Wave period, heading angle, wave force direction, response direction
    EXCITING_FORCE = periodToHeadingToFexcDim[WavePeriod][0][0][2] #Absolute value of excitation force [N] in surge direction for 0 heading angle. Multiplied by WA futher down (state space loop) to account for the wave steepness.
    print("Exciting force: ", EXCITING_FORCE)
    phase0 = np.deg2rad(periodToHeadingToFexcDim[WavePeriod][0][0][3]) #[rad]
    print("Phase: ",phase0)
    #From sealoads compendium
    lambd = (G/(2*np.pi))*(WavePeriod*WavePeriod)
    wf = (2*np.pi)/WavePeriod #Wave frequency [rad/s],
    k = (wf*wf)/G #wave number


    #5000-45000 N/m, stiffness range for one of the two springs, page 42 (Sølvberg,2023)
    d =periodToDamping[WavePeriod][0][0] #potential damping and viscous drag from Morison’s theory
    print("d: ",d)

    #setting system matrices
    A0, B, C, Mk = systemMatrices(n,m_inclAM,d)
    print("A0 :",A0)
    A = A0
    for i in range(len(Mk)):
        print("springVector: ",springVector)
        print("Mk: ",Mk)
        print("springVector[i]*Mk[i] :", springVector[i]*Mk[i])
        A += springVector[i]*Mk[i]

    #Initial conditions in the global coordinate system
    x0 = np.zeros(2*n)
    for i in range(2*n):
        if i>n-1:
            x0[i] = np.array([0])
        else:
            x0[i] =np.array([dx*(i+1)])
    print("x0: ",x0)

    #store global positions for each mass
    global_positions = np.zeros((n, num_steps))
    #x0Pos = x0[:n]
    for i in range(n):
        global_positions[i,0] = x0[i]

    print("A: ",A)
    print("B: ",B)

    #finding eigenvalues
    eigvalA, eigvecA = np.linalg.eig(A)
    eigvalA_array = np.array([x for x in eigvalA])
    print("eigvalA: ",eigvalA)

    """
    #plotting eigenvalues
    plt.plot(eigvalA_array.real,eigvalA_array.imag, marker='o', linestyle='', color='blue')
    plt.xlabel('Real part')
    plt.ylabel('Imaginary part')
    plt.axvline(0, color='black', linewidth=0.5)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.title('Eigenvalues of A')
    plt.grid(True)
    plt.plot()
    plt.show()
    """

    x = np.zeros((2*n, num_steps))
    w_time = np.zeros((n, num_steps))

    """
    #If other initial conditions is to be used for the states
    for i in range(states):
        x[i,0] = 0
        #x[:, 0] = x0.flatten()
        #print("x: ",x)
    """
    for i in tqdm(range(1, num_steps)):
        #Gradually increase the wave amplitude
        ramp_duration = 0.3*t_final
        t_ramp = np.arange(0, ramp_duration, dt)
        if i*dt < ramp_duration:
            ramp_value = np.linspace(0, 1, len(t_ramp))
            WA_run= ramp_value[i]*WA
        else:
            WA_run = WA

        #print("w(t[i]): i steg ",i, " ",w(t[i]))

        x_dot = np.matmul(A,x[:, i-1]) + np.matmul(B,waveForce.w(t[i],EXCITING_FORCE*WA_run,phase0,k,wf,x0,n,dx))

        #print("np.matmul(B1,w(t[i])): ",np.matmul(B1,w(t[i])))
        #print("x_dot i steg: ",i,": ",x_dot)
        #print("bidrag i steg: ",i,", Ax: ",np.matmul(A,x[:, i-1]),", B1w: ",np.matmul(B1,w(t[i])))

        x[:, i] = x[:, i-1] + x_dot * dt
        w_time[:,i] = waveForce.w(t[i],EXCITING_FORCE*WA_run,phase0,k,wf,x0,n,dx)

        # Update global positions for each mass
        for j in range(n):
            global_positions[j, i] = x[j, i] + x0[j]  #()+ j * dx) Add initial position to relative position
        #print("global_positions: ", global_positions)

        spring_Force_Vector = np.zeros(((n-1),num_steps))
        for i in range(len(spring_Force_Vector)):
            spring_Force_Vector[i,:] = springVector[i]*(x[i+1,:]-x[i,:])

    return(x,spring_Force_Vector,global_positions,w_time,WA,m_inclAM,Mk,A,B,C)

def plotTimeSeriers(WavePeriod,steepn,springVector):
    x, spring_Force_Vector,global_positions,w_time,WA,m_inclAM,Mk,A,B,C = main(WavePeriod,steepn,springVector,modules)
    plt.figure(figsize=(10, 6))
    for j in range(modules):
        plt.plot(t, x[j,:], label=f'Module {j + 1}')
        plt.plot(t, global_positions[j, :], label=f'Module {j + 1}')
        plt.plot(t,([global_positions[j,0]]*len(t)), label =f'Initial value for Module {j + 1}' )

    plt.xlabel('Time [s]')
    plt.ylabel('Surge [m]')
    plt.title('Surge response from the initial positions. k = '+ str(springVector[0]) + "[N/m]. T = " + str(WavePeriod) + "[s]")
    plt.legend()
    plt.grid(True)
    plt.show()

    figure, axis = plt.subplots(2, 2)
    for i in range(modules):
        axis[0, 0].plot(t,x[i,:],label="Module " + str(i+1))
    axis[0, 0].legend()
    axis[0, 0].set_title("Surge response")
    axis[0, 0].set_xlabel('Time [s]')
    axis[0, 0].set_ylabel('[m]')

    for i in range(modules-1):
        axis[1, 0].plot(t, (x[i+1,:]-x[i,:]),label="Between " + str(i+2)+" and "+ str(i+1))
    axis[1,0].legend()
    axis[1, 0].set_title("Relative motion between the modules")
    axis[1, 0].set_xlabel('Time [s]')
    axis[1, 0].set_ylabel('[m]')

    for i in range(modules-1):
        axis[1, 1].plot(t, spring_Force_Vector[i],"--",label="Spring "+str(i+1))
    axis[1,1].legend()
    axis[1, 1].set_title("Spring force between modules")
    axis[1, 1].set_xlabel('Time [s]',)
    axis[1, 1].set_ylabel('[N]')
    #Combine all the operations and display
    plt.suptitle('k = ' + str(springVector[0]) + "[N/m]. T = " + str(WavePeriod) + "[s]")
    plt.show()

    for j in range(modules):
        plt.plot(t, w_time[j, :], label=f'Wave force u on Module {j + 1}')

    plt.xlabel('Time [s]')
    plt.ylabel('Wave force [N]')
    plt.title('Wave excitation force. T = ' + str(WavePeriod) + "[s]")
    plt.legend()
    plt.grid(True)
    plt.show()



def runMultipleStiffnesses(modules,stiffnessVector,WavePeriod,WaveSteepness):
    #points for the k,response plots
    max_surge_response = np.zeros((modules,len(stiffnessVector)))
    max_abs_force_response = np.zeros((modules-1,len(stiffnessVector)))
    max_relative_response = np.zeros((modules-1,len(stiffnessVector)))

    after_transient = np.int(np.rint(0.5*num_steps))
    print("Ater transient: ",after_transient)

    for k in tqdm(range(len(stiffnessVector))):
        x,spring_Force_Vector,global_positions,w_time,WA,m_inclAM,Mk,A,B,C = main(WavePeriod,WaveSteepness,stiffnessVector[k],modules)
        print("x: ",x)
        print("spring_Force_Vector[after_transient:] :",spring_Force_Vector[after_transient:])
        #max_abs_force_response[k] = np.max(np.abs(spring_Force_Vector[after_transient:]))
        #max_relative_response[k] = np.max(x[1,after_transient:]-x[0,after_transient:])
        for i in range(modules):
            if i < modules-1:
                max_surge_response[i,k] = np.max(x[i,after_transient:]) #makes sure the max response is not taken from the transient response. Is 0.2 ok?
                print("x[i,after_transient:]: ",x[i,after_transient:])
                max_abs_force_response[i,k] = np.max(np.abs(spring_Force_Vector[i,after_transient:]))
                max_relative_response[i,k] = np.max(x[i+1,after_transient:]-x[i,after_transient:])
            else:
                max_surge_response[i,k] = np.max(x[i,after_transient:])
    print("max_surge_response: ",max_surge_response)
    print("max_abs_force_response: ",max_abs_force_response)
    print("max_relative_response: ",max_relative_response)
    #plotting
    for i in range(modules-1):
        #plt.plot(np.array([arr[0] for arr in stiffnessVector])/(2*1000),max_abs_force_response[i]/(1000),'--',label = "Force in spring "+ str(i+1)) #dividing spring force by 2 to get force per spring
        plt.plot(np.array([arr[0] for arr in stiffnessVector])/(2*1000),max_relative_response[i]*10,label= "Motion between mod. " + str(i+2) + "," + str(i+1) )
        #plt.plot(np.array([arr[0] for arr in stiffnessVector])/(1000),max_abs_force_response[i]/(1000),'--',label = "Force in spring "+ str(i+1)) #dividing spring force by 2 to get force per spring
        #plt.plot(np.array([arr[0] for arr in stiffnessVector])/(1000),max_relative_response[i]*10,'--',label= "Motion between mod. " + str(i+2) + "," + str(i+1) )
    plt.xlabel("Spring stiffness k [kN/m] per spring")
    plt.ylabel("Maximum relative motion [dm]") # "Max. total spring force Fs [kN], +
    plt.legend()
    plt.title("n = " + str(modules) + ", T = " + str(WavePeriod) + "[s]")
    plt.show()
    for j in range(modules):
        plt.plot(np.array([arr[0] for arr in stiffnessVector])/(2*1000),max_surge_response[j,:],label = "Module "+ str(j+1))
        #plt.plot(np.array([arr[0] for arr in stiffnessVector])/(1000),max_surge_response[j,:],'--',label = "Module "+ str(j+1))
    plt.xlabel("Spring stiffness k [kN/m] per spring")
    plt.ylabel("Max. Surge response [m]")
    plt.legend()
    plt.title("n = " + str(modules) + ", T = " + str(WavePeriod) + "[s]")
    plt.show()



def findEigenVal(Mk,stiffnvect,m_inclAm):
    mooringMatrix = np.zeros((modules,modules))
    mooringMatrix[0][0] = km
    mooringMatrix[modules-1][modules-1] = km
    mass = m_inclAm*np.identity(modules)

    K = mooringMatrix
    for i in range(modules-1):
     K+= -(stiffnvect[i]*m_inclAm*(Mk[i][modules:, :modules]))
    M_inv = np.linalg.inv(mass)

    eigenvalues, eigenvectors = np.linalg.eig(np.matmul(K,M_inv))
    eigenfrequencies = np.sqrt(eigenvalues) #[rad/s]
    eigen_periods =  2*np.pi/eigenfrequencies #[s]
    #eig_freq_rad = 2*np.pi/eigenfrequencies #[rad/s]
    return(eigenfrequencies,eigen_periods)

k_arr = np.linspace(5000,45000,35)
k_arr_2springs = 2*k_arr
k = []
for i in range(len(k_arr)):
    k.append([k_arr_2springs[i]]*modules)

#runMultipleStiffnesses(modules,k,10,1/60)

#optimal_k = np.array([5.0121*10**3,4.9870*10**3,5.0121*10**3])
plotTimeSeriers(4,1/80,[2*29332]*(modules-1)) #run this to update the system matrices
#plotTimeSeriers(10,1/60,optimal_k) #run this to update the system matrices

"""
resFrec = 2*10.981*1000
resfrecVect =resFrec*np.ones(modules-1)
x,spring_Force_Vector,global_positions,w_time,WA,m_inclAM,Mk,A,B,C = main(6,1/60,resfrecVect,modules)
ef,ep = findEigenVal(Mk,resfrecVect,m_inclAM)
print("Eigenfrequencies of k = "+ str(resfrecVect[0]) + "N/m,"+ str(modules)+ " modules: " + str(ef) + " [rad/s]")
print("ep: " + str(ep) + "[s]")

kl  = [5500]*(modules-1) #list of spring stiffnesses

x,spring_Force_Vector,global_positions,w_time,WA,m_inclAM,Mk,A,B,C = main(4,1/60,kl,modules)
ef,ep = findEigenVal(Mk,kl,m_inclAM)
print("Eigenfrequencies of k = "+ str(kl[0]) + "N/m, 2 modules: " + str(ef) + " [rad/s]")
print("ep: " + str(ep) + "[s]")
"""
