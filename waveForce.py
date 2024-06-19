import numpy as np

def w(time,ExcForce,phase,waveNumber,waveFreq,initPos,modules,dx): #send in like this w(time,EXCITING_FORCE,phase0,k,wf,x0)
    w = np.zeros(modules)
    for i in range(modules):
        w[i] = ExcForce*np.cos(waveFreq*time + phase -waveNumber*(initPos[i]-dx))
        #print("wf*time: ",wf*time)
        #print("phase0: ",phase0)
        #print("-k*(x0[i]-dx): ",-k*(x0[i]-dx))
        #print("wf:",wf)
        #print("np.cos(wf*time + phases[i]): ",np.cos(wf*time + phases[i]))
    return(w)
