import numpy as np
import matplotlib.pyplot as plt
import control as ct
class system: 
    def __init__(self, R, R2, L, C):
        self.R2 = R2 
        self.L = L 
        self.C = C
        self.R =R
        self.createSystem()

    def createSystem(self):
        self.A = np.array([[(-self.R2 -self.R)/(self.C*self.R*self.R2), 0], [0, 0]])
        self.B = np.array([[1/(self.R2*self.C)], [1/self.L]])
        self.Cm = np.array([1,0])
        self.D = np.array([0]) 

        
def wyliczSystem(system, wejscie):
    dT = 1e-9
    koniecSymulacji = 1e-6 
    podstawaCzasu = np.arange(0, koniecSymulacji, dT) 

    AdT = dT * system.A
    BdT  = dT * system.B 

    X = np.zeros((len(system.B), 1))
    Y = np.zeros(len(podstawaCzasu))

    wejscie = np.ones(len(podstawaCzasu))
    for i in range(len(podstawaCzasu)-1):
        u = wejscie[i]
        u_plus = wejscie[i+1]
        X_zgadywane = X + AdT @ X + BdT * u

        X = X + 0.5 * (AdT @ X + BdT * u + AdT @ X_zgadywane + BdT * u_plus)

        Y[i+1] =( system.Cm @ X + system.D * u_plus).item()
    return X, Y, podstawaCzasu
    




uklad = system(R=1e3, R2=10e3, L=125e-6, C=50e-12) 
 
X,Y,  podstawaCzasu = wyliczSystem(uklad, 1) 

w1= plt.figure(1)
plt.plot(podstawaCzasu,Y)
plt.title('wlasne')
w1.show()

sys = ct.ss(uklad.A, uklad.B, uklad.Cm, uklad.D)
czasControl, yControl = ct.step_response(sys)
w2= plt.figure(2)
plt.plot(czasControl, yControl) 
plt.title('Z Control')
w2.show()

input()