import numpy as np
import matplotlib.pyplot as plt
import control as ct
import os
from engfmt import Quantity #do ladnego pisania w tytulach


wspolczynnikDlugosci = 10 # decyduje o koncu symulacji i czasie trwania sygnalow, na podstawie ich okresu


class system: 
    def __init__(self, R, R2, L, C):
        self.R2 = R2 
        self.L = L 
        self.C = C
        self.R =R
        self.createSystem()
        self.wyliczBode()
        self.dobierzKrok()
        

    def createSystem(self):
        self.A = np.array([[(-self.R2 -self.R)/(self.C*self.R*self.R2), 0], [0, 0]])
        self.B = np.array([[1/(self.R2*self.C)], [1/self.L]])
        self.Cm = np.array([1,0])
        self.D = np.array([0]) 

    def wyliczBode(self):
        self.systemCT = ct.ss(self.A, self.B, self.Cm, self.D)
        self.amplituda, self.faza, self.omega = ct.frequency_response(self.systemCT)
    def dobierzKrok(self): 
        self.dT = 1/(100*max(abs(np.linalg.eigvals(self.A))))

        
def wyliczSystem(system, wejscie):
    koniecSymulacji = len(wejscie)*uklad.dT
    podstawaCzasu = np.arange(0, koniecSymulacji, uklad.dT) 
    
    AdT = uklad.dT * system.A
    BdT  = uklad.dT * system.B 

    X = np.zeros((len(system.B), 1))
    Y = np.zeros(len(podstawaCzasu))

    for i in range(len(podstawaCzasu)-1):
        u = wejscie[i]
        u_plus = wejscie[i+1]
        X_zgadywane = X + AdT @ X + BdT * u

        X = X + 0.5 * (AdT @ X + BdT * u + AdT @ X_zgadywane + BdT * u_plus)

        Y[i+1] =( system.Cm @ X + system.D * u_plus).item()
    return X, Y, podstawaCzasu
    

def trojkat(amplituda, okres, uklad):
    a = amplituda/(okres/2)
    dlugosc= wspolczynnikDlugosci*okres 
    czas = np.arange(0, dlugosc, uklad.dT)
    sygnal = np.zeros(len(czas))
    
    for i in range(len(czas)):
        t = czas[i] %okres
        if (t < okres/2):
            sygnal[i] = a*t
        elif (t > okres/2):
            sygnal[i] = a*(okres-t)
        else:
            sygnal[i] = amplituda
    return(sygnal)
def prostokat(amplituda, okres, wspolczynnikWypelnienia, dlugosc, uklad):
    dlugosc=    wspolczynnikDlugosci*okres
    czas = np.arange(0, dlugosc, uklad.dT)
    sygnal = np.zeros(len(czas))
    for i in range(len(czas)):
        t = czas[i] %okres
        if (t < okres/(1/wspolczynnikWypelnienia)):
            sygnal[i] = 1
        elif (t > okres/(1/wspolczynnikWypelnienia)):
            sygnal[i] = 0
        else:
            sygnal[i] = 0.5
    return(sygnal)
    
def sinus(amplituda, okres, deltaphi, uklad):
    czas = np.arange(0, wspolczynnikDlugosci*okres, uklad.dT)
    sygnal = amplituda*np.sin(2*np.pi*(1/okres)*czas + deltaphi)
    return(sygnal)


def zapiszWpng(uklad, wejscie): 
    X,Y,podstawaCzasu = wyliczSystem(uklad, wejscie)

    R2 =Quantity(uklad.R2, 'ohm')
    C = Quantity(uklad.C, 'farads')
    L = Quantity(uklad.L, 'henr')
    R = Quantity(uklad.R, 'ohms') 

    R2 = R2.to_eng()
    C = C.to_eng()
    L =L.to_eng()
    R = R.to_eng()


    if not(os.path.isdir('figures')):
        os.makedirs("./figures")

    plt.plot(podstawaCzasu, Y)
    plt.title(f"Symulacja układu R_2 = {R2}, R = {R}, L = {L}, C = {C}")
    plt.savefig(f"./figures/odpowiedz___{R2}___{R}___{L}___{C}.png")
    plt.clf()

    plt.semilogx(uklad.omega/(2*np.pi), 20*np.log10(np.squeeze(uklad.amplituda)))
    plt.title(f"Charakterystyka amplitudowa układu R_2 = {R2}, R = {R}, L = {L}, C = {C}")
    plt.xlabel('Częstotliwość [Hz]')
    plt.ylabel('Wzmocnienie [dB]')
    plt.savefig(f"./figures/bodeWzmocnienie___{R2}___{R}___{L}___{C}.png")
    plt.clf()

    plt.semilogx(uklad.omega/(2*np.pi), np.squeeze(uklad.faza)*180/np.pi)
    plt.title(f"Charakterystyka fazowa układu R_2 = {R2}, R = {R}, L = {L}, C = {C}")
    plt.xlabel('Częstotliwość [Hz]')
    plt.ylabel('Faza [stopnie]')
    plt.savefig(f"./figures/bodeFaza___{R2}___{R}___{L}___{C}.png")
    plt.clf()

uklad = system(R=2e3, R2=10e3, L=125e-6, C=6e-9) 
X,Y,  podstawaCzasu = wyliczSystem(uklad, sinus(2,20e-3,0,uklad)) 
#X,Y,  podstawaCzasu = wyliczSystem(uklad, trojkat(2,20e-3,uklad)) 
#plt.plot(podstawaCzasu, sinus(2,10e-3,0,uklad))
plt.plot(podstawaCzasu, Y)
plt.show()
