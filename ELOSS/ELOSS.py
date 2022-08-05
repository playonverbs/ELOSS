import numpy as np

Z = 18
A = 39.948   # g / mol
I = 188.0*(10**(-6)) # MeV
K = 0.307 # MeV * cm^2 / mol
Mmu = 105.658 # MeV for muon
Me  = 0.51 # MeV for electron
rho = 1.396 # g/cm3

def beta(gamma):
    return np.sqrt(1-(1./(gamma**2)))

def gamma(KE,mass):
    return (KE/mass)+1

beta = np.vectorize(beta)
gamma = np.vectorize(gamma)

def Wmax (KE,mass):
    g = gamma(KE,mass)
    b = beta(g)
    num = 2*Me*((b*g)**2)
    den = 1 + 2*g*Me/mass + (Me/mass)**2
    return num/den

# density correction for LAr
def density(bg):

    # constants and variable names obtained from :
    # PDG elos muons table [http://pdg.lbl.gov/2016/AtomicNuclearProperties/MUE/muE_liquid_argon.pdf]
    
    C  = -5.2146
    X0 = 0.2
    X1 = 3.0
    a  = 0.19559
    m  = 3.0
    N    = 2 * np.log(10)
    
    x = np.log10(bg)
    
    if (x < X0):
        return 0.
    if (x > X1):
        return N * x + C
    addition = a*((X1-x)**m)
    return N * x + C + addition


class Eloss:

    def __init__(self,I=188,rho=1.38):
        self.I   = I   # mean excitation energy [eV]
        self.rho = rho # argon density [g/cm3]

    # set density [input in g/cm3]
    def setDensity(rho):
        self.rho = rho

    # set mean excitation energy I [input in eV]
    def setI(I):
        self.I = I

    # print argon properties currently loaded
    def PrintArInfo():
        print
        print 'Argon density .............. : %.03f g/cm3'%self.rho
        print 'Argon mean excitation energy : %.03f eV'%self.I
        print
    

# KE in MeV
# x in cm
# mass in MeV
def dpdx(KE,x,mass,I=I):
    g = gamma(KE,mass)
    b = beta(g)
    epsilon = (K/2.)*(Z/A)*(x*rho/(b*b))
    A0 = (2*Me*(b*g)**2)/I
    A1 = epsilon/I
    return (1./x) * epsilon * (np.log(A0) + np.log(A1) + 0.2 - (b*b) - density(b*g))

# in MeV/cm
def dedx(KE,mass,I=I,dens=True):
    g = gamma(KE,mass)
    b = beta(g)
    F = K * (Z/A)*(1/b)**2
    wmax = Wmax(KE,mass)
    a0 = 0.5*np.log( 2*Me*(b*g)**2 * wmax / (I*I) )
    ret = a0 - b*b
    if (dens == True):
        ret -= density(b*g)/2.
    return F * ret


def Tfunc(g):
    return (g-1.)*Me
def tfunc(g):
    return (g-1.)

def Fminus(b,t):
    f = (1-b*b) * ( 1 + t*t/8. - (2*t+1)*np.log(2) )
    return f
def Fplus(b,t):
    f = 2*np.log(2) - (b*b/12.) * ( 23. + 14./(t+2.) + 10./(t+2.)**2 + 4./(t+2.)**3 )
    return f

def dedxelectrons(beta,dens=True):
    g = 1./np.sqrt(1-beta*beta)
    T = Tfunc(g)
    t = tfunc(g)
    F = (Z/A) * 0.153536 * ((1./beta)**2)
    ret =  (np.log((T/I)**2)) + np.log(1.+t/2.) + Fminus(beta,t)
    if (dens == True):
        ret -= density(beta*g)
    return ret * F

def dedxpositrons(beta,dens=True):
    g = 1./np.sqrt(1-beta*beta)
    T = Tfunc(g)
    t = tfunc(g)
    F = (Z/A) * 0.153536 * ((1./beta)**2)
    ret =   (np.log((T/I)**2)) + np.log(1+t/2.) + Fplus(beta,t)
    if (dens == True):
        ret -= density(beta*g)
    return ret * F
