from pylab import *
from qutip import *

def qubit_integrate(delta, eps0, A, gamma1, gamma2, psi0, tlist):
    
    # Hamiltonian
    sx = sigmax()
    sz = sigmaz()
    sm = destroy(2)
    
    H0 = - delta/2.0 * sx - eps0/2.0 * sz
    H1 = - A/2.0 * sz
 
    # collapse operators
    c_op_list = []
    
    n_th = 0.0 # zero temperature
    
    # relaxation
    rate = gamma1 * (1 + n_th)
    if rate > 0.0:
        c_op_list.append(sqrt(rate) * sm)
    
    # excitation
    rate = gamma1 * n_th
    if rate > 0.0:
        c_op_list.append(sqrt(rate) * sm.dag())
    
    # dephasing
    rate = gamma2
    if rate > 0.0:
        c_op_list.append(sqrt(rate) * sz)
    
    # evolve and calculate expectation values
    rho11 = basis(2,0)*basis(2,0).dag()
    rho12 = basis(2,1)*basis(2,0).dag()
    print rho12
    rho22 = basis(2,1)*basis(2,1).dag()

    # evolve and calculate expectation values
    H = [H0, [H1, 't']]
    W=H0+H1*tlist[-1]
    output = mesolve(H, psi0, tlist, c_op_list, [rho11,rho12,W*W,W], {})
    
    #return output.expect[0]
    return output.expect[1], output.expect[0], W #,output.expect[2]-output.expect[3]**2.

delta = 0.5 * 2 * np.pi   # qubit sigma_x coefficient
eps0  = 0.0 * 2 * np.pi   # qubit sigma_z coefficient
A     = 1. * 2 * np.pi   # sweep rate
gamma1 = 0.0          # relaxation rate
gamma2 = 0.1        # dephasing  rate
psi0 = basis(2,1)      # initial state

tlist = np.linspace(-50.0, 50.0, 5000)

aa=arange(0.1,15.2,0.2) * 2 * np.pi
output=open('ELZ_g1_%.3f_g1_%.3f.txt' %(gamma1,gamma2),'w')
for A in aa:
    rho12, rho11, H = qubit_integrate(delta, eps0, A, gamma1, gamma2, psi0, tlist)
    subplot(121)
    plot(A,abs(rho12[-1]),'ob')

    subplot(122)
    W=basis(2,0)*basis(2,0).dag()*rho11[-1]+basis(2,1)*basis(2,1).dag()*(1.-rho11[-1])+basis(2,0)*basis(2,1).dag()*(rho12[-1])+basis(2,1)*basis(2,0).dag()*(rho12[-1].conj())
    print H
    print W
    Es=W.eigenenergies()
    Psis=W.eigenstates()

    FI = 4*(Es[0]-Es[1])**2.*(Psis[1][0].dag()*H*Psis[1][1]).norm()**2.
    print FI, 4*(Psis[1][0].dag()*H*Psis[1][1]).norm()**2., 4*(Psis[1][0].dag()*H*Psis[1][1]*Psis[1][1].dag()*H*Psis[1][0])
    plot(A,FI/(4.*A*50)**2.,'ob')
    print >> output, A, 1-rho11[-1], FI, Es[0], Es[1], rho12[-1]


subplot(121)
plot(aa, 1 - np.exp(-np.pi * delta **2 / (2 * aa)) * np.ones(shape(aa)), '-r')
show()
