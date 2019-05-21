from pylab import *
from qutip import *

eec=arange(-1.5,1.,4.)
delta=0.
#ec=-1
for ec in eec:
    V=0.1
    
    # -3 -1.5 0 0
    
    H=basis(4)*basis(4).dag()*(-3.)+basis(4,1)*basis(4,1).dag()*(ec)+basis(4,2)*basis(4,2).dag()*(delta)
    H=H+basis(4)*basis(4,2).dag()*0.
    H=H+basis(4,2)*basis(4).dag()*0.

    H=H+basis(4,3)*basis(4,2).dag()*V
    H=H+basis(4,2)*basis(4,3).dag()*V

    H=H+basis(4,1)*basis(4,2).dag()*0.
    H=H+basis(4,2)*basis(4,1).dag()*0.


    W=H.eigenenergies()
    print(W)
    psis=H.eigenstates()

    H=basis(4)*basis(4).dag()*W[0]
    H=H+basis(4,1)*basis(4,1).dag()*W[1]
    H=H+basis(4,2)*basis(4,2).dag()*W[2]
    H=H+basis(4,3)*basis(4,3).dag()*W[3]

    print(H)
    
    T=1.
    B=0.1
    A=0.001
    C=0.0

    Gamma02=basis(4,2)*basis(4).dag()*sqrt(A)

    Gamma20=basis(4)*basis(4,2).dag()*sqrt(A)*sqrt(exp(-(W[0]-W[2])/T))

    Gamma03=basis(4,3)*basis(4).dag()*sqrt(A)
    
    Gamma30=basis(4)*basis(4,3).dag()*sqrt(A)*sqrt(exp(-(W[0]-W[3])/T))
    

    Gamma12=basis(4,2)*basis(4,1).dag()*sqrt(B)

    Gamma21=basis(4,1)*basis(4,2).dag()*sqrt(B)*sqrt(exp(-(W[1]-W[2])/T))

    Gamma13=basis(4,3)*basis(4,1).dag()*sqrt(B)
    
    Gamma31=basis(4,1)*basis(4,3).dag()*sqrt(B)*sqrt(exp(-(W[1]-W[3])/T))


    Gamma23=basis(4,3)*basis(4,2).dag()*sqrt(C)
    
    Gamma32=basis(4,2)*basis(4,3).dag()*sqrt(C)*sqrt(exp(-(W[2]-W[3])/T))

    #psi0=psis[1][2]*sqrt(1./3.)+psis[1][3]*sqrt(1./3.)+psis[1][0]*sqrt(1./3.)

    psi0=psis[1][2]*sqrt(1./2.)+psis[1][3]*sqrt(1./2.)



    tlist=arange(0.0,1000,.01)

    #rho0=basis(3,2)*basis(3,2).dag()
    #rho1=basis(3,1)*basis(3,1).dag()
    #rho2=basis(3,0)*basis(3,0).dag()

    rho0=psis[1][0]*psis[1][0].dag()
    rho1=psis[1][1]*psis[1][1].dag()
    rho2=psis[1][2]*psis[1][2].dag()
    rho3=psis[1][3]*psis[1][3].dag()

    figure(1)
    #output=essolve(H,rho0,tlist,Gamma12,rho1)
    output=mesolve(H, psi0, tlist, [Gamma12,Gamma21,Gamma02,Gamma20,Gamma03,Gamma30,Gamma23,Gamma32], [rho0,rho1,rho2,rho3])
    plot(tlist,output.expect[2],'-r',lw=2,alpha=1)
    plot(tlist,output.expect[1],'-b',lw=2,alpha=1)
    plot(tlist,output.expect[3],'-c',lw=2,alpha=1)
    plot(tlist,output.expect[0],'-k',lw=2,alpha=1)
    
    Q=exp(-W[0])+exp(-W[1])+exp(-W[2])+exp(-W[3])
    plot(tlist,output.expect[2]*0+exp(-W[2])/Q,'--r',lw=2,alpha=1)
    plot(tlist,output.expect[1]*0+exp(-W[1])/Q,'--b',lw=2,alpha=1)
    plot(tlist,output.expect[3]*0+exp(-W[3])/Q,'--c',lw=2,alpha=1)
    plot(tlist,output.expect[0]*0+exp(-W[0])/Q,'--k',lw=2,alpha=1)
    
    zz=zeros([len(tlist),5])
    zz[:,0]=tlist
    zz[:,1]=output.expect[0]
    zz[:,2]=output.expect[1]
    zz[:,3]=output.expect[2]
    zz[:,4]=output.expect[3]
    savetxt('pop4_0.5_0.5_en_1.5_del_0.0.txt',zz)
    
    #plot(tlist,output.expect[1]/output.expect[2],'-k',lw=2,alpha=1)

    #plot(tlist,tlist*0+exp(-1.5),'--')
    #xlabel(r'$t$')
    #ylabel(r'$\rho$')
    #plot(tlist,tlist*0+sqrt(B/A),'--')

    #plot(tlist,tlist*0+exp(1.5)/(1+exp(1.5)+exp(3)),'--')
    #plot(tlist,tlist*0+exp(3.)/(1+exp(1.5)+exp(3)),'--')

    #figure(2)
    #plot(ec,max(output.expect[1]),'ob')

    #plot(ec,exp(-ec)/(exp(-W[0])+exp(-W[1])+exp(-W[2])),'-or')
    #plot(ec,0.5*exp(-ec)/(exp(-W[2]))/(1+exp(-ec)/(exp(-W[2]))),'-ok')

#savefig('resource_lin.png')
show()
