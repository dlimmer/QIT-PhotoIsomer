from pylab import *
import matplotlib as mpl
mpl.style.use('classic')

from numba import jit
from scipy.interpolate import interp1d

#Function that constructs the Lorenze Curves
@jit
def majorization():
    qtrans=arange(0,1,.01)
    potenq=array([])
    for q in qtrans:
        for q1 in qtrans:
            for q2 in qtrans:
                pf=zeros(4)
                pf[0]=q1
                pf[1]=1.-q-q1-q2
                pf[2]=q
                pf[3]=q2
                
                wfs=exp(-energies)
                Rfs=pf*exp(energies)
                
                l=argsort(Rfs)
                l=l[::-1]
                
                Rfs=sort(Rfs)
                Rfs=Rfs[::-1]
                
                sumwf=zeros(5)
                sumwf[0]=0.
                sumwf[1]=wfs[l[0]]
                sumwf[2]=wfs[l[0]]+wfs[l[1]]
                sumwf[3]=wfs[l[0]]+wfs[l[1]]+wfs[l[2]]
                sumwf[4]=wfs[l[0]]+wfs[l[1]]+wfs[l[2]]+wfs[l[3]]
                
                sumpf=zeros(5)
                sumpf[0]=0.
                sumpf[1]=pf[l[0]]
                sumpf[2]=pf[l[0]]+pf[l[1]]
                sumpf[3]=pf[l[0]]+pf[l[1]]+pf[l[2]]
                sumpf[4]=pf[l[0]]+pf[l[1]]+pf[l[2]]+pf[l[3]]
                
                sumrf=zeros(5)
                sumrf[0]=0.
                sumrf[1]=Rfs[0]
                sumrf[2]=Rfs[0]+Rfs[1]
                sumrf[3]=Rfs[0]+Rfs[1]+Rfs[2]
                sumrf[4]=Rfs[0]+Rfs[1]+Rfs[2]+Rfs[3]
                
                # Now check
                yinterp = np.interp(sumwf, sumw,sump)
                #yinterp=fun(sumwf)
                
                if yinterp[1]>=sumpf[1] and yinterp[2]>=sumpf[2] and yinterp[3]>=sumpf[3] and yinterp[4]>=sumpf[4]: potenq=append(q,potenq)

    return max(potenq)



ec=arange(-3,0.01,.01)

thermo=array([])
resource=array([])
Yflag=1
lCflag=0

delta=-3.0

if lCflag: counter=1

for eec in ec:
    #for delta in deltas:
    
    # Initial state information
    pminus=0.9

    p=zeros(4)
    p[0]=0.
    p[1]=0.5 #1.-pminus
    p[2]=0.
    p[3]=0.5 #pminus

    energies=array([0.,-3.,eec,delta])

    ws=exp(-energies)
    Rs=p*exp(energies)

    l=argsort(Rs)
    l=l[::-1]

    Rs=sort(Rs)
    Rs=Rs[::-1]

    sumw=zeros(5)
    sumw[0]=0.
    sumw[1]=ws[l[0]]
    sumw[2]=ws[l[0]]+ws[l[1]]
    sumw[3]=ws[l[0]]+ws[l[1]]+ws[l[2]]
    sumw[4]=ws[l[0]]+ws[l[1]]+ws[l[2]]+ws[l[3]]

    sump=zeros(5)
    sump[0]=0.
    sump[1]=p[l[0]]
    sump[2]=p[l[0]]+p[l[1]]
    sump[3]=p[l[0]]+p[l[1]]+p[l[2]]
    sump[4]=p[l[0]]+p[l[1]]+p[l[2]]+p[l[3]]

    sumr=zeros(5)
    sumr[0]=0.
    sumr[1]=Rs[0]
    sumr[2]=Rs[0]+Rs[1]
    sumr[3]=Rs[0]+Rs[1]+Rs[2]
    sumr[4]=Rs[0]+Rs[1]+Rs[2]+Rs[3]

    #print Rs
    if lCflag: subplot(1,3,counter)
    if lCflag: plot(sumw,sump,'ob-',ms=10,lw=3)
    #print sumr, sumw
    # Final state

    fun=interp1d(sumw,sump)
    
    try: qtrans=majorization()
    except: qtrans=nan
    
    print eec,qtrans
    
    resource=append(resource,qtrans)
    thermo=append(thermo,exp(-eec)/sum(exp(-energies)))

    if lCflag:
        pf[0]=0.
        pf[1]=1.-qtrans
        pf[2]=qtrans
        pf[3]=0.

        wfs=exp(-energies)
        Rfs=pf*exp(energies)

        l=argsort(Rfs)
        l=l[::-1]

        Rfs=sort(Rfs)
        Rfs=Rfs[::-1]

        sumwf=zeros(5)
        sumwf[0]=0.
        sumwf[1]=wfs[l[0]]
        sumwf[2]=wfs[l[0]]+wfs[l[1]]
        sumwf[3]=wfs[l[0]]+wfs[l[1]]+wfs[l[2]]
        sumwf[4]=wfs[l[0]]+wfs[l[1]]+wfs[l[2]]+wfs[l[3]]

        sumpf=zeros(5)
        sumpf[0]=0.
        sumpf[1]=pf[l[0]]
        sumpf[2]=pf[l[0]]+pf[l[1]]
        sumpf[3]=pf[l[0]]+pf[l[1]]+pf[l[2]]
        sumpf[4]=pf[l[0]]+pf[l[1]]+pf[l[2]]+pf[l[3]]

        sumrf=zeros(5)
        sumrf[0]=0.
        sumrf[1]=Rfs[0]
        sumrf[2]=Rfs[0]+Rfs[1]
        sumrf[3]=Rfs[0]+Rfs[1]+Rfs[2]
        sumrf[4]=Rfs[0]+Rfs[1]+Rfs[2]+Rfs[3]
        subplot(1,3,counter)
        print sumrf, pf[l[0]]/wfs[l[0]],pf[l[0]]/wfs[l[0]]+pf[l[1]]/wfs[l[1]]
        plot(sumwf,sumpf,'or-',ms=10,lw=3)

        print "And the max value of qtrans is:", qtrans

        xlabel(r'$\sum_i \, e^{-\beta \epsilon_i}$',size=18)
        ylabel(r'$\sum_i \, p_i$',size=18)
    if lCflag: counter=counter+1

if Yflag:
    
    zz=zeros([len(ec),3])
    zz[:,0]=ec
    zz[:,1]=thermo
    zz[:,2]=resource
    savetxt('four_state_initial_%.2f_%.2f_%.2f_%.2f_energy_%.2f_%.2f_%.2f_%.2f.txt' %(p[0],p[1],p[2],p[3],energies[0],energies[1],energies[2],energies[3]),zz)
    plot(ec,resource,'-b',lw=3.)
    plot(ec,thermo,'-r',lw=3.)
    xlabel(r'$\beta \Delta \epsilon$',size=22)
    ylabel(r'$q_\mathrm{trans}$',size=22)

savefig('resource.png')
tight_layout()
show()
