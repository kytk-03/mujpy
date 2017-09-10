from mucomponents import muprompt # check
import numpy as np
import matplotlib.pyplot as P
from iminuit import Minuit as M # , describe, Struct
import probfit as PF

class MuSet(object):
    def __init__(self,run):
        '''
        checked with a notebook (YbTrensalAlphaCalibration, not yet as a class)       
        run = muload.musr2py()
        run.read('deltat_tdc_gps_0433.bin')
        mset = MuSet(run)
        initializes MuSet on a suitably high statistics run
        finds t0 from a mupropt fit
        finds self.firstbin, self.lastbin[detector] for background evaluation
        '''
        self.firstbin = 0
        self.prepeak = 7
        self.postpeak = 7
        self.second_plateau = 100
        self.peakheight = 100000.
        self.peakwidth = 1.   # broad guesses for default
        self.run = run



    def prepare(self, firstbin = 0, prepeak = 7, 
                      postpeak = 7, second_plateau = 100, 
                      peakheight = 100000., peakwidth = 1.):
        '''
        default is broad guesses from __init__,
        mut0 = MuSet(run)
        mut0.prepare(prepeak=10) # e.g. stores new prepeak value leaving others to default
        '''
        self.firstbin = firstbin
        self.prepeak = prepeak
        self.postpeak = postpeak
        self.second_plateau = second_plateau
        self.peakheight = peakheight
        self.peakwidth = peakwidth

    def promptfit(self, mplot = False, mprint = False):
        '''
        mset = MuSet(run)
        mset.fit(mplot=True, mprint=True)  
        fits peak positions 
        prints migrad results
        plots prompts and their fit (default no print, no plot)
        stores bins for background and t0        
        '''
        npeaks = np.array([np.argmax(self.run.get_histo_array_int(det)) for det in range(self.run.get_numberHisto_int())])
        # approximate positions of peaks
        nbin =  max(npeaks) + self.second_plateau # this sets a detector dependent second plateau bin interval
        x = np.arange(0,nbin,dtype=int) # nbin bins from 0 to nbin-1
        self.lastbin, np3s = npeaks - prepeak, npeaks + postpeak # final bin of first and initial bin of second plateaus
        mm = np.vectorize(muprompt)

        if mplot:
            fig,ax = P.subplots(2,3,figsize=(12,8)) 

        x0 = np.zeros(self.run.get_numberHisto_int())
        for detector in range(self.run.get_numberHisto_int()):
            # prepare for muprompt fit
            p = [ self.peakheight, float(npeaks[detector]), self.peakwidth, 
                  np.mean(self.run.get_histo_array_int(detector)[self.firstbin:self.lastbin[detector]]), 
                  np.mean(self.run.get_histo_array_int(detector)[np3s[detector]:nbin])]
            y = run.get_histo_array_int(detector)[:nbin]
            pars = dict(a=p[0],error_a=p[0]/100,x0=p[1]+0.1,error_x0=p[1]/100,dx=1.1,error_dx=0.01,
                        ak1=p[3],error_ak1=p[3]/100,ak2=p[4],error_ak2=p[4]/100)
            chi2 = PF.Chi2Regression(muprompt,x,y)
            level = 1 if mprint else 0
            m = M(chi2,pedantic=False,print_level=level,**pars)
            m.migrad()
            A,X0,Dx,Ak1,Ak2 = m.args
            x0[detector] = X0 # store float peak bin position (fractional) 
            if mplot:

                n1 = npeaks[detector]-50
                n2 = npeaks[detector]+50
                x3 = np.arange(n1,n2,1./10.)
                ax[divmod(detector,3)].plot(x[n1:n2],y[n1:n2],'.')
                ax[divmod(detector,3)].plot(x3,mm(x3,A,X0,Dx,Ak1,Ak2))
                ax[divmod(detector,3)].text(n1+5,0.8*max(y),'X0={:.2f}ns'.format(X0))
                ax[divmod(detector,3)].text(n1+5,0.6*max(y),'Dx={:.2f}ns'.format(Dx))
                P.show()
                ##################################################################################################
                # Simple cases: 
                # 1) Assume the prompt is entirely in bin nt0. (python convention, the bin index is 0,...,n,... 
                # The content of bin nt0 will be the t=0 value for this case and dt0 = 0.
                # The center of bin nt0 will correspond to time t = 0, time = (n-nt0 + mufit.offset + mufit.dt0)*mufit.binWidth_ns/1000.
                # 2) Assume the prompt is equally distributed between n and n+1. Then nt0 = n and dt0 = 0.5, the same formula applies
                # 3) Assume the prompt is 0.45 in n and 0.55 in n+1. Then nt0 = n+1 and dt0 = -0.45, the same formula applies.
                ##################################################################################################

        # these three are the sets of parameters used by other methods
        self.nt0 = x0.astype(int) # bin of peak, nd.array of shape run.get_numberHisto_int() 
        self.dt0 = x0-self.nt0 # fraction of bin, nd.array of shape run.get_numberHisto_int() 
        self.lastbin = self.nt0 - self.prepeak # nd.array of shape run.get_numberHisto_int() 
                                               # refresh, they may be slightly adjusted by the fit
        ##################################################################################################
