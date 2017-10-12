from numpy import cos, pi, exp, sqrt
import numpy as np
from scipy.special import dawsn,erf, j0

class mumodel(object):
    def _init_(self,x,y,e=1):
        ''' generic __init__ for all mufunction classes
        x, y, e are numpy arrays
        e is always defined, 
        either provided by the caller 
        or default to np.ones(x.shape[0])
        '''
        if x.shape[0]!=y.shape[0]:
            raise ValueError('x, y have different lengths')
        else:
            self.x = x
            self.y = y
        if e==1:
            self.e = np.ones(x.shape[0])
        else:
            if e.shape[0]!=x.shape[0]:
                raise ValueError('x, e have different lengths')           
            else:
                self.e = e
        self.radeg = pi/180.
        self.gamma_mu = 135.5
        self._help_={'bl':r'Lorentz decay: $\mbox{asymmetry}\exp(-\mbox{Lor_rate}t)$',
                     'bg':r'Gauss decay: $\mbox{asymmetry}\exp(-0.5(\mbox{Gau_rate}t)^2)$',
                     'bs':r'Gauss decay: $\mbox{asymmetry}\exp(-0.5(\mbox{rate}t)^\beta)$',
                     'mg':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-0.5(\mbox{Gau_rate}t)^2)$',
                     'ml':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-\mbox{Lor_rate}t)$',
                     'ms':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-(\mbox{rate}t)^\beta)$',
                     'jg':r'Gauss Bessel: $\mbox{asymmetry} j_0[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-0.5(\mbox{Lor_rate}t)^2)$',
                     'jl':r'Lorentz Bessel: $\mbox{asymmetry}j_0[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-0.5(\mbox{Lor_rate}t)^2)$',
                     'fm':r'FMuF: $\mbox{asymmetry}/6[3+\cos 2*\pi\gamma_\mu\mbox{dipfield}\sqrt{3} t +
               (1-1/\sqrt{3})\cos \pi\gamma_\mu\mbox{dipfield}(3-\sqrt{3})t +
               (1+1/\sqrt{3})\cos\pi\gamma_\mu\mbox{dipfield}(3+\sqrt{3})t ]\exp(-\mbox{Lor_rate}t)$',
                     'kg':r'Gauss Kubo-Toyabe: $\mbox{asymmetry}\exp(-0.5(\mbox{Gau_rate}t)^2)$',

    # ---- end generic __init__
    def _available_components_(self):
        from iminuit import describe
        '''
        Returns a template tuple of dictionaries (one per fit component):
        Each dictionary contains 'name' and 'pars', 
        the latter in turns is a list of dictionaries, one per parameter, 'name','error,'limits'
        ({'name':'bl','pars':[{'name':'asymmetry','error':0.01,'limits'[0,0]},
                              {'name':'Lor_rate','error':0.01,'limits'[0,0]}}, 
         ...)
        retreived magically from the class.
        Used in mufit.MuFit as 
           from mujpy.mucomponent.mucomponent import mumodel
           self.available_components = mumodel()._available_components_()         
        '''
        _available_components = [] # is a list, mutable
        # generates a template of available components.
        for name in [module for module in dir(self) if module[0]!='_']: # magical extraction of component names
            pars = describe(self.__class__.__dict__[name])[2:]            # the [2:] is because the first two arguments are self and x
            _pars = []
            # print('pars are {}'.format(pars))
            for parname in pars:
            # The style is like iminuit fitargs, but not exactly,
            # since the latter is a minuit instance:
            # it will contain parameter name: parname+str(k)[+'_'+str(nrun)]
            # error_parname, fix_parname (False/True), limits_parname, e.g.
            #   {'amplitude1_354':0.154,'error_amplitude1_354':0.01,'fix_amplitude1_354':False,'limits_amplitude1_354':[0, 0]
            # 
            # In this template only
            #   {'name':'amplitude','error':0.01,'limits':[0, 0]}
                error, limits = 0.01, [0, 0] # defaults
                if parname == 'field' or parname == 'phase' or parname == 'dipfield': error = 1.0
                if parname == 'beta': error,limits = 0.05, [1.e-2, 1.e2]
                # add here special cases for errors and limits, e.g. positive defined parameters
                _pars.append({'name':parname,'error':error,'limits':limits})
            _available_components.append({'name':name,'pars':_pars})
        return tuple(_available_components) # transformed in tuple, immutable

    def _add_(self,x,**kwargv):
        '''
        used by self._chisquare_ 
        kwargs is a dictionary of dictionaries
        use to plot e.g. a blmg model: 
          m=mumodel() 
          m._init_(x,y,e) 
          plt.errorbar(m.x,m.y,yerr=m.e)
          plt.plot(m.x,m._add_(m.x,**{'c1':{m.bl:[asy1,rate1]},'c2':{m.mg:[asy2,fld2,ph2,rate2]}})                    
        'c1','c2' are unique but dummy
        '''  
        f = np.zeros(x.shape[0])
        for num, dic in kwargv.items():
            for component, parameters in dic.items():
                f += component(x,*parameters)
        return f     

    def bl(self,x,asymmetry,Lor_rate): 
        '''
        fit component for a Lorentzian decay, 
        x [mus], asymmetry, Lor_rate [mus-1]
        '''
        return asymmetry*exp(-x*Lor_rate)

    def bg(self,x,asymmetry,Gau_rate): 
        '''
        fit component for a Gaussian decay, 
        x [mus], asymmetry, Gau_rate [mus-1]
        '''
        return asymmetry*exp(-0.5*(x*Gau_rate)**2)

    def bs(self,x,asymmetry,rate,beta): 
        '''
        fit component for a stretched decay, 
        x [mus], asymmetry, rate [mus-1], beta (>0)
        '''
        return asymmetry*exp(-(x*rate)**beta)

    def ml(self,x,asymmetry,field,phase,Lor_rate): 
        '''
        fit component for a precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        return asymmetry*cos(2*pi*self.gamma_mu*field*x+phase*self.radeg)*exp(-x*Lor_rate)

    def mg(self,x,asymmetry,field,phase,Gau_rate): 
        '''
        fit component for a precessing muon with Gaussian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Gau_rate [mus-1]
        '''
        return asymmetry*cos(2*pi*self.gamma_mu*field*x+phase*self.radeg)*exp(-0.5*(x*Gau_rate)**2)

    def ms(self,x,asymmetry,field,phase,rate,beta): 
        '''
        fit component for a precessing muon with stretched decay, 
        x [mus], asymmetry, field [T], phase [degrees], rate [mus-1], beta (>0)
        '''
        return asymmetry*cos(2*pi*self.gamma_mu*field*x+phase*self.radeg)*exp(-(x*rate)**beta)

    def fm(self,x,asymmetry,dipfield,Lor_rate):
        '''
        fit component for FmuF (powder average)
        '''
        return asymmetry/6.0*( 3.+cos(2*pi*self.gamma_mu*dipfield*sqrt(3.)*x)+
               (1.-1./sqrt(3.))*cos(pi*self.gamma_mu*dipfield*(3.-sqrt(3.))*x)+
               (1.+1./sqrt(3.))*cos(pi*self.gamma_mu*dipfield*(3.+sqrt(3.))*x) )*exp(-x*Lor_rate)

    def jl(self,x,asymmetry,field,phase,Lor_rate): 
        '''
        fit component for a Bessel j0 precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        return asymmetry*j0(2*pi*self.gamma_mu*field*x+phase*self.radeg)*exp(-x*Lor_rate)

    def jg(self,x,asymmetry,field,phase,Gau_rate): 
        '''
        fit component for a Bessel j0 precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        return asymmetry*j0(2*pi*self.gamma_mu*field*x+phase*self.radeg)*exp(-0.5*(x*Gau_rate)**2)

    def _kg(self,x,w,Gau_delta):
        '''
        auxiliary component for a static Gaussian Kubo Toyabe in longitudinal field, 
        x [mus], w [mus-1], Gau_delta [mus-1]
        w = 2*pi*gamma_mu*L_field
        '''
        Dt = Gau_delta*x
        DDtt = Dt**2
        DD = Gau_delta**2
        sqr2 = sqrt(2)
        argf = w/(sqr2*Gau_delta)
        fdc = dawsn(argf)
        wx = w*x
        if (w!=0): # non-vanishing Longitudinal Field
            Aa = np.real(exp(-0.5*DDtt + 1j*wx)*dawsn(-argf - 1j*Dt/sqr2) )
            Aa[Aa == np.inf] = 0 # bi-empirical fix
            np.nan_to_num(Aa,copy=False) # empirical fix 
            A=sqr2*(Aa + fdc)
            f = 1. - 2.*DD/w**2*(1-exp(-.5*DDtt)*cos(wx)) + 2.*(Gau_delta/w)**3*A
        else:
            f = (1. + 2.*(1-DDtt)*exp(-.5*DDtt))/3.
        return f

    def _kgdyn(self,x,w,Gau_delta,jump_rate,*argv):
        ''' 
        auxiliary dynamization of Gaussian Kubo Toyabe 
        by G. Allodi 
        N: number of sampling points;
        dt: time interval per bin [i.e. time base is t = dt*(0:N-1)]
        w [mus-1], Gau_delta [mus-1], jump_rate [MHz] 
        (longitudinal field freq, dGaussian distribution, scattering frequency 
        % alphaN: [optional argument] weighting coefficient alpha times N. Default=10 
        '''
        alphaN = 10. if not argv else argv[0] # default is 10.
        dt = x[1]-x[0]
        N = x.shape[0] + int(np.ceil(x[0]/dt)) # for function to include t=0
        Npad = N * 2 # number of total time points, includes as many zeros
        t = dt*np.linspace(0.,Npad-1,Npad)
        expwei = exp(-(alphaN/(N*dt))*t)

        gg = self._kg(t,w,Gau_delta)*(t < dt*N)  #  padded_KT
        # gg = 1/3*(1 + 2*(1 - s^2*tt.^2).*exp(-(.5*s^2)*tt.^2)) % 

        ff = np.fft.fft(gg*expwei*exp(-jump_rate*t)) # fft(padded_KT*exp(-jump_rate*t))
        FF = exp(-jump_rate*dt)*ff/(1.-(1.-exp(-jump_rate*dt))*ff) # (1-jump_rate*dt*ff)  

        dkt = np.real(np.fft.ifft(FF))/expwei  # ifft
        dkt = dkt[0:N] # /dkt(1) 

        #if (nargout > 1),
        #   t = t[0:intN-1]
        return dkt
         
    def kg(self,x,asymmetry,L_field,Gau_delta,jump_rate):
        '''
        Gaussian Kubo Toyabe in longitudinal field, static or dynamic
        x [mus], asymmetry, L_field [T], Gau_delta [mus-1], jump_rate (MHz)
        '''
        N = x.shape[0]
        w = 2*pi*L_field*self.gamma_mu
        if jump_rate==0: # static 
           f = self._kg(x,w,Gau_delta) # normalized to 1.
        else :            # dynamic
           # P=[w Gau_delta];
 
           f = self._kgdyn(x,w,Gau_delta,jump_rate)
# function generated from t=0, shift result nshift=data(1,1)/dt bins backward
           dt = x[1]-x[0]
           nshift = x[0]/dt
           Ns = N + np.ceil(nshift)
           if Ns%2: # odd
               Np = Ns//2
               Nm = -Np
           else: # even
               Np = Ns//2-1
               Nm = -Ns//2
           n = np.hstack((np.linspace(0,Np,Np+1),np.linspace(Nm,-1.,-Nm)))
           f = np.fft.ifft(np.fft.fft(f)*exp(nshift*1j*2*pi*n/Ns)); # shift back
        # multiply by amplitude
        f = asymmetry*np.real(f[0:N])
        return f

    def _chisquare_(self,**kwargv):
        '''
        kwargs is a dictionary of dictionaries
        use e.g. for a blmg model as: 
          m=mumodel() 
          m._init_(x,y,e) 
          iminuit.migrad(m._chisquare_(**{'c1':{m.bl:[asy1,rate1]},'c2':{m.mg:[asy2,fld2,ph2,rate2]}})                    
        'c1','c2' are unique but dummy
        chisquare is normalized if self.e was provided, unnormalized otherwise
        '''  
        return sum(((self._add_(self.x,**kwargv)-self.y)/self.e)**2)


# old.available_components =({'name':'da','npar':1,'par':[{'name':'dalpha','stepbounds':[0.,0.,0.]}],
#                                     'help':r'linear $f$ correction: $\frac{2\alpha-p[0](1-f)}{2\alpha+p[0](1-f)}$'},
#                                    {'name':'bl','npar':2,'par':[{'name':'blAsym','stepbounds':[0.,0.,0.]},
#                                                                {'name':'blDelo','stepbounds':[0.,0.,0.]}],
#                                     'help':r'Lorenz decay: $p[0]\exp(-p[1]t)$'},
#                                    {'name':'bg','npar':2,'par':[{'name':'bgAsym','stepbounds':[0.,0.,0.]},
#                                                                {'name':'muSigm','stepbounds':[0.,0.,0.]}],
#                                     'help':r'Gauss decay: $p[0]\exp(-(p[1]t)^2)/2$'},
#                                    {'name':'ml','npar':4,'par':[{'name':'mlAsym','stepbounds':[0.,0.,0.]},
#                                                                {'name':'mlBGau','stepbounds':[0.,0.,0.]},
#                                                                {'name':'mlPhiD','stepbounds':[0.,0.,0.]},
#                                                                {'name':'mlDelo','stepbounds':[0.,0.,0.]}],
#                                     'help':r'Lorenz decay cosine: $p[0]\exp(-p[3]t)\cos(2\pi(\gamma_\mu p[1] + p[2]/180.))$'},
#                                    {'name':'mg','npar':3,'par':[{'name':'mgAsym','stepbounds':[0.,0.,0.]},
#                                                                {'name':'mgBGau','stepbounds':[0.,0.,0.]},
#                                                                {'name':'mgPhiD','stepbounds':[0.,0.,0.]},
#                                                                {'name':'mgSigm','stepbounds':[0.,0.,0.]}],
#                                     'help':r'Gauss decay cosine: $p[0]\exp(-(p[3]t)^2/2)\cos(2\pi(\gamma_\mu p[1] + p[2]/180.))$'}
#                                    )

