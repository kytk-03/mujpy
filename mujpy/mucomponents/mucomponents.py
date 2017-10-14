from numpy import cos, pi, exp, sqrt
import numpy as np
from scipy.special import dawsn,erf, j0

class mumodel(object):
    def __init__(self):
        ''' 
        defines few constants and _help_ dictionary
        '''
        self._radeg_ = pi/180.
        self._gamma_mu_ = 135.5
        self._help_ = {'bl':r'Lorentz decay: $\mbox{asymmetry}\exp(-\mbox{Lor_rate}t)$',
                     'bg':r'Gauss decay: $\mbox{asymmetry}\exp(-0.5(\mbox{Gau_rate}t)^2)$',
                     'bs':r'Gauss decay: $\mbox{asymmetry}\exp(-0.5(\mbox{rate}t)^\beta)$',
                     'da':r'Linearized dalpha correction: $f = \frac{2f_0(1+\alpha/\mbox{dalpha})-1}{1-f_0+2\alpha/dalpha}$',
                     'mg':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-0.5(\mbox{Gau_rate}t)^2)$',
                     'ml':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-\mbox{Lor_rate}t)$',
                     'ms':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-(\mbox{rate}t)^\beta)$',
                     'jg':r'Gauss Bessel: $\mbox{asymmetry} j_0[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-0.5(\mbox{Lor_rate}t)^2)$',
                     'jl':r'Lorentz Bessel: $\mbox{asymmetry}j_0[2\pi(\gamma_\mu \mbox{field} t +\mbox{phase}/360)]\exp(-0.5(\mbox{Lor_rate}t)^2)$',
                     'fm':r'FMuF: $\mbox{asymmetry}/6[3+\cos 2*\pi\gamma_\mu\mbox{dipfield}\sqrt{3} t + \
               (1-1/\sqrt{3})\cos \pi\gamma_\mu\mbox{dipfield}(3-\sqrt{3})t + \
               (1+1/\sqrt{3})\cos\pi\gamma_\mu\mbox{dipfield}(3+\sqrt{3})t ]\exp(-\mbox{Lor_rate}t)$', 
                     'kg':r'Gauss Kubo-Toyabe: static and dynamic, in zero or longitudinal field by G. Allodi [Phys Scr 89, 115201]'}

    def _load_data_(self,x,y,_int,_alpha,e=1):
        ''' 
        Must be called before activating _chisquare_
        x, y, e are numpy arrays
        e is always defined, 
        either provided by the caller 
        or default to np.ones(x.shape[0])
        _int is a compact model list
        _alpha is ditto
        '''
        if x.shape[0]!=y.shape[0]:
            raise ValueError('x, y have different lengths')
        else:
            self._x_ = x
            self._y_ = y
            self._alpha_ = _alpha
            self._int = _int
        if e==1:
            self._e_ = np.ones(x.shape[0])
        else:
            if e.shape[0]!=x.shape[0]:
                raise ValueError('x, e have different lengths')           
            else:
                self._e_ = e

    # ---- end generic __init__

    def _add_(self,x,*argv):
        '''
        used by self._chisquare_ 
        e.g. a blmg model with 
        argv will be a tuple of parameter values (val1,val2.val3,val4,val5,val6) at this iteration 
        _add_ reconstructs how to distribute these parameter values
        use to plot : 
          plt.errorbar(x,y,yerr=e)
          plt.plot(x,mumodel()._add_(x,val1,val2.val3,val4,val5,val6)                    
        '''  
        _da_flag_ = False
        f = np.zeros(x.shape[0])   
        p = argv
        k = -1
        for component,parkeys in self._int:   # this allows only for single run fits       
            p_comp = []
            for key in parkeys:
                if key == '~':
                    k += 1
                    p_comp.append(p[k])
                else:
                    p_comp.append(eval(key[1:]))
            if component == self.da:
                _da_Flag = True
                da = p_comp[0]
            else:
                f += component(x,*p_comp)
        if _da_Flag:
            dada = da/self._alpha_
            f = ((2.+dada)*f-dada)/((2.+dada)-dada*f) # linearized correction 
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

    def da(self,x,dalpha):
        '''
        fit component for linearized alpha correction
        x [mus], dalpha
        '''
        # the returned value will not be used, correction in _add_
        return x*dalpha

    def ml(self,x,asymmetry,field,phase,Lor_rate): 
        '''
        fit component for a precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        return asymmetry*cos(2*pi*self._gamma_mu_*field*x+phase*self._radeg_)*exp(-x*Lor_rate)

    def mg(self,x,asymmetry,field,phase,Gau_rate): 
        '''
        fit component for a precessing muon with Gaussian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Gau_rate [mus-1]
        '''
        return asymmetry*cos(2*pi*self._gamma_mu_*field*x+phase*self._radeg_)*exp(-0.5*(x*Gau_rate)**2)

    def ms(self,x,asymmetry,field,phase,rate,beta): 
        '''
        fit component for a precessing muon with stretched decay, 
        x [mus], asymmetry, field [T], phase [degrees], rate [mus-1], beta (>0)
        '''
        return asymmetry*cos(2*pi*self._gamma_mu_*field*x+phase*self._radeg_)*exp(-(x*rate)**beta)

    def fm(self,x,asymmetry,dipfield,Lor_rate):
        '''
        fit component for FmuF (powder average)
        '''
        return asymmetry/6.0*( 3.+cos(2*pi*self._gamma_mu_*dipfield*sqrt(3.)*x)+
               (1.-1./sqrt(3.))*cos(pi*self._gamma_mu_*dipfield*(3.-sqrt(3.))*x)+
               (1.+1./sqrt(3.))*cos(pi*self._gamma_mu_*dipfield*(3.+sqrt(3.))*x) )*exp(-x*Lor_rate)

    def jl(self,x,asymmetry,field,phase,Lor_rate): 
        '''
        fit component for a Bessel j0 precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        return asymmetry*j0(2*pi*self._gamma_mu_*field*x+phase*self._radeg_)*exp(-x*Lor_rate)

    def jg(self,x,asymmetry,field,phase,Gau_rate): 
        '''
        fit component for a Bessel j0 precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        return asymmetry*j0(2*pi*self._gamma_mu_*field*x+phase*self._radeg_)*exp(-0.5*(x*Gau_rate)**2)

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
        w = 2*pi*L_field*self._gamma_mu_
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

    def _chisquare_(self,*argv):
        '''
        argv
        use e.g. for a blmg model as: 
          m=mumodel() 
          m._init_(x,y,e) 
          iminuit.migrad(m._chisquare_,**{'asy1':val1,'error_asy1:val,'fix_asy1':False,'limits_asy1':(val,val),
                                          'rate1':val2,'error_rate1:val,'fix_rate1':False,'limits_rate1':(val,val),
                                          'asy2':val3,'error_asy2:val,'fix_asy2':False,'limits_asy2':(val,val),
                                          'field2':val4,'error_field2:val,'fix_field2':False,'limits_field2':(val,val),
                                          'phase2':val5,'error_phase2:val,'fix_phase2':False,'limits_phase2':(val,val),
                                          'rate2':val6,'error_rate2:val,'fix_rate2':False,'limits_rate2':(val,val)}
        argv will be tuple of parameter values at this iteration (val1,val2.val3,val4,val5,val6)
        _add_ reconstructs how to distribute these parameter values
        chisquare is normalized if self.e was provided, unnormalized otherwise
        '''  
        return sum(((self._add_(self._x_,*argv)-self._y_)/self._e_)**2)


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

