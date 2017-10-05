import numpy as np
# import musuite
# import muset
from mujpy.mucomponents.mufcn import mufcn as fcn
from scipy.constants import physical_constants as C
import matplotlib.pyplot as P

class MuFit(object):                        # defines the python class
    def __init__(self,instance_name=None):
        """
        usage: 
        myfit = MuFit()  # this is model initialization
        """
        import traceback
        if instance_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        self.instance_name = instance_name  # a trick so that self.instance_name is the local name of the instance in the calling module
                                            # used by error messages in asymmetry (multi run compatibility check)
        self.TauMu_mus = 2.1969811 # numbers are from Particle Data Group 2017
        self.TauPi_ns = 2.6033 # numbers are from Particle Data Group 2017
        self.gamma_Mu_MHzperT = 3.183345142*C['proton gyromag. ratio over 2 pi'][0]  # numbers are from Particle Data Group 2017
        self.gamma_e_MHzperT = C['electron gyromag. ratio over 2 pi'][0]
        self.components = [] # becomes a list of directories
        self.model = '' # no model yet

        self.asymmetry = [] 
        self.asymerror = []
        self.numberHisto = []
        self.histoLength = []
        self.fig_asym = []
        self.available_components() # now use self.available_components as a list
        self.models() # now use self.model as a list

    def addcomponent(self, name):
        '''
        myfit = MuFit()
        myfit.addcomponent('ml') # adds e.g. a mu precessing, lorentzian decay, component
        this method adds a component selected from the available_component tuple of directories
        with zeroed values, stepbounds from available_components, flags set to '~' and empty functions
        '''
        if name in self.models:
            npar = self.available_components[self.models.index(name)]['npar']
            pars = self.available_components[self.models.index(name)]['par'] # list of dictionaries for the parameters
            for k in range(npar):
                pars[k].update({'value':0.0})
                pars[k].update({'flag':'~'})
                pars[k].update({'function':''}) # adds these three keys to each par dict
            self.components.append({'name':name,'npar':npar,'par':pars})
            return True # OK code
        else:
            print ('Warning: '+name+' is not a known component. Not added.\n'+
                   'With myfit = mufit(), type myfit.help to see the available components')
            return False # error code

    def addmodel(self, name):
        '''
        myfit = MuFit()       
        myfit.addmodel('daml') # adds e.g. the two component 'da' 'ml' model
        this method adds a model of components selected from the available_component tuple of directories
        with zeroed values, stepbounds from available_components, flags set to '~' and empty functions
        '''
        components = self.checkvalidmodel(name)
        if components: # exoloits the fact that [] is False and ['da'] is true
            self.model = name
            for component in components:
                self.addcomponent(component)
            return True
        else:
            return False

    def asymmetry(self,run):
        """
        myfit.setup(mset,msuite) # after initialization calls (see MuFit.setup)
        myfit.addmodel('daml') # e.g. for an alpha calibration run
        myfit.asymmetry(run) # run is a murs2py instance read from a data file
        #     mset, msuite are a muset and musuite instance, respectively check
        checks if this is the first run by previous_binWidth_ns
        adds aymmetry end error from a run over maximum available bins
        bin partial selection by interval is done elsewhere
        returns 0 for ok and -1 for error
        COMPLETE: 
        NEW:
        depending on self.global True or False
        creates time  or appends asym[metry&error]
        """ 
        question = lambda q: input(q+'? (y/n)').lower().strip()[0] == "y" or question(q)
        if self.binwidth_ns != run.get_binWidth_ns(): # includes first run and previous run with different binWidth
            if not self.binwidth_ns: # first run
                self.numberHisto = run.get_numberHisto_int()
                self.histoLength = run.get_histoLength_bin() - self.nt0.max() - self.offset # max available bins on all histos
                self.firstrun  = True
                self.binwidth_ns = run.get_binWidth_ns() 
                self.time = (np.arange(self.histoLength) + self.offset +
                             np.mean( self.dt0 [np.append(self.group['forward'] ,self.group['backward'])] ) 
                             )*self.binwidth_ns/1000. # in microseconds
##################################################################################################
# Simple cases: 
# 1) Assume the prompt is entirely in bin nt0. (python convention, the bin index is 0,...,n,... 
# The content of bin nt0 will be the t=0 value for this case and dt0 = 0.
# The center of bin nt0 will correspond to time t = 0, time = (n-nt0 + mufit.offset + mufit.dt0)*mufit.binWidth_ns/1000.
# 2) Assume the prompt is equally distributed between n and n+1. Then nt0 = n and dt0 = 0.5, the same formula applies
# 3) Assume the prompt is 0.45 in n and 0.55 in n+1. Then nt0 = n+1 and dt0 = -0.45, the same formula applies.
##################################################################################################
            else:
                print('WARNING! You are mixing runs with different resolutions. Not allowed.')
                return -1 # error code
        elif self.numberHisto != run.get_numberHisto_int() or self.histoLength != run.get_histoLength_bin():
            print('Mismatch in number and/or length of histograms')
            print('Analysing this run with the previous might make no sense')
            ans = question('Proceed anyway')
            if not ans:
                print('To delete previous run(s) type {:}.reset()'.format(self.instance_name))
                return -1 # error code
        # add run
        yforw = np.zeros(self.time.shape[0]) # counts with background substraction
        cforw = np.zeros(self.time.shape[0]) # pure counts for Poisson errors
        ybackw = np.zeros(self.time.shape[0]) # counts with background substraction
        cbackw = np.zeros(self.time.shape[0]) # pure counts for Poisson errors

        for detector in self.group['forward']:
            n1, n2 = self.nt0[detector]+self.offset, self.histoLength
            histo = run.get_histo_array_int(detector)
            background = np.mean(histo[self.firstbin:self.lastbin[detector]])
            yforw += histo[n1:n2]-background
            cforw += histo[n1:n2]

        for detector in self.group['backward']:
            n1, n2 = self.nt0[detector]+self.offset, self.histoLength
            histo = run.get_histo_array_int(detector)
            background = np.mean(histo[self.firstbin:self.lastbin[detector]])
            ybackw += histo[n1:n2]-background
            cbackw += histo[n1:n2]

        yplus = yforw + self.alpha*ybackw
        x = np.exp(-self.time/self.TauMu_mus)
        enn0 = np.polyfit(x,yplus,1)
        enn0 = enn0[0] # initial rate per ns
        y = (yforw-self.alpha*yback)/enn0*np.exp(self.time/self.TauMu_mus)  # since self.time is an np.arange, this is a numpy array
        ey = np.sqrt(cforw + self.alpha**2*cbackw)*np.exp(self.time/self.TauMu_mus)/enn0 # idem
        ey[ey==0] = 1 # substitute zero with one in ey
        if self.firstrun:
            self.asymm = y # np.array
            self.asyme = ey # idem
            self.firstrun = False
            self.nrun = [run.get_runNumber_int()]
        else:
            self.asymm = np.row_stack(self.asymmetry, y) # columns are times, rows are successive runs
            self.asyme = np.row_stack(self.asymerror, ey)
            self.nrun.append(run.get_runNumber_int())  # this is a list
        return 0 # no error

    def available_components(self):
        '''
        definition of the recognignized model components (invoked by __init__)
        the tuple contains one dictionary per component, each formed by
        'name', the 2-character name of the component
        'npar', its number of parameters, hereafter p[0]...p[npar-1]
        'par', a dictionary for each parameter,  with 
                'name', its name, 
                'stepbounds', a numpy array with step, lower bound, upper bound ([0.,0.,0.] means unbounded, step as 0.01*value)
        'help', a latex string, r'$...$', to be output by matplotlib in mufit.help
        '''
        self.available_components =({'name':'da','npar':1,'par':[{'name':'dalpha','stepbounds':[0.,0.,0.]}],
                                     'help':r'linear $f$ correction: $\frac{2\alpha-p[0](1-f)}{2\alpha+p[0](1-f)}$'},
                                    {'name':'bl','npar':2,'par':[{'name':'blAsym','stepbounds':[0.,0.,0.]},
                                                                {'name':'blDelo','stepbounds':[0.,0.,0.]}],
                                     'help':r'Lorenz decay: $p[0]\exp(-p[1]t)$'},
                                    {'name':'bg','npar':2,'par':[{'name':'bgAsym','stepbounds':[0.,0.,0.]},
                                                                {'name':'muSigm','stepbounds':[0.,0.,0.]}],
                                     'help':r'Gauss decay: $p[0]\exp(-(p[1]t)^2)/2$'},
                                    {'name':'ml','npar':4,'par':[{'name':'mlAsym','stepbounds':[0.,0.,0.]},
                                                                {'name':'mlBGau','stepbounds':[0.,0.,0.]},
                                                                {'name':'mlPhiD','stepbounds':[0.,0.,0.]},
                                                                {'name':'mlDelo','stepbounds':[0.,0.,0.]}],
                                     'help':r'Lorenz decay cosine: $p[0]\exp(-p[3]t)\cos(2\pi(\gamma_\mu p[1] + p[2]/180.))$'},
                                    {'name':'mg','npar':3,'par':[{'name':'mgAsym','stepbounds':[0.,0.,0.]},
                                                                {'name':'mgBGau','stepbounds':[0.,0.,0.]},
                                                                {'name':'mgPhiD','stepbounds':[0.,0.,0.]},
                                                                {'name':'mgSigm','stepbounds':[0.,0.,0.]}],
                                     'help':r'Gauss decay cosine: $p[0]\exp(-(p[3]t)^2/2)\cos(2\pi(\gamma_\mu p[1] + p[2]/180.))$'}
                                    )

    def checkvalidmodel(self,name):
        '''
        checkvalidmodel(name) checks that name is a 
                              2*component string of valid component names, e.g.
                              'daml' or 'mgmgbl'
        '''
        components = [name[i:i+2] for i in range(0, len(name), 2)]
        for component in components:            
            if component not in self.models:
                print ('Warning: '+component+' is not a known component. Not added.\n'+
                       'With myfit = mufit(), type myfit.help to see the available components')
                return [] # error code
            else:
                return components

    def clear_asymmetry(self):
        '''
        resets to a condition where no data were read
        perhap useless
        '''
        print('are you sure? anyways, too late, sorry')
        # implement as hidden accordion to clear, select one, get a message warning msg in output saying click again if you really want to delete
        del self.numberHisto, self.histoLength, self.firstrun, self.binwidth_ns, self.time, self.asymmetry, self.asymerror
          


    def deletemodel(self):
        '''
        myfit = MuFit.()
        myfit.deletemodel() 
        this method resets components to an empty list        
        '''
        self.components=[]
        self.model = ''

    def help(self):
        '''
        myfit = MuFit()
        myfit.help()
        print help strings of available_components
        (by matplotlib if no other Latex support is available)
        '''
        from matplotlib.pyplot import text, axis, show
        n = len(self.available_components)
        for k in range(n):
            text(0,n-k,available_components[k]['help'],fontsize=14)

    def load(self,path):
        """
        usage:
        import dill as pickle
        MuFit.load(path_and_filename) # 
        """
        with open(path,'rb') as f:
            self.__dict__.update(pickle.load(f).__dict__)

    def models(self):
        '''
        myfit = MuFit()
        myfit.models() # extracts the available model names 
        '''
        self.models = [self.available_components[i]['name'] for i in range(len(self.available_components))]

    def plot_limits(self, plot_init = 0, plot_last = 10000. ):
         self.plt_ini = plot_init
         self.plt_lst = plot_last

    def plot(self):
        ''' 
        to be completed
        this version is for the plot of a single run fit
        future version will deals with multi run plots 
        '''
        # for refresh the plot see http://nbviewer.jupyter.org/gist/branning/c8e63ce81be0391260b1
        if not self.fig_asym:
            self.fig_asym, self.ax_asym = P.suplots(2,1,figsize=(8,6),gridspec_kw = {'height_ratios':[3, 1]})
        x = self.time[self.plt_ini:self.plt_lst]
        y = self.asymm[self.plt_ini:self.plt_lst]
        ey = self.asyme[self.plt_ini:self.plt_lst]
        residues = y-fcn(x,p) # here the question is how does fcn take p
        # thsi is for plotting the best fit function
        x3 = np.linspace(self.time[self.plt_ini],self.time[self.plt_lst],3*(self.lst-self.ini))
        f3 = fcn(x3,p) # same question
        if binning>1:
            x,y,ey = rebin(x,y,ey,f)
        if (self.plt_lst-self.plt_ini)/self.binning > 200:
            # plot only markers
            self.ax_asym[0].plot(x,y,'o',c='g',mew=0.5,ms=2.5,lw=0.5,mfc='w',mec='g')
        else:
            # errorbars
            self.ax_asym[0].errorbar(x,y,yerr=ey,
                                     fmt='o',c='g',mew=0.5,ms=2.5,lw=0.5,mfc='w',mec='g')
        self.ax_asym[0].plot(x3,f3,'r-')  # best fit curve
        self.ax_asym[1].plot(x,residues,'b-')  # residues always unbinned, like the fit

    def rebin(self,x,y,*args):
        '''
        use either 
        xb,yb = rebin(x,y)
        or
        xb,yb,eyb = rebin(x,y,ey) # the 3rd is an error
        '''
        m = np.floor(len(x)/self.binning)
        mn = m*self.binning
        xx = x[:mn]
        xx.reshape(self.binning,m)
        yy = y[:mn]
        yy.reshape(self.binning,m)
        xb = xx.sum(0)/self.binning
        yb = yy.sum(0)/self.binning
        if args is not None:
            ey = args[:mn]
            ey.reshape(self.binning,m)
            eb = np.sqrt((ey**2).sum(0))/self.binning;
            return xb,yb,eb
        else:
            return xb,yb

    def reset(self):
        self.time = []
        self.asymm = []
        self.asyme = []
        self.nrun = []
        self.binwidth_ns =[]
        self.numberHisto = []
        self.histoLength = []

#    def setup(self,mset,msuite):
#        """
#        migrated to gui
#        run = musr2py()
#        run.read('deltat_tdc_gps_0433.bin')  # e.g. 
#        myf it = MuFit()
#        mset = muset(run)
#        mset.fit()
#        msuite = musuite()
#        myfit.setup(mset,msuite) # setup is a musuite instance containing alpha, group
#        """
#        self.alpha = mset.alpha
#        self.group = mset.group
#        self.offset = mset.offset
#        self.interval = mset.interval
#        self.plt_ini = self.interval[0]
#        self.plt_lst = self.interval[1]
#        self.binning = mset.binning
#        self.nt0 = msuite.nt0
#        self.dt0 = msuite.dt0
#        self.firstbin = msuite.firstbin
#        self.lastbin = msuite.lastbin
#        self.binwidth_ns = []

    def save(self,path):
        """
        usage:
        import dill as pickle
        MuFit.save(path_and_filename) # 
        """
        with open(path,'wb') as f:
            pickle.dump(self, f) 
        # read back also with: 
        # import dill as pickle
        # with open('musr2py/del.pkl','r') as f:
        #     b = pickle.load(f)

    def __getstate__(self):
        '''
        this method removes asymmetry and asymerror from the attributes that are pickled by the save command
        all the others end up in the pickled file
        '''
        state = self.__dict__
        del state['asymmetry']
        del state['asymerror']
        return state
