# write tlog
# write plot (inherit plot form mufit)
# change edit into an output
# write fft
# write suite suite [1.0]
######################################################
# Gui tabs correspond to distinct gui methods with independes scopes and additional local methods
# gui attributes: 
#     entities that must be shared between tabs 
#     including variables passed to functions  outside mugui
###############################################################
class mugui(object):

##########################
# INIT
##########################
    def __init__(self):
        '''
        initiates an instance and a few attributes, 
        launches the gui.
        Use as follows
         from mugui import mugui as MG
         MuJPy = MG() # instance is MuJPy

        '''
        import numpy as np
        import os
        from scipy.constants import physical_constants as C
        from mujpy import __file__ as MuJPyName
        from IPython.display import display

# check where constants are needed and how to define them
        self.TauMu_mus = 2.1969811 # numbers are from Particle Data Group 2017
        self.TauPi_ns = 2.6033 # numbers are from Particle Data Group 2017
        self.gamma_Mu_MHzperT = 3.183345142*C['proton gyromag. ratio over 2 pi'][0]  # numbers are from Particle Data Group 2017
        self.gamma_e_MHzperT = C['electron gyromag. ratio over 2 pi'][0]
# end check constants

        self.interval = np.array([0,7800], dtype=int)
        self.offset0 = 7 # initial value
        self.offset = [] # this way the first run load calculates get_totals with self.offset0
        self.firstbin = 0
        self.second_plateau = 100
        self.peakheight = 100000.
        self.peakwidth = 1.   # broad guesses for default

        self.binwidth_ns = [] # this way the first call to asymmetry(_the_run_) initializes self.time
        self.grouping = {'forward':np.array([1]),'backward':np.array([0])} # normal dict
        self._the_run_ = []  # if self._the_run_: is False, to check whether the handle is created
        self.first_t0plot = True
        self._global_fit_ = False # 
        # mujpy paths
        self.__path__ = os.path.dirname(MuJPyName)
        self.__logopath__ = os.path.join(self.__path__,"logo")
        self.__startuppath__ = os.getcwd() # working directory, in which, in case, to find mujpy_setup.pkl 

        self.binning = 100
        #####################################
        # actually produces the gui interface
        #####################################

        self.gui()
        self.setup()
        self.suite()
        self.fit()

        self.output()
        self.plot()
        # self fft
        # self.tlog
        self.about()
        try:
            whereamI = get_ipython().__class__
            if not str(whereamI).find('erminal')+1:
                display(self.gui)
            else:
                print(str(wheremI)) 
        except:
                print('Python test script')

##########################
# ABOUT
##########################
    def about(self):
        '''
        about tab:
        a few infos (version and authors)
        '''
        from ipywidgets import Textarea, Layout

        _version = 'MuJPy          version '+'0.1' # increment while progressing
        _authors = '\n\n  Authors: Roberto De Renzi, Pietro Bonfà, '
        _blahblah = ('\n\n  A Python MuSR data analysis graphical interface.'+
                     '\n  Based on classes, designed for jupyter.'+
                     '\n  Released under the MIT licence')
        _pronounce = ('\n  See docs in ReadTheDocs'+
                      '\n  Pronounce it as mug + pie')
        _about_text = _version+_blahblah+_pronounce+_authors
        _about_area = Textarea(value=_about_text,
                                   placeholder='Info on MuJPy',
                                   layout=Layout(width='100%',height='170px'),
                                   disabled=True)
        # now collect the handles of the three horizontal frames to the main fit window (see tabs_contents for index)
        self.mainwindow.children[7].children = [_about_area] # add the list of widget handles as the third tab, fit
       
##########################
# ASYMMETRY
##########################
    def asymmetry(self,run):   # lastbin):
        """
        the first run of the session must define self.time
        adds aymmetry end error from a run over maximum available bins
        bin partial selection by interval is done elsewhere
        # returns 0 for ok and -1 for error RETHINK
        MAYBE:
        depending on self.global True or False set by suite
        creates time or not and creates  or appends asymm asyme
        for single run if binwidth_ns does not exist, creates times 
        for _global_fit_ asymmetry must be called in suite sequence, starting from a master
            prior to which binwidth_ns must be set = []
        """ 
        import numpy as np

        question = lambda q: input(q+'? (y/n)').lower().strip()[0] == "y" or question(q)

        if not self.binwidth_ns or self._global_fit_: # first run
            self.numberHisto = run.get_numberHisto_int()
            self.histoLength = run.get_histoLength_bin() - self.nt0.max() - self.offset.value # max available bins on all histos
            self.firstrun  = True
            self.binwidth_ns = run.get_binWidth_ns() 
            self.time = (np.arange(self.histoLength) + self.offset.value +
                         np.mean(self.dt0 [np.append(self.grouping['forward'],self.grouping['backward'])] ) 
                         )*self.binwidth_ns/1000. # in microseconds
    ##################################################################################################
    # Simple cases: 
    # 1) Assume the prompt is entirely in bin self.nt0. (python convention, the bin index is 0,...,n,... 
    # The content of bin self.nt0 will be the t=0 value for this case and self.dt0 = 0.
    # The center of bin self.nt0 will correspond to time t = 0, time = (n-self.nt0 + self.offset.value + self.dt0)*mufit.binWidth_ns/1000.
    # 2) Assume the prompt is equally distributed between n and n+1. Then self.nt0 = n and self.dt0 = 0.5, the same formula applies
    # 3) Assume the prompt is 0.45 in n and 0.55 in n+1. Then self.nt0 = n+1 and self.dt0 = -0.45, the same formula applies.
    ##################################################################################################
        else:
            if self.binwidth_ns != run.get_binWidth_ns(): # includes first run and previous run with different binWidth
                print('WARNING! You are mixing runs with different resolutions. Not allowed.')
                return -1 # error code
            elif (self.numberHisto != run.get_numberHisto_int() or 
                 self.histoLength+self.nt0.max()+self.offset.value != run.get_histoLength_bin()):
                print('Mismatch in number and/or length of histograms')
                print('Analysing this run with the previous might make no sense')
                ans = question('Proceed anyway')
                if not ans:
                    print('To delete previous run(s) type {:}.reset()'.format(self.instance_name))
                    return -1 # error code
        # calculate asymmetry in y and error in ey
        yforw = np.zeros(self.time.shape[0]) # counts with background substraction
        cforw = np.zeros(self.time.shape[0]) # pure counts for Poisson errors
        ybackw = np.zeros(self.time.shape[0]) # counts with background substraction
        cbackw = np.zeros(self.time.shape[0]) # pure counts for Poisson errors

        for detector in self.grouping['forward']:
            n1, n2 = self.nt0[detector]+self.offset.value, self.nt0[detector]+self.offset.value+self.histoLength
            histo = run.get_histo_array_int(detector)
            background = np.mean(histo[self.firstbin:self.lastbin])
            yforw += histo[n1:n2]-background
            cforw += histo[n1:n2]

        for detector in self.grouping['backward']:
            n1, n2 = self.nt0[detector]+self.offset.value, self.nt0[detector]+self.offset.value+self.histoLength
            histo = run.get_histo_array_int(detector)
            background = np.mean(histo[self.firstbin:self.lastbin])
            ybackw += histo[n1:n2]-background
            cbackw += histo[n1:n2]

        yplus = yforw + self.alpha.value*ybackw
        x = np.exp(-self.time/self.TauMu_mus)
        enn0 = np.polyfit(x,yplus,1)
        enn0 = enn0[0] # initial rate per ns
        y = (yforw-self.alpha.value*ybackw)/enn0*np.exp(self.time/self.TauMu_mus)  # since self.time is an np.arange, this is a numpy array
        ey = np.sqrt(cforw + self.alpha.value**2*cbackw)*np.exp(self.time/self.TauMu_mus)/enn0 # idem
        ey[ey==0] = 1 # substitute zero with one in ey
        if self._global_fit_: # the first call of the suite the master, resets binwidth_ns, hence self.firstrun=True 
            if self.firstrun:
                self.asymm = y # np.array
                self.asyme = ey # idem
                self.firstrun = False
                self.nrun = [run.get_runNumber_int()]
            else:
                self.asymm = np.row_stack(self.asymm, y) # columns are times, rows are successive runs
                self.asyme = np.row_stack(self.asyme, ey)
                self.nrun.append(run.get_runNumber_int())  # this is a list
        else:
            self.asymm = y # np.array
            self.asyme = ey # idem
            self.nrun = [run.get_runNumber_int()]
            
        # return 0 # no error

##########################
# GUI
##########################
    def gui(self):
        '''
        gui layout
        designs external frame, 
                logo and title header
                tab structure
        '''
        from ipywidgets import Image, Text, Layout, HBox, Output, VBox, Tab
        import os
                        
        file = open(os.path.join(self.__logopath__,"logo.png"), "rb")
        image = file.read()
        logo = Image(value=image,format='png',width=132,height=132)
        self.title = Text(description='run title', value='none yet',layout=Layout(width='70%'),disable=True)
        self._the_run_display = Text(description='run number',value='no run',layout=Layout(width='30%'),Disable=True)
        title_content = [self._the_run_display, self.title]
        titlerow = HBox(description='Title')
        titlerow.children = title_content
        counts = ['Total counts', 'Group counts','ns/bin'] # needs an HBox with three Text blocks
        self.totalcounts = Text(value='0',description='Total counts',layout=Layout(width='40'),disabled=True)
        self.groupcounts = Text(value='0',description='Group counts',layout=Layout(width='40%'),disabled=True)
        self.nsbin = Text(description='ns/bin',layout=Layout(width='20%'),disabled=True)
        secondrow = HBox(description='counts',layout=Layout(width='100%'))
        secondrow.children = [self.totalcounts, self.groupcounts, self.nsbin]
        # self._output = Output(layout=Layout(width='100%'))
        # thirdrow = HBox([self._output],layout=Layout(height='60px',width='100%',overflow_y='scroll',overflow_x='scroll')) # x works y does scroll
        titlewindow = VBox()
        titlewindow_content = [titlerow, secondrow] # ,thirdrow (moved to 4th tab)
        titlewindow.children = titlewindow_content
        titlelogowindow = HBox()
        titlelogowindow_content = [logo, titlewindow]
        titlelogowindow.children = titlelogowindow_content

        # main layout: tabs
        tabs_contents = ['setup', 'suite', 'fit', 'output', 'plot','fft','tlog','about']
        tabs = [VBox(description=name,layout=Layout(border='solid')) for name in tabs_contents]
        self.mainwindow = Tab(children = tabs,layout=Layout(width='99.8%')) # '99.6%' works

        self.mainwindow.selected_index=0 # to stipulate that the first display is on tab 0, setup
        for i in range(len(tabs_contents)):
            self.mainwindow.set_title(i, tabs_contents[i])
        #
        self.gui = VBox(description='whole',layout=Layout(width='100%'))
        self.gui.children = [titlelogowindow, self.mainwindow]
        
##########################
# FIT
##########################
    def fit(self, model_in = 'daml'): # self.fit(model_in = 'mgmgbl') produces a different layout
        '''
        fit tab of mugui
        used to set: self.alpha.value, self.offset.value, forw and backw groups
                     fit and plot ranges, model version                    
             to display: model name 
             to activate: fit, plot and update buttons
             to select and load model (load from folder missing)
             to select parameters value, fix, function, fft subtract check
        '''
# rethink mufit inside mugui
# the calculation is performed in independent class mucomponents
# True inheritance is not viable because the mugui and mufit are interdependent 
# methods "inherited" by mugui: 
#     __init__ share initial attributes (constants) done
#     _available_components_ automagical list of mucomponents 
#     _eval_ evil but needed, integrate with muvalid and safetry
#     addcomponent to _the_model_ (rethink, it had these attributes, is it needed)
#     asymmetry (avoid passing attributes, same data load into mucomponents)
#     checkvalidmodel
#     clear_asymmetry: includes reset
#     create_model: lay out _the_model_
#     delete_model: for a clean start
#     help  
#     load 
#     rebin
#     save cannot dill a class, rather a selected list of attributes, including:
#       firstbin, lastbin, offset, nt0, dt0, _the_model_, _the_run_, _the_suite_, fitargs, ...
#       purpose: load will provide a status ready for a new fit


# write on_fit_request for fit_button.on_click(on_fit_request)
# test
# check that self attributes are all needed (some might just be local variables)

        # no need to observe parvalue, since their value is a perfect storage point for the latest value
        # validity check before calling fit
        from ipywidgets import FloatText, Text, IntText, Layout, Button, HBox, \
                               Checkbox, VBox, Dropdown
        from mujpy.mucomponents.mucomponents import mumodel
        import numpy as np

        def _available_components_():
            from iminuit import describe
            '''
            Method, returns a template tuple of dictionaries (one per fit component):
            Each dictionary contains 'name' and 'pars', 
            the latter in turns is a list of dictionaries, one per parameter, 'name','error,'limits'
            ({'name':'bl','pars':[{'name':'asymmetry','error':0.01,'limits'[0,0]},
                                  {'name':'Lor_rate','error':0.01,'limits'[0,0]}}, 
             ...)
            retreived magically from the mucomponents class.
            '''
            _available_components = [] # is a list, mutable
            # generates a template of available components.
            for name in [module for module in dir(mumodel()) if module[0]!='_']: # magical extraction of component names
                pars = describe(mumodel.__dict__[name])[2:]            # the [2:] is because the first two arguments are self and x
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
            self.available_components = (_available_components) # these are the mucomponents method directories
                                                                # transformed in tuple, immutable
            self.component_names = [self.available_components[i]['name'] 
                                    for i in range(len(self.available_components))] # list of just mucomponents method names

        def _eval(string):
            '''
            yes, I know eval is evil, but
            mujpy users with jupyter have already full control
            of the machine, hence I do not care!
            '''
            try: 
                return eval(string)
            except Exception as e: 
                print(e)

        def addcomponent(name):
            '''
            myfit = MuFit()
            addcomponent('ml') # adds e.g. a mu precessing, lorentzian decay, component
            this method adds a component selected from self.available_component, tuple of directories
            with zeroed values, stepbounds from available_components, flags set to '~' and empty functions
            '''
            if name in self.component_names:
                npar = len(self.available_components[self.component_names.index(name)]['pars']) # number of pars
                pars = self.available_components[self.component_names.index(name)]['pars'] # list of dictionaries for the parameters
                # e.g {'name':'asymmetry','error',0.01,'limits',[0, 0]}
                for k in range(npar):
                    pars[k].update({'value':0.0})
                    pars[k].update({'flag':'~'})
                    pars[k].update({'function':''}) # adds these three keys to each pars dict
                    # they serve to collect values in mugui
                self.model_components.append({'name':name,'pars':pars})
                return True # OK code
            else:
                self.mainwindow.selected_index=3
                with self._output:
                    print ('\nWarning: '+name+' is not a known component. Not added.\n'+
                           'With myfit = mufit(), type myfit.help to see the available components')
                return False # error code

        def arrange(bin_range,target = 'fit'):
                '''
                arrange(bin_range) 
                [arrange(bin_range, target='plot')]
                writes bin_range into self.fit_range 
                [self.plot_range] as csv integers
                '''
                return exec(target+'_range = str(bin_range[0])+", "+str(bin_range[1])')

        def create_model(name):
            '''
            myfit = MuFit()       
            myfit.create_model('daml') # adds e.g. the two component 'da' 'ml' model
            this method adds a model of components selected from the available_component tuple of directories
            with zeroed values, stepbounds from available_components, flags set to '~' and empty functions
            '''
    # name 2_0_mlml_blbl for 2 global parameters (A0 R), 0 kocal parameters (B end T) and two models 
    # e.g. alpha fit with a WTF and a ZF run, with two muon fractions of amplitude A0*R and A0*(1-R) respectively
    # find the three underscores in name by
    # [i for i in range(len(name)) if name.startswith('_', i)]

            components = checkvalidmodel(name)
            if components: # exploits the fact that [] is False and ['da'] is true
                self.model = name
                self.model_components = [] # start from empty model
                for component in components:
                    addcomponent(component)
                return True
            else:
                return False

        def checkvalidmodel(name):
            '''
            checkvalidmodel(name) checks that name is a 
                                  2*component string of valid component names, e.g.
                                  'daml' or 'mgmgbl'
            '''
            components = [name[i:i+2] for i in range(0, len(name), 2)]
            for component in components:            
                if component not in self.component_names:
                    print ('Warning: '+component+' is not a known component. Not added.\n'+
                           'With myfit = mufit(), type myfit.help to see the available components')
                    return [] # error code
            return components

        def derange(target = 'fit'):
            '''
            derange() 
            [derange(target='plot')]
            reads self.fit_range 
            [self.plot_range] assuming two csv positive integers
            if no comma returns -1,-1 as errors
            else if values are not numbers, returns -1,0 or 0,-1  
            '''
            string = eval('self.'+target+'_range.value') 
            try:
                comma = string.find(',')
            except:
                self.mainwindow.selected_index=3
                with self._output:
                    print("comma separated values, please")
                return -1,-1 
            try:
                first = int(string[:comma])
            except:
                self.mainwindow.selected_index=3
                with self._output:
                    print("Could not read first bin")
                return -1,0
            try:
                last = int(string[comma+1:])
            except:
                self.mainwindow.selected_index=3
                with self._output:
                    print("Could not read last bin")
                return 0,-1
            return first,last 

            # _alpha = self.alpha.value # alpha is float
            # _offset = self.offset.value # offset, integer, from t0 to first good bin for fit    
            # _grouping = self.grouping # {'forward':np.array,'backward':np.array}
            # _fit_range = [self.derange()[k] for k in range(2)] # changes to fit_range (derange('target='plot') for plot_range)
            # _plot_range = [self.derange(target='plot')[k] for k in range(2)] # changes to plot_range (derange('target='plot') for plot_range)

        def get_grouping(name):
            """
            name = 'forward' or 'backward'
            grouping(name) is an np.array wth detector indices
            group.value[k] for k=0,1 is a shorthand csv like '1:3,5' or '1,3,5' etc.
            """
            import numpy as np

            # two shorthands: either a list, comma separated, such as 1,3,5,6 
            # or a pair of integers, separated by a colon, such as 1:3 = 1,2,3 
            # only one column is allowed, but 1, 3, 5 , 7:9 = 1, 3, 5, 7, 8, 9 
            # or 1:3,5,7 = 1,2,3,5,7  are also valid
            #       get the shorthand from the gui Text 
            groups = ['forward','backward']
            groupcsv = self.group[groups.index(name)].value # self.group[0,1] are handles to the correspondimg Texts
            #       now parse groupcsv shorthand
            groupcsv = groupcsv.replace('.',',') # can only be a mistake: '.' means ','
            if groupcsv.find(':')==-1:
                # colon not found, csv only
                self.grouping[name] = np.array([int(s) for s in groupcsv.split(',')])
            else:
                # colon found
                try:
                    if groupcsv.find(',')+groupcsv.find(':')==-2:
                        self.grouping[name] = np.array([int(groupcsv)])
                    elif groupcsv.find(',')+1: # True if found, False if not found (not found yields -1)    
                        firstcomma = groupcsv.index(',')
                        lastcomma = groupcsv.rindex(',')
                        if firstcomma < groupcsv.find(':'): # first read csv then range
                            partial = np.array([int(s) for s in groupcsv[:lastcomma].split(',')])
                            fst = int(groupcsv[lastcomma:grouping.find(':')])
                            lst = int(groupcsv[groupcsv.find(':')+1:])
                            self.grouping[name] = np.concatenate((partial,arange(fst,lst+1,dtype=int)))
                        else: # first read range then csv
                            partial = np.array([int(s) for s in groupcsv[:lastcomma].split(',')])
                            fst = int(groupcsv[:groupcsv.find(':')])
                            lst = int(groupcsv[groupcsv.find(':')+1:firstcomma])
                            self.grouping[name] = np.concatenate((np.arange(fst,lst+1,dtype=int),partial))
                    else: # only range
                        fst = int(groupcsv[:groupcsv.find(':')])
                        lst = int(groupcsv[groupcsv.find(':')+1:])
                        self.grouping[name] = np.arange(fst,lst+1,dtype=int)
                    self.grouping[name] -=1 # python starts at 0
                except:
                    print('Wrong group syntax: {}'.format(self.group[groups.index(name)].value))
                    self.group[groups.index(name)].value = ''
    
        def int2_int():
            '''
            From internal parameters to the minimal representation 
            for the use of mucomponents._add_.
            Invoked just before submitting minuit 
            '''
            ntot = sum([len(self.model_components[k]['pars']) for k in range(len(self.model_components))])
            lmin = [-1]*ntot
            nint = -1 # initialize
            nmin = -1 # initialize
            _int = []
            for k in range(len(self.model_components)):  # scan the model
                name = self.model_components[k]['name']
                # print('name = {}, model = {}'.format(name,_the_model_))
                bndmthd = _the_model_.__getattribute__(name)
                component_dict = {name:bndmthd} # the method
                keys = []
                for j in range(len(self.model_components[k]['pars'])): # 
                    nint += 1  # internal parameter incremente always   
                    if flag[nint].value == '=': #  function is written in terms of nint
                        # must be translated into nmin 
                        str = function[nint].value
                        indices = [int(s) for s in str.split() if s.isdigit()]
                        for l in indices:
                            str.replace(l,lmin[l])
                        keys.append(function[nint].value)
                    else:
                        keys.append('~')
                        nmin += 1
                        lmin[nmin] = nint # number of minuit parameter
                _int.append([component_dict,keys])
           # for k in range(len(_int)):
           #     print(_int[k])      

            return _int


        def int2min():
            '''
            From internal parameters to minuit parameters.
            Invoked just before submitting minuit 
            Internal are numbered progressively according to the display:
               first global parameters not belonging to components - e.g. A0, R, 
                        such as for asymmetry1 = A0*R and asymmetry2= A0*(1.-R)
               then local parameters not belonging to components - e.g. B and T
                        from the data file headers
               then the list of components' parameters
            Minuit parameters are the same, including fixed ones, but 
               the ones defined by functions or sharing
            Each parameter requires name=value, error_name=value, fix_name=value, limits_name=value,value
            [plus
               the local replica of the non global component parameters
               to be implemented]
            '''
            ntot = sum([len(self.model_components[k]['pars']) for k in range(len(self.model_components))])
            ntot -= sum([1 for k in range(ntot) if flag[k]=='=']) # ntot minus number of functions 
            lmin = [-1]*ntot
            nint = -1 # initialize
            nmin = -1 # initialize
            fitargs = {}
            minuit_parameter_names = []
            for k in range(len(self.model_components)):  # scan the model
                component_name = self.model_components[k]['name'] # name of component
                keys = []
                for j, par in enumerate(self.model_components[k]['pars']): # list of dictionaries, par is a dictionary
                    nint += 1  # internal parameter incremented always   
                    if flag[nint].value != '=': #  skip functions, they are not new minuit parameter
                        keys.append('~')
                        nmin += 1
                        lmin[nmin] = nint # correspondence between nmin and nint, is it useful?
                    fitargs.update({par['name']:float(parvalue[nint].value)})
                    minuit_parameter_names.append(par['name'])
                    fitargs.update({'error_'+par['name']:float(par['error'])})
                    if flag[nint].value == '!':
                        fitargs.update({'fix_'+par['name']:True})
                    if not (par['limits'][0] == 0 and par['limits'][1] == 0):
                        fitargs.update({'limit_'+par['name']:par['limits']})

            # print('fitargs= {}'.format(fitargs))
            return fitargs, tuple(minuit_parameter_names)

        def path_file_dialog(path):
            import tkinter
            from tkinter import filedialog
            import os
            here = os.getcwd()
            oc.chdir(path)
            tkinter.Tk().withdraw() # Close the root window
            in_path = filedialog.askopenfilename()
            os.chdir(here)
            return in_path

        def load_fit():
            '''
            loads fit values such that the same fit can be reproduced on the same data
            '''
            import dill as pickle
            import os

            path_and_filename = path_file_dialog(self.paths[2].value)
            with open(path_and_filename,'rb') as f:
                fit_dict = pickle.load(f) 
            self.version.value = fit_dict['version']
            self.offset.value = fit_dict['self.offset.value']
            self.model_components = fict_dict['self.model_components']
            self.grouping = fit_dict['self.grouping']
            self.alpha.value = fit_dict['self.alpha.value']
            self.alpha.value = fit_dict['self.alpha.value']
            

        def min2int(fitargs):
            '''
            From minuit parameters to internal parameters,
            see int2min for a description   
            Invoked just after minuit convergence.
            '''
            nint = -1
            ntot = sum([len(self.model_components[k]['pars']) for k in range(len(self.model_components))])
            _parvalue =  []
            for k in range(len(self.model_components)):  # scan the model
                for j, par in enumerate(self.model_components[k]['pars']): # list of dictionaries, par is a dictionary
                    nint += 1  # internal parameter incremented always   
                    if flag[nint].value != '=': #  skip functions, they are not new minuit parameter
                        _parvalue.append(str(fitargs[par['name']]))
            return _parvalue
            
        def muvalid(string):
            '''
            parse function CHECK WITH MUCOMPONENT, THAT USES A DIFFERENT SCHEME
            accepted functions are RHS of agebraic expressions of parameters p[i], i=0...ntot  
            '''
            import re
            from mujpy.aux.safetry import safetry

            pattern = re.compile(r"\p\[(\d+)\]") # find all patterns p[*] where * is digits
            test = pattern.sub(r"a",string) # substitute "a" to "p[*]" in s
            #           strindices = pattern.findall(string)
            #           indices = [int(strindices[k]) for k in range(len(strindices))] # in internal parameter list
            #           mindices = ... # produce the equivalent minuit indices  
            valid = True
            message = ''
            try: 
                safetry(test) # should select only safe use (although such a thing does not exist!)
            except Exception as e:
                print('function: {}. Tested: {}. Wrong or not allowed syntax: {}'.format(string,test,e))
                valid = False
            return valid

        def on_fit_request(b):
            '''
            retrieve data from the gui:
            parameters values (parvalue[nint].value), flags (flag[nint].value), 
            errors, limits, functions (function[nint].value), self.alpha.value, self.fit_range.value
            construct _int, needed by mumodel._add_ to distribute minuit parameter
            costruct fitargs dictionary, needed by migrad (and provided by it at the end)
            pass them to minuit
            mumodel._load_data_
            call minuit(..., *fitargs)
            save values
            save 
            print summary
            '''
            from iminuit import Minuit as M
            import matplotlib.pyplot as P

            self.asymmetry(self._the_run_) # prepare asymmetry
            fitargs, minuit_parameter_names = int2min()
            # print('fitargs = {}\n minuit_parameter_names = {}'.format(fitargs,minuit_parameter_names))
            _the_model_._load_data_(self.time,self.asymm,
                                    int2_int(),self.alpha.value,
                                    e=self.asyme) # pass data to model
            level = 1
            # save values save fit print summary
            # self.mainwindow.selected_index=3
            with self._output:
                m = M(_the_model_._chisquare_,pedantic=False,forced_parameters=minuit_parameter_names,print_level=level,**fitargs)
                m.migrad()
            pars = [m.fitarg[name] for name in minuit_parameter_names]
            t,f = self.rebin(self.time,_the_model_._add_(self.time,*pars))
            t,y,ey = self.rebin(self.time,self.asymm,e=self.asyme)
            self.mainwindow.selected_index=4
            with self._plot:
                fig,ax = P.subplots(2,1,figsize=(6,4),sharex = True, gridspec_kw = {'height_ratios':[2, 1]})
                ax[0].errorbar(t,y,yerr=ey,fmt='ro',elinewidth=1.0,ms=2.0)
                ax[0].plot(self.time,_the_model_._add_(self.time,*pars),'g-',lw=1.0)
                ax[1].plot(t,y-f,'r-',lw=1.0)
                ax[0].set_ylabel('Asymmetry')
                ax[1].set_ylabel('Residues')
                ax[1].set_xlabel(r'Time [$\mu$s]')
                ax[0].set_title(self.title.value)
                P.subplots_adjust(wspace=0.07)
                P.show()
            save_fit()
         
        def on_flag_changed(change):
            '''
            observe response of fit tab widgets:
            set disabled on corresponding function (True if flag=='!' or '~', False if flag=='=') 
            '''
            dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            n = int(dscr[4:]) # description='flag'+str(nint), skip 'flag'
            function[n].disabled=False if change['new']=='=' else True

        def on_alpha_changed(change):
            '''
            observe response of fit tab widgets:
            validate float        
            '''
            string = change['owner'].value # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            try: 
                float(string)
            except:
                change['owner'].value = '{:.4f}'.format(alpha0)

        def on_function_changed(change):
            '''
            observe response of fit tab widgets:
            check for validity of function syntax
            '''
            dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            n = int(dscr[4:]) # description='func'+str(nint), skip 'func'
            if not muvalid(change['new']):
                function[n].value = ''     
  
        def on_group_changed(change):
            '''
            observe response of setup tab widgets:
            '''
            name = change['owner'].description
            get_grouping(name) # stores self.group shorthand in self.grouping dict

        def on_integer(change):
            name = change['owner'].description
            if name == 'offset':
                if self.offset.value<0: # must be positive
                   self.offset.value = self.offset0 # standard value

        def on_loadmodel_changed(change):
            '''
            observe response of fit tab widgets:
            check that change['new'] is a valid model
            relaunch MuJPy.fit(change['new'])
            '''
            if checkvalidmodel(change['new']): # empty list is False, non empty list is True
                self.fit(change['new']) # restart the gui with a new model
                self.mainwindow.selected_index=2
                self.gui
                del _the_model_
            else:
                loadmodel.value=''

        def on_parvalue_changed(change):
            '''
            observe response of fit tab widgets:
            check for validity of function syntax
            '''
            dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            n = int(dscr[5:]) # description='value'+str(nint), skip 'func'
            try:
                float(parvalue[n].value)
            except:
                parvalue[n].value = '' 
  
        def on_range(change):
            '''
            observe response of fit range widgets:
            check for validity of function syntax
            '''
            fit_or_plot = change['owner'].description[0] # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            name='fit' if fit_or_plot=='f' else 'plot'
            if sum(self.derange(target=name))<0: # errors return (-1,-1),(-1,0),(0,-1) 
               exec('self.'+name+'.value = ""') # reset to expty text

        def save_fit():
            '''
            saves fit values such that load_fit can reproduce the same fit
            '''
            import dill as pickle
            import os

            fitargs, minuit_parameter_names = int2min()
            _int = int2_int()
            version = self.version.value
            strgrp = self.group[0].value.replace(',','_')+'-'+self.group[1].value.replace(',','_')
            path = os.path.join(self.paths[2].value, model.value+'.'+str(self.nrun[0])+'.'+strgrp+'.fit')
            # create dictionary setup_dict to be pickled 
            names = ['self.alpha.value','self.offset.value',
                     'self.grouping','model.value',
                     'fitargs','self.model_components',
                     'version','nint'] # keys
            fit_dict = {}
            for k,key in enumerate(names):
               fit_dict[names[k]] = eval(key) # key:value
            _parvalue = min2int(fitargs)
            for k in range(len(self.model_components)):
               for k in range(nint):
                   fit_dict['_parvalue['+str(k)+']'] = _parvalue[k] # fit result
                   fit_dict['_flag['+str(k)+    ']'] = flag[k].value # from fit tab
                   fit_dict['_function['+str(k)+']'] = function[k].value # from fit tab
            with open(path,'wb') as f:
                pickle.dump(fit_dict, f) 
            self.mainwindow.selected_index=3
            with self._output:
                print('Saved {}'.format(path))

        def set_group():
            """
            return shorthand csv out of grouping
            name = 'forward' or 'backward'
            grouping[name] is an np.array wth detector indices
            group.value[k] for k=0,1 is a shorthand csv like '1:3,5' or '1,3,5' etc.
            """
            import numpy as np

            # two shorthands: either a list, comma separated, such as 1,3,5,6 
            # or a pair of integers, separated by a colon, such as 1:3 = 1,2,3 
            # only one column is allowed, but 1, 3, 5 , 7:9 = 1, 3, 5, 7, 8, 9 
            # or 1:3,5,7 = 1,2,3,5,7  are also valid
            #       get the shorthand from the gui Text 
            groups = ['forward','backward']
            for k, name in enumerate(groups):
                s = ''
                aux = np.split(self.grouping[name],np.where(np.diff(self.grouping[name]) != 1)[0]+1)
                for j in aux:                    
                    s += str(j[0]+1) # convention is from 1 python is from 
                    if len(j)>1:
                        s += ':'+str(j[-1]+1)
                    s += ','
                s = s[:-1]
                self.group[k].value = s

        ######### here starts the fit method of MuGui
        _available_components_() # creates tuple self.available_components automagically from mucomponents
        _the_model_ = mumodel() # local instance, need a new one each time a fit tab is reloaded (on_load_model)
        try:
            alpha0 = self.alpha.value
        except:
            alpha0 = 1.01 # generic initial value
        try:
            self.offset0 = self.offset.value
        except:
            self.offset0 = 7 # generic initial value 
        self.alpha = FloatText(description='alpha',value='{:.4f}'.format(alpha0),
                                layout=Layout(width='20%'),continuous_update=False) # self.alpha.value
        self.alpha.observe(on_alpha_changed,'value')
        self.offset = IntText(description='offset',value=self.offset0,
                                layout=Layout(width='20%'),continuous_update=False) # offset, is an integer
        # initialized to 7, only input is from an IntText, integer value, or saved and reloaded from mujpy_setup.pkl
        self.alpha.style.description_width='40%' 
        self.offset.style.description_width='50%' 
        self.offset.observe(on_integer,'value') # check validity, must be positive

        # group and grouping: csv shorthand 
        self.group = [Text(description='forward',layout=Layout(width='27%'),
                                    continuous_update=False),
                      Text(description='backward',layout=Layout(width='27%'),
                                    continuous_update=False)]
        set_group() # inserts shorthand from self.grouping into seld.group[k].value, k=0,1
        self.group[0].observe(on_group_changed,'value')
        self.group[1].observe(on_group_changed,'value')
        # end moved

        model = Text(description = '', layout=Layout(width='10%'), disabled = True) # this is static, empty description, next to loadmodel
        model.value = model_in
        loadmodel = Text(description='loadmodel',layout=Layout(width='19%'),continuous_update=False) # this is where one can input a new model name
        loadmodel.observe(on_loadmodel_changed,'value')
        loadmodel.style.description_width='35%'
        self.version = IntText(description='version',layout=Layout(width='11%',indent=False)) # version.value is an int
        self.version.style.description_width='43%'
        try:
            self.version.value = version
        except:
            self.version.value = 1
        fit_button = Button (description='Fit',layout=Layout(width='10%'))
        fit_button.style.button_color = 'lightgreen'
        fit_button.on_click(on_fit_request)
        self.fit_range = Text(description='ft range',value='0,10000',layout=Layout(width='15%'),continuous_update=False)
        self.fit_range.style.description_width='35%'
        self.fit_range.observe(on_range,'value')
        plot_button = Button (description='Plot',layout=Layout(width='10%'))
        plot_button.style.button_color = 'lightgreen'
        self.plot_range = Text(description='pt range',value='0,500',layout=Layout(width='15%'),continuous_update=False)
        self.plot_range.style.description_width='35%'
        self.plot_range.observe(on_range,'value')
        update_button = Button (description='Update',layout=Layout(width='10%'))
        update_button.style.button_color = 'lightgreen'
        
        topframe_handle = HBox(description = 'Model', children=[model, 
                                                                loadmodel,
                                                                self.version,
                                                                fit_button,
                                                                self.fit_range, 
                                                                plot_button, 
                                                                self.plot_range,
                                                                update_button])  #
        alphaframe_handle = HBox(description = 'Alpha', children=[self.alpha,
                                                                   self.offset,
                                                                   self.group[0],
                                                                   self.group[1]])  # 

        bottomframe_handle = HBox(description = 'Components', layout=Layout(width='100%',border='solid')) #

        try: 
            create_model(model.value) # this may be not a valid model, e.g. after fit('da#ò')
        except:
            self.fit()  # this starts over, producing model = 'daml', which is valid.

        leftframe_list, rightframe_list = [],[]
  
        words  = ['#','name','value','~!=','function']
        nint = -1 # internal parameter count, each widget its unique name
        ntot = np.array([len(self.model_components[k]['pars']) 
                         for k in range(len(self.model_components))]).sum()
        parvalue, flag, function = [], [], [] # lists, index runs according to internal parameter count nint
        self.compar = {} # dictionary: key nint corresponds to a list of two values, c (int index of component) and p (int index of parameter)
                # use: self.compar[nint] is a list of two integers, the component index k and its parameter index j

        for k in range(len(self.model_components)):  # scan the model            
            header = HBox([ Text(value=self.model_components[k]['name'],disabled=True,layout=Layout(width='8%')),
                           Checkbox(description='FFT',value=True) ]) # list of HBoxes, the first is the header for the component
                                                                     # composed of the name (e.g. 'da') and the FFT flag
                                                                     # fft will be applied to a 'residue' where only checked components
                                                                     # are subtracted
            componentframe_list = [header] # list of HBoxes, header and pars
            componentframe_handle = VBox()
            for j in range(len(self.model_components[k]['pars'])): # make a new par for each parameter 
                                                                                # and append it to component_frame_content
                nint += 1      # all parameters are internal parameters, first is pythonically zero 
                self.compar.update({nint:[k,j]}) # stores the correspondence between nint and component,parameter
                nintlabel_handle = Text(value=str(nint),layout=Layout(width='10%'),disabled=True)
                parname_handle = Text(value=self.model_components[k]['pars'][j]['name'],layout=Layout(width='15%'),disabled=True)
                # parname can be overwritten, not important to store

                parvalue.append(Text(value='{:.4}'.format(self.model_components[k]['pars'][j]['value']),
                                     layout=Layout(width='20%'),description='value'+str(nint),continuous_update=False))
                parvalue[nint].style.description_width='0%'
                try:
                    parvalue[nint].value = _parvalue[nint]
                except:
                    pass
                # parvalue handle must be unique and stored at position nint, it will provide the initial guess for the fit

                function.append(Text(value=self.model_components[k]['pars'][j]['function'],
                                     layout=Layout(width='38%'),description='func'+str(nint),continuous_update=False))
                function[nint].style.description_width='0%'
                try:
                    function[nint].value = _function[nint]
                except:
                    pass
                # function handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation 

                fdis = False if self.model_components[k]['pars'][j]['flag']=='=' else True 
                function[nint].disabled = fdis # enabled only if flag='='
                flag.append(Dropdown(options=['~','!','='], 
                                     value=self.model_components[k]['pars'][j]['flag'],
                                     layout=Layout(width='10%'),description='flag'+str(nint)))
                flag[nint].style.description_width='0%'
                try:
                    flag[nint].value = _flag[nint]
                except:
                    pass
                 # flag handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation to be evaluated
 
                # now put this set of parameter widgets for the new parameter inside an HBox
                par_handle = HBox([nintlabel_handle, parname_handle, parvalue[nint], flag[nint], function[nint]])
                           # handle to an HBox of a list of handles; notice that parvalue, flag and function are lists of handles
                
                # now make value flag and function active 
                parvalue[nint].observe(on_parvalue_changed,'value')
                flag[nint].observe(on_flag_changed,'value') # when flag[nint] is modified, function[nint] is z(de)activated
                function[nint].observe(on_function_changed,'value') # when function[nint] is modified, it is validated

                componentframe_list.append(par_handle) # add par widget to the frame list

            componentframe_handle.children = componentframe_list # add full component to the frame
            if k%2==0:                                         # and ...
                leftframe_list.append(componentframe_handle) # append it to the left if k even
            else:
                rightframe_list.append(componentframe_handle) # or to the right if k odd   

        # end of model scan, ad two vertical component boxes to the bottom frame
        bottomframe_handle.children = [VBox(leftframe_list),VBox(rightframe_list)]  # list of handles  

        # now collect the handles of the three horizontal frames to the main fit window 
        self.mainwindow.children[2].children = [alphaframe_handle, topframe_handle ,bottomframe_handle] 
                             # add the list of widget handles as the third tab, fit

##########################
# OUTPUT
##########################
    def output(self):
        '''
        create an Output widget in fourth tab       
        select by 
        self.mainwindow.selected_index=3
        '''
        from ipywidgets import Output, HBox, Layout 
                     # Output(layout={'height': '100px', 'overflow_y': 'auto', 'overflow_x': 'auto'})
        self._output = Output(layout={'height': '200px','width':'100%','overflow_y':'auto','overflow_x':'auto'})
        _output_box = HBox([self._output],layout=Layout(width='100%')) # x works y does scroll
        self.mainwindow.children[3].children = [_output_box] 
                             # add the list of widget handles as the fourth tab, output

##########################
# PLOT
##########################
    def plot(self):
        '''
        create an Output in fifth tab for plots     
        select by 
        self.mainwindow.selected_index=4
        '''
        from ipywidgets import Output, HBox, Layout 
                     # Output(layout={'height': '100px', 'overflow_y': 'auto', 'overflow_x': 'auto'})
        self._plot = Output(layout={'width':'100%'}) #,'height':'200px','overflow_y':'auto','overflow_x':'auto'})
        _plot_box = HBox([self._plot],layout=Layout(width='100%')) # x works y does scroll
        self.mainwindow.children[4].children = [_plot_box] 
                             # add the list of widget handles as the fourth tab, output

##########################
# REBIN
##########################
    def rebin(self,x,y,e=None):
        '''
        use either 
        xb,yb = rebin(x,y)
        or
        xb,yb,eyb = rebin(x,y,ey) # the 3rd is an error
        '''
        from numpy import floor, sqrt
        m = int(floor(len(x)/self.binning))
        mn = m*self.binning
        xx = x[:mn]
        xx = xx.reshape(m,self.binning)
        yy = y[:mn]
        yy = yy.reshape(m,self.binning)
        xb = xx.sum(1)/self.binning
        yb = yy.sum(1)/self.binning
        if e is not None:
            ey = e[:mn]
            ey = ey.reshape(m,self.binning)
            eb = sqrt((ey**2).sum(1))/self.binning
            return xb,yb,eb
        else:
            return xb,yb

##########################i
# SETUP
##########################
    def setup(self):
        '''
        setup tab of mugui
        used to set: paths, fileprefix and extension
                     prepeak, postpeak (for prompt peak fit)
                     prompt plot check, 
             to activate: fit, save and load setup buttons
        '''

        def load_setup(b):
            """
            when user presse this setup tab widget:
            loads mujpy_setup.pkl with saved attributes
            and replaces them in setup tab Text widgets
            """
            import dill as pickle
            import os
                            
            path = os.path.join(self.__startuppath__,'mujpy_setup.pkl')
            # print('loading {}, presently in {}'.format(path,os.getcwd()))
            try:
                with open(path,'rb') as f:
                    mujpy_setup = pickle.load(f) 
            except:
                print('File {} not found'.format(path))
            # _paths_content = [ self.paths[k].value for k in range(3) ] # should be 3 ('data','tlag','analysis')
            # _filespecs_content = [ self.filespecs[k].value for k in range(2) ] # should be 2 ('fileprefix','extension')
            # _prepostpk = [self.prepostpk[k].value for k in range(2)] # 'pre-prompt bin','post-prompt bin' len(bkg_content)
            # _nt0 = self.nt0 # numpy array
            # _dt0 = self.dt0 # numpy array
            try:
                for k in range(3):  # len(paths_contents)
                    self.paths[k].value =  mujpy_setup['_paths_content'][k] # should be 3 ('data','tlag','analysis')
                for k in range(2):  # len(filespecs.content)
                    self.filespecs[k].value = mujpy_setup['_filespecs_content'][k] # should be 2 ('fileprefix','extension')
                for k in range(2):  # len(bkg_content)
                    self.prepostpk[k].value = mujpy_setup['_prepostpk'][k] # 'pre-prompt bin','post-prompt bin' 
                self.nt0 = mujpy_setup['self.nt0'] # bin of peak, nd.array of shape run.get_numberHisto_int()
                self.dt0 = mujpy_setup['self.dt0'] # fraction of bin, nd.array of shape run.get_numberHisto_int()
                self.lastbin = mujpy_setup['self.lastbin'] # fraction of bin, nd.array of shape run.get_numberHisto_int()
            except Exception as e:
                print('Error in load_setup: {}'.format(e))

        def save_setup(b):
            """
            when user presses this setup tab widget:
            saves mujpy_setup.pkl with setup tab values
            """
            import dill as pickle
            import os

            path = os.path.join(self.__startuppath__, 'mujpy_setup.pkl')
            # create dictionary setup_dict to be pickled 
            _paths_content = [ self.paths[k].value for k in range(3) ] # should be 3 ('data','tlag','analysis')
            _filespecs_content = [ self.filespecs[k].value for k in range(2) ] # should be 2 ('fileprefix','extension')
            _prepostpk = [self.prepostpk[k].value for k in range(2)] # 'pre-prompt bin','post-prompt bin' len(bkg_content)
            names = ['_paths_content','_filespecs_content',
                      '_prepostpk','self.nt0','self.dt0','self.lastbin'] # keys
            setup_dict = {}
            for k,key in enumerate(names):
               setup_dict[names[k]] = eval(key) # key:value
            with open(path,'wb') as f:
                pickle.dump(setup_dict, f) # according to __getstate__()
            self.mainwindow.selected_index=3
            with self._output:
                print('Saved {}'.format(os.path.join(self.__startuppath__,'mujpy_setup.pkl')))

        def on_paths_changed(change):
            '''
            when user changes this setup tab widget:
            check that paths exist, in case creates analysis path
            '''
            import os

            path = change['owner'].description # description is paths[k] for k in range(len(paths)) () 
            k = paths_content.index(path) # paths_content.index(path) is 0,1,2 for paths_content = 'data','tlog','analysis'
            directory = self.paths[k].value # self.paths[k] = handles of the corresponding Text
            if not os.path.isdir(directory):
                print('Path {} is unreachable'.format(directory))  
                
                if k==2: # analysis, if it does not exixt mkdir
                    # eventualmente togli ultimo os.path.sep = '/' in directory
                    dire=directory
                    if dire.rindex(os.path.sep)==len(dire):
                        dire=dire[:-1]
                    # splitta all'ultimo os.path.sep = '/'
                    prepath=dire[:dire.rindex(os.path.sep)+1]
                    # controlla che prepath esista
                    # print('prepath for try = {}'.format(prepath))
                    try:
                        os.stat(prepath)
                        os.mkdir(dire+os.path.sep)
                        print ('Analysis path {} created'.format(directory))
                        # self.paths[k].value = dire+os.path.sep # not needed if path is made with os.path.join 
                    except:
                        self.paths[k].value = os.path.curdir
                else:
                    self.paths[k].value = os.path.curdir
 
        def on_prompt_fit_click(b):
            '''
            when user presses this setup tab widget:
            execute prompt fits
            '''
            promptfit(mplot=self.plot_check.value) # mprint we leave always False

        def promptfit(mplot = False, mprint = False):
            '''
            when user presses this setup tab widget:
            launches t0 prompts fit 
            fits peak positions 
            prints migrad results
            plots prompts and their fit (default no print, no plot)
            stores bins for background and t0        
            '''
            import numpy as np
            from iminuit import Minuit, describe
            import matplotlib.pyplot as P
            from mujpy.mucomponents.muprompt import muprompt

            font = {'family' : 'Ubuntu','size'   : 6}
            P.rc('font', **font)

            if not self._the_run_:
                self.mainwindow.selected_index=3
                with self._output:
                     print('No run loaded yet! Load one first (select suite tab).')
                return    
            else:
                npeaks = np.array([np.where(self._the_run_.get_histo_array_int(det) ==
                                   self._the_run_.get_histo_array_int(det).max())[0][0] 
                                   for det in range(self._the_run_.get_numberHisto_int())])
                # approximate positions of peaks
                nbin =  max(npeaks) + self.second_plateau # this sets a detector dependent second plateau bin interval
                x = np.arange(0,nbin,dtype=int) # nbin bins from 0 to nbin-1
                self.lastbin, np3s = npeaks.min() - self.prepostpk[0].value, npeaks.max() + self.prepostpk[1].value # final bin of first and 
                                                                                         # initial bin of second plateaus

                if mplot and self.first_t0plot:
                    with self.t0plot_container:
                        self.figt0,self.axt0 = P.subplots(2,3,figsize=(7.5,5),
                                                          num='Prompts fit') 
    
                x0 = np.zeros(self._the_run_.get_numberHisto_int())
                if self.first_t0plot:
                    self.prompt_fit_text = [None]*self._the_run_.get_numberHisto_int()
                    # print(describe(muprompt))
                for detector in range(self._the_run_.get_numberHisto_int()):
                    # prepare for muprompt fit
                    histo = self._the_run_.get_histo_array_int(detector)
                    p = [ self.peakheight, float(npeaks[detector]), self.peakwidth, 
                          np.mean(histo[self.firstbin:self.lastbin]), 
                          np.mean(histo[np3s:nbin])]
                    y = histo[:nbin]
                    pars = dict(a=p[0],error_a=p[0]/100,x0=p[1]+0.1,error_x0=p[1]/100,dx=1.1,error_dx=0.01,    
                         ak1=p[3],error_ak1=p[3]/100,ak2=p[4],error_ak2=p[4]/100)
                    level = 1 if mprint else 0
                    mm = muprompt()
                    mm._init_(x,y)
                    m = Minuit(mm,pedantic=False,print_level=level,**pars)
                    m.migrad()
                    A,X0,Dx,Ak1,Ak2 = m.args
                    x0[detector] = X0 # store float peak bin position (fractional)     
                    if mplot:    

                        n1 = npeaks[detector]-50
                        n2 = npeaks[detector]+50
                        x3 = np.arange(n1,n2,1./10.)
                        with self.t0plot_container:
                            if self.first_t0plot:
                                self.axt0[divmod(detector,3)].plot(x[n1:n2],y[n1:n2],'.')
                                self.axt0[divmod(detector,3)].plot(x3,mm.f(x3,A,X0,Dx,Ak1,Ak2))
                                self.prompt_fit_text[detector] = self.axt0[divmod(detector,3)].text(npeaks[detector]
                                                                 +10,0.8*max(y),'Det #{}'.format(detector+1))
                            else:
                                self.axt0[divmod(detector,3)].lines[0].set_ydata(y[n1:n2])
                                self.axt0[divmod(detector,3)].lines[1].set_ydata(mm.f(x3,A,X0,Dx,Ak1,Ak2))
                                x_text = self.prompt_fit_text[detector].get_position()[0]
                                self.prompt_fit_text[detector].set_position((x_text,0.8*max(y)))
                                self.axt0[divmod(detector,3)].relim() # find new dataLim
                                self.axt0[divmod(detector,3)].autoscale_view()

            if mplot:
                if self.first_t0plot:
                    P.show()
                    self.first_t0plot = False
                else:
                    P.draw() # TypeError: draw_wrapper() missing 1 required positional argument: 'renderer'

               ##################################################################################################
                # Simple cases: 
                # 1) Assume the prompt is entirely in bin nt0. (python convention, the bin index is 0,...,n,... 
                # The content of bin nt0 will be the t=0 value for this case and dt0 = 0.
                # The center of bin nt0 will correspond to time t = 0, time = (n-nt0 + mufit.offset + mufit.dt0)*mufit.binWidth_ns/1000.
                # 2) Assume the prompt is equally distributed between n and n+1. Then nt0 = n and dt0 = 0.5, the same formula applies
                # 3) Assume the prompt is 0.45 in n and 0.55 in n+1. Then nt0 = n+1 and dt0 = -0.45, the same formula applies.
                ##################################################################################################

            # these three are the sets of parameters used by other methods
            self.nt0 = x0.round().astype(int) # bin of peak, nd.array of shape run.get_numberHisto_int() 
            self.dt0 = x0-self.nt0 # fraction of bin, nd.array of shape run.get_numberHisto_int() 
            self.lastbin = self.nt0.min() - self.prepostpk[0].value # nd.array of shape run.get_numberHisto_int() 
                                                   # refresh, they may be slightly adjusted by the fit
            self.t0plot_results.clear_output()

            with self.t0plot_results:
                print('\n\n\n\nRun: {}'.format(self._the_run_.get_runNumber_int()))
                print(' Bin nt0')
                for detector in range(self._the_run_.get_numberHisto_int()):
                    print('#{}: {}'.format(detector,self.nt0[detector]))
                print('\n\n dt0 (bins)')
                for detector in range(self._the_run_.get_numberHisto_int()):
                    print('#{}: {:.2f}'.format(detector,self.dt0[detector]))
            ##################################################################################################

        from ipywidgets import HBox, Layout, VBox, Text, IntText, Checkbox, Button, Output, Accordion

        # first tab: setup for things that have to be set initially (paths, t0, etc.)
        # the tab is self.mainwindow.children[0], a VBox 
        # containing a setup_box of three HBoxes: path, and t0plot 
        # path is made of a firstcolumn, paths, and a secondcolumns, filespecs, children of setup_box[0]
        # agt0 is made of three 
        setup_contents = ['path','promptfit','nt0_dt0','t0plot'] # needs two VBoxes
        setup_hbox = [HBox(description=name,layout=Layout(border='solid',)) for name in setup_contents]
        self.mainwindow.children[0].children = setup_hbox # first tab (setup)
        # first path
        paths_content = ['data','tlog','analysis'] # needs a VBox with three Text blocks
        paths_box = VBox(description='paths',layout=Layout(width='60%'))
        self.paths = [Text(description=paths_content[k],layout=Layout(width='90%'),continuous_update=False) for k in range(len(paths_content))]
        # self.paths[k].value='.'+os.path.sep # initial value, 
        paths_box.children = self.paths
        for k in range(len(paths_content)):
            self.paths[k].observe(on_paths_changed,'value')
        filespecs_content = ['fileprefix','extension']
        filespecs_box = VBox(description='filespecs',layout=Layout(width='40%'))
        self.filespecs = [Text(description=filespecs_content[k],layout=Layout(width='90%'),continuous_update=False) 
                          for k in range(len(filespecs_content))]
        filespecs_box.children = self.filespecs
        #        for k in range(len(filespecs)):  # not needed, only check that data and tlog exixt
        #            self.filespecs_list[k].observe(on_filespecs_changed,'value')
        # paths finished
        # now agt0
        self.prepostpk = [IntText(description='prepeak',value = 7, layout=Layout(width='20%'),
                              continuous_update=False), 
                          IntText(description='postpeak',value = 7, layout=Layout(width='20%'),    
                              continuous_update=False)]
        self.prepostpk[0].style.description_width='60%'
        self.prepostpk[1].style.description_width='60%'
        self.plot_check = Checkbox(description='prompt plot',value=True,layout=Layout(width='15%'))
        self.plot_check.style.description_width='10%'
        fit_button = Button(description='prompt fit',layout=Layout(width='15%'))
        fit_button.on_click(on_prompt_fit_click)
        fit_button.style.button_color = 'lightgreen'
        save_button = Button(description='save setup',layout=Layout(width='15%'))
        save_button.style.button_color = 'lightgreen'
        load_button = Button(description='load setup',layout=Layout(width='15%'))
        load_button.style.button_color = 'lightgreen'
        prompt_fit = [self.prepostpk[0], self.prepostpk[1], self.plot_check, fit_button ,save_button, load_button] 
        # fit bin range is [self.binrange[0].value:self.binrange[1].value]
        save_button.on_click(save_setup)
        load_button.on_click(load_setup)
        nt0_dt0 = [Accordion(children=[HBox(children=[Text(description='t0 [bins]'), Text(description='dt0 [bins]')])])]
        nt0_dt0[0].set_title(0,'t0 bins and remainders')
        nt0_dt0[0].selected_index= None

        self.t0plot_container = Output(layout=Layout(width='85%'))     
        self.t0plot_results = Output(layout=Layout(width='15%')) 
        setup_hbox[0].children = [paths_box, filespecs_box]
        setup_hbox[1].children = prompt_fit
        setup_hbox[2].children = nt0_dt0
        setup_hbox[3].children = [self.t0plot_container,self.t0plot_results]
        load_setup([])
        nt0_dt0[0].children[0].children[0].value = ' '.join(map(str,self.nt0.astype(int)))
        nt0_dt0[0].children[0].children[1].value = ' '.join(map('{:.2f}'.format,self.dt0))

##########################
# SUITE
##########################
    def suite(self):
        '''
        suite tab of mugui
        used to select: run (single/suite)
                        load next previous, add next previous
             to print: run number, title, 
                       total counts, group counts, ns/bin
                       comment, start stop date, next run, last add                           
        '''
        def get_totals():
            '''
            calculates the grand totals and group totals after a single run is read
            '''
            import numpy as np
            # called only by self.suite after having loaded a run
            gr = set(np.concatenate((self.grouping['forward'],self.grouping['backward'])))
            tsum, gsum = 0, 0
            if self.offset:  # True if self.offste is already created by self.fit()
                offset_bin = self.offset.value # self.offset.value is int
            else:     # should be False if self.offset = [], as set in self.__init__()
                offset_bin = self.offset0  # temporary parking
           #       self.nt0 roughly set by suite model on_loads_changed
           #       with self._output:
           #            print('offset = {}, nt0 = {}'.format(offset_bin,self.nt0))
            for detector in range(self._the_run_.get_numberHisto_int()):
                n1 = offset_bin+self.nt0[detector] 
                # if self.nt0 not yet read it is False and totals include prompt
                histo = self._the_run_.get_histo_array_int(detector)[n1:].sum()
                tsum += histo
                if detector in gr:
                    gsum += histo
            #       print ('nt0.sum() = {}'.format(self.nt0.sum()))
            self.totalcounts.value = str(tsum)
            self.groupcounts.value = str(gsum)
            self.nsbin.value = '{:.3}'.format(self._the_run_.get_binWidth_ns())

        def muzeropad(runs):
            '''
            utility of the suite tab, not a method
            future:
            1) look for '+' and separate runs to be added
            2) determine how many leading zeros for padding
               read a file from data dir
               check number of digits before '.'
               count number of digits in run
               zero pad
            now:
            0) zeroth version pads a fixed number of zeros to 4 digits
            '''
            zeros='0000'
            if len(runs)<len(zeros):
                return zeros[:len(zeros)-len(runs)]+runs
            elif len(runs)==len(zeros):
                return runs
            else:
                print('Too long run number!')
                return []

        def on_loads_changed(change):
            '''
            observe response of suite tab widgets:
            try loading a run via musrpy 
            '''
            import mujpy.musr2py.musr2py as muload
            import numpy as np
            import os

            run_or_runs = change['owner'].description # description is either 'Single run' or 'Run  suite'
            if run_or_runs == loads[0]: # 'Single run'
                self._global_fit_ = False # check when introducing suites
                self._the_run_ = muload.musr2py() # *** fix a better integration between classes, mugui, mufit, muset, musuite ***
                path_and_filename = '' 
                filename = ''
                filename = filename.join([self.filespecs[0].value,muzeropad(self.loads_handles[0].value),'.',self.filespecs[1].value])
                path_and_filename = os.path.join(self.paths[0].value,filename)
                                    # data path + filespec + padded run rumber + extension)
                if self._the_run_.read(path_and_filename) == 1: # error condition, set by musr2py.cpp
                    self.mainwindow.selected_index=3
                    with self._output:
                        print ('\nFile {} not read. Check paths, filespecs and run rumber on setup tab'.
                                format(path_and_filename))
                else:
                    self.title.value = '{} {} {} {}'.format(self._the_run_.get_sample(),self._the_run_.get_field(),
                                                                self._the_run_.get_orient(),self._the_run_.get_temp())                    
                    self.comment_handles[0].value = self._the_run_.get_comment() 
                    self.temperatures = self._the_run_.get_temperatures_vector() # to be displaced to the .value of a widget in the tlog tab
                    self.temperaturedevs = self._the_run_.get_devTemperatures_vector()
                    self.comment_handles[1].value = self._the_run_.get_timeStart_vector() 
                    self.comment_handles[2].value = self._the_run_.get_timeStop_vector()
                    self._the_run_display.value = str(self._the_run_.get_runNumber_int()) # = self.loads_handles[0].value 
                    try:
                        dummy = self.nt0.sum() # fails if self.nt0 does not exist yet
                    except: # runs this only if self.nt0 does not exist
                        self.nt0 = np.zeros(self._the_run_.get_numberHisto_int(),dtype=int)
                        self.dt0 = np.zeros(self._the_run_.get_numberHisto_int(),dtype=float)
                        for k in range(self._the_run_.get_numberHisto_int()):
                            self.nt0[k] = np.where(self._the_run_.get_histo_array_int(k)==
                                                   self._the_run_.get_histo_array_int(k).max())[0][0]
                    get_totals() # sets totalcounts, groupcounts and nsbin                                
            else:
                # multi run
                self._global_fit_ = True
                print('to be implemented ...')
             
        from ipywidgets import HBox, Layout, Text, Button

        # second tab: select run or suite of runs (for sequential or global fits)
        # the tab is self.mainwindow.children[1], a VBox 
        # containing three HBoxes, loads_box, comment_box, speedloads_box
        # path is made of a firstcolumn, paths, and a secondcolumns, filespecs, children of setup_box[0]
        loads = ['Single run','Run suite']
        speedloads = ['Next run' 'Load next', 'Load previous', 'Add next', 'Add previous', 'Last added']
        loads_box = HBox(description='loads',layout=Layout(width='100%'))
        comment_box = HBox(description='comment',layout=Layout(width='100%'))
        speedloads_box = HBox(description='speedloads',layout=Layout(width='100%'))
        width = ['50%','150%']
        self.loads_handles = [Text(description=loads[k],layout=Layout(width=width[k]),continuous_update=False) 
                              for k in range(len(loads))]
        self.loads_handles[0].observe(on_loads_changed,'value')
        self.loads_handles[1].observe(on_loads_changed,'value')
        self.comment_handles = [Text(description='Comment',layout=Layout(width='30%'),disable=True),
                                Text(description='Start date',layout=Layout(width='30%'),disable=True),
                                Text(description='Stop date',layout=Layout(width='30%'),disable=True)]
        # the following doesn't work yet
        Ln_button = Button(description='Load nxt')
        Ln_button.style.button_color = 'lightgreen'
        Lp_button = Button(description='Load prv')
        Lp_button.style.button_color = 'lightgreen'
        An_button = Button(description='Add nxt')
        An_button.style.button_color = 'lightgreen'
        Ap_button = Button(description='Add prv')
        Ap_button.style.button_color = 'lightgreen'
        self.speedloads_handles = [Text(description='Next run',disabled=True),
                                   Ln_button, Lp_button, An_button, Ap_button, 
                                   Text(description='Last add',disabled=True)]
        loads_box.children = self.loads_handles
        comment_box.children = self.comment_handles
        speedloads_box.children = self.speedloads_handles

        self.mainwindow.children[1].children = [loads_box, comment_box, speedloads_box] # second tab (suite)


