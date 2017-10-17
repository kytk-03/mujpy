# link fit tab methods to self._the_fit_ 
# write on_fit_request for fit_button.on_click(on_fit_request)
# test
# check that self attributes are all needed (some might just be local variables)
# write tlog
# write plot
# write edit
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
        from mujpy import __file__ as MuJPyName
        from IPython.display import display

        self.interval = np.array([0,7800], dtype=int)

        self.firstbin = 0
        self.second_plateau = 100
        self.peakheight = 100000.
        self.peakwidth = 1.   # broad guesses for default

        self.nt0 = np.zeros((0,),dtype=int) # bin of peak, nd.array of shape run.get_numberHisto_int(), initialized as False        
        self.dt0 = np.zeros((0,),dtype=float) # fraction of bin, nd.array of shape run.get_numberHisto_int(), initialized as False 
        self.grouping = {'forward':np.array([2]),'backward':np.array([1])} # normal dict
        self._the_run_ = []  # if self._the_run_: is False, to check whether the handle is created
        self.group = [] # if self.group: is False, to check whether the handle is created
        self.alpha0 = 1.05
        self.alpha = [] # if self.alpha: is False, to check whether the handle is created
        self.offset0 = 7 
        self.offset = [] # if self.offset: is False, to check whether the handle is created
        self.first_t0plot = True
        # mujpy paths
        self.__path__ = os.path.dirname(MuJPyName)
        self.__logopath__ = os.path.join(self.__path__,"logo")
        self.__startuppath__ = os.getcwd() # working directory, in which, in case, to find mujpy_setup.pkl 

        #####################################
        # actually produces the gui interface
        #####################################

        self.gui()
        self.setup()
        self.suite()
        self.fit()
        # self.edit
        # self fft
        # self.plot
        # self.tlog
        self.about()
        display(self.gui)
        
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
        self.out = Output(layout=Layout(width='100%'))
        thirdrow = HBox([self.out],layout=Layout(height='60px',width='100%',overflow_y='scroll',overflow_x='scroll')) # x works y does scroll
        titlewindow = VBox()
        titlewindow_content = [titlerow, secondrow, thirdrow]
        titlewindow.children = titlewindow_content
        titlelogowindow = HBox()
        titlelogowindow_content = [logo, titlewindow]
        titlelogowindow.children = titlelogowindow_content

        # main layout: tabs
        tabs_contents = ['setup', 'suite', 'fit', 'edit', 'fft','plot','tlog','about']
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
    def fit(self, model_in='daml'):
        '''
        fit tab of mugui
        used to set: alpha, offset, forw and backw groups
                     fit and plot ranges, model version                    
             to display: model name 
             to activate: fit, plot and update buttons
             to select and load model (load from folder missing)
             to select parameters value, fix, function, fft subtract check
        '''
        # no need to observe parvalue, since their value is a perfect storage point for the latest value
        # validity check before calling fit

        def arrange(bin_range,target = 'fit'):
            '''
            arrange(bin_range) 
            [arrange(bin_range, target='plot')]
            writes bin_range into self.fit_range 
            [self.plot_range] as csv integers
            '''
            return exec(target+'_range = str(bin_range[0])+", "+str(bin_range[1])')

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
                with self.out:
                    print("comma separated values, please")
                return -1,-1 
            try:
                first = int(string[:comma])
            except:
                with self.out:
                    print("Could not read first bin")
                return -1,0
            try:
                last = int(string[comma+1:])
            except:
                with self.out:
                    print("Could not read last bin")
                return 0,-1
            return first,last 

            # _alpha = self.alpha.value # alpha is float
            # _offset = self.offset.value # offset, integer, from t0 to first good bin for fit    
            # _grouping = self.grouping # {'forward':np.array,'backward':np.array}
            # _fit_range = [self.derange()[k] for k in range(2)] # changes to fit_range (derange('target='plot') for plot_range)
            # _plot_range = [self.derange(target='plot')[k] for k in range(2)] # changes to plot_range (derange('target='plot') for plot_range)

        def get_grouping(self,name):
            """
            usage:

            MuJPy.get_grouping(name) # name = 'forward' or 'backward' groups for GPS &tc
            the numbers in the are parsed, 
            translated into histo indices (n-1, python style) 
            and stored in a dictionary of two np.arrays
            """
            import numpy as np

            # two shorthands: either a list, comma separated, such as 1,3,5,6 
            # or a pair of integers, separated by a colon, such as 1:3 = 1,2,3 
            # only one column is allowed, but 1, 3, 5 , 7:9 = 1, 3, 5, 7, 8, 9 
            # or 1:3,5,7 = 1,2,3,5,7  are also valid
            #       get the shorthand from the gui Text 
            groups = ['forward','backward']
            groupcsv = self.group[groups.index(name)].value # self.group[0,1] are handles to the correspomdimg Texts
            #       now parse groupcsv shorthand
            groupcsv = groupcsv.replace('.',',') # can only be a mistake: '.' means ','
            if groupcsv.find(':')==-1:
                # colon not found, csv only
                self.grouping[name] = np.array([int(s) for s in groupcsv.split(',')])
            else:
                # colon found
                try:
                    if groupcsv.find(',')+1: # True if found, False if not found (not found yields -1)    
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
                except:
                    print('Wrong group syntax: {}'.format(self.group[groups.index(name)].value))
                    self.group[groups.index(name)].value = ''
    
        def int2_int(self):
            '''
            From internal parameters to the minimal representation 
            for the use of mucomponents._add_.
            Invoked just before submitting minuit 
            '''
            from mujpy.mucomponents.mucomponents import mumodel

            ntot = sum([len(self._the_fit_.internal_components[k]['pars']) for k in range(len(self._the_fit_.internal_components))])
            lmin = [-1]*ntot
            nint = -1 # initialize
            nmin = -1 # initialize
            _int = []
            for k in range(len(self._the_fit_.internal_components)):  # scan the model
                component = mumodel().__class__.__dict__[self._the_fit_.internal_components[k]['name']] # the method
                keys = []
                for j in range(len(self._the_fit_.internal_components[k]['pars'])): # 
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
                _int.append([component,keys])
           # for k in range(len(_int)):
           #     print(_int[k])      

            return _int


        def int2min(self):
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
            [plus
               the local replica of the non global component parameters
               to be implemented]
            '''
            ntot = sum([len(self._the_fit_.internal_components[k]['pars']) for k in range(len(self._the_fit_.internal_components))])
            lmin = [-1]*ntot
            nint = -1 # initialize
            nmin = -1 # initialize
            fitarg = {}
            for k in range(len(self._the_fit_.internal_components)):  # scan the model
                component_name = self._the_fit_.internal_components[k]['name'] 
                keys = []
                for j in range(len(self._the_fit_.internal_components[k]['pars'])): # 
                    nint += 1  # internal parameter incremente always   
                    if flag[nint].value != '=': #  skip functions, they are not new minuit parameter
                        keys.append('~')
                        nmin += 1
                        lmin[nmin] = nint # number of minuit parameter
                _int.append([component,keys])

            fitarg.update()
            return fitarg

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

        def min2int(self):
            '''
            From minuit parameters to internal parameters,
            see int2min for a description   
            Invoked just after minuit convergence.
            '''
            
        def on_fit_request(b):
            '''
            retrieve data from the gui:
            parameters values (parvalue[nint].value), flags (flag[nint].value), 
            errors, limits, functions (function[nint].value), self.alpha.value, self.fit_range.value
            construct _int, needed by mumodel._add_ to distribute minuit parameter
            costruct fitar dictionary, needed by migrad (and provided by it at the end
            pass them to self._the_fit_
            mumodel._load_data_
            call minuit
            save values in self._the_fit_
            save 
            print summary
            '''
            from mujpy.mucomponents.mucomponents import mumodel
            from iminuit import Minuit as M 

            _the_model_ = mumodel() # local instance
            self._the_fit_.asymmetry(self._the_run_) # prepare asymmetry
            _the_model_._load_data_(self._the_fit_.time,self._the_fit_.asymm,
                                 self._int,self._alpha,e=self._the_fit_.asyme) # pass data to model
            pars = self.int2min
            m = M(_the_model_,pedantic=False,print_level=level,**pars)
      

        def on_flag_changed(change):
            '''
            observe response of fit tab widgets:
            set disabled on corresponding function (True if flag=='!' or '~', False if flag=='=') 
            '''
            dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            n = int(dscr[4:]) # description='flag'+str(nint), skip 'flag'
            function[n].disabled=False if change['new']=='=' else True

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
            self.get_grouping(name) # stores self.grouping shorthand in self.group dict

        def on_integer(change):
            name = change['owner'].description
            if name == 'offset':
                if self.offset.value<0: # must be positive
                   self.offset.value = 7 # standard value

        def on_loadmodel_changed(change):
            '''
            observe response of fit tab widgets:
            check that change['new'] is a valid model
            relaunch MuJPy.fit(change['new'])
            '''
            if self._the_fit_.checkvalidmodel(change['new']): # empty list is False, non empty list is True
                self.fit(change['new']) # restart the gui with a new model
                self.mainwindow.selected_index=2
                self.gui
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

        ######### here starts the fit method of MuGui
        from ipywidgets import FloatText, Text, IntText, Layout, Button, HBox, \
                               Checkbox, VBox, Dropdown
        from mujpy.mufit import MuFit as mufit
        import numpy as np

        self._the_fit_ = mufit()  # this is the official instance of mufit inside mugui
        # self._the_fit_.deletemodel() # this isn't really needed, because mufit is initialized bare of model

        self.alpha = FloatText(description='alpha',value='{:.4f}'.format(self.alpha0),
                                layout=Layout(width='20%'),continuous_update=False) # self.alpha.value
        self.offset = IntText(description='offset',value=self.offset0,
                                layout=Layout(width='20%'),continuous_update=False) # offset, is an integer
        # initialized to 7, only input is from an IntText, integer value, or saved and reloaded from mujpy_setup.pkl
        self.alpha.style.description_width='40%' 
        self.offset.style.description_width='50%' 
        self.offset.observe(on_integer,'value') # check validity, must be positive
        self.group = [Text(value=str(self.grouping['forward']).strip('[]'),description='forward',layout=Layout(width='27%'),
                                    continuous_update=False),
                      Text(value=str(self.grouping['backward']).strip('[]'),description='backward',layout=Layout(width='27%'),
                                    continuous_update=False)]
        self.group[0].observe(on_group_changed,'value')
        self.group[1].observe(on_group_changed,'value')
        # end moved

        model = Text(description = '', layout=Layout(width='10%'), disabled = True) # this is static, empty description, next to loadmodel
        model.value = model_in
        loadmodel = Text(description='loadmodel',layout=Layout(width='19%'),continuous_update=False) # this is where one can input a new model name
        loadmodel.observe(on_loadmodel_changed,'value')
        loadmodel.style.description_width='35%'
        version = IntText(description='version',value='1',layout=Layout(width='11%',indent=False)) # version.value is an int
        version.style.description_width='43%'
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
                                                                version,
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

        self._the_fit_.create_model(model.value)

        leftframe_list, rightframe_list = [],[]
  
        words  = ['#','name','value','~!=','function']
        nint = -1 # internal parameter count, each widget its unique name
        ntot = np.array([len(self._the_fit_.internal_components[k]['pars']) for k in range(len(self._the_fit_.internal_components))]).sum()
        parvalue, flag, function = [], [], [] # lists, index runs according to internal parameter count nint
        self.compar = {} # dictionary: key nint corresponds to a list of two values, c (int index of component) and p (int index of parameter)
                # use: self.compar[nint] is a list of two integers, the component index k and its parameter index j

        for k in range(len(self._the_fit_.internal_components)):  # scan the model            
            header = HBox([ Text(value=self._the_fit_.internal_components[k]['name'],disabled=True,layout=Layout(width='8%')),
                           Checkbox(description='FFT',value=True) ]) # list of HBoxes, the first is the header for the component
                                                                     # composed of the name (e.g. 'da') and the FFT flag
                                                                     # fft will be applied to a 'residue' where only checked components
                                                                     # are subtracted
            componentframe_list = [header] # list of HBoxes, header and pars
            componentframe_handle = VBox()
            for j in range(len(self._the_fit_.internal_components[k]['pars'])): # make a new par for each parameter 
                                                                                # and append it to component_frame_content
                nint += 1      # all parameters are internal parameters, first is pythonically zero 
                self.compar.update({nint:[k,j]}) # stores the correspondence between nint and component,parameter
                nintlabel_handle = Text(value=str(nint),layout=Layout(width='10%'),disabled=True)
                parname_handle = Text(value=self._the_fit_.internal_components[k]['pars'][j]['name'],layout=Layout(width='15%'),disabled=True)
                # parname can be overwritten, not important to store

                parvalue.append(Text(value='{:.4}'.format(self._the_fit_.internal_components[k]['pars'][j]['value']),
                                     layout=Layout(width='20%'),description='value'+str(nint),continuous_update=False))
                parvalue[nint].style.description_width='0%'
                # parvalue handle must be unique and stored at position nint, it will provide the initial guess for the fit

                function.append(Text(value=self._the_fit_.internal_components[k]['pars'][j]['function'],
                                     layout=Layout(width='38%'),description='func'+str(nint),continuous_update=False))
                function[nint].style.description_width='0%'
                # function handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation 

                fdis = False if self._the_fit_.internal_components[k]['pars'][j]['flag']=='=' else True 
                function[nint].disabled = fdis # enabled only if flag='='
                flag.append(Dropdown(options=['~','!','='], 
                                     value=self._the_fit_.internal_components[k]['pars'][j]['flag'],
                                     layout=Layout(width='10%'),description='flag'+str(nint)))
                flag[nint].style.description_width='0%'
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
                self.nt0 = mujpy_setup['_nt0'] # bin of peak, nd.array of shape run.get_numberHisto_int()
                self.dt0 = mujpy_setup['_dt0'] # fraction of bin, nd.array of shape run.get_numberHisto_int()
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
            _nt0 = self.nt0 # numpy array
            _dt0 = self.dt0 # numpy array 
            names = ['_paths_content','_filespecs_content',
                      '_prepostpk','_nt0','_dt0'] # keys
            setup_dict = {}
            for k,key in enumerate(names):
               setup_dict[names[k]] = eval(key) # key:value
            with open(path,'wb') as f:
                pickle.dump(setup_dict, f) # according to __getstate__()
            with self.out:
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
                print('No run loaded yet! load one first.')
                return    
            else:
                npeaks = np.array([np.argmax(self._the_run_.get_histo_array_int(det)) for det in range(self._the_run_.get_numberHisto_int())])
                # approximate positions of peaks
                nbin =  max(npeaks) + self.second_plateau # this sets a detector dependent second plateau bin interval
                x = np.arange(0,nbin,dtype=int) # nbin bins from 0 to nbin-1
                self.lastbin, np3s = npeaks - self.prepostpk[0].value, npeaks + self.prepostpk[1].value # final bin of first and 
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
                          np.mean(histo[self.firstbin:self.lastbin[detector]]), 
                          np.mean(histo[np3s[detector]:nbin])]
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
            self.lastbin = self.nt0 - self.prepostpk[0].value # nd.array of shape run.get_numberHisto_int() 
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

        def set_grouping(name):
            """
            usage:

            MuJPy.set_grouping(name) # name = 'forward' or 'backward' groups for GPS &tc
            reads from a dictionary of two np.arrays
            and produces the self.grouping[].value shorthand string that goes in the gui
            """
            import numpy as np

            # allows any succession of csv and 'n:m' formats
            # compatible with get_groupings, but not implemented there yet
            # get the shorthand from the gui Text 
            groups = ['forward','backward']
            group = self.grouping[name]
            groupsplit=np.split(group, np.where(np.diff(group) != 1)[0]+1)
            s = ','
            n=''
            a=''
            for k in range(len(groupsplit)):
                if len(groupsplit[k])>1: # two or more consecutive detectors
                    a = s.join([a,n.join([str(groupsplit[k][0]),':',str(groupsplit[k][-1])])])
                else: # singleton detectors
                    a = s.join([a,str(groupsplit[k]).strip('[]')])
            # unroll np.array group into shorthand
            groupcsv = a[1:]        
            self.group[groups.index(name)].value = groupcsv # +1 to account for Label before handles!

        from ipywidgets import HBox, Layout, VBox, Text, IntText, Checkbox, Button, Output

        # first tab: setup for things that have to be set initially (paths, t0, etc.)
        # the tab is self.mainwindow.children[0], a VBox 
        # containing a setup_box of three HBoxes: path, and t0plot 
        # path is made of a firstcolumn, paths, and a secondcolumns, filespecs, children of setup_box[0]
        # agt0 is made of three 
        setup_contents = ['path','promptfit','t0plot'] # needs two VBoxes
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

        self.t0plot_container = Output(layout=Layout(width='85%'))     
        self.t0plot_results = Output(layout=Layout(width='15%')) 
        setup_hbox[0].children = [paths_box, filespecs_box]
        setup_hbox[1].children = prompt_fit
        setup_hbox[2].children = [self.t0plot_container,self.t0plot_results]
        load_setup([])

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
           #       with self.out:
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
                self.globalfit = False
                self._the_run_ = muload.musr2py() # *** fix a better integration between classes, mugui, mufit, muset, musuite ***
                path_and_filename = '' 
                filename = ''
                filename = filename.join([self.filespecs[0].value,muzeropad(self.loads_handles[0].value),'.',self.filespecs[1].value])
                path_and_filename = os.path.join(self.paths[0].value,filename)
                                    # data path + filespec + padded run rumber + extension)
                if self._the_run_.read(path_and_filename) == 1: # error condition, set by musr2py.cpp
                    with self.out:
                        print ('File {} not read. Check paths, filespecs and run rumber'.format(path_and_filename))
                else:
                    self.title.value = '{} {} {} {}'.format(self._the_run_.get_sample(),self._the_run_.get_field(),
                                                                self._the_run_.get_orient(),self._the_run_.get_temp())                    
                    self.comment_handles[0].value = self._the_run_.get_comment() 
                    self.temperatures = self._the_run_.get_temperatures_vector() # to be displaced to the .value of a widget in the tlog tab
                    self.temperaturedevs = self._the_run_.get_devTemperatures_vector()
                    self.comment_handles[1].value = self._the_run_.get_timeStart_vector() 
                    self.comment_handles[2].value = self._the_run_.get_timeStop_vector()
                    self._the_run_display.value = self.loads_handles[0].value   # need to connect self._the_run_.get_runNumber_int, missing in musr2py  
                    self.nt0 = np.zeros(self._the_run_.get_numberHisto_int(),dtype=int)
                    self.dt0 = np.zeros(self._the_run_.get_numberHisto_int(),dtype=float)
                    for k in range(self._the_run_.get_numberHisto_int()):
                        self.nt0[k] = int(max(self._the_run_.get_histo_array_int(k)))
                    get_totals() # sets totalcounts, groupcounts and nsbin                                
            else:
                # multi run
                self.globalfit = True
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


