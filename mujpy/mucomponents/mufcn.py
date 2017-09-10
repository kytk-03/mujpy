import numpy as np
''' big common fcn
    inherits globals 
    MU_COMPONENT with list of all model components
    MU_MODEL with selected model and MU_MODEL.ALPHA
    par has the internal fminuit indices
    data(1,:) is TIME in microseconds
    if size(data,1)>1 data(2,:) is ASYMMETRY and data(3,:) its ERROR
 ROUGH DESCRIPTION
            rearrange parameters 
            calculate contributions
            and add corrections (e.g. dalpha)
'''
def mufcn(x,*kwargs):
    f=0.
    return f
# global MU_COMPONENT MU_MODEL MU_TABLE
# to be cythonized


#f  = np.zeros(mufit.asymmetry.shape); # start with zero f 
#for k in range(0,len(comp)):
#    kk = [d['name'] for d in available_components].index(comp[k])
#    % order number of model component k in the available_components tuple of directories
#    p=par; % restarts from MINUIT input values
#    command='[';
#    for kp in available_components[kk]['npar']: # here
#        command=[command ... 
#                 MU_MODEL.COMPONENT(k).PARAMETER(kp).MENU2MINUIT  ...
#                 ' '];
#        % MENU2MINUIT is a reallocated string where MUFIT
#        % recalculates fminuit parameter indices from Menu indices
#        % the string is either e.g. 'p(3)' (~!=)
#        % or e.g. 'p(2)*0.333+p(7)*0.666' (+) define it as pp[n]=lambda x: p[2]*0.333-p[7]*0.666 and use pp 
#    end
#    command=[command ']']; % command contains something like
#                             %‘p=[p(4) p(1)*0.333 p(2) ];'
#                             % or just 'p=[p(4) p(5) p(6) ];'
#    p=eval(command);% now p contains the right parameters for component k 
#    %disp(char(MU_COMPONENT.DEFINITION(kk)))
#    eval(['f=[f; ' char(MU_COMPONENT.DEFINITION(kk)) '];'])
#    % computes function from definition
#    % adding each component as a new row;
#    % pdf (polarization function) at center of time bin, data(:,1), 
#    % or, ideally, (cdf(data(1,:)+dt/2)-cdf(data(1,:)-dt/2))/dt, where 
#    % cdf is integral of pdf and dt is time bin width
#end

#% uncorrected sum
#f=sum(f,1); % sum along columns (components), calculates function

# % correction factors (like 'da')

#kda=find(strcmp(cellstr(strvcat(MU_MODEL.COMPONENT.NAME)),'da'));
#% is 'da' one of the ingredients? 
#% if it exists kda is the order of the 'da' component in MU_MODEL
#if ~isempty(kda)    
#    dalpha=par(MU_MODEL.COMPONENT(kda).SHARE);
#    dadalpha=dalpha/MU_MODEL.ALPHA;
#    f=((2.d0+dadalpha)*f-dadalpha)./...
#      ((2.d0+dadalpha)-dadalpha*f); 
#end

#if (size(data,1)>2), % error = 3rd row of data
#  f = sum(((data(2,:) - f)./data(3,:)).^2);
#end
