from mujpy.mugui import mugui

nogui = mugui() # instance
b = [] # dummy
nogui._load_fit(b) # _load_fit is a backdoor to open the load_fit window for choosing a model
nogui.loads_handles[0].value='433' # like writing 433 in the jupyter ipywidgets text, loads the data 
print('Total counts = {}'.format(nogui.totalcounts.value))
nogui._fit(b) # _fit is backdoor to on_fit_request
