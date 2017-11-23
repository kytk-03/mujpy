.. _reference:

+++++++++
Reference
+++++++++

This is the mujpy Reference Manual

------
Header
------
The top part of the mujpy gui contains a logo together with the information on the loaded run (in `single run mode`_) or the master loaded run (in the `suite mode`_). The information contains the run number, the run title, the total conuts, the total counts over the selected `grouping`_ of detectors.

----
Tabs
----
The lower part of the gui is divided in tabs: `setup`_, `suite`_, `fit`_, `output`_, `fft`_, `plots`_, `about`_.
Click on the tab to see its content and methods, activated by widget controls (buttons, checkboxes, input text field, etc.)

-----
setup
-----
The setup tab contains preliminary information and fields in three boxes. 

The first box contains paths on the right: the directory path (the folder where the data files are stored), the tlog path (the folder where the temperature log files are stored), the log path (the folder where log files will be saved). On the left the two strings that compose the data file name, the fileprefix and the file extension, the remaining string being the run number, possibly with a number of auto-generated leading zeros. 

The path from which the jupyter notebook is lauched, startup_path,  is saved automatically. 

The second box is related to the fit that identifies the implantation time, t=0, for the muon spin evolution inside the sample. At PSI this is identified with the center of the prompt peak, the sharp increase in count rate that corresponds to positrons detected erroneously as implanting muons by the start detector. The fit is performed by iminuit, minimizing the math:`\chi^2` of the initial data slice over the function  `muprompt`_. The prepeak and postpeak parameters represent the peak interval span (the number of bins) respectively before and after the maximum. The prompt plot checkbox determines whether a plot is produced. The three left buttons respectively launch this fit, load and save a setup file that stores the information of this tab in the file mujpy_setup.pkl

-----
suite
-----
Description

single run mode
---------------
Description

suite mode
----------
Description

---
fit
---
Description

grouping
--------
Description

muprompt
--------
Description

------
output
------
Description

---
fft
---
Description

-----
plots
-----
Description

-----
about
-----
Description

