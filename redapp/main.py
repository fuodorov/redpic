# Library
import redpic as rp
import kenv as kv
import numpy as np
import pandas as pd
import glob

# Bokeh basics
from bokeh.io import curdoc
from bokeh.models.widgets import Tabs
from bokeh.models import ColumnDataSource
from os.path import dirname, join

# Each tab is drawn by one script
from beam import beam_tab
from accelerator import accelerator_tab
from tuning import tuning_tab
from production import production_tab

print(rp.__doc__)
print('Version: ' + rp.__version__)
# Initial parameters
beam = rp.Beam(rp.electron)
beam.generate(rp.Distribution('KV', x=1, y=1, z=1, px=1, py=1, pz=1), n=1)
accelerator = rp.Accelerator(0, 1, 0.001)
accelerator.compile()

# Create each of the tabs
tab1 = beam_tab(beam)
tab2 = accelerator_tab(accelerator)
tab3 = tuning_tab(beam, accelerator)
tab4 = production_tab(beam, accelerator)

# Put all the tabs into one application
tabs = Tabs(tabs=[tab1, tab2, tab3, tab4])

# Put the tabs in the current document for display
curdoc().add_root(tabs)
curdoc().title = "REDPIC"
