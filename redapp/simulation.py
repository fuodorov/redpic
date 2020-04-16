import redpic as rp
import kenv as kv
import pandas as pd
import numpy as np
import glob
from multiprocessing import Process
from bokeh.models import (Panel, PreText, RadioButtonGroup, Select,
                           FileInput, Button, TextInput, ColumnDataSource,
                          RangeSlider, Div)
from bokeh.plotting import figure, curdoc
from bokeh.layouts import column, row
from bokeh.models.widgets import Slider
from os.path import dirname, join


def simulation_tab(beam, acc):
    '''Simulation

    '''
    pre = PreText(text='''Here you can run the beam simulation in the accelerator.''')

    start_button = Button(label='Start', button_type="success")

    refresh_button = Button(label='Refresh', button_type="warning")

    p = figure(x_axis_label='z [m]', y_axis_label='x [m]', x_range=(acc.z_start, acc.z_stop),
                        width=600, height=210)
    source = ColumnDataSource(data={'x': beam.df['x'], 'z': beam.df['z']})
    p.scatter('z', 'x', source=source, size=0.5, color="#3A5785", alpha=0.5)

    controls = row(start_button, refresh_button)
    tab = Panel(child=column(p, controls), title='Simulation')

    def start_handler(new, beam=beam, acc=acc):
        simulation = rp.Simulation(beam, acc)
        simulation.track(n_files=30, path=dirname(__file__) + '/data/')

    def play_handler(beam=beam, acc=acc):
        track_files = np.sort(glob.glob(dirname(__file__) + '/data/'+ '*.*[0-9].csv'))
        if len(track_files) > 0:
            df = pd.read_csv(track_files[-1], dtype='float32')
            source.data = {'x': df['x'], 'z': df['z']}

    start_button.on_click(start_handler)
    play_button.on_click(play_handler)

    return tab
