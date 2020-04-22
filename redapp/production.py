import redpic as rp
import kenv as kv
import pandas as pd
import numpy as np
import glob
from scipy import interpolate
from multiprocessing import Process
from bokeh.models import (Panel, PreText, RadioButtonGroup, RadioGroup, Select,
                          FileInput, Button, TextInput, ColumnDataSource,
                          RangeSlider, Div)
from bokeh.plotting import figure, curdoc
from bokeh.layouts import column, row
from bokeh.models.widgets import Slider, Tabs
from bokeh.models.ranges import Range1d
from os.path import dirname, join

def window_rms(a: np.array, window_size: int=1_000):
    '''windows RMS calculation

    '''
    a2 = np.power(a, 2)
    window = np.ones(int(window_size)) / float(window_size)
    return np.sqrt(np.convolve(a2, window, 'same'))

def window_mean(a: np.array, window_size: int=1_000):
    '''windows MEAN calculation

    '''
    a2 = a
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(a2, window, 'same')

def production_tab(beam, acc):
    '''Simulation

    '''
    pre = PreText(text='''Track beam in the accelerator.''')

    track_button = Button(label='Track', button_type="success", width=315)

    refresh_button = Button(label='Refresh', button_type="warning", width=315)

    mode_button = RadioButtonGroup(
        labels=['Phase', 'Line', 'Field'], active=0)

    phase_button = RadioButtonGroup(
        labels=['z-x', 'z-y', 'x-y', 'x-px', 'x-py', 'y-px', 'y-py', 'px-py'],
                active=0)
    field_button = RadioButtonGroup(
        labels=['z-Ez', 'z-Ex', 'z-Ey', 'z-Bz', 'z-Bx', 'z-By'], active=0)
    line_button = RadioButtonGroup(
        labels=['rms x', 'rms y', 'rms x Emitt', 'rms y Emitt', 'avg Energy',
                'x β-func', 'y β-func', 'avg x', 'avg y'], active=0)

    file_input = FileInput(accept='.csv')

    phase_source = ColumnDataSource(data={'x': beam.df['z'], 'y': beam.df['x']})
    phase_plot = figure(x_axis_label='z [m]', y_axis_label='x [m]',
                        width=650, height=250)
    phase_plot.scatter('x', 'y', source=phase_source, size=0.5, color="#3A5785",
                        alpha=0.5)

    beam.df['x2'] = beam.df['x']*beam.df['x']
    beam.df['rms x'] = beam.df['x2'].mean()**0.5
    line_source = ColumnDataSource(data={'x': beam.df['z'], 'y': beam.df['rms x']})
    line_plot = figure(x_axis_label='z [m]', y_axis_label='rms x [m]',
                       width=650, height=250)
    line_plot.line('x', 'y', source=line_source, line_color='#3A5785',
                   line_alpha=0.7, line_width=1)

    field_source = ColumnDataSource(data={'x': beam.df['z'], 'y': beam.df['rms x']})
    field_plot = figure(x_axis_label='z [m]', y_axis_label='Ez [MV/m]',
                        width=650, height=250)
    field_plot.line('x', 'y', source=field_source, line_color='#3A5785',
                    line_alpha=0.7, line_width=1)

    buttons = row(track_button, refresh_button)
    controls = row(pre, file_input)
    phase_tab = Panel(child=column(phase_button, phase_plot), title='Phase')
    line_tab = Panel(child=column(line_button, line_plot), title='Line')
    field_tab = Panel(child=column(field_button, field_plot), title='Field')
    tabs = Tabs(tabs=[phase_tab, line_tab, field_tab])
    tab = Panel(child=column(controls, tabs, buttons), title='Production')

    def calculate(df):
        df = df.sort_values('z')
        df['avg p'] = window_mean((df['px']*df['px'] + df['py']*df['py'] + df['pz']*df['pz'])**0.5)
        df['xp'] = df['px']/df['pz']
        df['yp'] = df['py']/df['pz']
        df['x xp'] = window_mean(df['x']*df['xp'])
        df['y yp'] = window_mean(df['y']*df['yp'])
        df['rms x'] = window_rms(df['x'])
        df['rms y'] = window_rms(df['y'])
        df['rms xp'] = window_rms(df['xp'])
        df['rms yp'] = window_rms(df['yp'])
        df['rms x Emittance'] = (df['rms x']**2 * df['rms xp']**2 - df['x xp']**2)**0.5
        df['rms y Emittance'] = (df['rms y']**2 * df['rms yp']**2 - df['y yp']**2)**0.5
        df['avg Energy'] = df['avg p'] - beam.type.mass*rp.c**2/rp.e/1e6
        df['x β-function'] = df['rms x']**2/df['rms x Emittance']
        df['y β-function'] = df['rms y']**2/df['rms y Emittance']
        df['centroid x'] = window_mean(df['x'])
        df['centroid y'] = window_mean(df['y'])
        df['avg Ez'] = window_mean(df['Ez'])
        df['avg Ex'] = window_mean(df['Ex'])
        df['avg Ey'] = window_mean(df['Ey'])
        df['avg Bz'] = window_mean(df['Bz'])
        df['avg Bx'] = window_mean(df['Bx'])
        df['avg By'] = window_mean(df['By'])
        return df

    def track_handler(new, beam=beam, acc=acc):
        simulation = rp.Simulation(beam, acc)
        simulation.track(n_files=30, path=dirname(__file__) + '/data/')

    def refresh_handler(beam=beam, acc=acc):
        fname = dirname(__file__) + '/data/' + file_input.filename
        df = pd.read_csv(fname, dtype='float32')
        df = calculate(df)
        df = df[df.z >= acc.z_start]
        df = df[df.z <= acc.z_stop]
        if tabs.active == 0:
            if phase_button.active == 0:
                phase_plot.xaxis.axis_label = 'z [m]'
                phase_plot.yaxis.axis_label = 'x [m]'
                phase_source.data = {'x': df['z'], 'y': df['x']}
            if phase_button.active == 1:
                phase_plot.xaxis.axis_label = 'z [m]'
                phase_plot.yaxis.axis_label = 'y [m]'
                phase_source.data = {'x': df['z'], 'y': df['y']}
            if phase_button.active == 2:
                phase_plot.xaxis.axis_label = 'x [m]'
                phase_plot.yaxis.axis_label = 'y [m]'
                phase_source.data = {'x': df['x'], 'y': df['y']}
            if phase_button.active == 3:
                phase_plot.xaxis.axis_label = 'x [m]'
                phase_plot.yaxis.axis_label = 'px [MeV/c]'
                phase_source.data = {'x': df['x'], 'y': df['px']}
            if phase_button.active == 4:
                phase_plot.xaxis.axis_label = 'x [m]'
                phase_plot.yaxis.axis_label = 'py [MeV/c]'
                phase_source.data = {'x': df['x'], 'y': df['py']}
            if phase_button.active == 5:
                phase_plot.xaxis.axis_label = 'y [m]'
                phase_plot.yaxis.axis_label = 'px [MeV/c]'
                phase_source.data = {'x': df['y'], 'y': df['px']}
            if phase_button.active == 6:
                phase_plot.xaxis.axis_label = 'y [m]'
                phase_plot.yaxis.axis_label = 'py [MeV/c]'
                phase_source.data = {'x': df['y'], 'y': df['py']}
            if phase_button.active == 7:
                phase_plot.xaxis.axis_label = 'px [MeV/c]'
                phase_plot.yaxis.axis_label = 'py [MeV/c]'
                phase_source.data = {'x': df['px'], 'y': df['py']}
        if tabs.active == 1:
            if line_button.active == 0:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'rms x [m]'
                line_source.data = {'x': df['z'], 'y': df['rms x']}
            if line_button.active == 1:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'rms y [m]'
                line_source.data = {'x': df['z'], 'y': df['rms y']}
            if line_button.active == 2:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'rms x Emittance [m rad]'
                line_source.data = {'x': df['z'], 'y': df['rms x Emittance']}
            if line_button.active == 3:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'rms y  Emittance [m rad]'
                line_source.data = {'x': df['z'], 'y': df['rms y Emittance']}
            if line_button.active == 4:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'avg Energy [MeV]'
                line_source.data = {'x': df['z'], 'y': df['avg Energy']}
            if line_button.active == 5:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'x β-function [m]'
                line_source.data = {'x': df['z'], 'y': df['x β-function']}
            if line_button.active == 6:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'y β-function [m]'
                line_source.data = {'x': df['z'], 'y': df['y β-function']}
            if line_button.active == 7:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'Centroid x [m]'
                line_source.data = {'x': df['z'], 'y': df['centroid x']}
            if line_button.active == 8:
                line_plot.xaxis.axis_label = 'z [m]'
                line_plot.yaxis.axis_label = 'Centroid y [m]'
                line_source.data = {'x': df['z'], 'y': df['centroid y']}
        if tabs.active == 2:
            if field_button.active == 0:
                field_plot.xaxis.axis_label = 'z [m]'
                field_plot.yaxis.axis_label = 'Ez [MV/m]'
                field_source.data = {'x': df['z'], 'y': df['avg Ez']}
            if field_button.active == 1:
                field_plot.xaxis.axis_label = 'z [m]'
                field_plot.yaxis.axis_label = 'Ex [MV/m]'
                field_source.data = {'x': df['z'], 'y': df['avg Ex']}
            if field_button.active == 2:
                field_plot.xaxis.axis_label = 'z [m]'
                field_plot.yaxis.axis_label = 'Ey [MV/m]'
                field_source.data = {'x': df['z'], 'y': df['avg Ey']}
            if field_button.active == 3:
                field_plot.xaxis.axis_label = 'z [m]'
                field_plot.yaxis.axis_label = 'Bz [T]'
                field_source.data = {'x': df['z'], 'y': df['avg Bz']}
            if field_button.active == 4:
                field_plot.xaxis.axis_label = 'z [m]'
                field_plot.yaxis.axis_label = 'Bx [T]'
                field_source.data = {'x': df['z'], 'y': df['avg Bx']}
            if field_button.active == 5:
                field_plot.xaxis.axis_label = 'z [m]'
                field_plot.yaxis.axis_label = 'By [T]'
                field_source.data = {'x': df['z'], 'y': df['avg By']}

    track_button.on_click(track_handler)
    refresh_button.on_click(refresh_handler)

    return tab
