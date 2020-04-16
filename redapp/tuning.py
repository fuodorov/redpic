
import redpic as rp
import kenv as kv
import pandas as pd
import numpy as np
from bokeh.models import (Panel, PreText, RadioButtonGroup, Select,
 						  FileInput, Button, TextInput, ColumnDataSource,
						  RangeSlider, Div)
from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.models.widgets import Slider
from os.path import dirname, join

def tuning_tab(beam, acc):
    '''Tuning accelerator

	'''
    pre = PreText(text='''You can configure accelerator elements.''')
    element_button = RadioButtonGroup(
        labels=['Accels', 'Solenoids', 'Quadrupoles'], active=0)
    tune_button = Button(label='Tune', button_type="success")
    select = Select(options=[itm.name for itm in acc.Ez_beamline.values()])
    select_maxfield = Slider(start=-5.0, end=5.0,
				    		 step=0.1, value=0.0,
							 title='MaxField [MV/m]', format='0[.]0')

    field_source = ColumnDataSource(data={'z': acc.z, 'Fz': acc.Ez(acc.z)})
    field_plot = figure(x_axis_label='z [m]', y_axis_label='Ez [MV/m]',
                        width=330, height=210)
    field_plot.line('z', 'Fz', source=field_source, line_color='#3A5785',
                    line_alpha=0.7, line_width=1)

    envelope_plot = figure(x_axis_label='z [m]', y_axis_label='Envelope, x [m]',
                        width=650, height=250)
    envelope_source = ColumnDataSource(data={'z': acc.z, 'x': 0*acc.z,
                                             '-x': 0*acc.z})
    envelope_plot.line('z', 'x', source=envelope_source, line_color='#3A5785',
                        line_alpha=0.7, line_width=1)
    envelope_plot.line('z', '-x', source=envelope_source, line_color='#3A5785',
                        line_alpha=0.7, line_width=1)

    controls = row(column(pre, element_button, select, select_maxfield, tune_button),
                   field_plot)
    tab = Panel(child=column(controls, envelope_plot), title='Tuning')

    def element_handler(new, acc=acc):
        if new == 0:
            select.options=[itm.name for itm in acc.Ez_beamline.values()]
            select_maxfield.start = -5.0
            select_maxfield.end = 5.0
            select_maxfield.step = 0.1
            if select.options == []:
                select_maxfield.value = 0.0
                select.value = ''
            else:
                select_maxfield.value = acc.Ez_beamline[select.options[0]].max_field
                select.value = [itm.name for itm in acc.Ez_beamline.values()][0]
            select_maxfield.format='0[.]0'
            select_maxfield.title = 'MaxField [MV/m]'
            field_source.data = {'z': acc.z, 'Fz': acc.Ez(acc.z)}
            field_plot.yaxis.axis_label = 'Ez [MV/m]'
        if new == 1:
            select.options=[itm.name for itm in acc.Bz_beamline.values()]
            select_maxfield.start = 0.000
            select_maxfield.end = 0.500
            select_maxfield.step = 0.001
            if select.options == []:
                select_maxfield.value = 0.000
                select.value = ''
            else:
                select_maxfield.value = acc.Bz_beamline[select.options[0]].max_field
                select.value = [itm.name for itm in acc.Bz_beamline.values()][0]
            select_maxfield.format='0[.]000'
            select_maxfield.title = 'MaxField [T]'
            field_source.data = {'z': acc.z, 'Fz': acc.Bz(acc.z)}
            field_plot.yaxis.axis_label = 'Bz [T]'
        if new == 2:
            select.options=[itm.name for itm in acc.Gz_beamline.values()]
            select_maxfield.start = -10.0
            select_maxfield.end = 10.0
            select_maxfield.step = 0.1
            if select.options == []:
                select_maxfield.value = 0.0
                select.value = ''
            else:
                select_maxfield.value = acc.Gz_beamline[select.options[0]].max_field
                select.value = [itm.name for itm in acc.Gz_beamline.values()][0]
            select_maxfield.format='0[.]0'
            select_maxfield.title = 'MaxField [T/m]'
            field_source.data = {'z': acc.z, 'Fz': acc.Gz(acc.z)}
            field_plot.yaxis.axis_label = 'Gz [T/m]'

    def update_select_maxfield(attrname, old, new, acc=acc):
        if new != '':
            if element_button.active == 0:
                select_maxfield.value = acc.Ez_beamline[new].max_field
            if element_button.active == 1:
                select_maxfield.value = acc.Bz_beamline[new].max_field
            if element_button.active == 2:
                select_maxfield.value = acc.Gz_beamline[new].max_field

    def tune_handler(new, beam=beam, acc=acc):
        if select.value != '':
            if element_button.active == 0:
                acc.Ez_beamline[select.value].max_field = select_maxfield.value
                acc.compile()
                field_source.data = {'z': acc.z, 'Fz': acc.Ez(acc.z)}
            if element_button.active == 1:
                acc.Bz_beamline[select.value].max_field = select_maxfield.value
                acc.compile()
                field_source.data = {'z': acc.z, 'Fz': acc.Bz(acc.z)}
            if element_button.active == 2:
                acc.Gz_beamline[select.value].max_field = select_maxfield.value
                acc.compile()
                field_source.data = {'z': acc.z, 'Fz': acc.Gz(acc.z)}
        beam.df['pz2'] = beam.df['pz']*beam.df['pz']
        kv_beam = kv.Beam(energy=(beam.df['pz2'].mean() +
                          beam.type.mass**2*rp.c**4/rp.e**2/1e6**2)**0.5,
                         current=beam.charge * rp.c,
                         radius_x=beam.df['x'].max(),
                         radius_y=beam.df['y'].max(),
                         radius_xp=beam.df['px'].max() / beam.df['pz'].max(),
                         radius_yp=beam.df['py'].max() / beam.df['pz'].max(),
                         normalized_emittance=1000e-6)
        kv_sim = kv.Simulation(kv_beam, acc)
        kv_sim.track()
        envelope_source.data = {'z': acc.z, 'x': kv_sim.envelope_x(acc.z),
                                '-x': -kv_sim.envelope_x(acc.z)}

    element_button.on_click(element_handler)
    select.on_change('value', update_select_maxfield)
    tune_button.on_click(tune_handler)

    return tab
