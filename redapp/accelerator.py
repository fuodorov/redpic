
import redpic as rp
import kenv as kv
import pandas as pd
import numpy as np
from bokeh.models import (Panel, PreText, RadioButtonGroup, Select,
                           FileInput, Button, TextInput, ColumnDataSource,
                          RangeSlider)
from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.models.widgets import Slider
from os.path import dirname, join


def accelerator_tab(acc):
    '''Creating an accelerator with parameters

    '''

    length = RangeSlider(start=0, end=25, value=(0,1), step=.1,
                         title='Length [m]', format='0[.]0')
    pre = PreText(text='''Here you can compile an accelerator.''')
    compile_button = Button(label='Compile', button_type="success", height=45)
    element_button = RadioButtonGroup(
        labels=['Accels', 'Solenoids', 'Quadrupoles'], active=0)

    name_input = TextInput(value="Acc. ", title="Name:", width=100)
    z_input = TextInput(value="0", title="z [m]:", width=100)
    field_input = TextInput(value="0", title="MaxField [MV/m]:", width=100)
    file_input = FileInput(accept='.csv, .dat, .txt')
    parameters_input = column(file_input,
                              row(name_input, z_input, field_input))

    source = ColumnDataSource(data={'z': acc.z, 'Fz': acc.Ez(acc.z)})
    field_plot = figure(x_axis_label='z [m]', y_axis_label='Ez [MV/m]',
                        width=650, height=250)
    field_plot.line('z', 'Fz', source=source, line_color='#3A5785',
                    line_alpha=0.7, line_width=1)

    controls = row(column(pre, length, compile_button),
                   column(element_button, parameters_input))

    # Create a row layout
    layout = column(controls, field_plot)
    tab = Panel(child=layout, title='Accelerator')

    def compile_handler(new, acc=acc):
        (acc.z_start, acc.z_stop) = length.value
        acc.z = acc.parameter = np.arange(acc.z_start, acc.z_stop, acc.dz)
        name = name_input.value
        file_name = dirname(__file__) + '/data/' + file_input.filename
        position = float(z_input.value)
        max_field = float(field_input.value)
        if element_button.active == 0:
            acc.add_accel(name, position, max_field, file_name)
            acc.compile()
            source.data = {'z': acc.z, 'Fz': acc.Ez(acc.z)}
            field_input.title = 'MaxField [MV/m]:'
            field_plot.yaxis.axis_label = 'Ez [MV/m]'
        if element_button.active == 1:
            acc.add_solenoid(name, position, max_field, file_name)
            acc.compile()
            source.data = {'z': acc.z, 'Fz': acc.Bz(acc.z)}
            field_input.title = 'MaxField [T]:'
            field_plot.yaxis.axis_label = 'Bz [T]'
        if element_button.active == 2:
            acc.add_quad(name, position, max_field, file_name)
            acc.compile()
            source.data = {'z': acc.z, 'Fz': acc.Gz(acc.z)}
            field_input.title = 'MaxField [T/m]:'
            field_plot.yaxis.axis_label = 'Gz [T/m]'
        print(acc)

    def element_handler(new):
        if new == 0:
            source.data = {'z': acc.z, 'Fz': acc.Ez(acc.z)}
            field_input.title = 'MaxField [MV/m]:'
            name_input.value = 'Acc. '
            field_plot.yaxis.axis_label = 'Ez [MV/m]'
        if new == 1:
            source.data = {'z': acc.z, 'Fz': acc.Bz(acc.z)}
            field_input.title = 'MaxField [T]:'
            name_input.value = 'Sol. '
            field_plot.yaxis.axis_label = 'Bz [T]'
        if new == 2:
            source.data = {'z': acc.z, 'Fz': acc.Gz(acc.z)}
            field_input.title = 'MaxField [T/m]:'
            name_input.value = 'Quad. '
            field_plot.yaxis.axis_label = 'Gz [T/m]'

    compile_button.on_click(compile_handler)
    element_button.on_click(element_handler)

    return tab
