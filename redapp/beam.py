import redpic as rp
import pandas as pd
import numpy as np
from bokeh.models import (Panel, PreText, RadioButtonGroup, Select,
                           FileInput, Button, TextInput, ColumnDataSource)
from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.models.widgets import Slider
from os.path import dirname, join

def beam_tab(beam):
    '''Creating a beam with parameters

    '''

    pre = PreText(text='''Here you can generate or upload a beam.''')

    species_button = RadioButtonGroup(
        labels=['electron', 'positron', 'proton', 'antiproton'], active=0)

    select_quantity = Slider(start=1_000, end=100_000,
                             step=1_000, value=1_000,
                             title='Number of particles')

    select_current = Slider(start=0, end=5_000,
                            step=100, value=0,
                            title='Curent [A]')
    # Phase ellipse
    x_input = TextInput(value="0", title="X [m]:", width=75)
    y_input = TextInput(value="0", title="Y [m]:", width=75)
    z_input = TextInput(value="0", title="Z [m]:", width=75)

    px_input = TextInput(value="0", title="Px [MeV/c]:", width=75)
    py_input = TextInput(value="0", title="Py [MeV/c]:", width=75)
    pz_input = TextInput(value="0", title="Pz [MeV/c]:", width=75)

    x_off_input = TextInput(value="0", title="X_off [m]:", width=75)
    y_off_input = TextInput(value="0", title="Y_off [m]:", width=75)
    px_off_input = TextInput(value="0", title="Px_off:", width=60)
    py_off_input = TextInput(value="0", title="Py_off:", width=60)

    phase_ell = column(row(x_input, y_input, z_input, x_off_input),
                       row(px_input, py_input, pz_input, y_off_input))

    select_distribution = RadioButtonGroup(labels=['Uniform', 'Gauss'],
                                           active=0)

    file_input = FileInput(accept='.csv, .ini')
    generate_button = Button(label='Generate', button_type="success",
                             width=155)
    upload_button = Button(label='Upload', button_type="warning",
                           width=155)

    source = ColumnDataSource(data={'x': beam.df['x'], 'y': beam.df['y']})
    p = figure(x_axis_label='x [m]', y_axis_label='y [m]',
               width=410, height=400)
    p.scatter('x', 'y', source=source, size=0.5, color='#3A5785', alpha=0.5)

    controls = column(pre, species_button, select_quantity,
                      select_distribution, select_current,
                      phase_ell, file_input,
                      row(upload_button, generate_button))
    # Create a row layout
    layout = row(controls, p)

    tab = Panel(child=layout, title='Beam')

    def generate_handler(new, beam=beam):
        I = select_current.value
        N = select_quantity.value
        X = float(x_input.value)
        Y = float(y_input.value)
        Z = float(z_input.value)
        Px = float(px_input.value)
        Py = float(py_input.value)
        Pz = float(pz_input.value)
        X_off = float(x_off_input.value)
        Y_off = float(y_off_input.value)
        Px_off = float(px_off_input.value)
        Py_off = float(py_off_input.value)
        if select_distribution.active == 0:
            distribution = rp.Distribution(name='KV', x=X, y=Y, z=Z,
                                           px=Px, py=Py, pz=Pz)
        if select_distribution.active == 1:
            distribution = rp.Distribution(name='GA', x=X, y=Y, z=Z,
                                           px=Px, py=Py, pz=Pz)
        if species_button.active == 0:
            species = rp.electron
        if species_button.active == 1:
            species = rp.positron
        if species_button.active == 2:
            species = rp.proton
        if species_button.active == 3:
            species = rp.antiproton
        beam.type = species
        Q = np.sign(species.charge) * I * Z / rp.c
        beam.generate(distribution, n=N, charge=Q, path=dirname(__file__) + '/data/',
                      x_off=X_off, y_off=Y_off, px_off=Px_off, py_off=Py_off)
        source.data = {'x': beam.df['x'], 'y': beam.df['y']}
        print(beam)

    def upload_handler(new, beam=beam):
        I = select_current.value
        if species_button.active == 0:
            species = rp.electron
        if species_button.active == 1:
            species = rp.positron
        if species_button.active == 2:
            species = rp.proton
        if species_button.active == 3:
            species = rp.antiproton
        beam.type = species
        Q = np.sign(species.charge)* I * (beam.df['z'].max()-beam.df['z'].min()) / rp.c
        beam.upload(dirname(__file__) + '/data/' + file_input.filename, charge=Q,
                    path=dirname(__file__) + '/data/')
        source.data = {'x': beam.df['x'], 'y': beam.df['y']}
        x_input.value = str(np.around(beam.df['x'].max(), 3))
        y_input.value = str(np.around(beam.df['y'].max(), 3))
        z_input.value = str(np.around(beam.df['z'].max(), 3) * 2)
        px_input.value = str(np.around(beam.df['px'].max(), 3))
        py_input.value = str(np.around(beam.df['py'].max(), 3))
        pz_input.value = str(np.around(beam.df['pz'].max(), 3))
        select_quantity.value = beam.n
        print(beam)

    generate_button.on_click(generate_handler)
    upload_button.on_click(upload_handler)

    return tab
