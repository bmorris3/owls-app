import os
from pathlib import Path
import solara
import re
import numpy as np
from echo import delay_callback
import s3fs

import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy.timeseries import LombScargle

import ipygoldenlayout
import ipysplitpanes
import ipyvue

import jdaviz
from jdaviz import Specviz
from jdaviz.app import custom_components

from lcviz import LCviz
from lightkurve import LightCurve

from aesop import EchelleSpectrum
from specutils import Spectrum1D
import pandas as pd


s3 = s3fs.S3FileSystem()
bucket = 'owls-spectra'
owls_bucket = list(filter(lambda x: '__' in x and '.fits' in x, s3.ls(bucket)))
spectra_paths = [f's3://{prefix}' for prefix in owls_bucket if 'fits.gz' in prefix]
standard_path = f's3://{bucket}/HZ44.0016.wfrmcpc.fits'
standard_spectrum = EchelleSpectrum.from_fits(standard_path)

mwo_path = f's3://{bucket}/s-indices/mwo_v1995.pkl'
mwo_v1995 = pd.read_pickle(mwo_path)

owls_path = f's3://{bucket}/s-indices/owls_paper1_sindices_2025-02-07.pkl'
owls = pd.read_pickle(owls_path)


hd_target_names_owls = {
    # Strip non-numeric characters, force to integer
    int(re.sub("[^0-9]", '', target)): target
    # for every target name in the owls measurements
    for target in owls.index.unique()
    # if the target name starts with HD
    if (target.startswith("HD") or target.startswith('hd')) and
    int(re.sub("[^0-9]", '', target)) in mwo_v1995.index.get_level_values("Star")
}


def to_spectrum1d(spec, meta=None):
    return Spectrum1D(
        flux=spec.flux[~spec.mask] / spec.flux[~spec.mask].max() * u.count,
        spectral_axis=spec.wavelength[~spec.mask],
        meta=meta,
    )


available_targets = sorted(
    set([
        os.path.basename(path).split('__')[0]
        for path in
        owls_bucket
    ])
)


def update_specviz(selected_paths, selected_orders):
    specviz = Specviz()

    # remove "import data" button
    specviz.app.state.tool_items = list(specviz.app.state.tool_items)[1:]

    # prevent saving files on the server
    specviz.plugins['Export']._obj.serverside_enabled = False

    # close plugin tray
    specviz.app.state.drawer = False

    for target_path in selected_paths:
        target_spectrum = EchelleSpectrum.from_fits(target_path)
        target_spectrum.continuum_normalize_from_standard(
            standard_spectrum, 5, only_orders=selected_orders
        )

        datetime = target_spectrum.header['DATE-OBS']
        time = Time(datetime, format='isot')
        skycoord = SkyCoord(
            ra=target_spectrum.header['RA'],
            dec=target_spectrum.header['DEC'],
            unit=(u.hourangle, u.deg)
        )
        target_spectrum.barycentric_correction(
            time=time, skycoord=skycoord, location=EarthLocation.of_site("APO")
        )
        for i in selected_orders:
            order_i = to_spectrum1d(target_spectrum[i], meta=target_spectrum.header)
            time_target_label = f"{os.path.basename(target_path)}"
            specviz.load_data(
                order_i,
                data_label=f'Order {i} ({time_target_label})'
            )

    viewer = specviz.app.get_viewer('spectrum-viewer')

    # skip these steps if no files selected
    if len(selected_paths):
        all_orders = list(range(len(target_spectrum)))
        order_labels = [
            f'{i} ({target_spectrum[i].wavelength.value.min():.0f} - {target_spectrum[i].wavelength.value.max():.0f} Å) '
            for i in all_orders
        ]

        with delay_callback(viewer.state, 'x_min', 'x_max', 'y_min', 'y_max'):
            viewer.state.x_min = min([
                layer.state.layer.get_component('World 0').data.min()
                for layer in viewer.layers
            ])
            viewer.state.x_max = max([
                layer.state.layer.get_component('World 0').data.max()
                for layer in viewer.layers
            ])
            viewer.state.y_min = 0
            viewer.state.y_max = max(
                [layer.state.v_max for layer in viewer.layers]
            )

        colors = [to_hex(plt.cm.viridis(x)) for x in np.linspace(0, 0.9, len(selected_paths))]

        for start in [0, 1]:
            for layer, color in zip(viewer.layers[start::2], colors):
                layer.state.color = color
    else:
        order_labels = []

    return specviz, order_labels


def update_lcviz(target_name):

    target_in_mwo = (
            target_name.startswith('HD') and
            int(target_name.split("_")[1]) in hd_target_names_owls.keys()
    )

    if target_in_mwo:
        mwo_measurements = mwo_v1995.loc[int(target_name.split("_")[1])]

        mwo_lc = LightCurve(
            time=Time(mwo_measurements.index.get_level_values('Date')).jd,
            flux=mwo_measurements['S'].astype(float)
        )

        ls = LombScargle(
            mwo_lc.time, mwo_lc.flux, mwo_lc.flux.std() / 2,
            normalization='psd',
            nterms=3
        )

        min_period = (5 * u.year).to(u.d)
        max_period = (15 * u.year).to(u.d)

        freq = np.geomspace(1 / max_period, 1 / min_period, 1000)
        power = ls.power(freq)
        freq_at_max_power = freq[np.argmax(power)]
        period_at_max_power = (1 / freq_at_max_power).value


    owls_measurements = owls.loc[target_name.replace("_", " ")]

    owls_lc = LightCurve(
        time=Time(np.atleast_1d(owls_measurements['owls_time']), format='jd'),
        flux=np.atleast_1d(owls_measurements['owls_s_mwo']).astype(float),
        flux_err=np.atleast_1d(owls_measurements['owls_s_mwo_err']).astype(float)
    )

    if target_in_mwo:
        time_model = np.linspace(mwo_lc.time.min(), owls_lc.time.max(), 500)

        model = ls.model(time_model, freq_at_max_power)
        model_lc = LightCurve(time=time_model, flux=model)

    lcviz = LCviz()

    # remove import data button
    lcviz.app.state.tool_items.pop(0)

    # prevent saving files on the server
    lcviz.plugins['Export']._obj.serverside_enabled = False

    # close plugin tray
    lcviz.app.state.drawer = False

    if target_in_mwo:
        freq_analysis = lcviz.plugins['Frequency Analysis']
        freq_analysis.auto_range = False
        freq_analysis.xunit = 'period'
        freq_analysis.minimum = min_period.value
        freq_analysis.maximum = max_period.value

    if target_in_mwo:
        lcviz.load_data(mwo_lc, data_label='MWO')

    lcviz.load_data(owls_lc, data_label='OWLS')

    if target_in_mwo:
        lcviz.load_data(model_lc, data_label='model')

        ephem = lcviz.plugins['Ephemeris']
        ephem.period = period_at_max_power

    plot_opts = lcviz.plugins['Plot Options']

    viewers = plot_opts.viewer.choices
    data_labels = ['MWO', 'OWLS', 'model']
    marker_size_scale = [1, 20, 0.5]
    marker_color = ['gray', 'r', 'b']
    marker_opacity = [1, 1, 0.2]

    for viewer in viewers:
        for layer, mss, mc, mo in zip(
                data_labels,
                marker_size_scale,
                marker_color,
                marker_opacity
        ):
            if layer in ['MWO', 'model'] and not target_in_mwo:
                continue

            plot_opts.viewer = viewer
            plot_opts.layer = layer
            plot_opts.marker_size_scale = mss
            plot_opts.marker_color = mc
            plot_opts.marker_opacity = mo

    return lcviz


@solara.component
def Page():
    ipysplitpanes.SplitPanes()
    ipygoldenlayout.GoldenLayout()
    for name, path in custom_components.items():
        ipyvue.register_component_from_file(None, name,
                                            os.path.join(os.path.dirname(jdaviz.__file__), path))

    ipyvue.register_component_from_file('g-viewer-tab', "container.vue", jdaviz.__file__)

    solara.Style(Path(__file__).parent / "solara.css")

    target_name, set_target_name = solara.use_state('HD_81809')

    def paths_for_target(target_name):
        return sorted(
            [path for path in spectra_paths if target_name in path]
        )

    paths = paths_for_target(target_name)
    maximum_files = 1
    selected_paths, set_selected_paths = solara.use_state(paths[:maximum_files])
    selected_orders, set_selected_orders = solara.use_state([16, 17])

    solara.Markdown("# OWLS – the Olin Wilson Legacy Survey")

    specviz, order_labels = update_specviz(selected_paths, selected_orders)

    with solara.Column(align='center'):
        with solara.Row():
            solara.display(specviz.app)
        with solara.Row():
            lcviz = update_lcviz(target_name)
            solara.display(lcviz.app)

        with solara.Row():
            with solara.Columns([1, 1]):
                with solara.Column(align='start'):
                    solara.Markdown("### Target")

                    def on_target_change(target_name):
                        set_target_name(target_name)
                        add_paths = paths_for_target(target_name)
                        set_selected_paths(selected_paths + add_paths[0:1])

                    solara.Select('Target', list(available_targets), target_name, on_target_change)
                    solara.SelectMultiple(
                        'Observations', selected_paths, paths, set_selected_paths, dense=True,
                    )

                    def on_all_obs():
                        set_selected_paths(paths)

                    solara.Button("Select all obs.", on_all_obs)

                    def on_orders_changed(labels, order_labels=order_labels):
                        set_selected_orders([order_labels.index(label) for label in labels])

                    def on_hk_orders_only():
                        set_selected_orders([16, 17])

                with solara.Column(align='start'):
                    solara.Markdown("### Spectral orders")
                    solara.SelectMultiple(
                        'Spectral orders',
                        [order_labels[i] for i in selected_orders],
                        order_labels,
                        on_value=on_orders_changed,
                        dense=True
                    )
                    solara.Button("Select H & K orders only", on_hk_orders_only)
