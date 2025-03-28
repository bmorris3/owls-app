import os
from pathlib import Path
import solara
import re
import numpy as np
from echo import delay_callback
import json

import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

import astropy.units as u
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
from astropy.utils.data import download_files_in_parallel


pkgname = 'owls-app'

urls = dict(
    manifest='https://stsci.box.com/shared/static/3duv3ywn7fzx75mybk4k7d6533f1sjj4.json',
    standard='https://stsci.box.com/shared/static/9cfsb23t644gbzh5b936wennypwzql2g.fits',
    mwo='https://stsci.box.com/shared/static/6brcqyfb0y1kx0phub1f6his9rcz7c26.pkl',
    owls='https://stsci.box.com/shared/static/72jp0im3qcs26etr8j4fiuq3cyxcfuxv.pkl',
)


def download(urls, **kwargs):
    if isinstance(urls, str):
        urls = [urls]

    return download_files_in_parallel(urls, pkgname=pkgname, cache=True, show_progress=False)


downloaded_paths = download(urls.values())
paths = {k: v for k, v in zip(urls.keys(), downloaded_paths)}

owls_latest = json.load(open(paths['manifest'], 'r'))
standard_spectrum = EchelleSpectrum.from_fits(paths['standard'])
mwo_v1995 = pd.read_pickle(paths['mwo'])
owls = pd.read_pickle(urls['owls'])
specviz = None
lcviz = None
order_labels = None

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
        obs['target'] for obs in owls_latest
    ])
)


def update_specviz(selected_paths, selected_times, selected_orders):
    global specviz, order_labels

    if specviz is None:
        specviz = Specviz()

        # close plugin tray
        specviz.app.state.drawer = False

    target_spectrum = None

    for target_path, spectrum_time in zip(selected_paths, selected_times):
        loaded_paths = [data.meta.get('path') for data in specviz.app.data_collection]
        loaded_orders = {k: list() for k in loaded_paths}

        for path in loaded_paths:
            for data in specviz.app.data_collection:
                loaded_orders[path].append(data.meta.get('order'))

        if target_path in loaded_paths and len(loaded_orders[path]) == selected_orders:
            continue

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
            data_label = f'{spectrum_time} (Order {i})'
            specviz.load_data(order_i, data_label=data_label)
            data = specviz.app.data_collection[data_label]
            data.meta['path'] = target_path
            data.meta['order'] = i

    # remove loaded datasets that have been deselected
    data_to_remove = []
    for data in specviz.app.data_collection:
        if data.meta.get('path') not in selected_paths or data.meta.get('order') not in selected_orders:
            data_to_remove.append(data)
    for data in data_to_remove:
        specviz.app.data_collection.remove(data)

    viewer = specviz.app.get_viewer('spectrum-viewer')

    # skip these steps if no files selected
    if len(selected_paths):
        if target_spectrum is not None:
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

        for layer in viewer.layers:
            layer.state.as_steps = False
            for i, path in enumerate(selected_paths):
                if layer.layer.meta['path'] == path:
                    layer.state.color = colors[i]


def update_lcviz(target_name, owls_measurements):
    global lcviz
    target_in_mwo = (
            target_name.startswith('HD') and
            int(target_name.split()[1]) in hd_target_names_owls.keys()
    )

    if target_in_mwo:
        mwo_measurements = mwo_v1995.loc[int(target_name.split()[1])]

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

    owls_lc = LightCurve(
        time=Time(np.atleast_1d(owls_measurements['owls_time']), format='jd'),
        flux=np.atleast_1d(owls_measurements['owls_s_mwo']).astype(float),
        flux_err=np.atleast_1d(owls_measurements['owls_s_mwo_err']).astype(float)
    )

    if target_in_mwo:
        time_model = np.linspace(mwo_lc.time.min(), owls_lc.time.max(), 500)

        model = ls.model(time_model, freq_at_max_power)
        model_lc = LightCurve(time=time_model, flux=model)

    if lcviz is None:
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

        lcviz.app.get_viewer(viewer).state.y_axislabel = 'S-index'

    return lcviz


@solara.component
def Page():
    global specviz, order_labels, lcviz
    ipysplitpanes.SplitPanes()
    ipygoldenlayout.GoldenLayout()
    for name, path in custom_components.items():
        ipyvue.register_component_from_file(None, name,
                                            os.path.join(os.path.dirname(jdaviz.__file__), path))

    ipyvue.register_component_from_file('g-viewer-tab', "container.vue", jdaviz.__file__)

    solara.Style(Path(__file__).parent / "solara.css")

    target_name, set_target_name = solara.use_state('HD 81809')

    def urls_for_target(target_name):
        return sorted(
            [obs['url'] for obs in owls_latest if obs['target'] == target_name]
        )

    def times_for_target(target_name):
        return sorted(
            [obs['datetime'] for obs in owls_latest if obs['target'] == target_name]
        )

    def paths_for_target(target_name):
        urls = urls_for_target(target_name)
        return download(urls)

    paths = paths_for_target(target_name)
    times = times_for_target(target_name)
    maximum_files = 1
    selected_paths, set_selected_paths = solara.use_state(paths[:maximum_files])
    selected_times, set_selected_times = solara.use_state(times[:maximum_files])
    selected_orders, set_selected_orders = solara.use_state([16, 17])

    target_name_owls = target_name.replace("_", " ")
    if target_name.startswith("TOI"):
        # for now, the owls table has the ".01" suffix, but it shouldn't
        target_name_owls = target_name_owls + ".01"
    elif target_name == 'GJ 29':
        # catch special case
        target_name_owls = 'GJ 29.1'

    owls_measurements = owls.loc[target_name_owls]

    if specviz is None:
        update_specviz(selected_paths, selected_times, selected_orders)
    if lcviz is None:
        update_lcviz(target_name, owls_measurements)

    solara.Markdown("# OWLS – The Olin Wilson Legacy Survey")

    with solara.Column(align='center'):
        with solara.Row():
            with solara.Columns([1, 1, 2]):
                with solara.Column(align='start'):
                    solara.Markdown("### Target")

                    def on_target_change(target_name):
                        global specviz, lcviz
                        urls = urls_for_target(target_name)
                        add_path = download(urls[0])
                        add_times = times_for_target(target_name)
                        specviz = lcviz = None  # force init of lcviz, specviz
                        update_specviz(add_path, add_times[:1], selected_orders)
                        set_target_name(target_name)
                        set_selected_paths(add_path)
                        set_selected_times(add_times[:1])

                    def on_times_change(times):
                        urls = [obs['url'] for obs in owls_latest if obs['datetime'] in times]
                        add_paths = download(urls)
                        update_specviz(add_paths, times, selected_orders)
                        set_selected_paths(add_paths)
                        set_selected_times(times)

                    solara.Select('Target', list(available_targets), target_name, on_target_change)
                    solara.SelectMultiple(
                        'Observations', selected_times, times, on_times_change, dense=True,
                    )

                    def on_all_obs():
                        all_times = [obs['datetime'] for obs in owls_latest if obs['target'] == target_name]
                        on_times_change(all_times)

                    solara.Button("Select all obs.", on_all_obs)

                    def on_orders_changed(labels, order_labels=order_labels):
                        new_orders = [order_labels.index(label) for label in labels]
                        update_specviz(selected_paths, selected_times, new_orders)
                        set_selected_orders(new_orders)

                    def on_hk_orders_only():
                        set_selected_orders([16, 17])
                        update_specviz(selected_paths, selected_times, selected_orders)

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

                if isinstance(owls_measurements, pd.DataFrame):
                    with solara.Column(align='start'):
                        solara.Markdown("### OWLS measurements")
                        solara.DataFrame(owls_measurements)


        with solara.Row():
            solara.display(specviz.app)
        with solara.Row():
            solara.display(lcviz.app)
