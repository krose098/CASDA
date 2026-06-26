from casatasks import imfit, imstat
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import click
import csv
import sys
from io import StringIO
from utils import unicoord, extract_fits_header

def flux_extractor_to_dict(path, coords, stokes, radius, rms_radius, scale, units='mJy'):
    unit_scale = 1e6 if units == 'uJy' else 1e3
    unit_label = 'uJy' if units == 'uJy' else 'mJy'

    header = extract_fits_header(path)

    result = {
        'Date of Observation (UTC)': header.get('DATE-OBS'),
        'Project': header.get('PROJECT'),
        'Freq (GHz)': header.get('CRVAL3'),
        f'S ({unit_label})': None,
        f'Err ({unit_label})': None,
        f'rms ({unit_label}/beam)': None,
        'SBID': header.get('SBID'),
        'Duration': header.get('DURATION'),
        'file': path,
        'ra_fit_deg': None,
        'ra_err_deg': None,
        'dec_fit_deg': None,
        'dec_err_deg': None,
        'fit_failed': False,
    }

    hmsdms = coords.to_string('hmsdms', sep=':', precision=1, pad=True)
    ra_str, dec_str = hmsdms.split()
    hh, mm, ss = ra_str.split(':')
    ra_casa = f"{hh}h{mm}m{ss}s"
    sign = '+' if dec_str.strip()[0] == '+' else '-'
    dd, dm, ds = dec_str.replace('+','').replace('-','').split(':')
    dec_casa = f"{sign}{dd}.{dm}.{ds[:4]}"
    region_temp = f"{ra_casa},{dec_casa}"
    region = 'circle[[{}], {}arcsec]'.format(region_temp, radius)
    annulus_region = 'annulus[[{}], [{}arcsec, {}arcsec]]'.format(region_temp, 2*radius, rms_radius)

    source_name = "J"+region_temp[0:2]+region_temp[3:5]
    name = f"{source_name}_{stokes}_{radius}_arcsec"
    summary = 'fit_summary_{}'.format(name)
    logfile = 'fit_log_{}'.format(name)

    try:
        fit = imfit(
            imagename=path,
            region=region,
            chans='',
            stokes="",
            summary=summary,
            logfile=logfile,
        )
        fit_failed = not fit.get('converged', [False])[0]
    except Exception as e:
        fit = None
        fit_failed = True

    stats = imstat(imagename=path, stokes=stokes, region=annulus_region)
    rms = stats["rms"] * unit_scale 
    result[f'rms ({unit_label}/beam)'] = round(rms[0], 4)
    result['fit_failed'] = fit_failed

    if not fit_failed and fit is not None:
        temp = fit['results']['component0']['peak']
        peak_flux = temp['value'] * unit_scale
        peak_err = temp['error'] * unit_scale
        combined_err = (peak_err**2 + rms**2 + (scale * peak_flux)**2)**(0.5)

        fit_coords = fit['results']['component0']['shape']['direction']
        ra = fit_coords['m0']['value'] * 180 / np.pi % 360
        ra_err = fit_coords['error']['longitude']['value'] / 3600
        dec = fit_coords['m1']['value'] * 180 / np.pi
        dec_err = fit_coords['error']['latitude']['value'] / 3600

        result.update({
            f'S ({unit_label})': round(peak_flux, 4),
            f'Err ({unit_label})': round(combined_err[0], 4),
            'ra_fit_deg': round(ra, 6),
            'ra_err_deg': round(ra_err, 6),
            'dec_fit_deg': round(dec, 6),
            'dec_err_deg': round(dec_err, 6),
        })

    return result


@click.group()
def cli():
    pass


@cli.command()
@click.option("-g", "--galactic", is_flag=True, default=False, help="If set, assumes galactic coordinates")
@click.option("-u", "--units", type=click.Choice(['mJy', 'uJy']), default='mJy', show_default=True, help="Flux density units")
@click.argument("path", nargs=1, type=str)
@click.argument("coords", nargs=1, type=str)
@click.argument("stokes", nargs=1, type=str, default="I")
@click.argument("radius", nargs=1, type=float, default=5)
@click.argument("rms_radius", nargs=1, type=float, default=30)
@click.argument("scale", nargs=1, type=float, default=0.1)
def single(path, coords, stokes, radius, rms_radius, scale, galactic,units):
    pos_eq, __ = unicoord(coords, galactic, display=False)
    flux_extractor(path, pos_eq, stokes, radius, rms_radius, scale)

@cli.command()
@click.option("-g", "--galactic", is_flag=True, default=False, help="If set, assumes galactic coordinates")
@click.option("-u", "--units", type=click.Choice(['mJy', 'uJy']), default='mJy', show_default=True, help="Flux density units")
@click.option("-o", "--output", default="results.csv", show_default=True, help="Output CSV file")
@click.argument("paths", nargs=-1, type=str, required=True)
@click.argument("coords", nargs=1, type=str)
@click.argument("stokes", nargs=1, type=str, default="I")
@click.argument("radius", nargs=1, type=float, default=5)
@click.argument("rms_radius", nargs=1, type=float, default=30)
@click.argument("scale", nargs=1, type=float, default=0.1)
def batch(paths, coords, stokes, radius, rms_radius, scale, galactic, output, units):
    pos_eq, __ = unicoord(coords, galactic, display=False)
    rows = []
    for path in paths:
        print(f"Processing {path}...")
        row = flux_extractor_to_dict(path, pos_eq, stokes, radius, rms_radius, scale, units=units)
        rows.append(row)

    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved {len(rows)} rows to {output}")


if __name__ == "__main__":
    cli()
