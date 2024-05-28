from .imports import *
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.table import QTable
from astropy.visualization import quantity_support
from gaiaxpy import calibrate
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
import lightkurve as lk
import os
import warnings 

# for information about available queries + columns, see 
# https://astroquery.readthedocs.io/en/latest/gaia/gaia.html

# set a high row limit to allow lots of stars in crowded fields
Gaia.ROW_LIMIT = 50000
Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"
gaia_epoch = 2016.0

# Gaia filter transformations from
# https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
terms = {}
terms["G-r_sloan"] = [-0.12879, 0.24662, -0.027464, -0.049465]
terms["G-i_sloan"] = [-0.29676, 0.64728, -0.10141]
terms["G-g_sloan"] = [0.13518, -0.46245, -0.25171, 0.021349]
terms["G-V_johnsoncousins"] = [-0.01760, -0.006860, -0.1732]
terms["G-R_johnsoncousins"] = [-0.003226, 0.3833, -0.1345]
terms["G-I_johnsoncousins"] = [0.02085, 0.7419, -0.09631]

uncertainties = {}
uncertainties["G-r_sloan"] = 0.066739
uncertainties["G-i_sloan"] = 0.98957
uncertainties["G-g_sloan"] = 0.16497
uncertainties["G-V_johnsoncousins"] = 0.0458
uncertainties["G-R_johnsoncousins"] = 0.048
uncertainties["G-I_johnsoncousins"] = 0.049


def estimate_other_filters(table):
    """
    Take a table of simple Gaia photometry
    from `get_gaia_data` and use
    color transformations to estimate the
    magnitudes in other filters.

    Data from:
    https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html

    Parameters
    ----------
    table : QTable
        Table with Gaia photometry (see `get_gaia_data`).

    Returns
    -------
    table : QTable
        Table updated with estimated filter magnitudes.
    """

    # define the inputs for color transformation functions
    G = table["G_gaia_mag"]
    BP_minus_RP = table["BP_gaia_mag"] - table["RP_gaia_mag"]

    # loop through filters
    for color, coefficients in terms.items():

        # calculate the polynomial
        G_minus_x = (
            np.sum([c * BP_minus_RP**i for i, c in enumerate(coefficients)], 0)
            * u.mag
        )

        # store in table
        filter_name = color.split("-")[-1]
        table[f"{filter_name}_mag"] = G - G_minus_x
    return table


def get_gaia(center, radius=6 * u.arcmin):
    """
    Get photometry and basic data based on Gaia EDR3.

    Use Gaia Early Data Release 3 to download data
    for all stars within a particular radius of a
    particular center. Gaia is an all-sky space survey;
    basically any star you can see with a moderate
    aperture ground-based telescope has been measured
    by Gaia, so it can provide a useful reference.
    Gaia positions, motions, and photometry will be
    included, along with magnitudes in other filters
    estimated with `estimate_other_filters` via
    Gaia color transformations.

    Parameters
    ----------
    center : SkyCoord
        An astropy SkyCoord object indicating the
        right ascension and declination of the center.
        This center can be created in a few ways:
        ```
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        center = SkyCoord(ra=123.45*u.deg, dec=67.89*u.deg)
        center = SkyCoord(ra=123.45*u.deg, dec=67.89*u.deg)
        other_center = SkyCoord.from_name('Vega')
        ```
    radius : astropy.units.Quantity
        The angular radius around which the star to
        search for photometry. Default is 6 arcminutes.

    Returns
    -------
    table : astropy.table.Table
        An astropy table containing the results,
        with columns for different coordinates
        or filters, and rows for different stars.
    """


    if isinstance(center, SkyCoord):
        center_skycoord = center
    elif isinstance(center, str):
        center_skycoord = SkyCoord.from_name(center)
    else:
        raise ValueError('Sorry! `center` must be a string or SkyCoord.')

    # get the data from the archive
    job = Gaia.cone_search_async(
        center_skycoord,
        radius=radius,
        columns=[
            "source_id",
            "ra",
            "dec",
            "phot_g_mean_mag",
            "phot_rp_mean_mag",
            "phot_bp_mean_mag",
            "parallax", "parallax_error",
            "pmra", "pmra_error",
            "pmdec", "pmdec_error"
        ],
    )
    results = job.get_results()

    # tidy up the table and convert to quantities
    table = QTable(results)
    table.rename_columns(
        ["phot_g_mean_mag", "phot_rp_mean_mag", "phot_bp_mean_mag", "dist"],
        ["G_gaia_mag", "RP_gaia_mag", "BP_gaia_mag", "distance_from_center"],
    )

    # add unit to the distance from the field center
    table["distance_from_center"].unit = u.deg

    # populate with other estimated filter magnitudes
    table = estimate_other_filters(table)

    # keep track of center and radius
    table.meta["center"] = center_skycoord
    table.meta["radius"] = radius
    table.meta['epoch'] = gaia_epoch

    # return the table
    return table


def plot_gaia(
    table,
    filter="G_gaia",
    faintest_magnitude_to_show=20,
    faintest_magnitude_to_label=16,
    size_of_zero_magnitude=100,
    unit=u.arcmin,
    show_spectra=False,
    **kwargs
):
    """
    Plot a finder chart using results from `get_gaia_data`.

    Use the table of positions and photometry returned by
    the `get_gaia_data` function to plot a finder chart
    with symbol sizes representing the brightness of the
    stars in a particular filter.

    Parameters
    ----------
    filter : str
        The filter to use for setting the size of the points.
        Options are "G_gaia", "RP_gaia", "BP_gaia", "g_sloan",
        "r_sloan", "i_sloan", "V_johnsoncousins", "R_johnsoncousins",
        "I_johnsoncousins". Default is "G_gaia".
    faintest_magnitude_to_show : float
        What's the faintest star to show? Default is 20.
    faintest_magnitude_to_label : float
        What's the faintest magnitude to which we should
        add a numerical label? Default is 16.
    size_of_zero_magnitude : float
        What should the size of a zeroth magnitude star be?
        Default is 100.
    unit : Unit
        What unit should be used for labels? Default is u.arcmin.
    **kwargs : dict 
        Additional keywords will be passed to `plt.scatter`.
    """

    # extract the center and size of the field
    center = table.meta["center"]
    radius = table.meta["radius"]

    # find offsets relative to the center
    dra = ((table["ra"] - center.ra) * np.cos(table["dec"])).to(unit)
    ddec = (table["dec"] - center.dec).to(unit)

    # set the sizes of the points
    mag = table[f"{filter}_mag"].to_value("mag")
    size_normalization = size_of_zero_magnitude / faintest_magnitude_to_show**2
    marker_size = (
        np.maximum(faintest_magnitude_to_show - mag, 0) ** 2 * size_normalization
    )

    # handle astropy units better
    with quantity_support():

        # plot the stars
        plt.scatter(dra, ddec, s=marker_size, **kwargs)
        plt.xlabel(
            rf"$\Delta$(Right Ascension) [{unit}] relative to {center.ra.to_string(u.hour, format='latex', precision=2)}"
        )
        plt.ylabel(
            rf"$\Delta$(Declination) [{unit}] relative to {center.dec.to_string(u.deg, format='latex', precision=2)}"
        )

        # add labels
        filter_label = filter.split("_")[0]
        to_label = np.nonzero(mag < faintest_magnitude_to_label)[0]
        for i in to_label:
            plt.text(
                dra[i],
                ddec[i],
                f"  {filter_label}={mag[i]:.2f}",
                ha="left",
                va="center",
                fontsize=5,
            )

        # add a grid
        plt.grid(color="gray", alpha=0.2)

        # plot a circle for the edge of the field
        circle = plt.Circle(
            [0, 0], radius, fill=False, color="gray", linewidth=2, alpha=0.2
        )
        plt.gca().add_patch(circle)
    
        # set the axis limits
        plt.xlim(radius, -radius)
        plt.ylim(-radius, radius)
        plt.axis("scaled")

        
def plot_star_spectra(center, 
                 starmap, 
                 ra_scale_fac=1, 
                 dec_scale_fac=1, 
                 ra_offset=-0.5, 
                 dec_offset=-0.5,
                 colormap="viridis"):
    
    """
    Plots the Gaia DR3 spectra of stars 
    within the given star map.
    
    Parameters
    ----------
    center : SkyCoord
        An astropy SkyCoord object indicating the
        right ascension and declination of the center.
    starmap : QTable
        An astropy table containing the results,
        with columns for different coordinates
        or filters, and rows for different stars.
    ra_scale_fac : float
        Scaling factor of overplotted spectra in the
        x-direction
    dec_scale_fac : float
        Scaling factor of overplotted spectra in the
        y-direction
    ra_offset : float
        Offset of spectra from their corresponding
        star in the x-direction
    dec_offset : float
        Offset of spectra from their corresponding
        star in the y-direction
    """

    # Loads spectroscopic data
    calibrated_spectra, sampling = calibrate(list(starmap['source_id']))
    
    # Unpacks starmap center ra and dec
    ra_center = center.ra
    dec_center = center.dec
    
    # Sets up color mapping for spectra
    cmap = get_cmap(colormap)
    norm = Normalize(vmin=0, vmax=1)

    # Iterates over every spectra found in Gaia
    for star_id, flux in zip(calibrated_spectra['source_id'], calibrated_spectra['flux']):

        # Corresponding (ra, dec) for every spectra
        idx = list(starmap['source_id']).index(star_id)
        ra = starmap['ra'][idx]
        dec = starmap['dec'][idx]

        # Calculates the offset used for each spectra
        star_ra = ((ra - ra_center).to(u.arcminute) * np.cos(dec)).value
        star_dec = (dec - dec_center).to(u.arcminute).value

        # Normalizes spectra data to range from [0, 1]
        normalized_sampling = (sampling - min(sampling)) / (max(sampling) - min(sampling))
        normalized_flux = (flux - min(flux)) / (max(flux) - min(flux))

        # Determines the color of the spectra
        c = normalized_sampling[np.argmax(normalized_flux)]
        color = cmap(norm(c))

        # Adds offsets and applies scales to spectra
        new_sampling = ra_scale_fac * ((sampling - min(sampling)) / (max(sampling) - min(sampling)) + ra_offset) + star_ra
        new_flux = dec_scale_fac * ((flux - min(flux)) / (max(flux) - min(flux)) + dec_offset) + star_dec

        # Plots spectra at the corresponding star
        plt.plot(new_sampling, new_flux, color=color)
        
    
def download_valid_data(save_address=os.getcwd(), 
                        target_name=None, 
                        ra=None, 
                        dec=None, 
                        missions=("TESS")):
    
    """
    Attempts to download all data products 
    for a given target star to a user-give 
    directory. Data products can be searched
    for using the star's name or its celestial 
    coordinates.
    
    Parameters
    -----------
    save_address : string
        The address at which downloaded 
        lightcurve data will be stored. Defaults 
        to present working directory
    target_name : string
        Name of a star
    ra : float
        Right ascension of a star
    dec : float
        Declination of a star
    missions : tuple
        List of which missions to download data
        for. Valid options are "Kepler," "K2," 
        or "TESS."
    
    Returns
    -------
    valid_data_tables : list
    
    """

    # Sets download location to folder in given directory
    save_address = save_address + "\lightcurve_data"
    
    # Throws an error if no target star info is given
    if all(i is None for i in [target_name, ra, dec]):
        raise ValueError("Must provide target name or position")
    
    # Handles downloads for a given star name
    elif type(target_name)==str:
        results = lk.search_targetpixelfile(target_name, mission=missions)
        
    # Handles downloads for a given star coordinate
    else:
        position = SkyCoord(ra=ra, dec=dec)
        results = lk.search_targetpixelfile(position, mission=missions)
        
    # Creates a list to hold data products
    valid_data_tables = [0 for i in range(len(results))]
    
    # Downloads every data product
    idx=0
    while idx < len(results):
        
        successful_download=False
        while successful_download==False:
            
            # Attempts to download a data product from the results table
            try:
                data_product = results[idx].download(download_dir=save_address)
                valid_data_tables[idx] = data_product
                successful_download=True
                idx += 1
                
            # Handles an error where data has been downloaded incorrectly
            except:
                print("WARNING: Failed to download data product...")
                print("This is typically caused by files that have only partially downloaded\n")
                response = input("Attempt to delete existing data and redownload (y/n): ")
                response = (response.upper() == "Y")
                
                # Deletes existing data
                if response:

                    idx=0
                    ids = np.unique(results.target_name)
                    
                    # Iterates over all data directories
                    for subdir in os.listdir(rf"{save_address}\mastDownload"):
                        
                        data_dir = rf"{save_address}\mastDownload\{subdir}"
                        
                        # Deletes each target star file in these directories
                        for file in os.listdir(data_dir):
                            for subfile in os.listdir(rf"{data_dir}\{file}"):
                                for target in ids:
                                    if target in file:
                                        os.remove(rf"{data_dir}\{file}\{subfile}")
                                        os.rmdir(rf"{data_dir}\{file}")

    # Drops leftover zeros in the data products list
    valid_data_tables = list(filter(lambda num: num != 0, valid_data_tables))
    
    return valid_data_tables

def generate_valid_lightcurves(starmap, 
                               save_address=os.getcwd(), 
                               missions=("TESS")):
    
    """
    Attempts to generate lightcurves for
    each star in a given starfield.
    
    Parameters
    -----------
    starmap : QTable
        An astropy table containing the results,
        with columns for different coordinates
        or filters, and rows for different stars.
    save_address : string
        The address at which downloaded 
        lightcurve data will be stored. Defaults 
        to present working directory
    missions : tuple
        List of which missions to download data
        for. Valid options are "Kepler," "K2," 
        or "TESS."
    """
    
    # Iterates over every star in the starmap
    for i in range(len(starmap)):
        
        # Extracts star location
        ra = starmap['ra'][i]
        dec = starmap['dec'][i]
        
        print(f"Looking for data products at ({ra}, {dec})...")
        
        # Downloads all valid data products for the given star
        valid_data = download_valid_data(save_address, ra=ra, dec=dec, missions=missions)
        
        print(f"Found {len(valid_data)} valid data products...\n")
        
        # Produces lightcurves for each data product
        for data_product in valid_data:
            data_product.to_lightcurve(method='pld').remove_outliers().flatten().scatter()
            plt.show()