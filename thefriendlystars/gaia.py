from .imports import *
from astropy.visualization import quantity_support
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.table import QTable
from gaiaxpy import calibrate
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from tqdm import tqdm
import lightkurve as lk
import pandas as pd
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

##############################################################
######################     NEW CODE     ######################
##############################################################

def search_result_to_df(search_result):
    
    """
    Converts lightkurve SearchResult object
    into a pandas dataframe
    
    Parameters
    ----------
    search_result : SearchResult
        A lightkurve SearchResult containing 
        information about data products for
        a given star
        
    Returns
    -------
    df : pd.DataFrame
        Pandas dataframe containing the same
        information as the input search result
    """

    # Initializes a dataframe to store search result data
    df = pd.DataFrame(columns=['mission', 'year', 'author', 'exptime', 'target_name', 'distance'])

    # Iterates over every row in the SearchResult
    for idx, row in enumerate(str(search_result).split('\n')[5:]):

        # Splits a row using spaces as a delimeter
        row = row.split()

        # Converts the split up data into the proper format
        mission = row[1] + " " + row[2] + " " + row[3]
        year = int(row[4])
        author = row[5]
        exptime = int(row[6])
        target_name = int(row[7])
        distance = float(row[8])

        # Adds row data to the dataframe
        df.loc[idx] = [mission, year, author, exptime, target_name, distance]
        
    return df

def delete_conflicting_data(save_address, search_results):
    
    """
    Deletes all local user data for a
    given star.
    
    Parameters
    ----------
    save_address : string
        Local save directory for lightkurve data
    search_results : SearchResult
        Lightkurve search result table with
        data product information
    """
    
    ids = np.unique(search_results.target_name)

    # Iterates over all data directories
    for subdir in os.listdir(rf"{save_address}\mastDownload"):

        data_dir = rf"{save_address}\mastDownload\{subdir}"

        # Deletes each target star file in these directories
        for file in os.listdir(data_dir):
            
            print(file)
            
            for subfile in os.listdir(rf"{data_dir}\{file}"):
                for target in ids:
                    if target in file:
                        os.remove(rf"{data_dir}\{file}\{subfile}")
                        os.rmdir(rf"{data_dir}\{file}")

def download_valid_data(save_address=os.getcwd(), 
                        target_name=None, 
                        ra=None, 
                        dec=None, 
                        missions=("TESS"),
                        exptime_to_download="longest"):
    
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
        search_results = lk.search_lightcurve(target_name, mission=missions)
        
    # Handles downloads for a given star coordinate
    else:
        position = SkyCoord(ra=ra, dec=dec)
        search_results = lk.search_lightcurve(position, mission=missions)
        
    if len(search_results) == 0:
        return
       
    # Extracts the exposures of each data product
    results_df = search_result_to_df(search_results)
    exp_list = list(results_df['exptime'])
    
    # Finds the index of a valid data product
    exp_dict = {"shortest":min(exp_list), "longest":max(exp_list)}
    exp_use = exp_dict[exptime_to_download]
    idx = exp_list.index(exp_use)
        
    # Runs until download works
    successful_download=False
    while successful_download==False:

        # Attempts to download a data product from the results table
        try:
            data_product = search_results[idx].download(download_dir=save_address)
            successful_download=True

        # Handles an error where data has been downloaded incorrectly
        except:
            print("WARNING: Failed to download data product...")
            print("This is typically caused by files that have only partially downloaded\n")
            do_delete = input("Attempt to delete existing data and redownload (y/n): ")
            do_delete = (do_delete.upper() == "Y")

            if do_delete:
                delete_conflicting_data(save_address, search_results)
    
    return data_product

def rescale(vals, old_min, old_max, new_min=0, new_max=1):

    """
    Rescales an array of values from
    one scale to another.
    
    Parameters
    ----------
    vals : ndarray
        Array of values to be rescaled
    old_min : float
        Original minimum of the passed
        vals array. Does not need to be
        the minimum array value.
    old_max : float
        Original maximum of the passed
        vals array. Does not need to be
        the maximum array value.
    new_min : float
        New minimum value of the rescaled
        vals array.
    new_max : float
        New maximum value of the rescaled
        vals array.
        
    Returns
    -------
    rescaled_vals : ndarray
        Rescaled array of values
    """
    
    return new_min + (new_max - new_min) * (vals - old_min)/(old_max - old_min)

def resolve_box_overlap(xs_orig, 
                        ys_orig, 
                        width, 
                        height, 
                        x_bound,
                        y_bound,
                        x_offset_fac=0.2, 
                        y_offset_fac=1.1):
    
    """
    Adjusts box locations so that they 
    no longer overlap.
    
    Parameters
    ----------
    xs_orig : ndarray
        X-coordinate positions of each
        box.
    ys_orig : ndarray
        Y-coordinate positions of each
        box.
    width : float
        Width of every box (must be
        uniform for every box)
    height : float
        Height of every box (must be
        uniform for every box)
    x_offset_fac : float
        Scaling factor used for displacing
        boxes in the x-direction.
    y_offset_fac : float
        Scaling factor used for displacing
        boxes in the y-direction.
        
    Returns
    -------
    xs : ndarray
        Modified x-coordinate positions
        of each box.
    ys : ndarray
        Modified y-coordinate positions
        of each box.
    """
    
    resolved = False     # Tracks whether or not any boxes are overlapping
    iterations = 0       # Number of iterations of attempted overlap resolution
    max_iterations = 25  # Add a maximum iteration limit to prevent infinite loops
    
    xs = xs_orig.copy()  # New x-array to prevent overwriting original x-data
    ys = ys_orig.copy()  # New y-array to prevent overwriting original y-data

    # Iterates until max iteration limit is hit
    for iteration in range(max_iterations):
        
        resolved = True

        # Iterates over every box
        for i in range(len(xs)):
            for j in range(i + 1, len(xs)):
                
                # Calculates x- and y-separations
                x_sep = xs[i] - xs[j]
                y_sep = ys[i] - ys[j]

                # Checks if two boxes are overlapping
                if abs(x_sep) < width and abs(y_sep) < height:
                    
                    # Calculates how much to offset these boxes
                    x_offset = x_offset_fac * (0.5 * (width - x_sep))
                    y_offset = y_offset_fac * (0.5 * (height - y_sep))

                    # Resolves y-direction overlaps
                    if y_sep > 0:
                        new_yi = min(ys[i] + y_offset, y_bound[1] - height)
                        new_yj = max(ys[j] - y_offset, y_bound[0])
                    else:
                        new_yi = max(ys[i] - y_offset, y_bound[0])
                        new_yj = min(ys[j] + y_offset, y_bound[1] - height)
                    ys[i], ys[j] = new_yi, new_yj

                    # Resolves x-direction overlaps
                    if x_sep > 0:
                        new_xi = min(xs[i] + x_offset, x_bound[1] - width)
                        new_xj = max(xs[j] - x_offset, x_bound[0])
                    else:
                        new_xi = max(xs[i] - x_offset, x_bound[0])
                        new_xj = min(xs[j] + x_offset, x_bound[1] - width)
                    xs[i], xs[j] = new_xi, new_xj

                    resolved = False
        
        # Ends the loop if no boxes overlap
        if resolved:
            break

    return xs, ys

def generate_valid_lightcurves(center,
                               starmap_radius,
                               starmap, 
                               save_address=os.getcwd(), 
                               missions=("TESS"),
                               x_factor=5,
                               y_factor=2,
			       exp_time_to_download="longest"):
    
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
    
    # Controls the size of the lightcurves
    x_factor = 5
    y_factor = 2
    
    # Unpacks starmap center coordinates
    ra_center = center.ra
    dec_center = center.dec
    
    # Generates alpha values to use for lightcurves based on magnitude
    mags = np.array(starmap['G_gaia_mag'].value)
    cmap = get_cmap('Greys')
    c = cmap(rescale(mags, old_min=max(mags), old_max=min(mags), new_min=0.3, new_max=1))
    
    # Initializes arrays for star locations
    star_ra = np.array([0. for i in range(len(starmap))])
    star_dec = np.array([0. for i in range(len(starmap))])
    
    # Calculates the offset of each lightcurve
    for idx, loc in enumerate(zip(starmap['ra'], starmap['dec'])):
        
        ra, dec = loc
        star_ra[idx] = ((ra - ra_center).to(u.arcminute) * np.cos(dec)).value
        star_dec[idx] = (dec - dec_center).to(u.arcminute).value
    
    # Resolves overlapping bounding boxes
    star_ra_updated, star_dec_updated = resolve_box_overlap(star_ra, star_dec, x_factor/2, y_factor/2,
                                                           x_bound=[-starmap_radius.value, starmap_radius.value],
                                                           y_bound=[-starmap_radius.value, starmap_radius.value])
    print("    Finished resolving overlapping spectra...")
    print("    Formatting lightcurve data...")
    
    # Iterates over every star in the starmap
    for i in tqdm(range(len(starmap))):
        
        # Downloads lightcurve for a given star
        valid_data = download_valid_data(save_address, 
                                         ra=starmap['ra'][i], 
                                         dec=starmap['dec'][i], 
                                         missions=missions,
                                         exptime_to_download='longest')
        
        # Runs when any data product is found
        if valid_data is not None:

            # Converts data products into lightcurve data arrays
            lk_data = valid_data.remove_outliers().normalize()
            times = lk_data['time'].value
            fluxes = lk_data['flux'].value
            
            # Rescales lightcurves
            norm_time = rescale(times, old_max=max(times), old_min=min(times))
            norm_flux = rescale(fluxes, old_max=1.01, old_min=0.99)
            
            # Moves and stretches lightcurves
            adjust_time = np.flip(x_factor*(norm_time-0.5) + star_ra_updated[i])
            adjust_flux = np.flip(y_factor*(norm_flux-0.5) + star_dec_updated[i])

            # Plots lightcurve
            plt.scatter(adjust_time, adjust_flux, s=0.01, color=c[i])
            
            # Plots a line connecting moved lightcurves to their host star
            #if np.abs(dy) != 0:
            plt.plot([star_ra[i], star_ra_updated[i]], [star_dec[i], star_dec_updated[i]], color='blue', ls='--')

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