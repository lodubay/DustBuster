# just copy-pasting coordinate transform functions developed in the dust_coordinates.ipynb notebook to access in other notebooks

import numpy as np
from astropy.coordinates import SkyCoord

def galactic_coords(x, y, z, sunpos=(300.5, 300.5, 40.5), gridstep=10):
    """
    Convert X, Y, Z position in the data array into Galactic coordinates.
    
    Parameters
    ----------
    x : float or array of floats
        Data x-coordinate
    y : float or array of floats
        Data y-coordinate
    z : float or array of floats
        Data z-coordinate
    sunpos : 3-tuple, optional
        Position of the sun in data coordinates. The default is (300.5, 300.5, 40.5).
    gridstep : float, optional
        Data grid step size in pc. The default is 10.
        
    Returns
    -------
    l, b, distance : float or array of floats
        Galactic (longitude, latitude, distance) in (deg, deg, pc)
    """
    datapoint = np.array((x, y, z))
    # convert data coordinates to physical cartesian coordinates in pc
    cartesian = (datapoint - np.array(sunpos)) * gridstep
    sc = SkyCoord(w=cartesian[2], u=cartesian[0], v=cartesian[1],
                  unit='pc', frame='galactic', representation_type='cartesian')
    sc.representation_type = 'spherical'
    return sc

def radec_coords(x, y, z, sunpos=(300.5, 300.5, 40.5), gridstep=10):
    """
    Convert X, Y, Z position in the data array into equatorial coordinates (ICRS).
    
    Parameters
    ----------
    x : float or array of floats
        Data x-coordinate
    y : float or array of floats
        Data y-coordinate
    z : float or array of floats
        Data z-coordinate
    sunpos : 3-tuple, optional
        Position of the sun in data coordinates. The default is (300.5, 300.5, 40.5).
    gridstep : float, optional
        Data grid step size in pc. The default is 10.
        
    Returns
    -------
    l, b, distance : float or array of floats
        ICRS (right ascension, declination, distance) in (deg, deg, pc)
    """
    datapoint = np.array((x, y, z))
    # convert data coordinates to physical cartesian coordinates in pc
    cartesian = (datapoint - np.array(sunpos)) * gridstep
    sc = SkyCoord(w=cartesian[2], u=cartesian[0], v=cartesian[1],
                  unit='pc', frame='galactic', representation_type='cartesian')
    return sc.icrs