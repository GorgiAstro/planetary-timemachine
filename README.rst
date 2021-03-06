planetary-timemachine
=====================

Explore the solar system planets (and the Earth's Moon) from 1957 to 2057. This can be useful to visualize the lighting conditions at the time a picture was taken by an interplanetary probe.

The VTS Timeloop software is used for visualization and runs on Windows and Linux. It can be downloaded for free at https://timeloop.fr/vts/.

It is different than NASA's eyes on the solar system: it does not show historical interplanetary missions, but it offers more 2D and 3D visualization options.

VTS installation
################

Just unpack the archive. Some extra packages might be needed in Linux, as described in https://timeloop.fr/vts/download/ :

* libpng12: https://packages.ubuntu.com/xenial/libpng12-0
* libjpeg62: https://packages.ubuntu.com/bionic/libjpeg62
* libssl1.0: https://packages.ubuntu.com/bionic/libssl1.0.0

WMS Layers
##########

Web Map Service (WMS) layers are configured in the VTS project for Mercury, Venus, the Earth, and Mars. The WMS URLs were found here: https://astrowebmaps.wr.usgs.gov/webmapatlas/Layers/maps.html

In addition, the Sentinel cloudless WMS from EOX is used from https://s2maps.eu/.

Ephemerides
###########

By default, VTS has ephemerides from 1999 to 2101 for the planets. Celestia has its own ephemerides, in this Celestia is configured to use the VSOP87 propagation model.

In the VTS project, the default ephemerides were overriden by the ones generated in this repository, which are valid from 1957 to 2057.

Running the VTS visualization
#############################

In VTS, open the project ``visualization.vts``. When running the project, three windows should appear:

* Broker: controls the time and visualization parameters
* SurfaceView: 2D view of the planet
* Celestia: 3D view

**********************************
Advanced: generate new ephemerides
**********************************

New ephemerides can be generated by using the Orekit/Python script ``generate_ephemerides.py``.

Two modified JPL BSP ephemerides are used in Orekit to retrieve the planets' position and attitude, located in ``orekit-data/DE-430-ephemerides/``. For dates prior to 1990, the DE405 model is used, otherwise DE430 is used.

The Orekit Python wrapper is only available via conda-forge. Create a new conda environment:

``conda env create -f environment.yml``

If the environment needs to be updated (this also forces the pip requirements to update):

``conda env update -f environment.yml``

Enter the conda environment:

``conda activate planetary``

Edit the script ``generate_ephemerides.py`` as you like, then run it. The script can take several hours to complete as Python loops are slow.

The CIC ephemerides files are generated using the `odmadmpy library <https://github.com/GorgiAstro/ccsds-cic-odm-adm-py>`_.
