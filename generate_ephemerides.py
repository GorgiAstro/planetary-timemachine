import orekit
orekit.initVM()

# Modified from https://gitlab.orekit.org/orekit-labs/python-wrapper/blob/master/python_files/pyhelpers.py
from java.io import File
from org.orekit.data import DataProvidersManager, DirectoryCrawler
from orekit import JArray

orekit_data_dir = 'orekit-data'
DM = DataProvidersManager.getInstance()
datafile = File(orekit_data_dir)
if not datafile.exists():
    print('Directory :', datafile.absolutePath, ' not found')

crawler = DirectoryCrawler(datafile)
DM.clearProviders()
DM.addProvider(crawler)

from org.orekit.frames import FramesFactory
icrf = FramesFactory.getICRF()
eme2000 = FramesFactory.getEME2000()

from org.orekit.time import TimeScalesFactory
tai = TimeScalesFactory.getTAI()
utc = TimeScalesFactory.getUTC()

from org.orekit.time import AbsoluteDate
mjd_epoch_tai = AbsoluteDate(1858, 11, 17, 0, 0, 0.0, tai)

planets = [
    'Mercury',
    'Venus',
    'Earth',
    'Moon',
    'Mars',
    'Jupiter',
    'Saturn',
    'Uranus',
    'Neptune'
]

originator = 'GorgiAstro'
from odmadmpy.core import Oem, Aem
oem = Oem(originator, standard='CIC')
aem = Aem(originator, standard='CIC')
comment = 'Generated from Orekit using DE405 ephemerides'

cic_output_folder = 'CIC-data'

meta_mandat_oem_cic_sample = {
    'OBJECT_NAME': 'MARS',
    'OBJECT_ID': 'MARS',
    'CENTER_NAME': 'SUN',
    'REF_FRAME': 'ICRF',
    'TIME_SYSTEM': 'TAI'
}

meta_mandat_oem_cic_moon = {
    'OBJECT_NAME': 'MOON',
    'OBJECT_ID': 'MOON',
    'CENTER_NAME': 'EARTH',
    'REF_FRAME': 'EME2000',
    'TIME_SYSTEM': 'TAI'
}

meta_mandat_aem_cic_sample = {
    'OBJECT_NAME': 'MARS',
    'OBJECT_ID': 'MARS',
    'REF_FRAME_A': 'ICRF',
    'REF_FRAME_B': 'BODY',
    'ATTITUDE_DIR': 'A2B',
    'TIME_SYSTEM': 'TAI',
    'ATTITUDE_TYPE': 'QUATERNION'
}
meta_opt_aem_cic = {
    'QUATERNION_TYPE': 'FIRST'
}

meta_mandat_aem_cic_moon = {
    'OBJECT_NAME': 'MOON',
    'OBJECT_ID': 'MOON',
    'REF_FRAME_A': 'EME2000',
    'REF_FRAME_B': 'BODY',
    'ATTITUDE_DIR': 'A2B',
    'TIME_SYSTEM': 'TAI',
    'ATTITUDE_TYPE': 'QUATERNION'
}

from org.orekit.bodies import CelestialBodyFactory
from org.orekit.utils import PVCoordinatesProvider
from orekit.pyhelpers import datetime_to_absolutedate, absolutedate_to_datetime
import pandas as pd
import numpy as np
import os
from tqdm import tqdm

date_init = AbsoluteDate(1957, 1, 1, tai)
date_end = AbsoluteDate(2057, 1, 1, tai)
dt = 86400.0 / 5

for planet_name in planets:
    print(f'Started generating ephemerides for {planet_name}')

    if planet_name == 'Moon':
        meta_mandat_oem_cic = meta_mandat_oem_cic_moon
        meta_mandat_aem_cic = meta_mandat_aem_cic_moon
        ref_frame = eme2000
    else:
        meta_mandat_oem_cic = meta_mandat_oem_cic_sample
        meta_mandat_aem_cic = meta_mandat_aem_cic_sample
        meta_mandat_oem_cic['OBJECT_NAME'] = planet_name.upper()
        meta_mandat_oem_cic['OBJECT_ID'] = planet_name.upper()
        meta_mandat_aem_cic['OBJECT_NAME'] = planet_name.upper()
        meta_mandat_aem_cic['OBJECT_ID'] = planet_name.upper()
        ref_frame = icrf

    df = pd.DataFrame(columns=['MJD', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'qs', 'qx', 'qy', 'qz'])
    cbody = CelestialBodyFactory.getBody(planet_name.upper())
    cbody_pv_provider = PVCoordinatesProvider.cast_(cbody)
    cbody_frame = cbody.getBodyOrientedFrame()

    date_current = date_init
    years_from_start = 0.0
    years_from_start_previous = 0.0
    with tqdm(total=100) as pbar:
        while date_current.compareTo(date_end) < 0:
            pv = cbody_pv_provider.getPVCoordinates(date_current, ref_frame)

            ref_frame_to_body = ref_frame.getTransformTo(cbody_frame, date_current)
            rot = ref_frame_to_body.getRotation()
            quat = [rot.getQ0(), rot.getQ1(), rot.getQ2(), rot.getQ3()]

            mjd = date_current.offsetFrom(mjd_epoch_tai, tai) / 86400

            df.loc[absolutedate_to_datetime(date_current)] = np.concatenate((
                [mjd], 
                np.array(pv.getPosition().toArray()),
                np.array(pv.getVelocity().toArray()),
                quat))
            date_current = date_current.shiftedBy(dt)

            date_comp = date_current.getComponents(tai)
            is_new_year = (date_comp.getDate().getDayOfYear() == 1) and (date_comp.getTime().getSecondsInLocalDay() < dt / 2)
            if is_new_year:
                pbar.update(1)

    df['datetime'] = df.index

    oem_segments = []
    oem_segments += oem.format_segment(df, meta_mandat_oem_cic)
    oem.write_file(oem_segments, os.path.join(cic_output_folder, f'{planet_name}_OEM_EME2000.txt'), comments=[comment])

    aem_segments = []
    aem_segments += aem.format_segment(df, meta_mandat_aem_cic, meta_opt_aem_cic)
    aem.write_file(aem_segments, os.path.join(cic_output_folder, f'{planet_name}_AEM_EME2000.txt'), comments=[comment])

    print(f'Finished generating ephemerides for {planet_name}')
