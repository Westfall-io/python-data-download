from astropy import units as u
from astropy.coordinates import EarthLocation

def make_targets(config, debug):
    if debug:
        print('-'*40)
        print('Making targets.')
    # Make Earth Locations
    els = []
    a_f = []
    for k,v in enumerate(config.scenario.targets.lats):
        lat = v
        lon = config.scenario.targets.lons[k]
        els.append(
            EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=0*u.m)
        )
        a_f.append(open("access_output/access_{}_{}_{}.csv".format(k,int(lat*100),int(lon*100)),'w'))
        a_f[-1].write('sat,access_start,access_end')

    return els, a_f
