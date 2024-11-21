import math

from astropy import units as u
from astropy.constants import k_B, c
from hapsira.bodies import Earth

def calc_graze_snr_floor(config, debug):
    if debug:
        print('-'*40)
        print('Building comms constraints.')

    # Transmit power
    transmit_pwr = 10*math.log10(config.payload.comms.power)
    gain_pwr = config.payload.comms.gain
    loss = config.payload.comms.cable_loss
    eirp = transmit_pwr + gain_pwr - loss #dB

    # Receive power
    recv_pwr = config.scenario.targets.recvr_gain #dB

    # Noise
    B = config.payload.comms.bandwidth * 1000 * u.Hz # Hz
    T = 290 * u.K # Room temp
    noise_pwr = 10*math.log10((k_B * B * T).to(u.W).value)

    # Max range
    snr = transmit_pwr + recv_pwr - noise_pwr
    max_path_loss = 10**((snr-config.payload.comms.min_snr)/10)
    p1 = 4 * math.pi * (config.payload.comms.frequency*1e9)
    p2 = math.sqrt(max_path_loss) * c.to(u.m/u.s).value
    max_distance = (p2/p1)*u.m
    print('Maximum Sensor Distance: {}'.format(max_distance.to(u.km)))

    ER = Earth.R.to(u.km)
    o_a = config.constellation.orbit.h * u.km + ER

    # Law of cosines
    a_side = max_distance.to(u.km).value
    b_side = ER.value
    c_side = o_a.to(u.km).value

    p1 = a_side**2+b_side**2-c_side**2
    p2 = 2*a_side*b_side

    graze = math.acos(p1/p2)-math.pi/2
    if graze < 0:
        # We can see further than the horizon allows
        graze = 0

    graze = graze * u.rad

    print('Graze Limit: {} deg.'.format(graze.to(u.deg).value))
    return graze

def calc_look_angle(config, debug):
    if debug:
        print('-'*40)
        print('Building sensor constraints.')


    # Calculate look angle
    ppm = config.payload.sensor.pixel_pitch_micron * u.micron
    fl = config.payload.sensor.focal_length * u.cm
    mr = config.payload.sensor.max_resolution * u.m
    ER = Earth.R.to(u.km)
    o_a = config.constellation.orbit.h * u.km + ER

    # Calculate instantaneous field of view
    ifov = ppm.to(u.m)/fl.to(u.m) #rad

    # Calculate maximum distance based on resolution
    # Projecting one pixel at the max resolution, how far away
    # does the angle that represent one pixel reach the desired size
    #max_distance = mr/math.tan(ifov.value)
    # Law of sines
    s1 = mr/math.sin(ifov.value)
    other_angle = (math.pi-ifov.value)/2 * u.rad
    s2 = 1/math.sin(other_angle.to(u.rad).value)
    max_distance = s1/s2
    print('Maximum Sensor Distance: {}'.format(max_distance.to(u.km)))

    # Law of cosines
    a = max_distance.to(u.km).value
    b = ER.value
    c = o_a.to(u.km).value

    p1 = a**2+b**2-c**2
    p2 = 2*a*b

    graze = math.acos(p1/p2)-math.pi/2
    if graze < 0:
        # We can see further than the horizon allows
        graze = 0

    graze = graze * u.rad

    print('Graze Limit: {} deg.'.format(graze.to(u.deg).value))
    return graze

def make_orbits(config, debug):
    import copy

    from hapsira.twobody import Orbit
    from astropy.time import Time
    epoch = Time(config.constellation.orbit.epoch)

    # Create first satellite
    raan = 0 * u.deg
    nu = 0 * u.deg
    ER = Earth.R.to(u.km)
    o_a = config.constellation.orbit.h * u.km + ER

    if debug:
        print('-'*40)
        print('Building satellites in constellation.')
        print("Sat 1: 0,0")

    orb = [copy.deepcopy(
        Orbit.from_classical(
            Earth,
            o_a,
            config.constellation.orbit.ecc * u.one,
            config.constellation.orbit.inc * u.deg,
            raan,
            config.constellation.orbit.argp * u.deg,
            nu,
            epoch=epoch)
        )
    ]

    if config.constellation.sats == 1:
        return orb

    plane = 1
    stp = 1
    plane_draan = 360/config.constellation.planes #36 deg/plane in raan
    sats_per_plane = config.constellation.sats/config.constellation.planes
    plane_dnu = 360/sats_per_plane #360 deg/sat in nu (only 1 sat)
    for sat in range(2,config.constellation.sats+1):
        if stp + 1 > sats_per_plane:
            plane = plane + 1
            raan += plane_draan*u.deg
            nu = 0 * u.deg
        else:
            nu += plane_dnu * u.deg
        if debug:
            print("Sat {}: {},{}".format(sat,raan.value,nu.value))
        orb.append(
            copy.deepcopy(
                Orbit.from_classical(
                    Earth,
                    o_a,
                    config.constellation.orbit.ecc * u.one,
                    config.constellation.orbit.inc * u.deg,
                    raan,
                    config.constellation.orbit.argp * u.deg,
                    nu,
                    epoch=epoch
                )
            )
        )

    return orb

def propagate_orbit(sat_now, ss):
    #import math
    from astropy.coordinates import GCRS

    sat_now = sat_now.propagate(ss)
    #t_now = Time((time.to_datetime()+datetime.timedelta(sim_step.to(u.s).value)).isoformat())
    #print(sat_now.r, math.sqrt(sum([x.to(u.km).value**2 for x in sat_now.r])))
    sat_eci = GCRS(
        sat_now.r[0].to(u.m),
        sat_now.r[1].to(u.m),
        sat_now.r[2].to(u.m),
        obstime=sat_now.epoch,
        representation_type='cartesian'
    )
    return sat_now, sat_eci

def check_targets(cnt_accesses, sat_eci, epoch, sat_num, ss, step, els, graze,
                  accesses, a_f, a_tf, debug):
    from astropy.coordinates import AltAz
    for k, el in enumerate(els):
        # Find the satellite in the sky above this target
        sat_altaz = sat_eci.transform_to(
            AltAz(obstime=epoch, location=el)
        )

        if sat_altaz.alt.to(u.deg) >= graze[k].to(u.deg) and accesses[k]==0:
            if debug:
                print('Access start at {} on target {}.'.format(epoch.to_datetime(),k))
            accesses[k]=1
            a_tf[k]=step
            #print(sat_altaz.distance.to(u.km))
        elif sat_altaz.alt.to(u.deg) < graze[k].to(u.deg) and accesses[k]==1:
            if debug:
                print('Access end at {} on target {}.'.format(epoch.to_datetime(),k))
            cnt_accesses[k] += 1
            accesses[k]=0
            a_f[k].write('\n{},{},{}'.format(sat_num,(a_tf[k]*ss).to(u.min).value,(step*ss).to(u.min).value))
            a_tf[k]=-1
        else:
            continue

    # Return a binary for each target whether it is in access or not [accesses]
    # Return the time at which this accesses started
    return accesses, a_tf

def make_accesses(config, orb, els, a_f, graze, debug):
    # import datetime
    if debug:
        print('-'*40)
        print('Creating accesses.')

    if len(graze) == 1:
        # Make the graze constraint the same for all targets
        graze = graze*len(els)

    d = config.scenario.duration * u.day
    ss = config.scenario.sim_step * u.min
    steps = (d.to(u.min)/ss.to(u.min)).to(u.one).value

    cnt_accesses = [0]*len(els)
    for sat_num in range(len(orb)):
        sat_now = orb[sat_num]
        if debug:
            print(
                'Checking this satellite: {} -- {},{}'.format(
                    sat_num, sat_now.raan, sat_now.nu
                )
            )

        sat_now = orb[sat_num]

        step = 0
        accesses = [0]*len(els)
        a_tf = [-1]*len(els)
        while step < steps:
            #print('Sim Step: {}'.format(i))
            sat_now, sat_eci = propagate_orbit(sat_now, ss)
            accesses, a_tf = check_targets(
                cnt_accesses, sat_eci, sat_now.epoch, sat_num, ss, step, els,
                graze, accesses, a_f, a_tf, debug
            )
            step += 1

        # Check if any had open accesses at sim end
        for k,v in enumerate(a_tf):
            if v != -1:
                a_f[k].write(
                    '\n{},{},{}'.format(
                        sat_num,
                        (a_tf[k]*ss).to(u.min).value,
                        (step*ss).to(u.min).value
                    )
                )

    # Close all files
    for k,v in enumerate(a_f):
        v.close()
