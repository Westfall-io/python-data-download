import math
from os import listdir, path
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u

def find_csv_filenames( path_to_dir, suffix=".csv" ):
    filenames = listdir(path_to_dir)
    return [ filename for filename in filenames if filename.endswith( suffix ) ]

def collect_access_dataframes(config):
    o = []
    filenames = find_csv_filenames("access_output")
    for name in filenames:
        # For each target add a dataframe
         o.append(pd.read_csv(path.join('access_output',name)))
    return o

def generate_ground_schedule(config, o):
    s = 0
    ## TODO: More work to be done to schedule multiple sites and sats
    #sm = [0]*duration*1440

    v = pd.DataFrame()
    for k,df in enumerate(o):
        dfs = df[df['sat']==s]
        dfs['target']=k
        v = pd.concat([v, dfs])

    v = v.reset_index().sort_values('access_start')
    v['access_time'] = v['access_end']-v['access_start']
    v = v[v['access_time']>=config.constellation.satellite.comms.min_comms_access]

    schedule = pd.DataFrame()
    t = 0
    while t < len(v):
        # For each entry
        if len(schedule.index) == 0:
            schedule = v.iloc[[t]]
        else:
            if schedule.iloc[-1]['access_end'] < v.iloc[t]['access_start']:
                schedule = pd.concat([schedule, v.iloc[[t]]])
            else:
                pass
            #pass
        t += 1

    schedule = schedule.reset_index()

    return schedule

def calc_data_storage(config, debug):
    if debug:
        print('-'*40)
        print('Build data storage.')

    pixel_array = config.payload.sensor.pixels_x*config.payload.sensor.pixels_y
    bits_per_image = config.payload.sensor.bpp*pixel_array
    compressed_image_size = config.payload.sensor.comp_ratio*bits_per_image
    image_separation = config.payload.sensor.integration_time_s*u.s
    payload_data_rate = compressed_image_size / image_separation.to(u.min).value
    telemetry_data_rate = config.payload.telemetry.data_rate/(u.s)
    telemetry_data_rate = telemetry_data_rate.to(1/u.min).value*1024

    if config.scenario.comms.comms_while_collecting:
        payload_data_rate_in_comms = payload_data_rate
        payload_data_rate_comm_setup = payload_data_rate
    else:
        payload_data_rate_in_comms = 0
        payload_data_rate_comm_setup = 0

    # stateDiagram-v2
    #   [*] --> Collect
    #   Collect --> Comm_Attitude_Change
    #   Comm_Attitude_Change --> Comm
    #   Comm --> Collect

    duration = config.scenario.duration * u.day

    accesses = collect_access_dataframes(config)
    schedule = generate_ground_schedule(config, accesses)

    # Let's make a discrete event sim
    import queue
    images = queue.Queue()
    ground = queue.Queue()

    t = 0
    last_image = 0
    a_bool = False
    ssdr_list = []
    ssdr = 0
    while t <= duration.to(u.s).value:
        # Check if in access
        a = schedule[schedule['access_start'] <= t/60]
        a = a[a['access_end'] > t/60]
        if len(a)>1:
            # Too many accesses were returned
            raise NotImplementedError

        if len(a)==0:
            # Not in an access
            if a_bool:
                # Capture the ssdr size at the end of the comms window
                ssdr_list.append([t/60, ssdr/1024/1024/1024])
                a_bool = False

            pca = config.constellation.satellite.attitude.pre_comms_attitude
            b = schedule[schedule['access_start'] <= t/60 + pca]
            b = b[b['access_end'] > t/60 + pca]
            if len(b) > 1:
                # Too many accesses were returned
                raise NotImplementedError

            if len(b)==0:
                # Attitude maneuver not about to start
                images.put((t, compressed_image_size))
                t += image_separation.to(u.s).value
                ssdr += compressed_image_size
                #print('Adding image at t: {}, SSDR: {}'.format(t, ssdr/1024/1024/1024))
            else:
                # Check if image while comms
                if config.scenario.comms.comms_while_collecting:
                    # TODO
                    raise NotImplementedError
                else:
                    # Don't do anything until comms starts
                    t = schedule[schedule['access_start'] >= t/60].iloc[0]['access_start']*60
                    print('Waiting for attitude at t: {}'.format(t))
        else:
            # In an access -- downlink files
            if config.scenario.comms.comms_while_collecting:
                # TODO: Add this
                raise NotImplementedError

            if not a_bool:
                # Capture the current SSDR before we start downlinking
                ssdr_list.append([t/60, ssdr/1024/1024/1024])
                a_bool = True

            if not images.empty():
                i = images.get()
                ct = i[1]/(config.payload.comms.dl_rate*1024*1024) # b/(b/s) == s
                t += ct
                ssdr += -i[1]
                #print('Downloading image at t: {}, SSDR: {}'.format(t, ssdr/1024/1024/1024))
                ground.put((i[0], t))
            else:
                t = schedule[schedule['access_end'] >= t/60].iloc[0]['access_end']*60
                #print('No more images to downlink moving to t: {}'.format(t))

    ssdr_list.append([t/60, ssdr/1024/1024/1024])
    msft = pd.DataFrame(ssdr_list, columns=['Time [min]', 'SSDR Usage [Gb]'])
    plt.boxplot(msft['SSDR Usage [Gb]'], whis=(0,100))
    plt.title("SSDR Sizing -- Summary")
    plt.xlabel("Satellite")
    plt.ylabel("SSDR Usage [Gb]")
    plt.savefig('ssdr.png')
    plt.clf()

    plt.plot(msft['Time [min]'], msft['SSDR Usage [Gb]'])
    plt.title("SSDR Sizing -- Time Series")
    plt.xlabel("Time [min]")
    plt.ylabel("SSDR Usage [Gb]")
    plt.savefig('ssdr_timeseries.png')
    plt.clf()

    space_latency = []
    while not ground.empty():
        gi = ground.get()
        space_latency.append((gi[1]-gi[0])/60/60)

    msft = pd.DataFrame(space_latency, columns=['Latency'])
    plt.boxplot(msft.Latency, whis=(0,100))
    plt.title("Space Latency")
    plt.xlabel("Satellite")
    plt.ylabel("Latency [hr]")
    plt.savefig('ssdr_box.png')

    # File output
    with open('latency.csv', 'w') as file:
        file.write('latency')
        for i in space_latency:
            file.write("\n"+i)
