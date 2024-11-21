import time
start_time = time.time()

from functions.filesystem import get_config, make_output_folder
from functions.orbit import calc_look_angle, make_orbits, make_accesses, calc_graze_snr_floor
from functions.errors import check_targets_error, check_sat_num
from functions.targets import make_targets
from functions.dataoperations import calc_data_storage

def main(debug=False):
    config = get_config()

    if config.post_process.run_access:
        # Check targets for errors
        check_targets_error(config)
        # Check sats for errors
        check_sat_num(config)

        # Make a output folder if needed
        make_output_folder()

        # Make orbits
        orb = make_orbits(config, debug)

        # Make Targets
        els, a_f = make_targets(config, debug) # Returns earth locations and file pointers

        # Get payload graze Limit
        graze = []
        for k,v in enumerate(els):
            if config.scenario.targets.target_type[k] == 'optical':
                graze.append(calc_look_angle(config, debug))
            elif config.scenario.targets.target_type[k] == 'comms':
                graze.append(calc_graze_snr_floor(config, debug))
            else:
                raise NotImplementedError('Could not find this target type')

        # Generate accesses
        make_accesses(config, orb, els, a_f, graze, debug)

    if config.post_process.generate_data_curve:
        calc_data_storage(config, debug)

    return None

if __name__ == '__main__':
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
