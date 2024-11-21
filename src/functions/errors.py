def check_targets_error(config):
    if len(config.scenario.targets.lats) != len(config.scenario.targets.lons):
        raise NotImplementedError('Error with targets.')
    return None

def check_sat_num(config):
    if config.constellation.sats < 0:
        raise NotImplementedError('Not enough satellites.')
    return None
