import yaml

from box import Box

def get_config():
    config = Box.from_yaml(filename="./config.yml", Loader=yaml.FullLoader) 
    return config

def make_output_folder():
    import os
    if not os.path.exists('access_output'):
        os.makedirs('access_output')
    else:
        import shutil
        shutil.rmtree('access_output')
        os.makedirs('access_output')
    return None
