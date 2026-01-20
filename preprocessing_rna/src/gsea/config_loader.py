# %%
import yaml
import os

def load_config(config_name: str = 'config.yaml') -> dict:
    """
    Load configuration from YAML file.
    
    Parameters:
    ----------
    config_name : str
        Name of the config file (default: 'config.yaml')
    
    Returns:
    -------
    dict
        Configuration dictionary
    """
    config_path = os.path.join(
        os.path.dirname(__file__),
        config_name
    )
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config
