import os

def get_package_path():
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
