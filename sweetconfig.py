__author__ = 'pschuprikov'

import os
import configparser

if not os.path.exists('sweet.ini'):
    initial_cfg = configparser.ConfigParser()
    initial_cfg['general'] = {}
    initial_cfg['genblasta'] = {}
    with open('sweet.ini', 'w') as configfile:
        initial_cfg.write(configfile)

cfg = configparser.ConfigParser()
cfg.read('sweet.ini')

general_cfg = cfg['general']
genblasta_cfg = cfg['genblasta']