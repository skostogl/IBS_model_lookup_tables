# %%
import tree_maker
import pandas as pd
from tree_maker import NodeJob
from tree_maker import initialize
import time
import os
from pathlib import Path
import itertools
import numpy as np
import yaml
from user_defined_functions import generate_run_sh
from user_defined_functions import generate_run_sh_htc

# Import the configuration
config=yaml.safe_load(open('config_sequences.yaml'))

En      = 7000.0
vrf     = 12.0
npart   = 1.8e11
on_disp = 1
on_x1   = 160
# Find all Run3 optics and store in dictionary
optics_path_2023 = "/afs/cern.ch/eng/lhc/optics/runIII/RunIII_dev/2022_V5/PROTON/"

path = optics_path_2023
files = os.listdir(path)
files = [i for i in files if 'forLSA' not in i and 'READ' not in i]
files = [path + '/' + i for i in files]

phi_tot         = []
betax_tot       = []
file_name_tot   = []
optics_name_tot = []

# Find beta* in all optics files
for counter, file in enumerate(files):
    flag=True
    with open(file) as fp:
       line = fp.readline()
       cnt = 1
       while flag:
           line = fp.readline()
           cnt += 1
           if 'betx_IP1' in line:
            flag=False
            betax_IP1 = float(line.split('=')[-1].split(';')[0])
            if (betax_IP1<=1.5) & (betax_IP1>=0.2) & (not file.split('/')[-1]=='opticsfile.22'):
                betax_tot.append(betax_IP1)
                file_name_tot.append(file)
                optics_name_tot.append(file.split('/')[-1])

df_optics = pd.DataFrame({'beta':betax_tot,
                          'optics':optics_name_tot,
                          'path': optics_path_2023})
df_optics['on_x1'] = on_x1

sequence_folder = os.getcwd() + '/sequences'
sequence_dir    = 'create_sequences'
os.makedirs(sequence_folder, exist_ok=True)

# Create tree
children={}
for optics_job, (idx, row) in enumerate(df_optics.iterrows()):
    current_job = f"Run3_beta{row['beta']}_nrj{En}_phi{row['on_x1']}_VRF{vrf}_ondisp{on_disp}"
    children[f'{sequence_dir}/{current_job}'] = {
                                    'optics_file':f"{row['path']}{row['optics']}",
                                    'beam_npart':npart,
                                    'beam_energy_tot':En,
                                    'vrf_total': vrf,
                                    'knob_settings':{'on_x1':row['on_x1'], 'on_x5':row['on_x1'], 'on_disp':on_disp},
                                    'save_to':f"{sequence_folder}/{current_job}.seq",
                                    }

if config['root']['use_yaml_children']== False:
    config['root']['children'] = children
config['root']['setup_env_script'] = os.getcwd() + '/../miniconda/bin/activate'

# Create tree object
start_time = time.time()
root = initialize(config)
print('Done with the tree creation.')
print("--- %s seconds ---" % (time.time() - start_time))

# From python objects we move the nodes to the file-system.
start_time = time.time()
#root.make_folders(generate_run_sh)
root.make_folders(generate_run_sh_htc)
print('The tree folders are ready.')
print("--- %s seconds ---" % (time.time() - start_time))

os.rename("tree_maker.json", "tree_maker_sequences.json")
os.rename("tree_maker.log", "tree_maker_sequences.log")
