# %%
import tree_maker
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
config=yaml.safe_load(open('config.yaml'))

sequence    = ["/home/HPC/skostogl/workspace_Jul22/IBS_model_new/master_study/sequences/Run3_beta0.6_nrj7000.0_phi160_VRF12.0_ondisp1.seq"]
En          = [6800.0]
V0max       = [12.0]
blns        = np.arange(0.85, 1.4, 0.05)
blns        = [round(i,6) for i in blns]
emitx       = np.arange(1., 4.1, 0.2)
emity       = np.arange(1., 4.1, 0.2)
b_intensity = np.arange(1e10, 2.4e11, 2e10)
save_to     = 'IBS_output.parquet'

chunks_per_job = 50 # number of scan points in each job

all_jobs = itertools.product(sequence, En, V0max, blns, emitx, emity, b_intensity)
all_jobs = np.array([i for i in all_jobs])
all_jobs_split =np.array([all_jobs[i:i + chunks_per_job, :] for i in range(0, len(all_jobs), chunks_per_job)])

study_name = "IBS_scan_6p8TeV_12MV_log"

children={}
for child in range(len(all_jobs_split)):
    children[f"{study_name}/{child:03}"] = {
                                    'sequence':all_jobs_split[child].T[0].tolist(),
                                    'En': all_jobs_split[child].T[1].tolist(),
                                    'V0max':all_jobs_split[child].T[2].tolist(),
                                    'blns':all_jobs_split[child].T[3].tolist(),
                                    'bunch_intensity':all_jobs_split[child].T[6].tolist(),
                                    'emit_x':all_jobs_split[child].T[4].tolist(),
                                    'emit_y':all_jobs_split[child].T[5].tolist(),
                                    'save_to': f"{os.getcwd()}/{study_name}/{child:03}/{save_to}",
                                    'log_file': f"{os.getcwd()}/{study_name}/{child:03}/tree_maker.log"
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
