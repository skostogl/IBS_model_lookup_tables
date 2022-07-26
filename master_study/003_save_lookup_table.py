import pandas as pd
import tree_maker
from tree_maker import NodeJob
import time
import numpy as np
import matplotlib.pyplot as plt

def load_tree(filename):
    try:
        root=tree_maker.tree_from_json(filename)
        return root
    except Exception as e:
        print(e)
        print('Probably you forgot to edit the address of you json file...')

root = load_tree(f'tree_maker.json')
save_to = "lookup_table_En_6800_Vrf_12.0_NEW.parquet"

start_time = time.time()

appended_data = []
for counter, node in enumerate(root.descendants[:]):
    try:
        #if True:
        current_df = pd.read_parquet(f'{root.get_abs_path()}/{node}/IBS_output.parquet')
        current_df["exin"] = current_df["exin"]*1e6
        current_df["eyin"] = current_df["eyin"]*1e6

        if counter%500 == 0:
            print("Current study ", counter, node)
        appended_data.append(current_df)
    except:
            print(f"Missing output from {node}")
        
appended_data = pd.concat(appended_data, axis=0)
appended_data.to_parquet(f"{save_to}")
end_time = time.time()
print( f"Time needed to append IBS output (mins): ",(end_time-start_time)/60.)
