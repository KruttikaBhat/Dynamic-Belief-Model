import scipy.io as sio
import numpy as np
import pandas as pd
import json

matdat = sio.loadmat('E:/Kruttika/LabWork/Code_python/Behavior-Modelling-all/Behavior-Modelling/data/tACS_40Hz_woETrej.mat')
stim_locations = ['leftPPC', 'rightPPC']
session_types = ['Sham', 'Stim', 'Post']


to_save={'subject_num':[], 'block_num':[], 'trial_num':[], 'session_label':[], 'session_idx':[], 'cue_validity':[], 'cue_location':[], 'stimuli':[],
            'left_change':[],'right_change':[], 'choice':[],'rt':[]}

for subject_num in range(26):
    for j_index, j in enumerate(stim_locations):
        data = matdat[j]['Trials_Info']
        for session_type_index, session_type in enumerate(session_types):
            session_data = np.array(data[0][0][subject_num][0][session_type][0][0])
            session_idx=j_index*len(session_types) + session_type_index
            for block_idx in range(5):
                block_data=session_data[:,:,block_idx]
                to_save['subject_num'].extend([subject_num]*block_data.shape[0])
                to_save['block_num'].extend([block_idx]*block_data.shape[0])
                to_save['trial_num'].extend(block_data[:,0].astype(int).tolist())
                to_save['session_label'].extend([j+'-'+session_type]*block_data.shape[0])
                to_save['session_idx'].extend([session_idx]*block_data.shape[0])
                to_save['cue_validity'].extend(block_data[:,1].astype(int).tolist())
                to_save['cue_location'].extend(block_data[:,2].astype(int).tolist())
                to_save['stimuli'].extend(block_data[:,3].astype(int).tolist())
                to_save['left_change'].extend(block_data[:,4].astype(int).tolist())
                to_save['right_change'].extend(block_data[:,5].astype(int).tolist())
                to_save['choice'].extend(block_data[:,8].astype(int).tolist())
                to_save['rt'].extend(block_data[:,9].tolist())


for key in to_save.keys():
	print(key,len(to_save[key]))

json.dump(to_save,open('E:/Kruttika/LabWork/Code_python/Behavior-Modelling-all/data/Ankita_data.json','w'))

df=pd.DataFrame(to_save)
print(df)
