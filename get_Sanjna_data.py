import os
import scipy.io as sio
import numpy as np
import pandas as pd
import json

folder='E:/Kruttika/LabWork/Code_python/Behavior-Modelling-all/ctbs_data_Sanjna/WorkData/'

print(os.listdir(folder))

to_save={'subject_num':[], 'block_num':[], 'trial_num':[], 'session_label':[], 'session_idx':[], 'cue_validity':[], 'cue_location':[], 'stimuli':[],
            'left_change':[],'right_change':[], 'choice':[],'rt':[]}

subjects=[3,6,7,8,9,10,12,13,16,17,19,20,24,25,26,27,28,31,32,34,36,39,40,47,54,55,58,59]
sess_label=['Sham','Stim']

print(enumerate(subjects))
num_trials=50

for sub_idx,sub in enumerate(subjects):
    for sess_idx,sess in enumerate(sess_label):
        path=folder+'WorkData - S'+str(sub)+sess+'.mat'
        matdata=sio.loadmat(path)
        for block_num in range(5):
            start_tr=block_num*num_trials
            end_tr=(block_num+1)*num_trials
            sess_data=matdata['CData']
            print(path, start_tr,end_tr,sess_data[start_tr:end_tr,:].shape[0])
            to_save['subject_num'].extend([sub_idx]*num_trials)
            to_save['block_num'].extend([block_num]*num_trials)
            to_save['trial_num'].extend(sess_data[start_tr:end_tr,0].astype(int).tolist())
            to_save['session_label'].extend([sess]*num_trials)
            to_save['session_idx'].extend([sess_idx]*num_trials)
            to_save['cue_validity'].extend(sess_data[start_tr:end_tr,2].astype(int).tolist())
            for j in sess_data[start_tr:end_tr,1]:
                if j==0:
                    to_save['cue_location'].append(-1)
                else:
                    to_save['cue_location'].append(1)

            to_save['stimuli'].extend(sess_data[start_tr:end_tr,8].astype(int).tolist())
            to_save['left_change'].extend(sess_data[start_tr:end_tr,5].astype(int).tolist())
            to_save['right_change'].extend(sess_data[start_tr:end_tr,6].astype(int).tolist())
            to_save['choice'].extend(sess_data[start_tr:end_tr,9].astype(int).tolist())
            to_save['rt'].extend(sess_data[start_tr:end_tr,10].tolist())

for key in to_save.keys():
	print(key,len(to_save[key]))

json.dump(to_save,open('E:/Kruttika/LabWork/Code_python/Behavior-Modelling-all/data/Sanjna_data.json','w'))

df=pd.DataFrame(to_save)
print(df)
