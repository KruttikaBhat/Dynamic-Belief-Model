import os
import scipy.io as sio
import numpy as np
import pandas as pd
import json

folder='E:/Kruttika/LabWork/Code_python/Behavior-Modelling-all/behavdata_Varsha/'

print(os.listdir(folder))

to_save={'subject_num':[], 'block_num':[], 'orig_block_num':[],'trial_num':[], 'session_label':[], 'session_idx':[], 'cue_validity':[], 'cue_location':[], 'stimuli':[],
            'left_change':[],'right_change':[], 'choice':[],'rt':[]}

subjects=[24,27,28,29,32,34,35,36,38,39,40,41,42,43,44,45,46,47,48,49,50,51]
sessions=['Results']
sess_label=['Sham']

dashed=[24,29,32,34,35,36,38,39,40,41,42,43,44,45,46,47,48,49,50,51]
not_dashed=[27,28]

for sub_idx,sub in enumerate(subjects):
    for sess_idx,sess in enumerate(sessions):
        path=folder+'S'+str(sub)+'/'+sess+'/'
        for block_count in range(1,11):
            for i in os.listdir(path):
                block_label='block-'+str(block_count)+'_'
                if sub in not_dashed:
                    block_label='block'+str(block_count)+'_'
                if os.path.isfile(os.path.join(path,i)) and block_label in i:
                    data_path=os.path.join(path,i)
                    matdata=sio.loadmat(data_path)
                    block_data=matdata['Data']
                    print(data_path, block_data.shape[0])
                    to_save['subject_num'].extend([sub_idx]*block_data.shape[0])

                    # if block_count in [1,3,5,7,9]:
                    #     to_save['session_label'].extend([sess_label[sess_idx]+'-1']*block_data.shape[0])
                    #     to_save['session_idx'].extend([sess_idx]*block_data.shape[0])
                    #     to_save['block_num'].extend([int((block_count+1)/2-1)]*block_data.shape[0])
                    # elif block_count in [2,4,6,8,10]:
                    #     to_save['session_label'].extend([sess_label[sess_idx]+'-2']*block_data.shape[0])
                    #     to_save['session_idx'].extend([sess_idx+1]*block_data.shape[0])
                    #     to_save['block_num'].extend([int(block_count/2-1)]*block_data.shape[0])

                    to_save['session_label'].extend([sess_label[sess_idx]]*block_data.shape[0])
                    to_save['session_idx'].extend([sess_idx]*block_data.shape[0])
                    to_save['block_num'].extend([block_count-1]*block_data.shape[0])
                    to_save['orig_block_num'].extend([block_count]*block_data.shape[0])
                    to_save['trial_num'].extend(block_data[:,0].astype(int).tolist())
                    to_save['cue_validity'].extend(block_data[:,2].astype(int).tolist())

                    for j in block_data[:,1]:
                        if j==0:
                            to_save['cue_location'].append(-1)
                        else:
                            to_save['cue_location'].append(1)

                    to_save['stimuli'].extend(block_data[:,8].astype(int).tolist())
                    to_save['left_change'].extend(block_data[:,5].astype(int).tolist())
                    to_save['right_change'].extend(block_data[:,6].astype(int).tolist())
                    to_save['choice'].extend(matdata['ResponseKey'][:,1].astype(int).tolist())
                    to_save['rt'].extend(matdata['ResponseTime'][:,0].tolist())

for key in to_save.keys():
	print(key,len(to_save[key]))

json.dump(to_save,open('E:/Kruttika/LabWork/Code_python/Behavior-Modelling-all/data/Varsha_data_10Blocks.json','w'))

df=pd.DataFrame(to_save)
print(df)
print(df['session_idx'].unique())
