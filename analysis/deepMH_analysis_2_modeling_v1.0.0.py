### deepMH - identify methylation heterogeneity from single-cell bisulfite sequencing using deep learning
# Copyright (C) 2022  Dong, et. al.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

### updates
# v1.0.0 - 2022.10.06 - initial version.

# in python3

# source /data/x001/xdong/apps/venv_gpu/bin/activate
# https://datascience.stackexchange.com/questions/21955/tensorflow-regression-model-giving-same-prediction-every-time
# z-score normalization: https://en.wikipedia.org/wiki/Feature_scaling

from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import sys
import string
import numpy as np
import datetime
import dill

import collections
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow.keras import layers
from tensorflow import keras

version="v1.0.0"
date="2022.10.06"

###
parser=argparse.ArgumentParser(description="deepMH_analysis_2_modeling_" + version + ".py\nModeling methylation with deepMH.\nDate: " + date + "; Version: " + version)
parser.add_argument("-i","--input",type=str, required=True, help="Input header of bed file")
parser.add_argument("-b","--num_bp",type=int, default=100, help="No. upstream and downstream basepairs to consider. Default=100")
args=parser.parse_args()

### notes ###

print(args.input)

dillfile = args.input + ".pkl"
dill.load_session(dillfile)

fileout = open(args.input + ".predict2", "w")

nepochs=5
bsize=512

### import data ###
tf.keras.backend.clear_session()

print(data_label[0:20])

## layer DNA sequence ##
data_sequence=[k[1] for k in data_feature]
data_sequence=np.array(data_sequence, dtype="float32")

## layer methylation ##
data_meth=np.array(data_meth, dtype="float32")

# z-score normalization
tmp1=np.array([])
tmp2=np.array([])
for e in data_meth:
    tmp1=np.append(tmp1, e[0])
    tmp2=np.append(tmp2, e[1])

tmp1_mean=tmp1.mean()
tmp1_sd=tmp1.std()
tmp2_mean=tmp2.mean()
tmp2_sd=tmp2.std()

data_meth_z=[]
for e in data_meth:
    data_meth_z.append(np.array([np.array([(e1-tmp1_mean)/tmp1_sd for e1 in e[0]], dtype="float32"), np.array([(e1-tmp2_mean)/tmp2_sd for e1 in e[1]], dtype="float32")], dtype="float32"))

data_meth_z=np.array(data_meth_z, dtype="float32")

data_meth = data_meth_z

### model 4.3.1.3 ###
tf.keras.backend.clear_session()
dna_input=keras.Input(shape=(args.num_bp*2+1,), name="dna")
dna_model=layers.Embedding(input_dim=args.num_bp*2+1, output_dim=64)(dna_input)
dna_model=layers.LSTM(64)(dna_model)
dna_model=layers.Dense(4, activation='relu')(dna_model)

meth_input=keras.Input(shape=(2,min_cpg-1), name="meth")
meth_model=layers.GRU((min_cpg-1)*2)(meth_input)
meth_model=layers.Dense(4, activation='relu')(meth_model)

merged_model=layers.concatenate([dna_model, meth_model])
merged_model=layers.Dense(8, activation='relu')(merged_model)
merged_model=layers.Dense(8, activation='relu')(merged_model)
merged_model=layers.Dense(1)(merged_model)

model = keras.Model(inputs=[dna_input, meth_input], outputs=merged_model, name="combined")
model.summary()

optimizer = tf.keras.optimizers.RMSprop(0.001)
model.compile(loss='mse',
        optimizer=optimizer,
        metrics=['mae', 'mse'])

optimizer = tf.keras.optimizers.RMSprop(0.001)
model.compile(loss='mse',
        optimizer=optimizer,
        metrics=['mae', 'mse'])

model.fit({"dna":data_sequence,"meth":data_meth}, data_label, epochs=nepochs, batch_size=bsize)

pre3=model.predict({"dna":data_sequence,"meth":data_meth})
print(pre3[0:20])

for e in range(0, len(data_label)):
    lineout = str(data_label_chr[e]) + "\t" + str(data_label_pos[e]) + "\t" + str(data_label[e]) + "\t" + str(pre3[e][0]) + "\n"
    fileout.write(lineout)

fileout.close()

