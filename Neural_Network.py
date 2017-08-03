#!/usr/bin/env python
# coding: utf-8

import os
import sys
from configparser import ConfigParser
import socket

print('Running on Hostcomputer {}'.format(socket.gethostname()))

parser = ConfigParser()
try:
	parser.read('config.cfg')
except:
	raise Exception('Config File is missing!!!!')  
backend = parser.get('Basics', 'keras_backend')
cuda_path = parser.get('Basics', 'cuda_installation')
os.environ["PATH"] += os.pathsep + cuda_path

if not os.path.exists(cuda_path):
	raise Exception('Given Cuda installation does not exist!')

if backend == 'tensorflow':
	print('Run with backend Tensorflow')
	import tensorflow as tf
elif backend == 'theano':
	print('Run with backend Theano')
	import theano
	os.environ["THEANO_FLAGS"] = "mode=FAST_RUN,device=gpu,floatX=float32" 
else:
	raise NameError('Choose tensorflow or theano as keras backend')

import numpy as np
import keras
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Activation, Flatten, Convolution2D,\
 BatchNormalization, MaxPooling2D,Convolution3D,MaxPooling3D
from keras import regularizers
import h5py
import datetime
import argparse
import math
import time
import resource
import shelve
import shutil
from keras_exp.multigpu import get_available_gpus
from keras_exp.multigpu import make_parallel

################# Function Definitions ####################################################################

def parseArguments():

  parser = argparse.ArgumentParser()
  parser.add_argument("--project", help="The name for the Project", type=str ,default='some_NN')
  parser.add_argument("--input", help="Name of the input files seperated by :", type=str ,default='all')
  parser.add_argument("--model", help="Name of the File containing th qe model", type=str, default='simple_CNN.cfg')
  parser.add_argument("--virtual_len", help="Use an artifical array length (for debugging only!)", type=int , default=-1)
  parser.add_argument("--continue", help="Give a folder to continue the training of the network", type=str, default = 'None')
  parser.add_argument("--date", help="Give current date to identify safe folder", type=str, default = 'None')
  parser.add_argument("--ngpus", help="Number of GPUs for parallel training", type=int, default = 1)
  parser.add_argument("--version", action="version", version='%(prog)s - Version 1.0')
  args = parser.parse_args()
  return args

def read_files(input_files, virtual_len=-1):

  input_data = []
  out_data = []
  file_len = []
  print input_files
  for run, input_file in enumerate(input_files):
    data_file = os.path.join(file_location, 'training_data/{}'.format(input_file))

    if virtual_len == -1:
      data_len = len(h5py.File(data_file)['charge'])
    else:
      data_len = virtual_len
      print('Only use the first {} Monte Carlo Events'.format(data_len))

    input_data.append(h5py.File(data_file, 'r')['charge'])
    out_data.append(h5py.File(data_file, 'r')['reco_vals'])
    file_len.append(data_len)

  return input_data, out_data, file_len


def add_layer(model, layer, args, kwargs):
    eval('model.add({}(*args,**kwargs))'.format(layer))
    return


def base_model(conf_model_file):
  model = Sequential()
  with open(conf_model_file) as f:
      args = []
      kwargs = dict()
      layer = ''
      mode = 'args'
      for line in f:
          cur_line = line.strip()
          if cur_line == '' and layer != '':
              add_layer(model, layer, args,kwargs)
              args = []
              kwargs = dict()
              mode = 'args'
              layer = ''
          elif cur_line[0] == '#':
              continue
          elif cur_line == '[kwargs]':
              mode = 'kwargs'
          elif layer == '':
              layer = cur_line[1:-1]
          elif mode == 'args':
              try:
                  args.append(eval(cur_line.split('=')[1].strip()))
              except:
                  args.append(cur_line.split('=')[1].strip())
          elif mode == 'kwargs':
              split_line = cur_line.split('=')
              try:
                  kwargs[split_line[0].strip()] = eval(split_line[1].strip())
              except:
                  kwargs[split_line[0].strip()] = split_line[1].strip()
      if layer != '':
          add_layer(model, layer, args,kwargs)
  
  print(model.summary())
  return model

class MemoryCallback(keras.callbacks.Callback):
    def on_epoch_end(self, epoch, log={}):
        print(' \n RAM Usage {:.2f} GB \n \n'.format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6))
        os.system("nvidia-smi")

def generator(batch_size, input_data, out_data, inds):
  batch_input = np.zeros((batch_size, 1, 21, 21, 51))
  batch_out = np.zeros((batch_size,1))
  cur_file = 0
  cur_event_id = inds[cur_file][0]
  cur_len = 0
  up_to = inds[cur_file][1]
  while True:
    temp_in = []
    temp_out = []
    while cur_len<batch_size:
      fill_batch = batch_size-cur_len
      if fill_batch < (up_to-cur_event_id):
        temp_in.extend(input_data[cur_file][cur_event_id:cur_event_id+fill_batch])
        temp_out.extend(out_data[cur_file][cur_event_id:cur_event_id+fill_batch])
        cur_len += fill_batch
        cur_event_id += fill_batch
      else:
        temp_in.extend(input_data[cur_file][cur_event_id:up_to])
        temp_out.extend(out_data[cur_file][cur_event_id:up_to])
        cur_len += up_to-cur_event_id
        cur_file+=1
        if cur_file == len(inds):
          cur_file = 0
          cur_event_id = inds[cur_file][0]
          cur_len = 0
          up_to = inds[cur_file][1]
          break
        else:
          cur_event_id = inds[cur_file][0]
          up_to = inds[cur_file][1]
    for i in range(len(temp_in)):
      batch_input[i] = temp_in[i]
      batch_out[i] = np.log10(temp_out[i][0])
    cur_len = 0 
    yield (batch_input, batch_out)


if __name__ == "__main__":

#################### Process Command Line Arguments ######################################

  file_location = parser.get('Basics', 'thisfolder')

  args = parseArguments()
  print("\n ---------")
  print("You are running the script with arguments: ")
  for a in args.__dict__:
      print('{} : {}'.format(a, args.__dict__[a]))
  print("--------- \n")

######################## Setup the training variables ########################################################

  if args.__dict__['continue'] != 'None':
    shelf = shelve.open(os.path.join(file_location, 
      args.__dict__['continue'], 
      'run_info.shlf'))

    project_name = shelf['Project']
    input_files = shelf['Files']
    train_inds = shelf['Train_Inds'] 
    valid_inds = shelf['Valid_Inds']
    test_inds = shelf['Test_Inds']
    model = load_model(os.path.join(file_location, args.__dict__['continue'], 'best_val_loss.npy'))
    today = args.__dict__['continue'].split('/')[1]
    print(today)
    shelf.close()

    input_data, out_data, file_len = read_files(input_files.split(':'))

  else:
    project_name = args.__dict__['project']

    if args.__dict__['input'] =='all':
      input_files = os.listdir(os.path.join(file_location, 'training_data/'))
    else:
      input_files = (args.__dict__['input']).split(':')

    tvt_ratio=[float(parser.get('Training_Parameters', 'training_fraction')),
    float(parser.get('Training_Parameters', 'validation_fraction')),
    float(parser.get('Training_Parameters', 'test_fraction'))] 

    ## Create Folders
    if args.__dict__['date'] != 'None':
      today = args.__dict__['date']
    else:
      today = datetime.date.today()
      folders=['train_hist/',
       'train_hist/{}'.format(today),
       'train_hist/{}/{}'.format(today, project_name)]

    input_data, out_data, file_len = read_files(input_files,
     virtual_len = args.__dict__['virtual_len'])

    train_frac  = float(tvt_ratio[0])/np.sum(tvt_ratio)
    valid_frac = float(tvt_ratio[1])/np.sum(tvt_ratio)
    train_inds = [(0, int(tot_len*train_frac)) for tot_len in file_len] 
    valid_inds = [(int(tot_len*train_frac), int(tot_len*(train_frac+valid_frac))) for tot_len in file_len] 
    test_inds = [(int(tot_len*(train_frac+valid_frac)), tot_len-1) for tot_len in file_len] 

    # print(train_inds)
    # print(valid_inds)
    # print(test_inds)

    ### Create the Model
    conf_model_file = os.path.join('Networks', args.__dict__['model'])
    ngpus = args.__dict__['ngpus']

    adam = keras.optimizers.Adam(lr=float(parser.get('Training_Parameters', 'learning_rate')))
    if ngpus > 1 :
      with tf.device('/cpu:0'):
        # define the serial model.
        model_serial = base_model(conf_model_file)

      gdev_list = get_available_gpus()
      print('Using GPUs: {}'.format(gdev_list))
      model = make_parallel(model_serial, gdev_list)
    else:
      model = base_model(conf_model_file)

    model.compile(loss='mean_squared_error', optimizer=adam, metrics=['accuracy'])
    os.system("nvidia-smi")  

    ## Save Run Information
    shelf = shelve.open(os.path.join(file_location,'train_hist/{}/{}/run_info.shlf'.format(today, project_name)))
    shelf['Project'] = project_name
    shelf['Files'] = args.__dict__['input']
    shelf['Train_Inds'] = train_inds
    shelf['Valid_Inds'] = valid_inds
    shelf['Test_Inds'] = test_inds
    shelf.close()

    shutil.copy(conf_model_file, os.path.join(file_location,'train_hist/{}/{}/model.cfg'.format(today, project_name)))

#################### Train the Model #########################################################

  CSV_log = keras.callbacks.CSVLogger( \
    os.path.join(file_location,'train_hist/{}/{}/loss_logger.csv'.format(today, project_name)), 
    append=True)
  
  early_stop = keras.callbacks.EarlyStopping(\
    monitor='val_loss',
    min_delta = int(parser.get('Training_Parameters', 'delta')), 
    patience = int(parser.get('Training_Parameters', 'patience')), 
    verbose = int(parser.get('Training_Parameters', 'verbose')), 
    mode = 'auto')

  best_model = keras.callbacks.ModelCheckpoint(\
    os.path.join(file_location,'train_hist/{}/{}/best_val_loss.npy'.format(today, project_name)), 
    monitor = 'val_loss', 
    verbose = int(parser.get('Training_Parameters', 'verbose')), 
    save_best_only = True, 
    mode='auto', 
    period=1)

  batch_size = ngpus*int(parser.get('Training_Parameters', 'single_gpu_batch_size'))
  model.fit_generator(generator(batch_size, input_data, out_data, train_inds), 
                steps_per_epoch = math.ceil(np.sum([k[1]-k[0] for k in train_inds])/batch_size),
                validation_data = generator(batch_size, input_data, out_data, valid_inds),
                validation_steps = math.ceil(np.sum([k[1]-k[0] for k in valid_inds])/batch_size),
                callbacks = [CSV_log, early_stop, best_model, MemoryCallback()], 
                epochs = int(parser.get('Training_Parameters', 'epochs')), 
                verbose = int(parser.get('Training_Parameters', 'verbose')),
                max_q_size = int(parser.get('Training_Parameters', 'max_queue_size'))
                )


#################### Saving and Calculation of Result for Test Dataset ######################

  print('\n Save the Model \n')
  model.save(os.path.join(\
  file_location,'train_hist/{}/{}/final_network.h5'.format(today,project_name)))  # save trained network

  print('\n Calculate Results... \n')
  res = []
  test_out = []

  for i in range(len(input_data)):
    print('Predict Values for File {}/{}'.format(i, len(input_data)))
    test  = input_data[i][test_inds[i][0]:test_inds[i][1]]
    test_out_chunk = np.log10(out_data[i][test_inds[i][0]:test_inds[i][1],0:1])
    res_chunk= model.predict(test, verbose=int(parser.get('Training_Parameters', 'verbose')))
    res.extend(list(res_chunk))
    test_out.extend(list(test_out_chunk))


  np.save(os.path.join(file_location,'train_hist/{}/{}/test_results.npy'.format(today, project_name)), 
    [res, np.squeeze(test_out)])

  print(' \n Finished .... ')
