#!/usr/bin/env python

from ROOT import TMVA, TFile, TTree, TCut
from subprocess import call
from os.path import isfile
import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'

import numpy as np
import keras
import h5py

from keras.models import Sequential
from keras.layers import Dense, Activation, Conv2D, MaxPooling2D, Flatten, Dropout
from keras.regularizers import l2
from keras.optimizers import SGD
from keras.callbacks import ModelCheckpoint

# Functions
def getKerasModel(inputDim, modelName, layerSize = 100, nLayers = 5, dropValue = 0.5):

    model = Sequential()

    model.add(Dense(layerSize, activation='relu', kernel_initializer='normal', input_dim=inputDim))
    model.add(Dropout(dropValue))

    for i in range(1, nLayers):

        model.add(Dense(layerSize, activation='relu', kernel_initializer='normal'))
        model.add(Dropout(dropValue))

    model.add(Dense(2, activation='softmax'))

    # Set loss and optimizer
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])

    # Store model to file
    model.save(modelName)
    model.summary()

    return modelName


# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

DNNFLAG= True
BDTFLAG = True

year = sys.argv[1]

if year != '2016' and year != '2017':
    print 'Invalid argument year'
    sys.exit()

# Load data
file = './BsMC/ntuBsMC' + year +'.root'
#file ='ntu.root'

data = TFile.Open(file)

tree = data.Get('PDsecondTree')

cut = 'osMuon==1'
cut += '&&!TMath::IsNaN(muoDxy)'

cutSgn = cut + '&&osMuonTag==1'
cutBkg = cut + '&&osMuonTag==0'

# Prepare factory

name = 'OsMuon' + year

outputName = 'TMVA' + name + '.root'

output = TFile.Open(outputName, 'RECREATE')

factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification')

dataloader = TMVA.DataLoader('dataset')

varList = [
    ('muoPt', 'F')
    ,('muoEta', 'F')
    ,('muoDxy', 'F')
    ,('muoDz', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDR_B', 'F')
    ,('muoPFIso', 'F')
    ,('muoJetPt', 'F')
    ,('muoPtRel', 'F')
    ,('muoDrJet', 'F')
    ,('muoEnergyRatio', 'F')
    ,('muoJetCSV', 'F')
    ,('muoJetDFprob', 'F')
    ,('muoJetSize', 'I')
    ,('muoQCone', 'F')
    ]

for var in varList:
    if year == '2016' and var[0] == 'muoJetDFprob':
        continue
    dataloader.AddVariable( var[0], var[1] )
    nVars += 1

# prepare dataloader

dataloader.AddSignalTree(tree, 1.0)
dataloader.AddBackgroundTree(tree, 1.0)

dataloaderOpt = 'SplitMode=Random:NormMode=NumEvents:V'

dataloader.PrepareTrainingAndTestTree(TCut(cutSgn), TCut(cutBkg), dataloaderOpt)

# Define Keras Model
if DNNFLAG:
    modelName = getKerasModel(nVars, 'model' + year + '.h5')

    # Book methods
    dnnOptions = '!H:!V:NumEpochs=1000:TriesEarlyStopping=50:BatchSize=32:FilenameModel='

    dnnOptions = dnnOptions + modelName

    preprocessingOptions = ':VarTransform=N'
    preprocessingOptions += ',G('

    for var in varList:
        if year == '2016' and var[0] == 'muoJetDFprob':
            continue
        if var[1] == 'F':
            preprocessingOptions += var[0] + ','

    preprocessingOptions = preprocessingOptions[:-1]
    preprocessingOptions +=  '),N'

    dnnName = 'DNN' + name

    factory.BookMethod(dataloader, TMVA.Types.kPyKeras, dnnName, dnnOptions + preprocessingOptions)

# BDT
if BDTFLAG:
    bdtOptions = '!H:!V:UseBaggedBoost:BaggedSampleFraction=0.2:NTrees=600:MaxDepth=6:nCuts=50:MinNodeSize=1.5%:BoostType=RealAdaBoost:AdaBoostBeta=0.3:DoBoostMonitor=True'
    bdtName = 'BDT' + name
    factory.BookMethod(dataloader, TMVA.Types.kBDT, bdtName, bdtOptions)


# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
