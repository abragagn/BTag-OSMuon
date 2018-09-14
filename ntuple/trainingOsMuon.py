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
def getKerasModel(inputDim, modelName, layerSize = 100, , nLayers = 5, dropValue = 0.5):

    model = Sequential()

    model.add(Dense(layerSize, activation='relu', kernel_initializer='normal', input_dim=inputDim))
    model.add(Dropout(dropValue))

    for i in range(1, nLayers)

        model.add(Dense(layerSize, activation='relu', kernel_initializer='normal'))
        model.add(Dropout(dropValue))

    model.add(Dense(2, activation='softmax'))

    # Set loss and optimizer
    sgd = SGD(lr=0.1, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])

    # Store model to file
    model.save(modelName)
    model.summary()

    return modelName



# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

year = sys.argv[1]

if year != '2016' and year != '2017':
    print 'Invalid argument year'
    sys.exit()

name = 'OsMuon' + year

outputName = 'TMVA' + name + '.root'

output = TFile.Open(outputName, 'RECREATE')

factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification')

# Load data
file = './BsMC/ntuBsMC' + year +'.root'
#file ='ntu.root'

data = TFile.Open(file)

tree = data.Get('PDsecondTree')

dataloader = TMVA.DataLoader('dataset')

# add variables (name, type, useForNoJet)
varList = [
    ('muoPt', 'F', True)
    ,('muoEta', 'F', True)
    ,('muoDxy', 'F', True)
    ,('muoDz', 'F', True)
    ,('muoSoftMvaValue', 'F', True)
    ,('muoDR_B', 'F', True)
    ,('muoPFIso', 'F', True)
    ,('muoJetPt', 'F', False)
    ,('muoPtRel', 'F', False)
    ,('muoDrJet', 'F', False)
    ,('muoEnergyRatio', 'F', False)
    ,('muoQJet', 'F', False)
    ,('muoJetCSV', 'F', False)
    ,('muoJetDFprob', 'F', False)
    ,('muoJetSize', 'I', False)
    ,('muoQCone', 'F', False)
    ]

nVars = 0
nVarsJet = 0
nVarsNoJet = 0
JetVars = ''
noJetVars = ''

for var in varList:
    if year == '2016' and var[0] == 'muoJetDFprob':
        continue
    dataloader.AddVariable( var[0], var[1] )
    JetVars += var[0] + ':'
    nVars += 1
    nVarsJet += 1
    if var[2] == True:
        noJetVars += var[0] + ':'
        nVarsNoJet += 1

JetVars = JetVars[:-1]
noJetVars= noJetVars[:-1]

# prepare dataloader

mycutgen = TCut('osMuon==1&&!TMath::IsNaN(muoDxy)')

mycuts = mycutgen + TCut('osMuonTag==1')
mycutb = mycutgen + TCut('osMuonTag==0')

dataloader.AddSignalTree(tree, 1.0)
dataloader.AddBackgroundTree(tree, 1.0)

dataloaderOpt = 'nTrain_Signal=-1:nTrain_Background=-1:nTest_Signal=-1:nTest_Background=-1:SplitMode=Random:NormMode=NumEvents:!V'

dataloader.PrepareTrainingAndTestTree(mycuts, mycutb, dataloaderOpt)

# Define Keras Model
modelName = getKerasModel(nVars, 'model' + year + '.h5')
modelNameJet = getKerasModel(nVarsJet, 'modelJet' + year + '.h5')
modelNameNoJet = getKerasModel(nVarsNoJet, 'modelnoJet' + year + '.h5')

# Book methods
dnnOptions = '!H:!V:NumEpochs=1000:TriesEarlyStopping=50:BatchSize=32:LearningRateSchedule=[5,0.01;100,0.001]:FilenameModel='

dnnOptionsJet = dnnOptions + modelNameJet
dnnOptionsNoJet = dnnOptions + modelNameNoJet
dnnOptions = dnnOptions + modelName

transOptJet = ':VarTransform=N'
transOptJet += ',G('
transOptnoJet = transOptJet

for var in varList:
    if year == '2016' and var[0] == 'muoJetDFprob':
        continue
    if var[1] == 'F':
        transOptJet += var[0] + ','
        if var[2] == True:
            transOptnoJet += var[0] + ','

transOptJet = transOptJet[:-1]
transOptJet +=  '),N'

transOptnoJet = transOptnoJet[:-1]
transOptnoJet +=  '),N'

dnnName = 'DNN' + name

factory.BookMethod(dataloader, TMVA.Types.kPyKeras, dnnName, dnnOptions)

dnnCat = factory.BookMethod(dataloader, TMVA.Types.kCategory, 'DnnCategory', '')
dnnCat.AddMethod(TCut('muoJetSize>=0'), JetVars, TMVA.Types.kPyKeras, dnnName + 'Jet', dnnOptionsJet + transOptJet)
dnnCat.AddMethod(TCut('muoJetSize<0'), noJetVars, TMVA.Types.kPyKeras, dnnName + 'noJet', dnnOptionsNoJet +  transOptnoJet)

# BDT
bdtOptions = '!H:!V:UseBaggedBoost:BaggedSampleFraction=0.6:NTrees=600:MaxDepth=4:nCuts=50:MinNodeSize=1.5%:BoostType=RealAdaBoost:AdaBoostBeta=0.3:DoBoostMonitor=True'
bdtName = 'BDT' + name
factory.BookMethod(dataloader, TMVA.Types.kBDT, bdtName, bdtOptions)

bdtCat = factory.BookMethod(dataloader, TMVA.Types.kCategory, 'BdtCategory', '')
bdtCat.AddMethod(TCut('muoJetSize>=0'), JetVars, TMVA.Types.kBDT, bdtName + 'Jet', bdtOptions)
bdtCat.AddMethod(TCut('muoJetSize<0'), noJetVars, TMVA.Types.kBDT, bdtName + 'noJet', bdtOptions)


# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
