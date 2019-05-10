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
from keras.optimizers import SGD, Adam
from keras.callbacks import ModelCheckpoint

##### FUNCTIONS

def getKerasModel(inputDim, modelName, nLayers = 5, layerSize = 100, dropValue = 0.25):
    model = Sequential()
    model.add(Dense(layerSize, activation='relu', kernel_initializer='normal', input_dim=inputDim))
    if dropValue != 0:
        model.add(Dropout(dropValue))

    for i in range(1, nLayers):
        model.add(Dense(layerSize, activation='relu', kernel_initializer='normal'))
        if dropValue != 0:
            model.add(Dropout(dropValue))

    # Anything below this point should not be changed in order for the network to work with TMVA
    # Exception: the optimizer
    model.add(Dense(2, activation='softmax'))

    #opt = SGD(lr=0.05, decay=1e-5, momentum=0.9, nesterov=True)
    opt = Adam(lr=0.001)
    model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])
    model.save(modelName)
    model.summary()
    return modelName


# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

# Load data
file = '../ntuBsMC2017_skim.root'

data = TFile.Open(file)

tree = data.Get('PDsecondTree')

# Event wise selection
cut = 'osMuon==1&&hltJpsiMu==1'
#cut += '&&!isnan(muoDxy)&&!isnan(muoJetDFprob)&&!isinf(muoJetEnergyRatio)&&!isinf(muoConeEnergyRatio)'

cutSgn = cut + '&&osMuonTag==1' #Correcly tagged events selection e.g if sign(charge lepton) -> correct flavour 
cutBkg = cut + '&&osMuonTag==0' #Uncorrecly tagged events

# Prepare factory

nTest = sys.argv[1] ## just a label
nLayers = sys.argv[2]
layerSize = sys.argv[3]
dropValue = sys.argv[4]

nLayers = int(nLayers)
layerSize = int(layerSize)
dropValue = float(dropValue)

name = 'OsMuonHLTJpsiMu_test' + nTest

outputName = 'TMVA' + name + '.root'

output = TFile.Open(outputName, 'RECREATE')

factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification')

dataloader = TMVA.DataLoader('dataset')

# variable list
varList = [
    ('muoPt', 'F')
    ,('abs_muoEta := fabs(muoEta)', 'F')
    ,('muoDxy', 'F')
    ,('muoExy', 'F')
    ,('muoDz', 'F')
    ,('muoEz', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDrB', 'F')
    ,('muoPFIso', 'F')
    ,('muoConePt', 'F')
    ,('muoConePtRel', 'F')
    ,('muoConeDr', 'F')
    ,('muoConeEnergyRatio', 'F')
    ,('muoConeQ', 'F')
    ,('muoConeCF', 'F')
    ,('muoConeNCH', 'I')
    ,('muoJetDFprob', 'F')
    ]

varListClean = [
    ('muoPt', 'F')
    ,('muoEta', 'F')
    ,('muoDxy', 'F')
    ,('muoExy', 'F')
    ,('muoDz', 'F')
    ,('muoEz', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDrB', 'F')
    ,('muoPFIso', 'F')
    ,('muoConeCleanPt', 'F')
    ,('muoConeCleanPtRel', 'F')
    ,('muoConeCleanDr', 'F')
    ,('muoConeCleanEnergyRatio', 'F')
    ,('muoConeCleanQ', 'F')
    #,('muoJetDFprob', 'F')
    ]

varList = varListClean

# automatic variable counting and adding
nVars = 0
for var in varList:
    dataloader.AddVariable( var[0], var[1] )
    nVars += 1

# prepare dataloader
# same tree, add selection later
dataloader.AddSignalTree(tree)
dataloader.AddBackgroundTree(tree)

# evtWeight variable in the ntuple, address simulation bias
dataloader.SetWeightExpression( 'evtWeight' );

dataloaderOpt = 'SplitMode=Random:NormMode=NumEvents:V:nTrain_Signal=155331:nTrain_Background=66249'

dataloader.PrepareTrainingAndTestTree(TCut(cutSgn), TCut(cutBkg), dataloaderOpt)

# Create Keras Model
modelName = getKerasModel(nVars, 'model' + name + '.h5', nLayers, layerSize, dropValue)

# Book methods
dnnOptions = '!H:!V:NumEpochs=50:TriesEarlyStopping=10:BatchSize=512:ValidationSize=30%'
dnnOptions = dnnOptions + ':Tensorboard=./logs:FilenameModel=' + modelName
# 

# Preprocessing string creator, loop was for selection of which variable to apply gaussianification
preprocessingOptions = ':VarTransform=N,G,N'
# preprocessingOptions += ',G('

# iVar = 0
# for var in varList:
#     preprocessingOptions += '_V' + str(iVar) + '_' + ','
#     iVar += 1

# preprocessingOptions = preprocessingOptions[:-1]
# preprocessingOptions +=  '),N'

dnnName = 'DNN' + name

factory.BookMethod(dataloader, TMVA.Types.kPyKeras, dnnName, dnnOptions + preprocessingOptions)

# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
