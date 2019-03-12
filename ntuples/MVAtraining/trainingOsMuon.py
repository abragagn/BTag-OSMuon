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

def getKerasModel(inputDim, modelName, layerSize = 100, nLayers = 5, dropValue = 0.25):
    model = Sequential()
    model.add(Dense(layerSize, activation='relu', kernel_initializer='normal', input_dim=inputDim))
    if dropValue != 0:
        model.add(Dropout(dropValue))

    for i in range(1, nLayers):
        model.add(Dense(layerSize, activation='relu', kernel_initializer='normal'))
        if dropValue != 0:
            model.add(Dropout(dropValue))

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
#file ='ntu.root'

data = TFile.Open(file)

tree = data.Get('PDsecondTree')

cut = 'osMuon==1&&hltJpsiMu==1&&muoJetPt != -1'
#cut += ''
cut += '&&!isnan(muoDxy)&&!isnan(muoJetDFprob)&&!isinf(muoJetEnergyRatio)&&!isinf(muoConeEnergyRatio)'

cutSgn = cut + '&&osMuonTag==1'
cutBkg = cut + '&&osMuonTag==0'

# Prepare factory

nTest = sys.argv[1]
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

varListJetCone = [
    ('muoPt', 'F')
    ,('abs_muoEta := fabs(muoEta)', 'F')
    ,('muoDxy', 'F')
    ,('abs_muoDz := fabs(muoDz)', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDrB', 'F')
    ,('muoPFIso', 'F')
    ,('muoJetConePt := muoJetPt != -1 ? muoJetPt : muoConePt', 'F')
    ,('muoJetConePtRel := muoJetPt != -1 ? muoJetPtRel : muoConePtRel', 'F')
    ,('muoJetConeDr := muoJetPt != -1 ? muoJetDr : muoConeDr', 'F')
    ,('muoJetConeEnergyRatio := muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio', 'F')
    ,('muoJetDFprob', 'F')
    ,('muoJetConeSize := muoJetPt != -1 ? muoJetSize : muoConeSize', 'I')
    ,('muoJetConeQ := muoJetPt != -1 ? muoJetQ : muoConeQ', 'F')
    ]

varListJet = [
    ('muoPt', 'F')
    ,('abs_muoEta := fabs(muoEta)', 'F')
    ,('muoDxy', 'F')
    ,('abs_muoDz := fabs(muoDz)', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDrB', 'F')
    ,('muoPFIso', 'F')
    ,('muoJetPt', 'F')
    ,('muoJetPtRel', 'F')
    ,('muoJetDr', 'F')
    ,('muoJetEnergyRatio', 'F')
    ,('muoJetDFprob', 'F')
    ,('muoJetSize', 'I')
    ,('muoJetNF', 'F')
    ,('muoJetCF', 'F')
    ,('muoJetNCH', 'I')
    ]

varList = varListJet


nVars = 0

for var in varList:
    dataloader.AddVariable( var[0], var[1] )
    nVars += 1

# prepare dataloader

dataloader.AddSignalTree(tree, 1.0)
dataloader.AddBackgroundTree(tree, 1.0)

dataloaderOpt = 'SplitMode=Random:NormMode=NumEvents:V'

dataloader.PrepareTrainingAndTestTree(TCut(cutSgn), TCut(cutBkg), dataloaderOpt)

# Define Keras Model
modelName = getKerasModel(nVars, 'model' + name + '.h5', layerSize, nLayers, dropValue)

# Book methods
dnnOptions = '!H:!V:NumEpochs=100:TriesEarlyStopping=20:BatchSize=512:FilenameModel='

dnnOptions = dnnOptions + modelName

preprocessingOptions = ':VarTransform=N'
preprocessingOptions += ',G('

iVar = 0
for var in varList:
    preprocessingOptions += '_V' + str(iVar) + '_' + ','
    iVar += 1

preprocessingOptions = preprocessingOptions[:-1]
preprocessingOptions +=  '),N'

dnnName = 'DNN' + name

factory.BookMethod(dataloader, TMVA.Types.kPyKeras, dnnName, dnnOptions + preprocessingOptions)

# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
