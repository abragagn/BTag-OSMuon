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


# Load data
file = '../BsMC/ntuBsMC2017.root'
#file ='ntu.root'

data = TFile.Open(file)

tree = data.Get('PDsecondTree')

cut = 'osMuon==1'
cut += '&&!isnan(muoDxy)&&!isnan(muoQCone)&&!isnan(muoJetDFprob)'

cutSgn = cut + '&&osMuonTag==1'
cutBkg = cut + '&&osMuonTag==0'

# Prepare factory

nTest = sys.argv[1]

name = 'OsMuon2017test' + nTest

outputName = 'TMVA' + name + '.root'

output = TFile.Open(outputName, 'RECREATE')

factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification')

dataloader = TMVA.DataLoader('dataset')

varListTest1 = [
    ('muoPt', 'F')
    ,('abs_muoEta := abs(muoEta)', 'F')
    ,('muoDxy', 'F')
    ,('abs_muoDz := abs(muoDz)', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDrB', 'F')
    ,('muoPFIso', 'F')
    ,('muoJetConePt := muoJetPt != -1 ? muoJetPt : muoConePt', 'F')
    ,('muoJetConePtRel := muoJetPt != -1 ? muoJetPtRel : muoConePtRel', 'F')
    ,('muoJetConeDr := muoJetPt != -1 ? muoJetDr : muoConeDr', 'F')
    ,('muoJetConeEnergyRatio := muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio', 'F')
    ,('muoJetCSV', 'F')
    ,('muoJetDFprob', 'F')
    ,('muoJetConeSize := muoJetPt != -1 ? muoJetSize : muoConeSize', 'I')
    ,('muoJetConeQ := muoJetPt != -1 ? muoJetQ : muoConeQ', 'F')
    ]

varListTest2 = [
    ('muoPt', 'F')
    ,('abs_muoEta := abs(muoEta)', 'F')
    ,('muoDxy', 'F')
    ,('abs_muoDz := abs(muoDz)', 'F')
    ,('muoSoftMvaValue', 'F')
    ,('muoDrB', 'F')
    ,('muoPFIso', 'F')
    ,('muoJetPt := muoJetPt != -1 ? muoJetPt : muoConePt', 'F')
    ,('muoJetPtRel := muoJetPt != -1 ? muoJetPtRel : muoConePtRel', 'F')
    ,('muoJetDr := muoJetPt != -1 ? muoJetDr : muoConeDr', 'F')
    ,('muoJetEnergyRatio := muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio', 'F')
    ,('muoJetSize := muoJetPt != -1 ? muoJetSize : muoConeSize', 'I')
    ,('muoJetQ := muoJetPt != -1 ? muoJetQ : muoConeQ', 'F')
    ]


if nTest == '1':
    varList = varListTest1

if nTest == '2':
    varList = varListTest2

if nTest == '3':
    varList = varListTest3

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
if DNNFLAG:
    modelName = getKerasModel(nVars, 'model2017.h5')

    # Book methods
    dnnOptions = '!H:!V:NumEpochs=1000:TriesEarlyStopping=50:BatchSize=128:FilenameModel='

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

# BDT
if BDTFLAG:
    bdtOptions = '!H:!V:UseBaggedBoost:BaggedSampleFraction=0.6:NTrees=600:MaxDepth=6:nCuts=50:MinNodeSize=1.5%:BoostType=RealAdaBoost:AdaBoostBeta=0.3:DoBoostMonitor=True:VarTransform=N'
    bdtName = 'BDT' + name
    factory.BookMethod(dataloader, TMVA.Types.kBDT, bdtName, bdtOptions)


# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
