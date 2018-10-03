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

CATFLAG = True
FIXFLAG = False

year = sys.argv[1]

if year != '2016' and year != '2017':
    print 'Invalid argument year'
    sys.exit()

# Load data
file = './BsMC/ntuBsMC' + year +'.root'
#file ='ntu.root'

data = TFile.Open(file)
if FIXFLAG:
    temp = TFile.Open("temp.root", 'RECREATE')

tree = data.Get('PDsecondTree')

cut = 'osMuon==1'
cut += '&&!TMath::IsNaN(muoDxy)&&muoEnergyRatio<100'

cutSgn = cut + '&&osMuonTag==1'
cutBkg = cut + '&&osMuonTag==0'

if FIXFLAG:
    treeSgnTrain = tree.CopyTree(cutSgn + '&&(evtNumber%2)==0')
    treeBkgTrain = tree.CopyTree(cutBkg + '&&(evtNumber%2)==0')

    treeSgnTest = tree.CopyTree(cutSgn + '&&(evtNumber%2)==1')
    treeBkgTest = tree.CopyTree(cutBkg + '&&(evtNumber%2)==1')

# Prepare factory

name = 'OsMuon' + year

outputName = 'TMVA' + name + '.root'

output = TFile.Open(outputName, 'RECREATE')

factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification')

dataloader = TMVA.DataLoader('dataset')

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

if FIXFLAG:
    dataloader.AddSignalTree(treeSgnTrain, 1.0, TMVA.Types.kTraining)
    dataloader.AddBackgroundTree(treeBkgTrain, 1.0, TMVA.Types.kTraining)

    dataloader.AddSignalTree(treeSgnTest, 1.0, TMVA.Types.kTesting)
    dataloader.AddBackgroundTree(treeBkgTest, 1.0, TMVA.Types.kTesting)

    cutSgn = ''
    cutBkg = ''

else:
    dataloader.AddSignalTree(tree, 1.0)
    dataloader.AddBackgroundTree(tree, 1.0)

dataloaderOpt = 'SplitMode=Random:NormMode=NumEvents:V'

dataloader.PrepareTrainingAndTestTree(TCut(cutSgn), TCut(cutBkg), dataloaderOpt)

# Define Keras Model
if DNNFLAG:
    modelName = getKerasModel(nVars, 'model' + year + '.h5')
    modelNameJet = getKerasModel(nVarsJet, 'modelJet' + year + '.h5')
    modelNameNoJet = getKerasModel(nVarsNoJet, 'modelnoJet' + year + '.h5')

    # Book methods
    dnnOptions = '!H:!V:NumEpochs=1000:TriesEarlyStopping=50:BatchSize=32:FilenameModel='

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

    if CATFLAG:
        dnnCat = factory.BookMethod(dataloader, TMVA.Types.kCategory, dnnName + 'Category', '')
        dnnCat.AddMethod(TCut('muoJetSize>=0'), JetVars, TMVA.Types.kPyKeras, dnnName + 'Jet', dnnOptionsJet + transOptJet)
        dnnCat.AddMethod(TCut('muoJetSize<0'), noJetVars, TMVA.Types.kPyKeras, dnnName + 'noJet', dnnOptionsNoJet +  transOptnoJet)

# BDT
if BDTFLAG:
    bdtOptions = '!H:!V:UseBaggedBoost:BaggedSampleFraction=0.2:NTrees=600:MaxDepth=6:nCuts=50:MinNodeSize=1.51%:BoostType=RealAdaBoost:AdaBoostBeta=0.3:DoBoostMonitor=True'
    bdtName = 'BDT' + name
    factory.BookMethod(dataloader, TMVA.Types.kBDT, bdtName, bdtOptions)

    if CATFLAG:
        bdtCat = factory.BookMethod(dataloader, TMVA.Types.kCategory, bdtName + 'Category', '')
        bdtCat.AddMethod(TCut('muoJetSize>=0'), JetVars, TMVA.Types.kBDT, bdtName + 'Jet', bdtOptions)
        bdtCat.AddMethod(TCut('muoJetSize<0'), noJetVars, TMVA.Types.kBDT, bdtName + 'noJet', bdtOptions)


# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

if FIXFLAG:
    os.remove("./temp.root")
