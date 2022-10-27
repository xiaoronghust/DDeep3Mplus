# DDeep3Mplus

A docker-powered weakly supervised learning framwork for neuron segmentation of fMOST and public BigNeuron datasets, based on DDeep3M network. 

Thanks for the great work from DDeep3M!

1. build the Docker image for cuda-9 and cudnn-7
   cd cuda-9.0-cudnn7-devel 
   docker build -t cuda9-cudnn7 .
  
2. build the DDeep3M image
   cd ../ddeep3m
   docker build -t ddeep3m .
   
3. run the DDeep3M image and obtain an interactive bash prompt in the container
   nvidia-docker run -it ddeep3m:latest bash

4. pre-process the training data of fMOST with augmentation
   ./PreprocessTrainingData.m ./MOST/train/images ./MOST/train/labels ./preparation/
   
5. run the model with training data
   ./runtraining.sh --numiterations 1000 ./preparation/ ./MOST_trainout
   
6. predict the test data with trained model
   ./runprediction.sh ./MOST_trainout ./MOST/test/images ./MOST_predictout/ 
   
#License
See LICENSE for DDeep3M
