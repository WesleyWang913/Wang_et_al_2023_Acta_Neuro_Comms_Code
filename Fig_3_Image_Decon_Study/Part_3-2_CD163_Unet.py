#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 13:47:27 2021
@author: wesley
"""
#Read in libraries
import tensorflow as tf
import os
import random
import numpy as np
from tqdm import tqdm 
from skimage.io import imread, imshow
import cv2 as cv

import matplotlib.pyplot as plt

#Set seed for replication
seed = 42
np.random.seed = seed

IMG_WIDTH = 512
IMG_HEIGHT = 512
IMG_CHANNELS = 3

##Note: 2 neural nets need to be formed, 1 for the processes and 1 for the bodies
#For original model performance studies, please see the CD163-GIMP-2 Folder
#Set paths
IMAGE_PATH = 'Image_Tiles'
PROCESS_PATH = 'Process_Tiles'
BODY_PATH = 'Body_Tiles'

#Read in file lists
image_ids = os.listdir(IMAGE_PATH)
image_ids = np.sort(image_ids)
pro_ids = os.listdir(PROCESS_PATH)
pro_ids = np.sort(pro_ids)
bod_ids = os.listdir(BODY_PATH)
bod_ids = np.sort(bod_ids)

#Read in images
images = np.zeros((len(image_ids), IMG_HEIGHT, IMG_WIDTH, IMG_CHANNELS), dtype=np.uint8)
pro_images = np.zeros((len(image_ids), IMG_HEIGHT, IMG_WIDTH, 1), dtype=np.bool)
bod_images = np.zeros((len(image_ids), IMG_HEIGHT, IMG_WIDTH, 1), dtype=np.bool)

print('Resizing imageing images and masks')
for n, id_ in tqdm(enumerate(image_ids), total=len(image_ids)):   
    img = imread(IMAGE_PATH + '/' +  image_ids[n])
    images[n] = img[...,:3]  #Fill empty X_image with values from img
    pro = imread(PROCESS_PATH + '/' + pro_ids[n], as_gray=True)
    pro_images[n]  = np.expand_dims(pro, axis=-1)
    bod = imread(BODY_PATH + '/' + bod_ids[n], as_gray=True)
    bod_images[n] =  np.expand_dims(bod, axis=-1)


image_x = random.randint(0, len(image_ids))
imshow(images[image_x])
plt.show()
imshow(np.squeeze(pro_images[image_x]))
plt.show()
imshow(np.squeeze(bod_images[image_x]))
plt.show()


################################
#Build the model
inputs = tf.keras.layers.Input((IMG_HEIGHT, IMG_WIDTH, IMG_CHANNELS))
s = tf.keras.layers.Lambda(lambda x: x / 255)(inputs)

#Contraction path
c1 = tf.keras.layers.Conv2D(16, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(s)
c1 = tf.keras.layers.Dropout(0.1)(c1)
c1 = tf.keras.layers.Conv2D(16, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c1)
p1 = tf.keras.layers.MaxPooling2D((2, 2))(c1)

c2 = tf.keras.layers.Conv2D(32, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(p1)
c2 = tf.keras.layers.Dropout(0.1)(c2)
c2 = tf.keras.layers.Conv2D(32, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c2)
p2 = tf.keras.layers.MaxPooling2D((2, 2))(c2)
 
c3 = tf.keras.layers.Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(p2)
c3 = tf.keras.layers.Dropout(0.2)(c3)
c3 = tf.keras.layers.Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c3)
p3 = tf.keras.layers.MaxPooling2D((2, 2))(c3)
 
c4 = tf.keras.layers.Conv2D(128, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(p3)
c4 = tf.keras.layers.Dropout(0.2)(c4)
c4 = tf.keras.layers.Conv2D(128, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c4)
p4 = tf.keras.layers.MaxPooling2D(pool_size=(2, 2))(c4)
 
c5 = tf.keras.layers.Conv2D(256, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(p4)
c5 = tf.keras.layers.Dropout(0.3)(c5)
c5 = tf.keras.layers.Conv2D(256, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c5)

#Expansive path 
u6 = tf.keras.layers.Conv2DTranspose(128, (2, 2), strides=(2, 2), padding='same')(c5)
u6 = tf.keras.layers.concatenate([u6, c4])
c6 = tf.keras.layers.Conv2D(128, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(u6)
c6 = tf.keras.layers.Dropout(0.2)(c6)
c6 = tf.keras.layers.Conv2D(128, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c6)
 
u7 = tf.keras.layers.Conv2DTranspose(64, (2, 2), strides=(2, 2), padding='same')(c6)
u7 = tf.keras.layers.concatenate([u7, c3])
c7 = tf.keras.layers.Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(u7)
c7 = tf.keras.layers.Dropout(0.2)(c7)
c7 = tf.keras.layers.Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c7)
 
u8 = tf.keras.layers.Conv2DTranspose(32, (2, 2), strides=(2, 2), padding='same')(c7)
u8 = tf.keras.layers.concatenate([u8, c2])
c8 = tf.keras.layers.Conv2D(32, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(u8)
c8 = tf.keras.layers.Dropout(0.1)(c8)
c8 = tf.keras.layers.Conv2D(32, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c8)
 
u9 = tf.keras.layers.Conv2DTranspose(16, (2, 2), strides=(2, 2), padding='same')(c8)
u9 = tf.keras.layers.concatenate([u9, c1], axis=3)
c9 = tf.keras.layers.Conv2D(16, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(u9)
c9 = tf.keras.layers.Dropout(0.1)(c9)
c9 = tf.keras.layers.Conv2D(16, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c9)
 
outputs = tf.keras.layers.Conv2D(1, (1, 1), activation='sigmoid')(c9)
 
model = tf.keras.Model(inputs=[inputs], outputs=[outputs])
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.summary()

model2 = tf.keras.Model(inputs=[inputs], outputs=[outputs])
model2.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model2.summary()

################################
#Modelcheckpoint
checkpointer = tf.keras.callbacks.ModelCheckpoint('model_for_nuclei.h5', verbose=1, save_best_only=True)

callbacks = [
        tf.keras.callbacks.EarlyStopping(patience=2, monitor='val_loss'),
        tf.keras.callbacks.TensorBoard(log_dir='logs')]

#Process Model
results = model.fit(images, pro_images, validation_split=0.1, batch_size=16, epochs=50, callbacks=callbacks)

#Bodies Model
results_2 = model2.fit(images, bod_images, validation_split=0.1, batch_size=16, epochs=50, callbacks=callbacks)

####################################

#Set path to images to predict
Pred_PATH = 'Tiles'
pred_ids = os.listdir(Pred_PATH)
pred_ids = np.sort(pred_ids)

#Read in images for prediction
test = np.zeros((len(pred_ids), IMG_HEIGHT, IMG_WIDTH, IMG_CHANNELS), dtype=np.uint8)

print('Resizing imageing images and masks')
for n, id_ in tqdm(enumerate(pred_ids), total=len(pred_ids)):   
    img = imread(Pred_PATH + '/' +  pred_ids[n])
    test[n] = img[...,:3]  #Fill empty X_image with values from img

#Body Prediction
preds_test = model.predict(test, verbose=1)

#Save Prediction Maps
directory = "UNET_Maps_Bodies"
parent_dir = "/Users/sketc/Desktop/CD163_Analysis"
path = os.path.join(parent_dir, directory)
os.mkdir(path)
for n, id_ in tqdm(enumerate(pred_ids), total=len(pred_ids)):
    file_name = path + '/' + pred_ids[n] + '_bmap.txt'
    np.savetxt(file_name, np.reshape(preds_test[n], (512,512)))

#Threshold Masks
preds_test_t = (preds_test > 0.1).astype(np.uint8)
pred_mask = preds_test_t.astype(bool)

image_x = random.randint(0, len(pred_ids))
imshow(preds_test[image_x])
plt.show()
imshow(np.squeeze(pred_mask[image_x]))
plt.show()

#Save Masks
directory = "UNET_Masks_Bodies"
parent_dir = "/Users/sketc/Desktop/CD163_Analysis"
path = os.path.join(parent_dir, directory)
os.mkdir(path)
for n, id_ in tqdm(enumerate(pred_ids), total=len(pred_ids)):
    file_name = path + '/' + pred_ids[n] + '_bmask.txt'
    np.savetxt(file_name, np.reshape(pred_mask[n], (512,512)))

#Process Prediction
preds_test2 = model2.predict(test, verbose=1)

#Save prediction maps
directory = "UNET_Maps_Processes"
parent_dir = "/Users/sketc/Desktop/CD163_Analysis"
path = os.path.join(parent_dir, directory)
os.mkdir(path)
for n, id_ in tqdm(enumerate(pred_ids), total=len(pred_ids)):
    file_name = path + '/' + pred_ids[n] + '_pmap.txt'
    np.savetxt(file_name, np.reshape(preds_test2[n], (512,512)))

#Threshold Masks
preds_test2_t = (preds_test2 > 0.1).astype(np.uint8)
pred_mask2 = preds_test2_t.astype(bool)

image_x = random.randint(0, len(pred_ids))
imshow(preds_test2[image_x])
plt.show()
imshow(np.squeeze(pred_mask2[image_x]))
plt.show()

directory = "UNET_Masks_Processes"
parent_dir = "/Users/sketc/Desktop/CD163_Analysis"
path = os.path.join(parent_dir, directory)
os.mkdir(path)
for n, id_ in tqdm(enumerate(pred_ids), total=len(pred_ids)):
    file_name = path + '/' + pred_ids[n] + '_pmask.txt'
    np.savetxt(file_name, np.reshape(pred_mask2[n], (512,512)))


