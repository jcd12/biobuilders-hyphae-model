#!/usr/bin/env python3
# coding=utf-8
import sys
import time
import argparse
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import copy
import pdb
from skimage.morphology import skeletonize
from skimage import img_as_bool, io, color, morphology, data
from skimage.util import invert


#image = img_as_bool(color.rgb2gray(io.imread('../../image_processing/test_image_colored.jpg')))
image = img_as_bool(color.rgb2gray(io.imread('../../image_processing/2020_fluorescent.jpg')))
skeleton = skeletonize(image)
skeleton_lee = skeletonize(image, method='lee')

print('Pixels:', len(skeleton), '*', len(skeleton[0]))

fig, axes = plt.subplots(1, 3, figsize=(8, 4), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(image, cmap=plt.cm.gray)
ax[0].set_title('original')
ax[0].axis('off')

ax[1].imshow(skeleton, cmap=plt.cm.gray)
ax[1].set_title('skeletonize')
ax[1].axis('off')

ax[2].imshow(skeleton_lee, cmap=plt.cm.gray)
ax[2].set_title('skeletonize (Lee 94)')
ax[2].axis('off')

fig.tight_layout()
plt.show()