from pandas import DataFrame
from scipy import ndimage
from array import array
import math
from numpy.linalg import inv
import numpy
import scipy.ndimage.filters as fi
import os
from util_functions import read_config
import argparse
import pandas
import scipy


# Python implementation of Dan Larson's Local_Background IDL function.
#
# 02/20/2017
# Coded by:
# Prabhakar R. Gudla, Ph.D.
# HiTIF, LRBGE, CCR, NCI

def local_background_python ( pic ):
    #print pic
    # Get the dimensions of the local patch
    
    pic_dim = pic.shape
    # pic is a numpy array, 1st-dim = rows, 2nd-dim = cols
    x_dim = pic_dim[1]
    y_dim = pic_dim[0]
    
    x_border = numpy.zeros(2*x_dim, float)

    x_border[0:x_dim] = pic[0,]
    x_border[x_dim:(2*x_dim)] = pic[y_dim-1,]
    
    x = numpy.zeros(2*x_dim, float)
    x[0:x_dim] = numpy.arange(x_dim)
    x[x_dim:(2*x_dim)] = numpy.arange(x_dim)

    y_border = numpy.zeros(2*y_dim, float)
    y_border[0:y_dim] = pic[:,0]
    y_border[y_dim:(2*y_dim)] = pic[:,x_dim-1]
    
    y = numpy.zeros(2*y_dim, float)
    y[0:y_dim] = numpy.arange(y_dim)
    y[y_dim:(2*y_dim)] = numpy.arange(y_dim)

    # Following the method of Bevington, p. 96
    delta_x = 2*x_dim*(x**2).sum()-(x.sum())**2
    a=(1./delta_x)*( ((x**2).sum()) * (x_border.sum())-(x.sum())*((x*x_border).sum()) )
    b=(1./delta_x)*( 2*x_dim*((x*x_border).sum())-(x.sum())*(x_border.sum()) )

    delta_y = 2*y_dim*(y**2).sum()-(y.sum())**2
    c=(1./delta_y)*( ((y**2).sum()) * (y_border.sum())-(y.sum())*((y*y_border).sum()) )
    d=(1./delta_y)*( 2*y_dim*((y*y_border).sum())-(y.sum())*(y_border.sum()) )

    # The offset which is returned is averaged over each edge in x, and each edge in y.
    # The actual offset needs to be corrected for the tilt of the plane.
    # Then, the 2 offsets are averaged together to give a single offset.
    offset=(a- d*(y_dim-1)/2.0+c - b*(x_dim-1)/2.0 )/2.0

    # Print some values
    #print 'slope:', b, d
    #print 'offset:', offset
    #print 'adjusted: ', a-d*(y_dim-1)/2.0, c-b*(x_dim-1)/2.0
    
    # now define the background plane in terms of the fit parameters
    plane=numpy.zeros((y_dim, x_dim), float)
    for i in range(0, x_dim):
    	for j in range(0,y_dim):
    		plane[j, i]=offset+b*float(i)+d*float(j)
    		
    return plane

# Function for calcuating Gaussian weights
def gkern2(kernlen=21, nsig=3):
    """Returns a 2D Gaussian kernel array."""

    # create nxn zeros
    inp = numpy.zeros((kernlen, kernlen))
    # set element at the middle to one, a dirac delta
    inp[kernlen//2, kernlen//2] = 1
    # gaussian-smooth the dirac, resulting in a gaussian filter mask
    return fi.gaussian_filter(inp, nsig)

def getBkgCorrectedGaussianWeightedIntensity(input_img, imagej_point, windowSize=5, sigma=1.0):
    """Returns a background corrected, Gaussian Weighted, Normalized Intensity for a RNA spot in a 2D-t (3D numpy array, imput_img) at:
    ARGS:
      imput_img (ndarray or image file name): The 2D input image
      imagej_point (tuple float, float):      X-Y coordinates in imagej format (switch x-y axis when we go to numpy) 
      windowSize (int):                       odd patch size  (X,Y will be at the center of the patch)
      sigma (float):                          sigma of the Gaussian Kernel
    """

    image_array = None
    if isinstance(input_img, basestring):
      image_array = scipy.ndimage.imread(input_img)
    else:
      image_array = input_img 

    # Round the coordinates JUST incase they are @ sub-pixel resolution
    ptR = (int(round(imagej_point[0])), int(round(imagej_point[1])))

    # Get half-width of the window size to extract pixels
    pixelWindow = (windowSize/2)

    # Get the gaussian weights
    gWts = gkern2(pixelWindow*2 + 1, sigma)
    numpy.set_printoptions(precision=4)
    #print gWts
    #print numpy.sum(gWts)

    #print ""
    # Dan's Normalization factor
    gWtsSqSum = numpy.sum(numpy.dot(gWts,gWts))

    #Get the image dimensions:
    ims = image_array.shape

    # Let's get the range of pixels in x and y (image space) and also check for 
    # out-of-bounds
    yRange = list(range(ptR[1]-pixelWindow, ptR[1]+pixelWindow+1))
    xRange = list(range(ptR[0]-pixelWindow, ptR[0]+pixelWindow+1))

    # Initialize a NUMPY array (patch) of size windwosSize by windowSize
    pixelsR = numpy.zeros((windowSize, windowSize))    
    # 
    # Populate intensity values into the patch.
    # Not that most efficient way but IS EASY TO READ
    # Note that the 
    for yInd,yval in enumerate(yRange):
        for xInd,xval in enumerate(xRange):
            if ((xval >= 0) and (xval < ims[1])) and ((yval >= 0) and (yval < ims[0])):
                pixelsR[yInd][xInd] = image_array[yval][xval]


                
    #print pixelsR
    # Get MEAN and SUM of pixels in the patch
    pixelsSum = numpy.sum(pixelsR)
    pixelsMean = numpy.mean(pixelsR)

    # Get Gaussian weighted pixel intensity (NO BACKGROUND correction)
    pixelsSumgWt = numpy.sum(numpy.dot(pixelsR,gWts))

    # Get the Background plane for the patch
    pixelsRBkg = local_background_python(pic=pixelsR)
    # subtract the background plane from patch
    pixelsRBkgSub = pixelsR - pixelsRBkg
    # set any negative background pixels within the patch 0
    pixelsRBkgCorr = (pixelsRBkgSub > 0)*pixelsRBkgSub
    # Calculate the background corrected and gaussian weighted intensity
    pixelsRBkgCorrgWt = numpy.sum(numpy.dot(pixelsRBkgCorr,gWts))
    # Normalize background corrected, gaussian weighted intensity usinf Dan's normalization factor
    #print gWtsSqSum
    pixelsRBkgCorrgWtNorm = pixelsRBkgCorrgWt/gWtsSqSum

    return pixelsRBkgCorrgWtNorm



def calculate_track_intensity(image_prefix, image_suffix, tracked_particles, box_size=5, sigma_in=1.0):
  """
  Get the intensity of overy point defined in the tracked_particle dataframe. Append a new colum "intensity_values" with the calculated intensity. 
  The image file will be <image_prefix><frame-id-3digits><image_suffix> 
  Args: 
    image_prefix (str): The prefix of every image.
    image_suffix (str): The image file suffix
    tracked_particles (data_frame): The data frame that contains the partiles. It must have the columns "x_index", "y_index" and "frame". "frame: will be converted to 3 digits.
    box_size(int): odd patch size  (X,Y will be at the center of the patch)
    sigma (float): sigma of the Gaussian Kernel
  """

  if not isinstance(image_prefix, basestring):
    raise Exception("ERROR: The image_prefix {0} should be a string".format(image_prefix))

  if not isinstance(image_suffix, basestring):
    raise Exception("ERROR: The image_suffix {0} should be a string".format(image_suffix))


  tracked_particles['intensity_values']=tracked_particles.apply(lambda row:\
    getBkgCorrectedGaussianWeightedIntensity(image_prefix+str(int(row['frame'])).zfill(3)+image_suffix,\
      (float(row['x_index']),float(row['y_index'])),\
      windowSize=box_size, \
      sigma=sigma_in), \
    axis=1)

  return tracked_particles


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description="A script to calculate the adjusted intesity of the particles, and write the intensities out in a new data frame files ")
  parser.add_argument('config_file', help='The config file for tracking')
  parser.add_argument('tracked_particles', help='The tracked particles csv file. This file should contain the columns "frame", "x_index", and "y_index"')
  parser.add_argument('output_file', help='The output file name.')
  parser.add_argument('--box', dest='box', default=5, type = int, help="The size of one side of the square to be considered around every particle ")
  parser.add_argument('--sigma', default=1, type = float, dest='sigma', help="The standard deviation to be used while calculating the intensity ")


  args = parser.parse_args()

  #Check that the tracked file exists
  tracked_particles_file = args.tracked_particles
  if not os.path.isfile(tracked_particles_file):
    raise Exception ('Error the tracked particles file "{0}" does not exist.'.format(tracked_particles_file))

  #Parse config to get the box size, sigma, prefix and suffix 
  options_dic = read_config(args.config_file, mandatory_attr = ['stack-dir', 'images-prefix', 'images-suffix'])
  images_stack = os.path.abspath(options_dic['stack-dir'])
  file_prefix = options_dic['images-prefix'] 
  file_suffix= options_dic['images-suffix']

  image_path_prefix = os.path.join(images_stack,file_prefix)


  #If box size is passed use it, otherwise look in the dictionary
  box = None
  if args.box != None:
    box = args.box
  elif "box-size" in options_dic.keys():
    box = int(options_dic['box-size'])

  sigma = None
  if args.sigma != None:
    sigma = args.sigma
  elif "sigma" in options_dic.keys():
    sigma = float(options_dic['sigma'])

  #Parse the tracked_particles
  particles = pandas.read_csv(tracked_particles_file)

  #Calculate the new intensities 
  augmented_particles = calculate_track_intensity(image_path_prefix, file_suffix, particles, box_size=box, sigma_in=sigma)

  #Write the results 
  augmented_particles.to_csv(args.output_file, index=False)
  
  print "Wrote:" + os.path.abspath(args.output_file)

