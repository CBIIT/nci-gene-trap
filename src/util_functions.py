import sys
import os
import ConfigParser
import skimage.io
import numpy as np
from skimage import exposure
from matplotlib import cm
from PIL import Image, ImageDraw, ImageFont
import math
import logging




#Author: George Zaki
#Common utility functions

def error_exit(error_message):
  """Print an error message and exit the program with the error code 1.
    """
  error_message = error_message + "\n"
  sys.stderr.write(error_message)
  exit(1)



def read_config(config_file, mandatory_attr=None):
  """
  Parse and read a the config file. Return a dictionary with the arguments
    ARGS: 
      config_file (str): The file name of the configuration file
      mandatory_attr (list of str): The mandatory attributes that should exist in the config file. 
        If it is set to None, check all manadatory attributes.

    RETURNS:
      options_dic (dictionary): Contrains the attribute/value of the configuration parameters
  """
  if not os.path.isfile(config_file):
    raise Exception ('Error the configuration file "{0}" does not exist.'.format(config_file))

  config = ConfigParser.ConfigParser()
  #config.optionxform = str

  try:
    config.read(config_file)
  except ConfigParser.Error as error:
    error_exit('ERROR: unable to parse "config_file"\nERROR details: {0}'.format(error))

  if mandatory_attr == None:
    mandatory_attr = ['stack-dir', 'images-prefix', 'images-suffix', 'processing-dir', 'override-registration',\
      'segmented-rnas', 'point-or-index', 'registration-file', 'destination', 'stride', 'resolution']

  #Check the input file contains the mandatory sections:
  if not "arguments" in config.sections():
    error_exit('ERROR: The config file does not contain the mandatory section {0}\n'.format(section))


  arguments = config.items("arguments")

  options_dic = dict(arguments)
  option_names=options_dic.keys()

  #Verify all options are present
  for item in mandatory_attr: 
    if not item in option_names: 
      error_exit('ERROR: Section "{0}" does not contain the mandatory attribute "{1}"\n'\
        .format("arguments", item))

  #Parse override-registration as a boolean argument
  if 'override-registration'  in option_names:
    options_dic['override-registration'] = config.getboolean('arguments', 'override-registration') 

  return options_dic


class HighlightPoints:
  """
  A class to highlight different objects in an image by drawing a shape around it and giving it an ID.
  """

  def __init__(self, nobjects=1, alpha_in=0.9):
    """Initialize the colour map that will be used to give a different colour for every object in the image.
       Every object will have an ID starting from 0 to 254.

    Args:
      nobjects (int): The number of objects in the image
      alpha_in (float): the transparance of the circles. It should be between 0 and 1. 
    """


    if not isinstance(nobjects, int):
      raise Exception('nobjects should be a positive integer less than 255. Received:"{0}"'.format(nobjects))

    if not (nobjects > 0 and nobjects < 255):
      raise Exception('nobjects should be a positive integer less than 255. Received:"{0}"'.format(nobjects))

    if alpha_in < 0 or alpha_in > 1: 
      raise Exception ('ERROR: alpha_in should be between 0 and 1. Received:"{0}"'.format(alpha_in))


    #Distribute the colour over the 255 range for the color map
    color_step =  255  / float(nobjects) 
    streched_ids = [ int(math.floor(float(set_id) * color_step)) for set_id in range (0,nobjects)] 
    scalar_colors = [  cm.Accent(s_id, bytes=True, alpha = alpha_in) for s_id in streched_ids] 

    rgb_color_map = zip(range(0, nobjects), scalar_colors)

    #rgb_color_map = [ (particle, scalar_color[set_id]) \
    #  for particle, set_id in zip(range(0, nobjects), range(0, nobjects))] 

    #Create the dictionary that maps an object to a colour. 
    self.color_dict = dict(rgb_color_map)

    #define a font
    self.fnt = ImageFont.load_default()




  def highlight_points(self,image_in, image_out, points):
    """ Highlight the points provided in image_in and generate the file image_out.
       
    Args:
      image_in (str): The path to the input image.
      image_out (str): The path to the output image.
      points  (list or tuple): A list of tuples in the format 
        (float, float, int) for the values of (x_pos,y_pos,object_id)i 
        for the points to be highlighted.
    """

    if not os.path.exists(image_in):
      raise Exception('The input image:"{0}" does not exist.'.format(image_in))


    points_list = points
    if  not isinstance (points, list):
      points_list  = [points]
    

    temp_image_name='temp8bit.tif' 

    #First, read the image using skimage as it supports all types then convert it to unit8.
    original_image = skimage.io.imread(image_in)
    

    #Enhance the contrast of the image.
    v_min, v_max = np.percentile(original_image, (0.2, 99.8))
    better_contrast = exposure.rescale_intensity(original_image, in_range=(v_min, v_max))
    new_image = better_contrast * 255.0 / better_contrast.max()
    skimage.io.imsave(temp_image_name, new_image.astype('uint8'))

    #Convert to frame to  RGBA
    orig_frame_image = Image.open(temp_image_name)
    rgba_frame = Image.new("RGBA", orig_frame_image.size)
    rgba_frame.paste(orig_frame_image)
    
    #At this point I do not need the temp image.
    os.remove(temp_image_name)
    contours_frame = Image.new('RGBA', rgba_frame.size, (255, 255, 255, 0))
    draw= ImageDraw.Draw(contours_frame)

    
    # Create a colored circle around every point  
    for x_pos, y_pos, obj_id in points_list:
      if x_pos  < 0 or x_pos >  rgba_frame.size[0] or y_pos  < 0 or y_pos >  rgba_frame.size[1]:
        logging.warning("The point ({0}, {1}) falls outside the image boundry:{2}".format(x_pos,y_pos, rgba_frame.size)) 
      top_x = x_pos - 10
      top_y = y_pos -10
      bottom_x = x_pos + 10 
      bottom_y = y_pos + 10 
      color = self.color_dict[obj_id]
      draw.ellipse((top_x, top_y, bottom_x , bottom_y), outline=color)

      # typethe particle id at the lower corner
      draw.text((bottom_x,bottom_y), str(int(obj_id)), font=self.fnt, fill=color)

    #Blend the original image with the contours. 
    out = Image.alpha_composite(rgba_frame, contours_frame)

    # Save the new file in the output directory
    out.save(image_out)
     





