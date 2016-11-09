
from pandas import DataFrame as df
import pandas
import numpy as np
import math
from PIL import Image
import trackpy
import os


def check_file(filename): 
  """
    Checks if a file exists and raise and exception with the appropriate error
    Argument:
      filename: the file to be checked.
  """

  if not os.path.exists(filename):
    raise Exception("The file {0} does not exist.".format(os.path.abspath(filename)))

def process_segmented_stack(particles_analysis, output_dir, image_stack_file, segmented_stack_file, resolution=None):
  """
    Track the segmented cells. Generate a floder for every tracked cell where
    the masked and cropped cell frames will reside. Only particles that appear
    in 500 frames will be considered. For these particles, the missing frames
    be filled with the segmentation information of the nearest identified
    frame. 


    Arguments:
      particles_analysis: A csv file generated from imageJ containing metadata 
        about the segmented paticles. Columns must include the series 'X', 'Y',
        "Slice", "BX", "BY", "Width", and "Height"

      output_dir: The directory that will contain the cropped segmented cell.
        Every cell will have its own directory with the name /output_dir/cell_<id>

      image_stack_file: A tif stack (512 frames) of 8-bit precision. 
      
      segmented_stack_file: A tif stack (512 frames) for the 8-bit segmented cell. This
        is a binary image.

      resolution: A tuple specifiying the image resolution in X and Y dimensions 
        If not passed the resolution will be infered from the image_stack metadata.

  """

  #some constants
  # Number of frames in the input image.
  n_frames = 512

  # The extension ratio for the cropped cells 
  extention = 1.10
  
  #Verify the inputs: 
  check_file(particles_analysis)
  check_file(image_stack_file)
  check_file(segmented_stack_file)

  mask_stack = Image.open(segmented_stack_file)

  #Make sure the image_stack is '8b'
  image_stack = Image.open(image_stack_file)
  if image_stack.mode != 'P':
    raise Exception('The image stack {0} does not have the correct 8-bit mode. Mode = {1}'.\
      format(image_stack_file,image_stack.mode))

  # Make sure the image stack have 512 frames.
  #This is called multi-image tiff
  if image_stack.n_frames != 512:
    raise Exception("Image stack {0} only contains {1} frames out of 512."\
      .format(image_stack_file, image_stack.n_frames))

  if mask_stack.n_frames != 512:
    raise Exception("Mask stack {0} only contains {1} frames out of 512."\
      .format(mask_stack_file, mask_stack.n_frames))

  #Make sure the image resolution is passed
  if resolution == None:
    if not "resolution" in  image_stack.info:
      raise Exception("Can not retrive image resolution from {0}".format(image_stack_file))
    else:
      x_voxel_size = 1 / image_stack.info['resolution'][0]
      y_voxel_size = 1 / image_stack.info['resolution'][1]
  else: 
    x_voxel_size = resolution[0] 
    y_voxel_size = resolution[1] 


  # Create the output directory
  if not os.path.exists(output_dir):
    os.mkdir(output_dir) 

  #Track the particles
  tracks_info = pandas.DataFrame.from_csv(particles_analysis)
  tracks_info.rename(columns={'X': 'x', 'Y': 'y', 'Slice' : 'frame'}, inplace=True)
  tracks = trackpy.link_df(tracks_info, search_range =  5, memory = 3)
  particles_ids = tracks.particle.unique()
  
 
  #Pick the best tracked cells.
  sizes=[]
  good_ids = []
  for x in np.nditer(particles_ids):
      n_frames = tracks[tracks.particle == int(x)].shape[0]
      sizes.append(n_frames)
      if n_frames > 500:
          good_ids.append(int(x))
      
  sizes.sort()
  print sizes
  print good_ids
      
  good_ids = [int(i) for i in good_ids] 

  #Loop over all the cells
  for cell_id in good_ids:
 
    #Create the cell directory 
    cell_dir = os.path.join(output_dir, "cell-{0}".format(cell_id))
    if not os.path.exists(cell_dir):
      os.mkdir(cell_dir)

    # To get the bounding box will be in image spacing between (BX, BY) and (BY + Width, BY + Height) 
    # Get this information for particle 1
    particle_info = tracks[tracks.particle == cell_id].loc[:,['x','y', 'frame', 'particle', 'BX', 'BY', 'Width', 'Height']]
  
    #LX and LY are teh right lower border of the segmented cell. 
    particle_info['LX'] = particle_info.BX + particle_info.Width
    particle_info['LY'] = particle_info.BY + particle_info.Height

    # Convert the physical indexes to pixel coordinates for cropping
    particle_info['Left'] = particle_info.BX  / x_voxel_size
    particle_info.Left = particle_info.Left.apply(math.floor).apply(int)
    particle_info['Upper'] = particle_info.BY  / y_voxel_size
    particle_info.Upper = particle_info.Upper.apply(math.floor).apply(int)
    particle_info['Right'] = particle_info.LX  / x_voxel_size
    particle_info.Right = particle_info.Right.apply(math.ceil).apply(int)
    particle_info['Lower'] = particle_info.LY  / y_voxel_size
    particle_info.Lower = particle_info.Lower.apply(math.ceil).apply(int)
  
    new_max_width = int(math.ceil(particle_info.Width.max() / x_voxel_size) * extention)
    new_max_height = int(math.ceil(particle_info.Height.max() / y_voxel_size) * extention)
    
    #Loop over all the frames
    for frame_id in range(1,513):


      #Get the frame if it is part of the track, otherwise get the information
      # from the closest frame. 
      frame_info = particle_info[particle_info.frame == frame_id]
      if len(frame_info.index) == 0:

        closest_row_id = (particle_info.frame - frame_id).apply(abs).argmin()  
        frame_info = particle_info.loc[closest_row_id]

        print("Warning: cell:{0} frame_id:{1} not found, using frame_id:{2} instead."\
          .format(cell_id, frame_id, int(frame_info.frame)))

      try: 
        bounding_box = (frame_info.Left, frame_info.Upper, \
          frame_info.Right, frame_info.Lower)
      except:
        print "frame_id=" , frame_id
        raise

      # To move the image object to the frame_id different frame.
      image_stack.seek(frame_id-1)

      #Move the mask to the closest frame id (which can be the frame_id)
      closest_frame_id = int(frame_info.frame)
      mask_stack.seek(closest_frame_id-1)


      # Crop the bounding box from the image and the mask 
      cropped_frame = image_stack.crop(bounding_box)
      mask_frame = mask_stack.crop(bounding_box).convert('1')

      #Create the masked image
      zero_image = Image.new(image_stack.mode, cropped_frame.size)
      masked = Image.composite(cropped_frame.convert('P'),zero_image.convert('P'), mask_frame)

      #cropped_frame.save("cropped.tif")
      #mask_frame.save("mask.tif")
      #masked.save('masked.tif')
     
      #Expand the masked image with the added zero boundry
      expanded_image = Image.new(image_stack.mode,(new_max_width, new_max_height))
      old_width = bounding_box[2]-bounding_box[0]
      old_height = bounding_box[3] - bounding_box[1]
      expanded_image.paste(masked, (int((new_max_width-old_width)/2), int((new_max_height-old_height)/2)))

      output_file_name = os.path.join(cell_dir,"fragment-{0}.tif".format(str(frame_id).zfill(3)  ))
      expanded_image.save(output_file_name)


process_segmented_stack("Results.csv", "segmented-stacks","MAX_Image11-8bit.tif", "MAX_Image11_segmented_filtered.tif")

