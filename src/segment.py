from pandas import DataFrame as df
import pandas
import numpy as np
import math
from PIL import Image
import trackpy
import os
import argparse
from ast import literal_eval
import pdb


def check_file(filename): 
  """
    Checks if a file exists and raise and exception with the appropriate error
    Argument:
      filename: the file to be checked.
  """

  if not os.path.exists(filename):
    raise Exception("The file {0} does not exist.".format(os.path.abspath(filename)))

def read_Dataframe_file(file_name, mandatory_columns):
  """
    Read and Verifies that the Dataframe object contains the mandatory columns
      
    Arguments:
      file_name: The csv file 
      columns:    List of mandatory columns

    Returns:
      pandas.Dataframe object

  """
  pass


def process_segmented_stack(particles_analysis, output_dir, image_stack_file, segmented_stack_file, resolution=None, segmented_rna_points = None):
  """
    Track the segmented cells. Generate a floder for every tracked cell where
    the masked and cropped cell frames will reside. Only particles that appear
    in 98% frames will be considered. For these particles, the missing frames
    be filled with the segmentation information of the nearest identified
    frame. 


    Arguments:
      particles_analysis: A csv file generated from imageJ containing metadata 
        about the segmented paticles. Columns must include the series 'X', 'Y',
        "Slice", "BX", "BY", "Width", and "Height"

      output_dir: The directory that will contain the cropped segmented cell.
        Every cell will have its own directory with the name /output_dir/cell_<id>

      image_stack_file: A tif stack  of 8-bit precision. 
      
      segmented_stack_file: A tif stack (512 frames) for the 8-bit segmented cell. This
        is a binary image.

      resolution: A tuple specifiying the image resolution in X and Y dimensions 
        If not passed the resolution will be infered from the image_stack metadata.

      segmented_rna_pionts: A csv file that contains the information about the
        segmented RNA transcription points. The format should be x, y, intensity,
        frame.  This script will collect the points that falls within a given
        cell bounding box (per frame) and saves a csv file with these points in in
        the corresponding folder of the cell. The points will saved in image
        indexes.

  """

  #some constants

  # The extension ratio for the cropped cells 
  extention = 1.10

  
  #Verify the inputs: 
  check_file(particles_analysis)
  check_file(image_stack_file)
  check_file(segmented_stack_file)
  if segmented_rna_points != None:
    check_file(segmented_rna_points)

  mask_stack = Image.open(segmented_stack_file)

  #Make sure the image_stack is '8b'
  image_stack = Image.open(image_stack_file)
  if image_stack.mode != 'P' and image_stack.mode != 'L':
    raise Exception('The image stack {0} does not have the correct 8-bit mode. Mode = {1}'.\
      format(image_stack_file,image_stack.mode))


  # Number of frames in the input image.
  n_frames = image_stack.n_frames 


  #percentage of frames where the particle should exist to be processed
  minimum_frames_per_particle = int(math.ceil(0.98 * n_frames)) 
  print "Minimum number of frames where the particle should exist is {0}".format(minimum_frames_per_particle)

  # Make sure the mask stack have n_frames.
  #This is called multi-image tiff

  if mask_stack.n_frames != n_frames:
    raise Exception("Mask stack {0} contains {1} frames while Image stack {2} contains {3} frames."\
      .format(mask_stack_file, mask_stack.n_frames, image_stack_file, n_frames))

  #Make sure the image resolution is passed
  if resolution == None:
    if not "resolution" in  image_stack.info:
      raise Exception("Can not retrive image resolution from {0}".format(image_stack_file))
    else:
      x_voxel_size = 1 / image_stack.info['resolution'][0]
      y_voxel_size = 1 / image_stack.info['resolution'][1]
  else: 
    x_voxel_size = float(resolution[0])
    y_voxel_size = float(resolution[1])

  print "Image resolution in X:{0}, Y:{1}".format(x_voxel_size, y_voxel_size)
  # Create the output directory
  if not os.path.exists(output_dir):
    os.mkdir(output_dir) 

  #Parse the RNA file and rename its columns.
  if segmented_rna_points != None:
    RNAs = pandas.read_csv(segmented_rna_points)
    RNAs.rename(columns={'frame':'rna_frame', 'x':'org_x', 'y':'org_y', 'intensity':'rna_intensity'}, inplace = True)
    RNAs.rna_frame = RNAs.rna_frame.apply(int)


  #Track the particles
  tracks_info = pandas.DataFrame.from_csv(particles_analysis)
  tracks_info.rename(columns={'X': 'x', 'Y': 'y', 'Slice' : 'frame'}, inplace=True)

  #Make the search range proportional to the resolution 
  physical_search_range = 6 * math.sqrt(x_voxel_size ** 2 + y_voxel_size ** 2)
  tracks = trackpy.link_df(tracks_info, search_range =  physical_search_range, memory = 3)
  particles_ids = tracks.particle.unique()
  
 
  #Pick the best tracked cells.
  sizes=[]
  good_ids = []
  for x in np.nditer(particles_ids):
      n_frames_for_particle = tracks[tracks.particle == int(x)].shape[0]
      sizes.append(n_frames_for_particle)
      if n_frames_for_particle > minimum_frames_per_particle:
          good_ids.append(int(x))
      
  #sizes.sort()
  #print sizes
  print "Found {0} cells to process".format(len(good_ids))

  good_ids = [int(i) for i in good_ids] 

  #Loop over all the cells
  iteration = 1
  n_cells = len(good_ids)
  for cell_id in good_ids:

    print "processing cell {0} out of {1}".format(iteration, n_cells) 
    #Create the cell directory 
    cell_dir = os.path.join(output_dir, "cell-{0}".format(cell_id))
    if not os.path.exists(cell_dir):
      os.mkdir(cell_dir)

    # To get the bounding box will be in image spacing between (BX, BY) and (BY + Width, BY + Height) 
    # Get this information for particle 1
    particle_info = tracks[tracks.particle == cell_id].loc[:,['x','y', 'frame', 'particle', 'BX', 'BY', 'Width', 'Height']]
    
    #LX and LY are the right lower borders of the segmented cell. 
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

    #Create a csv file to put the bounding box corresponding to every frame from the original image. 
    frames_info_file = os.path.join(cell_dir,"frames-info.csv")
    particle_info.to_csv(frames_info_file, index = False)

    new_max_width = int(math.ceil(particle_info.Width.max() / x_voxel_size) * extention)
    new_max_height = int(math.ceil(particle_info.Height.max() / y_voxel_size) * extention)

    #Initialize the per cell RNA Dataframe
    if RNAs is not  None:
      new_x = 'x_new'
      new_y = 'y_new'
      rna_spots_info = pandas.DataFrame(columns=[new_x, new_y, 'rna_frame', 'org_x', 'org_y', 'rna_intensity'])  

    #Loop over all the frames
    for frame_id in range(1,n_frames + 1):


      #Get the frame if it is part of the track, otherwise get the information
      # from the closest frame. 
      frame_info = particle_info[particle_info.frame == frame_id].iloc[0]
      if len(frame_info.index) == 0:

        closest_row_id = (particle_info.frame - frame_id).apply(abs).argmin()  
        frame_info = particle_info.loc[closest_row_id]

        print("Warning: cell:{0} frame_id:{1} not found, using frame_id:{2} instead."\
          .format(cell_id, frame_id, int(frame_info.frame)))

      try: 
        bounding_box = (frame_info.Left, frame_info.Upper, \
          frame_info.Right, frame_info.Lower)
      except:
        print "Can not calculate the bounding box for frame_id=" , frame_id
        raise

      #Get the RNA points that fall within the bounding box for this frame    
      if RNAs is not None:

        #Get all the RNA points that satisfy: frame = frame_id, and  frame_info.Left < X < frame_info.Right, andframe_info.Upper < Y < frame_info. Lower 
        frame_cell_points = RNAs[(RNAs.rna_frame == frame_id) &\
          (RNAs.org_x > frame_info.Left) & ( RNAs.org_x < frame_info.Right) &\
          (RNAs.org_y > frame_info.Upper) & (RNAs.org_y < frame_info.Lower)] 

        if len(frame_cell_points.index) != 0:
          #Calculate the new index for the RNA transcription points using the formulat x_new = (x_old - frame_info.x) / x_voxel_size  + new_max_width / 2.0 
          #Note that everything is in image dimention at this point
          frame_cell_points[new_x] = (frame_cell_points.org_x - frame_info.x / x_voxel_size + new_max_width / 2.0)
          frame_cell_points[new_y] = (frame_cell_points.org_y - frame_info.y / y_voxel_size + new_max_height / 2.0)

          rna_spots_info = pandas.concat([rna_spots_info, frame_cell_points])   

      # Move the image object to the frame_id.
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
   
    #Write the csv file for the RNAs spots 
    if RNAs is not None: 
      rna_info_file = os.path.join(cell_dir,"rna-info.csv")
      rna_spots_info.to_csv(rna_info_file, index = False, columns = ['rna_frame', new_x, new_y, 'org_x', 'org_y', 'rna_intensity'] )
      print "Found {0} RNA spots".format(len(rna_spots_info.index))
  
    exit()    
    iteration = iteration + 1

if __name__ == '__main__':
  
  parser = argparse.ArgumentParser(description="A script to track the cells and generate a directory that includes the crropped cells")
  parser.add_argument('measurments_file', help='A csv file that contains the measurments about the segmented cell using imagej format')
  parser.add_argument('image_stack', help='The 8 bit image stack')
  parser.add_argument('binary_stack', help="The binary stack to crop the region of interest for every cell.")
  parser.add_argument('--output_dir', dest='output_dir', default="./segmented-stacks", help="Directory to put the cells")
  parser.add_argument('--resolution', dest='resolution', help="Tuple representing the resolution of the stack in X and Y. Should use the format (x,y)")
  parser.add_argument('--rna_segmentation', dest='rnas', help="A csv file containing the x,y,frame,intensity for the RNA spots")


  args = parser.parse_args()

  stack_resolution = None
  if args.resolution != None:
    try:
      stack_resolution = literal_eval(args.resolution)
    except Exception as e: 
      print "Can not parse the resolution {0} as tuple. Error{1}".format(args.resolution, e)
      raise
  process_segmented_stack(args.measurments_file, args.output_dir, args.image_stack,args.binary_stack, resolution = stack_resolution, segmented_rna_points = args.rnas)


