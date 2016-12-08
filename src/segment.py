from pandas import DataFrame as df
import pandas
import numpy as np
import math
import trackpy
import os
import argparse
from ast import literal_eval
import pdb
import skimage.io
import skimage.external.tifffile
import glob
import re
from skimage import measure


def check_file(filename, file_list = None): 
  """
    Checks if a file exists and raise and exception with the appropriate error
    Argument:
      filename: the file to be checked.
      file_list: Check if the filename is part of the list
  """

  file_exists= False
  if file_list != None:
    if filename in file_list:
      file_exists= True 
  else:  
    if os.path.exists(filename):
      file_exists= True 

  if not file_exists:
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


def process_segmented_stack(cells_analysis, output_dir, image_stack_in, segmented_stack_in, resolution, segmented_rna_points = None, search_distance = None, maximum_memory = None, track_length = None):
  """
    Track the segmented cells. Generate a floder for every tracked cell where
    the masked and cropped cell frames will reside. Only cells that appear
    in 98% frames will be considered. For these cells, the missing frames
    be filled with the segmentation information of the nearest identified
    frame. 


    Arguments:
      cells_analysis: A csv file generated from imageJ containing metadata 
        about the segmented paticles. Columns must include the series 'X', 'Y',
        "Slice", "BX", "BY", "Width", and "Height"

      output_dir: The directory that will contain the cropped segmented cell.
        Every cell will have its own directory with the name /output_dir/cell_<id>

      image_stack_in: a prefix to the image files 
      
      segmented_stack_in: a prefix to to the segmentation files. All values greater than one are foreground.

      resolution: A tuple specifiying the image resolution in X and Y dimensions 
        If not passed the resolution will be infered from the image_stack metadata.

      segmented_rna_points: A csv file that contains the information about the
        segmented RNA transcription points. The format should be x, y, intensity,
        frame.  This script will collect the points that falls within a given
        cell bounding box (per frame) and saves a csv file with these points in in
        the corresponding folder of the cell. The points will saved in image
        indexes.

  """

  #some constants

  # The extension ratio for the cropped cells 
  extention = 1.10
  if maximum_memory == None:
    maximum_memory = 3
  
  if search_distance == None:
    search_distance = 25

  #Verify the inputs: 
  check_file(cells_analysis)

  #Get information about the segmented cells
  tracks_info = pandas.DataFrame.from_csv(cells_analysis)
  tracks_info.rename(columns={'X': 'x', 'Y': 'y', 'Slice' : 'frame'}, inplace=True)

  minimum_frame = int(tracks_info.frame.min()) 
  maximum_frame = int(tracks_info.frame.max())

  image_files = glob.glob(image_stack_in + "*")
  image_files.sort()
  first_file = image_files[0]
  last_file = image_files[-1]

  
  m_first = re.match(r"{0}(?P<id_extension>.*)".format(image_stack_in), first_file)
  first_id, extension = os.path.splitext(m_first.group('id_extension'))


  m_last = re.match(r"{0}(?P<id_extension>.*)".format(image_stack_in), last_file)
  last_id, extension = os.path.splitext(m_last.group('id_extension'))

  if len(first_id) != len(last_id):
    raise Exception("The number of digits representing the id should match. Got {0} and {1}".format(first_file, last_file))

  ndigits = len(first_id) 

  first_frame = int(first_id)
  last_frame = int(last_id)
  n_frames = last_frame - first_frame + 1


  #Make sure all the image files and the segmentation files between first_id and last_id exists

  for frame_id in range (first_frame, last_frame + 1):
    image_file = image_stack_in + str(frame_id).zfill(ndigits) + extension
    mask_file = segmented_stack_in + str(frame_id).zfill(ndigits) + extension
    check_file(image_file, file_list = image_files)
    try:
      check_file(mask_file)
    except Exception as e:
      raise Exception("The image file {0} does not have a matching mask file {1}. Exception:{2} ".format(image_file, mask_file, e))


  if minimum_frame < first_frame:  
    raise Exception("The analysis file contains a frame {0} that has an ID less than the first input id:{1}".format(minimum_frame, first_frame))

  if maximum_frame > first_frame + n_frames - 1:  
    raise Exception("The analysis file contains a frame {0} that has an ID greater than the last input id:{1}".format(maximum_frame, first_frame+n_frames))

  if segmented_rna_points != None:
    check_file(segmented_rna_points)

  #percentage of frames where the tracked cell should exist to be processed


  if track_length == None:
    in_percentage =   0.98
    minimum_frames_per_cell = int(math.ceil(in_percentage * n_frames)) 
  else: 
    if not isinstance(track_length, int):
      raise Exception("The minimum track length should be a positive interger. Got{0}".format(track_length)) 
    elif track_length < 1:
      raise Exception("The minimum track length should be a positive interger. Got{0}".format(track_length)) 
    else:
      minimum_frames_per_cell = track_length

  print "Minimum number of frames where the cell should exist is {0}".format(minimum_frames_per_cell)
  if minimum_frames_per_cell >  n_frames:
    raise Exception ("ERROR: The minimum number of frames:{0} is greater than the number of frames:{1}."\
      .format(minimum_frames_per_cell, n_frames))

  #Make sure the image resolution is passed
  x_voxel_size = float(resolution[0])
  y_voxel_size = float(resolution[1])

  print "Image resolution is X:{0}, Y:{1}".format(x_voxel_size, y_voxel_size)
  # Create the output directory
  if not os.path.exists(output_dir):
    os.mkdir(output_dir) 

  #Parse the RNA file and rename its columns.
  RNAs = None
  if segmented_rna_points != None:
    RNAs = pandas.read_csv(segmented_rna_points)
    RNAs.rename(columns={'frame':'rna_frame', 'x':'org_x', 'y':'org_y', 'intensity':'rna_intensity'}, inplace = True)
    RNAs.rna_frame = RNAs.rna_frame.apply(int)


  #Make the search range proportional to the resolution 
  physical_search_range = search_distance * math.sqrt(x_voxel_size ** 2 + y_voxel_size ** 2)
  tracks = trackpy.link_df(tracks_info, search_range =  physical_search_range, memory = maximum_memory)
  cell_ids = tracks.particle.unique()
 
  #Pick the best tracked cells.
  sizes=[]
  good_ids = []
  for x in np.nditer(cell_ids):
      n_frames_for_cell = tracks[tracks.particle == int(x)].shape[0]
      sizes.append(n_frames_for_cell)
      if n_frames_for_cell > minimum_frames_per_cell:
          good_ids.append(int(x))
      
  print "Found {0} cells to process".format(len(good_ids))

  good_ids = [int(i) for i in good_ids] 

  #Loop over all the cells
  iteration = 1
  n_cells = len(good_ids)



  for cell_id in good_ids:

    print "processing cell {0} out of {1}, cell-id:{2}".format(iteration, n_cells, cell_id) 
    #Create the cell directory 
    cell_dir = os.path.join(output_dir, "cell-{0}".format(cell_id))
    if not os.path.exists(cell_dir):
      os.mkdir(cell_dir)
    else:    
      iteration += 1
      continue

    # To get the bounding box will be in image spacing between (BX, BY) and (BY + Width, BY + Height) 
    # Get this information for cell 1
    cell_info = tracks[tracks.particle == cell_id].loc[:,['x','y', 'frame', 'particle', 'BX', 'BY', 'Width', 'Height']]
    
    #LX and LY are the right lower borders of the segmented cell. 
    cell_info['LX'] = cell_info.BX + cell_info.Width
    cell_info['LY'] = cell_info.BY + cell_info.Height

    # Convert the physical indexes to pixel coordinates for cropping
    cell_info['Left'] = cell_info.BX  / x_voxel_size
    cell_info.Left = cell_info.Left.apply(math.floor).apply(int)
    cell_info['Upper'] = cell_info.BY  / y_voxel_size
    cell_info.Upper = cell_info.Upper.apply(math.floor).apply(int)
    cell_info['Right'] = cell_info.LX  / x_voxel_size
    cell_info.Right = cell_info.Right.apply(math.ceil).apply(int)
    cell_info['Lower'] = cell_info.LY  / y_voxel_size
    cell_info.Lower = cell_info.Lower.apply(math.ceil).apply(int)

    #Create a csv file to put the bounding box corresponding to every frame from the original image. 
    frames_info_file = os.path.join(cell_dir,"frames-info.csv")
    cell_info.to_csv(frames_info_file, index = False)

    new_max_width = int(math.ceil(cell_info.Width.max() / x_voxel_size) * extention)
    new_max_height = int(math.ceil(cell_info.Height.max() / y_voxel_size) * extention)

    #Initialize the per cell RNA Dataframe
    if RNAs is not  None:
      new_x = 'x_new'
      new_y = 'y_new'
      rna_spots_info = pandas.DataFrame(columns=[new_x, new_y, 'rna_frame', 'org_x', 'org_y', 'rna_intensity'])  
      #print "Before analyzing the cell, number of spots info = {0}".format(len(rna_spots_info.index))

    #Loop over all the frames
    #The first frame in the analysis file start with 1. However, this might not be the case for  
    # images


    for frame_id in range(first_frame,n_frames + first_frame):

      #Get the frame if it is part of the track, otherwise get the information
      # from the closest frame. 

      frame_info_df =  cell_info[cell_info.frame == frame_id]
      if len(frame_info_df.index) > 0:
        frame_info = cell_info[cell_info.frame == frame_id].iloc[0]
      else: 
        closest_row_id = (cell_info.frame - frame_id).apply(abs).argmin()  
        frame_info = cell_info.loc[closest_row_id]

        print("Warning: cell:{0} frame_id:{1} not found, using frame_id:{2} instead."\
          .format(cell_id, frame_id, int(frame_info.frame)))

      try: 
        left =  int(frame_info.Left)
        right = int(frame_info.Right)
        upper = int(frame_info.Upper)
        lower = int(frame_info.Lower)
        bounding_box = (int(frame_info.Left), int(frame_info.Upper), \
          int(frame_info.Right), int(frame_info.Lower))
      except:
        print "Can not calculate the bounding box for frame_id=" , frame_id
        raise
      
      #Get the RNA points that fall within the bounding box for this frame    
      if RNAs is not None:

        #Get all the RNA points that satisfy: frame = frame_id, and
        #frame_info.Left < X < frame_info.Right,and
        #frame_info.Upper < Y < frame_info.Lower 

        frame_filter = RNAs[(RNAs.rna_frame == frame_id)]
        x_filter = frame_filter[(frame_filter.org_x > left) & ( frame_filter.org_x < right)]
        frame_cell_points = x_filter[(x_filter.org_y > upper) & (x_filter.org_y < lower)] 

        if len(frame_cell_points.index) != 0:
          #Calculate the new index for the RNA transcription points using the formulat x_new = (x_old - frame_info.x) / x_voxel_size  + new_max_width / 2.0 
          #Note that everything is in image dimention at this point

          x_df = (frame_cell_points.loc[:,'org_x'] - frame_info.x / x_voxel_size + new_max_width / 2.0).apply(int)
          frame_cell_points.loc[:,new_x] = x_df 
          frame_cell_points.loc[:,new_y] = \
            (frame_cell_points.loc[:,'org_y'] - frame_info.y / y_voxel_size + new_max_height / 2.0).apply(int)

          rna_spots_info = pandas.concat([rna_spots_info, frame_cell_points])   

      if 1:
        # Move the image object to the frame_id.
        image_array = skimage.io.imread(image_stack_in + str(frame_id).zfill(ndigits) + extension)

  
        #Move the mask to the closest frame id (which can be the frame_id)
        mask_array = skimage.io.imread(segmented_stack_in + str(int(frame_info.frame)).zfill(ndigits) + extension)
  
        # Crop the bounding box from the image and the mask 
        cropped_image = image_array[upper:lower,left:right]
        cropped_mask = mask_array[upper:lower,left:right]

        #Make sure cropped_mask includes one area only
        all_labels, num= measure.label(cropped_mask, return_num = True)

        #print "Warning, found:{0} regions in the mask corresponding to frame:{1}. Filtering out the small regions"\
        #  .format(num,frame_id)
        properties = measure.regionprops(all_labels) 
        details = [(p.area, p.label) for p in properties]
        largest_region_tuple = max(details, key = lambda t: t[0]) 
        largest_region_id = largest_region_tuple[1]
        #pdb.set_trace()
        cropped_image[all_labels != largest_region_id] = 0

  
        #Create the masked image
        #mask_ids = cropped_mask < 1
        #cropped_image[mask_ids] = 0

        #Extend the image
        expanded_image = np.zeros((new_max_height,new_max_width), dtype=cropped_image.dtype)

        #Paste the cropped image in the extended image]
        height_offset =  (new_max_height - cropped_image.shape[0]) / 2
        width_offset = (new_max_width - cropped_image.shape[1]) /2  
        expanded_image[height_offset:height_offset+cropped_image.shape[0],\
          width_offset:width_offset+cropped_image.shape[1]]  = cropped_image

        output_file_name = os.path.join(cell_dir,"fragment-{0}.tif".format(str(frame_id).zfill(ndigits)))
        skimage.external.tifffile.imsave(output_file_name, expanded_image)
   
    #Write the csv file for the RNAs spots 
    if RNAs is not None: 
      rna_info_file = os.path.join(cell_dir,"rna-info.csv")
      rna_spots_info.to_csv(rna_info_file, index = False, columns = ['rna_frame', new_x, new_y, 'org_x', 'org_y', 'rna_intensity'] )
      print "Found {0} RNA spots".format(len(rna_spots_info.index))
  
    iteration = iteration + 1
    #if iteration > 200:
    #  exit()

if __name__ == '__main__':
  
  parser = argparse.ArgumentParser(description="A script to track the cells and generate a directory that includes the crropped cells")
  parser.add_argument('measurments_file', help='A csv file that contains the measurments about the segmented cells using imageJ format. Note the Slice value should match the file name')
  parser.add_argument('image_stack', help='The image stack prefix')
  parser.add_argument('binary_stack', help="The binary stack prefix to crop the region of interest for every cell.")

  parser.add_argument('--output_dir', dest='output_dir', default="./segmented-stacks", help="Directory to put the cells")
  parser.add_argument('--resolution', dest='resolution', default = "(1,1)", help="Tuple representing the resolution of the stack in X and Y. Should use the format (x,y)")
  parser.add_argument('--rna_segmentation', dest='rnas', help="A csv file containing the x,y,frame,intensity for the RNA spots")
  parser.add_argument('--minimum_track_length', type = int, dest='track_length', help="The minimun number of frames a cell should exist.")
  parser.add_argument('--memory', type = int, dest='memory', help="The maximum number of frames a cell can disappear")
  parser.add_argument('--search_range', dest='search_range', type = int,  help="Number of pixels a cell can move between two appearances")


  args = parser.parse_args()

  stack_resolution = None
  if args.resolution != None:
    try:
      stack_resolution = literal_eval(args.resolution)
    except Exception as e: 
      print "Can not parse the resolution {0} as tuple. Error{1}".format(args.resolution, e)
      raise

  process_segmented_stack(args.measurments_file, args.output_dir,\
    args.image_stack, args.binary_stack, resolution = stack_resolution, \
    segmented_rna_points = args.rnas, search_distance = args.search_range, \
    maximum_memory = args.memory, track_length = args.track_length )


