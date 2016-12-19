import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import pandas
import pims
import trackpy as tp
import re
import shutil
import threading
import multiprocessing
import os
import csv 
import math
import glob
import skimage.io
from skimage import exposure

class Track:

  def __init__(self, output_dir, override_results=False):

    #work directory 
    self.output_directory = os.path.abspath(output_dir) 
    self.create_dir(self.output_directory, override=override_results)

    self.elastix_dir = self.output_directory + '/elastix/' 
    self.transformix_dir = self.output_directory + '/transformix/' 

    self.create_dir(self.elastix_dir, override=override_results)
    self.create_dir(self.transformix_dir, override=override_results)

    self.registration_results = os.path.join(self.elastix_dir, 'registration_results') 
    self.create_dir(self.registration_results)

    #Elastics convensions.
    self.elastix_result_prefix="transform" 
    #This can be retrived from the registration params file.
    self.registered_image_name = "result.0.tif"

    """ 
    frames_dictionary has the following structure:

                                                                                      |->result_tag(string):value #The group id provided by the tracking method.
                                                                                      |
                                                                                      |->original_group_id(key):value(int) #The ground truth group id.
                                                                                      |
                                          |-->'original_points'(dict)-->point_id(dict)-->index(key):index[x,y](value)
                                          | 
    frames_dictionary-->frame_id(dict)--->|-->'transformed_points'(dict)-->original_frame_id(dict)-->point_id(key):index[x,y](value)
                                          |
                                          |-->'n_points'(int): Number of original point in this frame_id. Used to generate successive points_id. 
    """
    self.frames_dictionary = {} 
    
    #Initialize the images information
    self.images_dir = None 
    self.frame_prefix = None
    self.frame_suffix = None


    #Transcription point can be either index of point
    #It should be set when points are added
    self.point_index = None 

    #A dataframe to store the segmented particles
    self.particles_df = None

#####PUBLIC CLASS FUNCTIONS

  def preprocess_track_file(self, track, index, original_group_id = None): 
    """Create a numpy array that contains the frame numer and location of the object in only the frames where the object was physically present

        Args: 
          track: file that indicates the spacial location of the object in each frame.
          index: file that indicates which objects are physically present in each frame.
          original_group_id: Add a group id for the points that will be added. This will be useful to evaluate the tracking results 


        Returns:
        [length(file),3] np array.
    """
    self.check_file(index)
    index_file = open(index)
    a = index_file.readlines()
    length = len(a)
    indices = []
    count = 0
    while count < length-1:
      if int(a[count]) == 1:
        indices.append(count)
      count+=1
    locations = np.zeros((length,3))
    location_array = self.find_points_from_trk_file(track)
    x = location_array[:,0].tolist()
    y = location_array[:,1].tolist()
    x.pop(0)
    y.pop(0)
    locations[:,1] = x
    locations[:,2] = y
    locations[:,0] = list(range(1,length+1))
  
    positions = locations[indices]

    #George's Edit 
    for row in positions:
      frame_id, x, y = row   
      frame_id = str(int(frame_id)).zfill(3)
      self.add_original_point_to_dictionary(frame_id, [x, y], original_group_id)

   
    #if self.particles_df is not None:
    #self.particles_df =  DataFrame(data = positions, columns = ['frame', 'x', 'y'])
    #self.particles_df.to_csv('particles-segmentation.csv')
    return positions


  def populate_rna_segmentation(self, segmentations, voxel_size, image_dimensions): 
    """Create a dataframe and populate it using the segmented spots in a csv file.
       The first row of csv file should contain the column tags.  
       The column tags should contain the 'frame', and  'x', 'y' in physical coordinates.
        

       Parameters: 
        segmentations: file that contains the segmented spots in the format mentiond above. 
        voxel_size: A tuple (x,y) containing the spacing between two pixels in a frame.
        image_dimensions: A tuple (x,y) containing the dimensions of a frame.
       
       Retruns:
        The read dataframe with an extra column called unique_id.
    """

    #Check that the input file can be read   
    frames_list = []

    if (not isinstance(voxel_size[0], float)) and (not isinstance(voxel_size[0], int)):
      if (not isinstance(voxel_size[1], float)) and (not isinstance(voxel_size[1], int)):
        raise Exception("The voxel size should be int or float")

    if voxel_size[0] <= 0 or voxel_size[1] <= 0:
      raise Exception ("The voxel size should be a positive number")


    if (not isinstance(image_dimensions[0], int)) or (not isinstance(image_dimensions[1], int)):
      raise Exception ("The image dimensions should be integers")

    if image_dimensions[0] <= 0 or image_dimensions[1] <= 0:
      raise Exception ("The image dimensions should be positive integers")

   
    self.particles_df = pandas.read_csv(segmentations)
    
    #Convert thunderstorm column names: 
    if 'x [nm]' in self.particles_df.columns.tolist():
      self.particles_df.rename(columns={'x [nm]': 'x'}  , inplace =True)

    if 'y [nm]' in self.particles_df.columns.tolist():
      self.particles_df.rename(columns={'y [nm]': 'y'}  , inplace =True)

    #Make sure the required columns are present:
    if not 'x' in self.particles_df.columns.tolist():
      raise Exception ("The csv file {0} does not contain the column 'x'".format(segmentations))

    if not 'y' in self.particles_df.columns.tolist():
      raise Exception ("The csv file {0} does not contain the column 'y'".format(segmentations))

    if not 'frame' in self.particles_df.columns.tolist():
      raise Exception ("The csv file {0} does not contain the column 'frame'".format(segmentations))

    #Convert the frames to int:
    self.particles_df.frame = self.particles_df.frame.astype(int)
    
    #Convert x physical coordinates to index
    self.particles_df['x_index'] = self.particles_df.apply(lambda row: \
      self.convert_coordinates(row['x'], voxel_size[0],image_dimensions[0]), axis = 1) 

    #Convert y physical coordinates to index
    self.particles_df['y_index'] = self.particles_df.apply(lambda row: \
      self.convert_coordinates(row['y'], voxel_size[1],image_dimensions[1]), axis = 1) 
  
    
    self.particles_df['unique_id'] = self.particles_df.apply(lambda row: self.add_original_point_to_dictionary(int(row['frame']),[row['x_index'], row['y_index']]), axis=1) 

    return self.particles_df['frame'].tolist()


  def convert_coordinates(self, physical_value, voxel_size, image_size):
    """
      Converts the physical coordinates to index values within the image and verifies the converted result lies within the image boundries  
      Arguments:
        physical_value: The coordinate in physical domain  
        voxel_size: The voxel size in that dimension
        image_size: The image size in that dimension

      Returns:
        index_value: The index value 
    """

    index_value = physical_value / voxel_size
    if index_value > image_size:
      raise Exception("ERROR: The physical coordinate:{0} is converted to index: {1} using voxel size:{2} and it can no exceed the image size {3}.".format(physical_value, index_value, voxel_size, image_size))

    return index_value


  def preprocess_rna_segmentation(self, segmentations, voxel_size, image_dimensions): 
    """Create a numpy array that contains the frame number and location of the RNA spot from a csv file.
       The first column of the CSV file should be the frame number.
       The second column is the X position in physical coordinates.  
       The third column should be the Y position in physical coordinates.
        
       All non numeric rows will be ignored.

       Parameters: 
        track: file that contains the segmented spots in the format mentiond above. 
        voxel_size: A tuple containing the spacing between two pixels in the image.
        image_dimensions: A tuple containing the dimensions of the frames.
       
       Retruns:
        A list of frame_ids where the spots are segmented
    """

    #Check that the input file can be read   
    frames_list = []

    if (not isinstance(voxel_size[0], float)) and (not isinstance(voxel_size[0], int)):
      if (not isinstance(voxel_size[1], float)) and (not isinstance(voxel_size[1], int)):
        raise Exception("The voxel size should be int or float")

    if voxel_size[0] <= 0 or voxel_size[1] <= 0:
      raise Exception ("The voxel size should be a positive number")


    if (not isinstance(image_dimensions[0], int)) or (not isinstance(image_dimensions[1], int)):
      raise Exception ("The image dimensions should be integers")

    if image_dimensions[0] <= 0 or image_dimensions[1] <= 0:
      raise Exception ("The image dimensions should be positive integers")

    csv_file = open(segmentations, 'r')
    segmented_spots = csv.reader(csv_file)  
    

    row_id = 1 
    for row in segmented_spots:
      if len(row) != 0:      
        
        if str.isdigit(row[0].replace('.',"")):
          frame_id = int(float(row[0]))    
          x = int(float(row[1])/ voxel_size[0])
          y = int(float(row[2])/ voxel_size[1])
  
          if x > image_dimensions[0]:
            raise Exception("Error in {0} row {1}, the physical coordinate x={2} can not exceed the image size {3}.".format(segmentations, row_id, x, image_dimensions[0]))
  
          if y > image_dimensions[1]:
            raise Exception("Error in {0} row {1}, the physical coordinate y={2} can not exceed the image size {3}.".format(segmentations, row_id, y, image_dimensions[1]))
          
          frames_list.append([frame_id, x, y])
          self.add_original_point_to_dictionary(frame_id, [x, y])
        else: 
          print "Warning: ignoring row:{0}".format(row)
      else:
        print "Warning: ignoring row:{0}".format(row)
      row_id = row_id + 1

    return frames_list

  def transform_points(self, float_frame):
    """Transform the points using transfromix 
      Parameters:
        float_frames: A dictionary where the keys are the reference frames where points resides,
        and the values are the float frame to which the every point in the reference frame will be transformed.

    """

    #Create the directory where the transformix input files reside  
    input_files_dir = os.path.join(self.transformix_dir, "transformix_inputs") 
    self.create_dir(input_files_dir)

    for ref_id in  float_frame:
      #Get the registration result
      float_id = float_frame[ref_id]
      registration_result = self.get_registration_results(ref_id, float_id)
      if not os.path.isfile(registration_result):
        raise Exception ("The registration result file {0} does not exists".format(registration_result))

      original_points = self.get_original_points(ref_id)
      transformed_points = self.get_transformed_points(ref_id)
      points_to_be_transformed = original_points + transformed_points
   
      if points_to_be_transformed != None:
        n_points = len(points_to_be_transformed)
        points_indexes_list = [ point_info[2] for point_info in points_to_be_transformed ]

        #Create the input file
        frame_input_file = os.path.join(input_files_dir,"input_points" + ref_id + "-" + float_frame[ref_id] + ".txt")
        self.points_to_transformix(points_indexes_list, frame_input_file)
  
        #Create the transformix directory
        transform_directory = self.get_transform_dir(ref_id, float_id) 
        self.create_dir(transform_directory, override=True)
  
        #Call transformix
        output_points_file = self.transform(frame_input_file, registration_result, transform_directory)
  
  
        #Parse the transformix results
        output_points = self.find_points_from_transformix_output(output_points_file)
        n_transformed_points = len(output_points)
        
        if n_transformed_points != n_points:
          raise Exception ("ERROR: The number of transformed points {0} does not equal the number of input points {1}.\n Output points file:{2}"\
            .format(n_transformed_points, n_points, output_points_file))
 

        #Write the results back in the dictionary 

        for original_point, transformed_indexes in zip(points_to_be_transformed,  output_points):
          original_frame_id = original_point[0]
          point_id = original_point[1]
          self.add_transformed_point_to_dictionary(float_id, original_frame_id, point_id, transformed_indexes)

        #for original_frame_id, point_id, transformed_indexes in points_to_be_transformed:#zip(original_frame_id_list, point_id_list, output_points):
        #  self.add_transformed_point_to_dictionary(float_id, original_frame_id, point_id, transformed_indexes)

      else:
        print "WARNING: there is no points to be transformed for frame_id {0}.".format(ref_id)

  def points_to_transformix(self, points,  destination_file): 
    """Create a transformix input file for points in the transform directory
      Parameters:  
        points: a list of points  
        destination_file: The path to the file where the transformix input files will reside

    """

    #self.create_dir(transform_directory)
    #input_points_file_name = os.path.join(transform_directory, "input_points.txt")

    if self.point_index == None:
      raise Exception ("ERROR: point_index should be set before calling transformix")

    input_points_file = open(destination_file, 'w')
    input_points_file.write(self.point_index + '\n')
    input_points_file.write(str(len(points)) + '\n')

    for point in points:
      input_points_file.write(str(point[0]) + ' ' + str(point[1]) + '\n')

    input_points_file.close()

  def transform(self, input_points_file, registration_result, transform_directory=None):
    """Transform the input_points  using transformix
      Parameters:
        input_points_file: File containing the input points
        registration_result: The file containing the registraton result to use in the transformation 
        transform_directory: The directory where the transformation will take place 


      Returns:
        trasnform_result: Path to the transformix output file
    """

    #Make sure ELASTIX_PATH is exported
    try:
      os.environ['ELASTIX_PATH'] 
    except KeyError as error:
      raise Exception('ERROR: The environment variable "ELASTIX_PATH" is not exported')


    #Make sure GENE_TRAP_ROOT  is exported
    try:
      os.environ['GENE_TRAP_ROOT'] 
    except KeyError as error:
      raise Exception ('ERROR: The environment variable "GENE_TRAP_ROOT" is not exported')

    transformix_script = os.environ['GENE_TRAP_ROOT'] + '/bin/transformix.sh'

    #Verify the input arguments
    self.check_file(input_points_file)
    self.check_file(registration_result)
        
    if transform_directory != None:  
      self.create_dir(transform_directory, override=True)
      output_dir = transform_directory
    else:
      output_dir = self.transformix_dir
    output_dir = os.path.abspath(output_dir)



    input_points = ' -def ' + os.path.abspath(input_points_file)
    transform_dir = ' -out ' +  output_dir
    transform_parameters = ' -tp ' +    os.path.abspath(registration_result)
    output_arguments = ' "' + output_dir + '" '
    transformix_call =  transformix_script + ' "' + input_points + transform_dir + transform_parameters +  ' " ' + output_arguments

    print transformix_call
    status = os.system(transformix_call)
   
    if status != 0:
      raise Exception ('ERROR: The transformix call "{0}" returned with code {1}'.format(transformix_call, status))
    else: 
      transformed_points =  os.path.join(output_dir, "outputpoints.txt" ) 
         
    if not os.path.isfile(transformed_points):
      transformed_points = None

    return transformed_points


  def collect_registered_images_to_single_dir(self, reference_frame_id, moving_frames_ids, destination_dir=None):
    """Group the transformed image frames into the destination directory
      Parameters:
        reference_frame: The fixed image id to which the images were registered. 
        moving_frames_ids: List of the ids for the moving frames.
        destination_dir: The directory where the results should be accumulated.

      Returns:
        destination_directory: The folder where the results were grouped.
    """
    
    if destination_dir != None:
      destination_name  = destination_dir     
    else:
      destination_name  = os.path.join(self.elastix_dir, str(reference_frame_id) + "-results/") 
    destination_name = os.path.abspath(destination_name)
    self.create_dir(destination_name, override = True)
    reference_list = [str(reference_frame_id)] * len (moving_frames_ids)
    
    results_files = [os.path.join(self.get_registration_dir(reference_id,moving_id), self.registered_image_name)  for reference_id, moving_id in zip(reference_list, moving_frames_ids)]
    for file_name, file_id in zip(results_files, moving_frames_ids) :
      if not os.path.exists(file_name):
        print 'WARNING: The file "{0}" does not exist.'.format(file_name)
      else:
        shutil.copy(file_name, os.path.join(destination_name, "frame" + str(file_id)))
    return destination_name

  def set_image_information(self, images_dir, frame_prefix, frame_suffix):
    """Sets the information about the cell frames before registration  
       Every frame is located at: /path/to/images_dir/<image_prefix><frame_id><image_suffix>
  
       Parameters:
          images_dir: The directory that contains the images files 
          frame_prefix: The prefix of the frame file 
          imagee_suffix: The suffix for the frame file 
    """
      
    if images_dir == None:
      raise Exception ("ERROR: Null value passed to images_dir")

    if not os.path.isdir(images_dir): 
      raise Exception ("The direcotry {0} does not exist".format(images_dir)) 
    else:
      self.images_dir = os.path.abspath(images_dir)

 
    if frame_prefix == None or not isinstance(frame_prefix, basestring):
      raise Exception ("ERROR: frame_prefix should be a null string")
    self.frame_prefix = frame_prefix 


    if frame_suffix == None or not isinstance(frame_suffix, basestring):
      raise Exception ("ERROR: frame_suffix should be a null string")
    self.frame_suffix = frame_suffix 


  def register_series_to_single_frame(self, reference_frame, moving_frames, registration_params, reverse=False):
    """Register the list of moving_frames to the reference_frame. 
        The images file name follows the format <frame_prefix><frame_id><frame_suffix>
      Parameters:
        reference_frame: The ID of the reference image
        moving_frames: The list of IDs of the float images
        registration_params: Path to the elastix parameters file
        reverse: Switch the registration direction
    """
 

    moving_list = [str(frame) for frame in moving_frames]
    reference_list = [str(reference_frame)] * len(moving_frames)

    if reverse:
      temp = moving_list
      moving_list = reference_list
      reference_list = temp

    self.register_frames(reference_list, moving_list, registration_params)

  def register_frames(self, reference_list, moving_list, registration_params):
    """Register the moving_list to the reference_list. 
      Parameters:
        reference_list: List containing the reference's frame_id(s). 
        moving_list: List containing the moving frame_id(s).
        registration_parameters: The registration parameters to be passed to elastix
 
    """

    if self.images_dir == None: 
      raise Exception ("ERROR:the image_dir is not set ")

    if self.frame_prefix == None: 
      raise Exception ("ERROR: the frame_prefix is not set ")

    if self.frame_suffix == None:
      raise Exception ("ERROR: the frame_suffix is not set ")

    if reference_list != None  and moving_list != None:
      if len(reference_list) != len(moving_list):
        raise Exception ("ERROR: The length of the reference and moving list should match")
   
    image_prefix = self.frame_prefix 
    image_suffix = self.frame_suffix

    fixed_images = [os.path.join(self.images_dir,image_prefix + str(frame_id) + image_suffix) for frame_id in reference_list]
    moving_images = [os.path.join(self.images_dir,image_prefix + str(frame_id) + image_suffix) for frame_id in moving_list]
    params = [registration_params] * len(moving_list)
    registration_directories = [self.get_registration_dir(reference_id,moving_id) for reference_id, moving_id in zip(reference_list, moving_list)]
    registration_results_filenames = [self.get_registration_results(fixed, moving) for fixed, moving in zip(reference_list, moving_list)]
    #transformed_images_names = [self.get_transformed_image(fixed, moving) for fixed, moving in zip (reference_list, moving_list)] 

    register_arguments = zip(fixed_images, moving_images, params, registration_directories,registration_results_filenames)

  
    #Check if registration is already done  
    work_list = []
    for fixed,moving, arguments in zip(reference_list, moving_list, register_arguments):
      if not os.path.isfile(self.get_registration_results(fixed, moving)):
        work_list.append(arguments)

    print "Registering {0} pairs.".format(len(work_list))
    from multiprocessing.pool import ThreadPool
    nthreads = multiprocessing.cpu_count()
    #TODO This needs to be implemented using processess instead of threads. Setting nthreads to 1 for now.
    nthreads = 1
    print "nthreads = " + str(nthreads)

    pool = ThreadPool(processes=nthreads)
    async_results = []
    for registration in range(len(work_list)):
      async_result = pool.apply_async(self.register, work_list[registration])
      async_results.append(async_result)

    for registration in range(len(work_list)):
      transformed_image, registration_solution = async_results[registration].get()
      if registration_solution == None: 
        raise Exception ("ERROR: The {0} file does not exist".format(registration_solution))


  def register(self, fixed_image, moving_image, registration_params, temp_directory=None, registration_result_name=None):
    """ Register two images using elastix
        Parameters:
          fixed_image: Path to the fixed image
          moving_image: Path to the moving image
          registration_params: Path to the elastix parameters file
          registration_result_name: The complete path to the file where the registration result should be copied

        Returns:
          [transformed_image, registration_results] 
          transformed_image: Path to the transfomed image. "None" if the file does not exist.
          registration_result: The results of elastix regsitration. "None" if the file does not exist. 
    """

    #Make sure ELASTIX_PATH is exported
    try:
      os.environ['ELASTIX_PATH'] 
    except KeyError as error:
      raise Exception('ERROR: The environment variable "ELASTIX_PATH" is not exported')


    #Make sure GENE_TRAP_ROOT  is exported
    try:
      os.environ['GENE_TRAP_ROOT'] 
    except KeyError as error:
      raise Exception ('ERROR: The environment variable "GENE_TRAP_ROOT" is not exported')

    elastix_script = os.environ['GENE_TRAP_ROOT'] + '/bin/elastix.sh'

    self.check_file(fixed_image)
    self.check_file(moving_image)
    self.check_file(registration_params)
        
    if temp_directory != None:  
      self.create_dir(temp_directory, override=True)
      output_dir = temp_directory
    else:
      output_dir = self.elastix_dir
      
    output_dir = os.path.abspath(output_dir) 

  
    moving = ' -m ' + moving_image
    fixed = ' -f ' + fixed_image 
    parameter = ' -p ' + registration_params 
    output = ' -out ' + output_dir 

    elastix_arguments = ' " ' + moving + fixed + parameter + output + '"'
    output_arguments = ' "' + output_dir + '" '

    elastix_call =  elastix_script + elastix_arguments + output_arguments
    print elastix_call

    status = os.system(elastix_call)
   
    if status != 0:
      raise Exception ('ERROR: The elastix call "{0}" returned with code {1}'.format(elastix_call, status))
    else: 
      transformed_image =  os.path.join(output_dir, "result.0.tif") 
      registration_result =  os.path.join(output_dir, "TransformParameters.0.txt")
         
    if not os.path.isfile(transformed_image):
      transformed_image = None
    if not os.path.isfile(registration_result):
      registration_result = None


    #Rename the results if the new name is provided
#    if transfomed_image != None and transformed_image_name != None:
#      destitation = os.path.join(output_dir, transformed_image_name)
#      shutil.copy(transformed_image, destination)
#      transformed_image = destination

    if registration_result != None and registration_result_name != None:
      #destination = os.path.join(output_dir, registration_result_name)
      shutil.copy(registration_result,registration_result_name)
      registration_result = registration_result_name 

       
    return [transformed_image, registration_result]

  def track(self, np_array, memory, distance):
    """Tracks the particles and populate the particles_df with the correct id.
        Parameters:
          np_array: generated from transformix_results_to_array
          memory: maximum number of frames an object can be absent before it is declared a separate object
          distance: search range between frames for the movement of one object
    """

    t = tp.link_df(np_array,search_range = distance, memory = memory)

    t['unique_id'] = t.apply(lambda row: self.get_unique_id(int(row['frame']),int(row['particle_id']) ), axis = 1)

    #Merge  partile id with the original particles dataframe
    if self.particles_df is not  None:
      subset_df = t[['unique_id', 'particle']]      
      merged_df = self.particles_df.merge(subset_df,  on='unique_id')
      merged_df.to_csv("tracked-particles.csv", index=False)
      print "Merged the dataframes"
      
    #tp.plot_traj(t)


###PRIVATE CLASS FUNCTIONS   

  def find_points_from_trk_file(self, file_name):
    """Finds input points from trk file

      Parameters:
       .trk file


      Returns:
        np array. array[:,0] is the x coordinate and array[:,1] is the y coordinate
    """

    self.check_file(file_name)
    points = open(file_name).readlines()
    count = 0
    x = []
    y = []

    for a in range(len(points)):
      p = re.findall("[0-9]*\.[0-9]*", points[a])
      x.append(float(p[0]))
      y.append(float(p[1]))
      count += 1
    array = np.zeros((len(points),2))
    array[:,0] = x
    array[:,1] = y
    return array


  def create_dir(self, dir_name, override=False):
    """Creates a directory  
      Prameters:
        dir_name: The directory path.
        override: Remove the directory if it exists before creating it.

    """

    if os.path.isdir(dir_name):
      if override:
        shutil.rmtree(dir_name)
        os.makedirs(dir_name)
    else:
        os.makedirs(dir_name)


  def check_file(self, file_path):
    """Checks the existance of a file
      Parameters:
        file_path: The path of the file
    """

    try:
      if not os.path.exists(file_path):
        raise Exception(file_path)  
    except Exception as inst:
     file_name = inst.args[0]
     raise Exception ('The input file:"{0}" does not exists'.format(file_name))

  def get_registration_results(self, fixed_id, float_id):
    """returns the file name that should have the registration results 
      Parameters:
        fixed_id: The id of the fixed image
        float_id: the id of the float image

    """

    file_path = os.path.join(self.registration_results,self.elastix_result_prefix + str(fixed_id) + '-' + str(float_id) + ".txt")
    return file_path

  def get_registration_dir(self, fixed_id, float_id): 
    """returns the folder name that should have the registration results 
      Parameters:
        fixed_id: The id of the fixed image
        float_id: the id of the float image
    """

    return os.path.join(self.elastix_dir,str(fixed_id) + "-" + str(float_id))



  def get_transform_dir(self, fixed_id, float_id): 
    """returns the folder name that should have the registration results 
      Parameters:
        fixed_id: The id of the fixed image
        float_id: the id of the float image
    """
    return os.path.join(self.transformix_dir,  str(fixed_id) + "-" + str(float_id))


  def set_point_index(self, point_index): 
    """Set the type of points for the object
      Parameters:
        point_index: A string that can be either "point" for physical points or "index" for image indexes
    """

    if point_index != "index" and point_index != "point":
      raise Exception ("ERROR: The point_index variable should be either 'point' or 'index'  ")
    else:
      self.point_index = point_index

  def find_points_from_transformix_output(self,file_name):
    """Finds output points in a transformix outputpoints.txt file
      Parameters:
        file_name: Transformix outputpoints.txt file
        point_index: Match the either the point or the index. Possible values are "index" or "point"  

      Returns:
        ouptut_points: A 2D list of points    size 2 list
    """
    
    self.check_file(file_name)
    if self.point_index == None:
      raise Exception ("ERROR: point_index should be set before calling transformix")

    point_index = self.point_index
    points_file = open(file_name)
    if point_index ==  "index":
      regex = re.compile(".*OutputIndexFixed\s*=\s*\[\s*(?P<x>[-+]?[0-9]*)\s*(?P<y>[-+]?[0-9]*)\s*\]")
    else:
      regex = re.compile(".*OutputPoint\s*=\s*\[\s*(?P<x>[-+]?\d+\.\d+)\s*(?P<y>[-+]?\d+\.\d+)\s*\]") 

    results = regex.findall(points_file.read()) 

    #Parse the results
    output_points = []
    for result in results: 
      if point_index ==  "index":
        output_points.append([int(result[0]) , int(result[1])])
      else:
        output_points.append([float(result[0]) ,  float(result[1])])
      
    points_file.close() 
    return output_points


  def add_transformed_point_to_dictionary(self, dest_frame_id, original_frame_id, point_id, indexes):
    """Adds a transformed point to the frames_dictionary.
      Prameters:
        dest_frame_id: The frame where the point will be added
        original_frame_id: The original frame where the point belongs. 
        point_id: The point id in the origninal frame.
        indexes: the indexes of the point in the destination frame id.
    """

    #Make sure indexes contains two number for X and Y.
    if len(indexes) != 2:
      print indexes
      raise Exception ("ERROR: the lenght of indexes does not equal 2.")

    if dest_frame_id == original_frame_id: 
      raise Exception ("ERROR: Destination and original frame_id are equivalent{0}.".format(dest_frame_id))
    
    if not dest_frame_id in self.frames_dictionary:
      self.frames_dictionary[dest_frame_id] = {}

    if not 'transformed_points' in self.frames_dictionary[dest_frame_id]: 
      self.frames_dictionary[dest_frame_id]['transformed_points']={}

    if not original_frame_id in self.frames_dictionary[dest_frame_id]['transformed_points']:
      self.frames_dictionary[dest_frame_id]['transformed_points'][original_frame_id]  = {}

    self.frames_dictionary[dest_frame_id]['transformed_points'][original_frame_id][point_id] = indexes  


  def add_original_point_to_dictionary(self, original_frame_id_raw, indexes, original_group_id = None):
    """Adds an original point to the frames_dictionary.
      Parameters:
        original_frame_id_raw: The original frame where the point belongs. 
        indexes: the indexes of the point in the destination frame id.
        original_group_id: The original group id that the point belongs too

      Returns:
        unique_id: A unique id for every original point in the dictionary. The first part will be the frame_id, 
          and the second part will be the point_id.
    """

    #Make sure indexes contains two number for X and Y.
    if len(indexes) != 2:
      print indexes
      raise Exception ("ERROR: the lenght of indexes does not equal 2.")
 

    #Get the integer and consistantly add the zero fill 
    original_frame_id = self.id2str(original_frame_id_raw)

    if not original_frame_id in self.frames_dictionary:
      self.frames_dictionary[original_frame_id] = {}

    if not 'original_points' in self.frames_dictionary[original_frame_id]:
      self.frames_dictionary[original_frame_id]['original_points']={}
      self.frames_dictionary[original_frame_id]['n_points']=0
  
    #Increament 
    self.frames_dictionary[original_frame_id]['n_points'] += 1
    point_id = self.frames_dictionary[original_frame_id]['n_points']

    if not point_id in self.frames_dictionary[original_frame_id]['original_points']:
      self.frames_dictionary[original_frame_id]['original_points'][point_id] = {}

    self.frames_dictionary[original_frame_id]['original_points'][point_id]['indexes'] = indexes  

    if original_group_id != None:
      self.frames_dictionary[original_frame_id]['original_points'][point_id]['original_group_id'] = original_group_id

    return self.get_unique_id(int(original_frame_id_raw), point_id) 
    

  def get_unique_id(self,frame_id, point_id):
    """
      Returns: a unique string id for every point in the dictionary. 
      Parameters: 
        frame_id: A positive  integer for the frame_id.
        point_id: A positive integer for the point id

      Returns: 
        a unique id for every point in the dictionary. 
    """

    #Verify the inputs
    if not isinstance(frame_id, int):
      raise Exception ("frame_id should be an integer. Got:{0}".format(frame_id))

    if not isinstance(point_id, int):
      raise Exception ("point_id should be an integer. Got:{0}".format(point_id))

    if frame_id < 0:
      raise Exception ("frame_id should be a positive integer. Got:{0}".format(frame_id))

    #Assume the maximum point id is 999
    if point_id < 0 or point_id > 999 :
      raise Exception ("point_id should be a positive integer between 0 and 999. Got:{0}".format(point_id))

    return  frame_id * 1000  + point_id


  def frame_points_to_dataframe(self, frame_id):
    """Collect the original and transformed points from a frame and returns a data frame
      Parameters:
        frame_id: The ID of the frame where the points will be collected 

      Returns:
        data_frame: The DataFrame object, None if there are points in the frame.
    """

    if frame_id == None: 
      raise Exception ("Invalid frame_id")

    if not frame_id in self.frames_dictionary:
      print "WARNING: the frame_id {0} is not part of the frames dictionary".format(frame_id)
      return None

    temp_np = np.empty((0,4))

    original_points = self.get_original_points(frame_id)
    transformed_points = self.get_transformed_points(frame_id)
    frame_points = original_points + transformed_points

    for frame_id, point_id, indexes in frame_points:
        temp_np = np.append(temp_np,[[int(frame_id),indexes[0], indexes[1], point_id]], axis=0)
  
    if temp_np.shape[0] == 0:
      return None
    else:    
      return DataFrame(data = temp_np, columns = ['frame', 'x', 'y', 'particle_id'])

    return temp_np


  def get_original_points(self, frame_id):
    """Generate a list of tuples in the form (frame_id, point_id, indexes) for the original points in a given frame_id
       Parameters:
          frame_id: The ID of the frame where the points will be colleceted

       Returns:
          original_points_list: List of tuples in the form (frame_id, point_id, indexes)

    """

    if frame_id == None: 
      raise Exception ("Invalid frame_id")

    original_points_list = []

    if not frame_id in self.frames_dictionary:
      #print "WARNING: the frame_id {0} is not part of the frames dictionary"
      return original_points_list

    if 'original_points' in  self.frames_dictionary[frame_id]:
      for point_id in self.frames_dictionary[frame_id]['original_points'].keys():
        indexes = self.frames_dictionary[frame_id]['original_points'][point_id]['indexes']
        original_points_list.append((frame_id, point_id, indexes))

    return original_points_list
          
  def get_transformed_points(self, frame_id):
    """Generate a list of tuples in the form (frame_id, point_id, indexes) for the transformed points in a given frame_id
       Parameters:
          frame_id: The ID of the frame where the points will be colleceted

       Returns:
          transformed_points_list: List of tuples in the form (frame_id, point_id, indexes)

    """

    if frame_id == None: 
      raise Exception ("Invalid frame_id")

    transformed_points_list = []

    if not frame_id in self.frames_dictionary:
      #print "WARNING: the frame_id {0} is not part of the frames dictionary"
      return transformed_points_list 

    #Add transformed points
    if 'transformed_points' in  self.frames_dictionary[frame_id]:
      for original_frame_id in self.frames_dictionary[frame_id]['transformed_points'].keys():
        for point_id , indexes in self.frames_dictionary[frame_id]['transformed_points'][original_frame_id].items():
          transformed_points_list.append((original_frame_id,point_id, indexes))
          
    return transformed_points_list 

  def add_track_result(self, frame_id, point_id, result_tag, result_id):
    """ Add to the dictionary the results information
        Parametes:
          frame_id: The original frame where the point belongs
          point_id: The original id of the point
          result_tag(string): The tag for the tracking operation
          result_id(int): The object_id the tracking operation gave to the point
    """

    if not frame_id in self.frames_dictionary:
      raise Exception("ERROR: Adding a result to non exisiting frame")

    if not 'original_points' in  self.frames_dictionary[frame_id]:
      raise Exception("ERROR: frame_id contains no points")

    if not point_id in self.frames_dictionary[frame_id]['original_points']:
      raise Exception("ERROR: Adding a result to non exisiting point")

    if not result_tag or not isinstance(result_tag, basestring):
      raise Exception("ERROR: results_tag should be a non empty string")

    if not isinstance(result_id, int):
      raise Exception("ERROR: result_id should be an interger")

    self.frames_dictionary[frame_id]['original_points'][point_id][result_tag] = result_id

  def track_results_to_dictionary(self, tracking_results, result_tag): 
    """ Extrack the tracking results and add them to the dictionary
        Parameters:
          tracking_results: The resultant data_frame from the tracking algorithm
          result_tag: A string to represent that result
    """

    #Make sure the result_tag is not an empty  
    if not result_tag or not isinstance(result_tag, basestring):
      raise Exception("ERROR: results_tag should be a non empty string")

    #Make sure the data_frame has the correct type:
    if not isinstance(tracking_results,pd.DataFrame):
      raise Exception("ERROR: tracking_results should be of type DataFrame")


    if not "particle" in tracking_results:
      raise Exception ("ERROR: the tracking_results do not contain the 'particle' column")
 

    for index, row in tracking_results.iterrows(): 
      frame_id = str(int(row['frame'])).zfill(3)
      point_id = int(row['particle_id'])
      result_particle_id = int(row['particle'])
      self.add_track_result(frame_id, point_id, result_tag, result_particle_id)

  def assess_tracking_result(self, original_id, result_tag, verbose=False):     
    """ Check the tracking result a given particle across the frames
        The function prints if the tracking algorithm returned the same id for all the 
        particle activation and where it missed.
        Parameters:
          original_id(int): The original id of the particle
          result_tag(string): The tag given for the tracking method 
          verbose: Print extra information 
    """

    #Make sure the result_tag is not an empty  
    if not result_tag or not isinstance(result_tag, basestring):
      raise Exception("ERROR: results_tag should be a non empty string")

    if not isinstance(original_id, int):
      raise Exception("ERROR: the original_id should be of type integer")
    
    #Track list for the points that have the same original_group_id 
    original_tracking_list = []  

    #Track list for the points that have different id than the original_group_id 
    other_points_list = []  

    dic = self.frames_dictionary
    #Loop over all the point to retrieve the tracking results
    for  frame_id in dic:
      if 'original_points' in dic[frame_id]:
        for point_id in dic[frame_id]['original_points'].keys():
          if 'original_group_id' in dic[frame_id]['original_points'][point_id]:
              original_particle_id = dic[frame_id]['original_points'][point_id]['original_group_id']
              if original_particle_id == original_id:
                if result_tag in dic[frame_id]['original_points'][point_id]:
                  new_track_id = dic[frame_id]['original_points'][point_id][result_tag]
                  original_tracking_list.append((frame_id, point_id,new_track_id)) 
                else:    
                  raise Exception("ERROR: No {0} found for frame:{1}, particle{2}".format(result_tag, frame_id, point_id))
              elif result_tag in dic[frame_id]['original_points'][point_id]:
                new_track_id = dic[frame_id]['original_points'][point_id][result_tag]
                other_points_list.append((frame_id, point_id, new_track_id))


    #Get the list and frequency of new_track_id
    frequency= {}
    if len(original_tracking_list) == 0: 
      print "ERROR: No {0} points found".format(original_id)
      return 

    for point in original_tracking_list: 
      if not point[2] in frequency:
        frequency[point[2]] = 1
      else:        
        frequency[point[2]] += 1

    #Get the new particle_id with the maximum frequency
    maximum = 0
    track_id = -1 
    for particle, rate in frequency.items():
      if rate > maximum:
        track_id = particle
        maximum = rate

    print "The track result id for the original particle {0} is {1}".format(original_id, track_id)
    nb_false_neg = 0  
    
    #Get the points that was not added to the same track_id 
    for entry in original_tracking_list: 
      if entry[2] != track_id:  
        nb_false_neg +=1 
        if verbose:
          print "WARNING:Frame{0} point{1} got a different group id ".format(entry[0], entry[1])


    nb_false_pos = 0
    #Get the points form a different group that got the same track_id 
    for entry in other_points_list: 
      if entry[2] == track_id:  
        nb_false_pos +=1
        if verbose:
          print "WARNING:Frame{0} point{1} got wrongly added to the group ".format(entry[0], entry[1])

    print "There exists {0} points for original_id:{1}. Tracker miscalculated {2} points, \
    and added to the group {3} points".format(len(original_tracking_list), original_id, nb_false_neg, nb_false_pos)  



  def iterative_tracking(self, destination, stride, frames_ids, params):
    """ Iteratevely transform the frames the destination frame using strides    
        A list of pilot frames which ID are distanted  n * stride from the destination will be constructed.  
        Every frame is registrared to the nearest pilot frame in the direction of the final destination. 
        Then the pilots will be registred to each other to move all the points to the final destination frame space.
        The function will update the frames dictionary with the transformed points 
        Parametes:
          destination(int): The id of the destination frame
          stride(positive integer): The stride to take in transforming the points in every frame to the destination frame. 
          frames_ids(list): List of frames that contains trascription points
          params: The registration parameters
    """ 
    
    #Get the registration pairs
    frames2pilots, upper_p2p, lower_p2p = self.get_iterative_registration_pairs(destination, stride, frames_ids) 

    reg_list = frames2pilots + upper_p2p + lower_p2p
    
    fixed_list = [self.id2str(pair[0]) for pair in reg_list]
    moving_list = [self.id2str(pair[1]) for pair in reg_list]
 
    #Run the registrations
    self.register_frames(fixed_list, moving_list, params)

    #Run the frames to pilots transformations
    for fixed, moving in frames2pilots:
      trans_dic = {}
      trans_dic[self.id2str(fixed)] = self.id2str(moving)
      self.transform_points(trans_dic)
      
    #Run the pilots to pilots transformations
    for fixed, moving in upper_p2p + lower_p2p:
      trans_dic = {}
      trans_dic[self.id2str(fixed)] = self.id2str(moving)
      self.transform_points(trans_dic)

  def get_iterative_registration_pairs(self, destination, stride, frames_ids):
    """ Returns three list of tuples showing the reference and float frames ids that should be registered and transformed.
        Parametes:
          Similar to iterative tracking. 

        Returns:
          [frames_to_pilots_list, upper_pilot_to_pilot_list, lower_pilot_to_pilot_list]
          frames_to_pilots_list: list of frame ids and the closest pilot frame it should register to
          upper_pilot_to_pilot_list:  Ordered list for the pilot to pilot registration for pilot frame_id > destination
          lower_pilot_to_pilot_list:  Ordered list for the pilot to pilot registration for pilot frame_id < destination
    """

    if not isinstance(destination, int):
      raise Exception("ERROR: The destination frame id should be an integer.")

    if destination  < 0 :
      raise Exception("ERROR: The destination frame id should be a non negative integer.")

    if not isinstance(stride, int):        
      raise Exception("ERROR: The stride should be an integer.")

    if stride <= 0:
      raise Exception("ERROR: The stride should be a positive integer.")
           
    if not all(isinstance(x, int) and x >= 0 for x in frames_ids):
      raise Exception("ERROR: All the frames_ids should be  positive integer.")


    #Get the ids of the pilot frames
    frames_array = np.asarray(frames_ids)
    lower_ids = frames_array[frames_array < destination]       
    upper_ids = frames_array[frames_array > destination]
    
    lower_steps = (destination  - lower_ids) / stride
    lower_float = destination - lower_steps * stride
    
    upper_steps = (upper_ids - destination) / stride
    upper_float = destination +  upper_steps * stride


    frames_to_pilots_list = [] 
    pairs = zip(lower_ids, lower_float) + zip(upper_ids, upper_float) 

    #Register the pairs that have different fixed and float ids
    for fixed_frame, float_frame in pairs:
      if fixed_frame != float_frame:
        frames_to_pilots_list.append((fixed_frame, float_frame))


    min_frame = min(frames_array)
    max_frame = max(frames_array)

    upper_pilot_to_pilot_list = [] 
    #Register the upper pilots iteratively
    if max_frame > destination:
      max_pilot = destination + ((max_frame - destination) / stride) * stride
      upper_ref_pilots = range(max_pilot, destination, -stride)             
      upper_flt_pilots = range(max_pilot - stride, destination-1, -stride)
      for fixed_frame, float_frame in zip(upper_ref_pilots, upper_flt_pilots):
        upper_pilot_to_pilot_list.append((fixed_frame, float_frame))


    lower_pilot_to_pilot_list = [] 
    if min_frame < destination:
      min_pilot =  destination - ((destination - min_frame) / stride) * stride
      lower_ref_pilots = range(min_pilot, destination, stride )
      lower_flt_pilots = range(min_pilot + stride, destination + 1, stride)
      for fixed_frame, float_frame in zip(lower_ref_pilots, lower_flt_pilots):
        lower_pilot_to_pilot_list.append((fixed_frame, float_frame))

    return [frames_to_pilots_list, upper_pilot_to_pilot_list, lower_pilot_to_pilot_list]


  def id2str(self, int_id, zero_fill=None):
    """ Convert the frame id to a string with the appropriate zero fill digits
        Parameters: 
          int_id: The integer frame id (can be string too)

        Return:
          str_id: The id converted to string with the approprite number of zero fills
    """
    
    if isinstance(int_id, str):
      int_id = int(int_id)

    if (not isinstance(int_id, int)) or int_id < 0:
      raise Exception("ERROR the frame id should be positive integer. Got {0}".format(int_id)) 

    if zero_fill != None:
      return str(int_id).zfill(zero_fill)
    else:
      return  str(int_id).zfill(3)

  def overlay_tracking_results(self, destination_dir, tracking_results, alpha_in = 0.40, source_dir=None):
    """ Assign a color for every set of tracked transcription points and overlay 
        a circle around the points in the corresponding original frame. 


        Parameters:
          destination_dir: The directory to save overlayed image series 
          tracking_results: The dataframe that includes the tracking results. 
            It should has the columns: frame, x, y, and particle.
          alpha_in: the transparance of the circles. It should be between 0 and 1. 
          source_dir: Optiona directory for the stacked images to be overlayed. 
            The images should have the same sizes, and names as the images used during the tracking tracking. 
  """

    from matplotlib import cm
    from PIL import Image, ImageDraw, ImageFont

    #Make sure the data_frame has the correct type:
    if not isinstance(tracking_results,pd.DataFrame):
      raise Exception("ERROR: tracking_results should be of type DataFrame")

    if not "particle" in tracking_results:
      raise Exception ("ERROR: the tracking_results do not contain the 'particle' column")

    if not "frame" in tracking_results:
      raise Exception ("ERROR: the tracking_results do not contain the 'frame' column")

    if not "x" in tracking_results:
      raise Exception ("ERROR: the tracking_results do not contain the 'x' column")

    if not "y" in tracking_results:
      raise Exception ("ERROR: the tracking_results do not contain the 'y' column")

    if alpha_in < 0 or alpha_in > 1: 
      raise Exception ("ERROR: alpha_in should be between 0 and 1")
      

    self.create_dir(destination_dir)

    tr = tracking_results
    # The number of particles 
    #n_particles =  len(set(tr['particle'].tolist()))
    particles_ids =  tr.particle.unique().tolist()
    n_particles = len(particles_ids)

    color_step =  255  / float(n_particles) 

    #Distribute the colour over the 255 range for the color map
    streched_ids = [ int(math.floor(float(set_id) * color_step)) for set_id in range (0,n_particles)] 
    scalar_color = [  cm.Accent(s_id, bytes=True, alpha = alpha_in) for s_id in streched_ids] 
    rgb_color_map = [ (particle, scalar_color[set_id]) \
      for particle, set_id in zip(particles_ids, range(0, n_particles))] 
    color_dict = dict(rgb_color_map)
   
    #define a font
    fnt = ImageFont.load_default()

    #Get the images: 
    if source_dir !=None:
      if not os.path.exists (source_dir):
        raise Exeption("The directory {0} does not exist.".format(source_dir))
      image_dir = source_dir
    else:
      image_dir = self.images_dir

    images_regex= self.frame_prefix + "*" + self.frame_suffix
    image_files = glob.glob(os.path.join(image_dir, images_regex))

    temp_image_name='temp8bit.tif' 
    file_regex = re.compile(".*" + self.frame_prefix + "(?P<id>[0-9]*)" + self.frame_suffix)
    for frame in image_files:
      #Get the frame id:
      frame_id = file_regex.search(frame).groupdict()['id'] 
      original_frame_points = self.get_original_points(self.id2str(frame_id))

      #A list of tuples that will contains the information about the particles in a given frame
      points_color_information = []

      #Get the particles in that frame and get the corresponding color information
      for original_frame_id, point_id, indexes in original_frame_points:

        mapped_particle_row = tr[(tr.frame == int(frame_id)) & (tr.particle_id == point_id)]

        #Double check that we have only one particle with a given point_id in a frame
        if mapped_particle_row.shape[0] != 1:  
          raise Exception ("Error got more than one particle with the same id in one frame")

        mapped_particle_id = mapped_particle_row.iloc[0].particle
        points_color_information.append((color_dict[mapped_particle_id], indexes, mapped_particle_id))

            
    
      #First, read the image using skimage as it supports all types then convert it to unit8
      original_image = skimage.io.imread(frame)
      v_min, v_max = np.percentile(original_image, (0.2, 99.8))
      better_contrast = exposure.rescale_intensity(original_image, in_range=(v_min, v_max))
      #better_contrast = exposure.rescale_intensity(original_image)   
      new_image = better_contrast * 255.0 / better_contrast.max()
      skimage.io.imsave(temp_image_name, new_image.astype('uint8'))

      #Convert to frame to  RGBA
      orig_frame_image = Image.open(temp_image_name)
      rgba_frame = Image.new("RGBA", orig_frame_image.size)
      rgba_frame.paste(orig_frame_image)
      contours_frame = Image.new('RGBA', rgba_frame.size, (255, 255, 255, 0))
      draw= ImageDraw.Draw(contours_frame)

      # Create a colored circle around every point  
      for color, indexes, mapped_id in  points_color_information: 
        top_x = indexes[0] - 10
        top_y = indexes[1] -10
        bottom_x = indexes[0] + 10 
        bottom_y = indexes[1] + 10 
        draw.ellipse((top_x, top_y, bottom_x , bottom_y), outline=color)

        # typethe particle id at the lower corner
        draw.text((bottom_x,bottom_y), str(int(mapped_id)), font=fnt, fill=color)

      #Blend the original image with the contours
      out = Image.alpha_composite(rgba_frame, contours_frame)

      # Save the new file in the output directory
      out.save(os.path.join(destination_dir,os.path.basename(frame)))
     
    os.remove(temp_image_name)
    print "saved images in {0}".format(os.path.abspath(destination_dir))

  
