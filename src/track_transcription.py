import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import pims
import trackpy as tp
import re
import shutil
import threading
import multiprocessing
import os

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

                                          |-->'original_points'(dict)-->point_id(key):index[x,y](value)
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


#####PUBLIC CLASS FUNCTIONS

  def preprocess_track_file(self, track, index): 
    """Create a numpy array that contains the frame numer and location of the object in only the frames where the object was physically present

        Parameters: 
        track: file that indicates the spacial location of the object in each frame
        index: file that indicates which objects are physically present in each frame


        Returns:
        [length(file),3] np array
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
      self.add_original_point_to_dictionary(frame_id, [x, y])

    return positions


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


      #Get the points in the frame (original and transformed)
      original_points = {}
      transformed_points = {}
      if ref_id in self.frames_dictionary: 
        if "original_points" in self.frames_dictionary[ref_id]:
          original_points =  self.frames_dictionary[ref_id]['original_points']
        if 'transformed_points' in self.frames_dictionary[ref_id]:
          transformed_points = self.frames_dictionary[ref_id]['transformed_points']


      #Loop over original and transformed to create the list of points
      points_indexes_list = []
      original_frame_id_list = []
      point_id_list = []
      n_points = len(original_points)

      for  point_id, indexes in original_points.items(): 
        points_indexes_list.append(indexes)
        original_frame_id_list.append(ref_id)
        point_id_list.append(point_id) 

      for original_frame_id, points_id in transformed_points.items():
        n_points += len(points_id)
        for point, indexes in points_id.items(): 
          points_indexes_list.append(indexes)
          original_frame_id_list.append(original_frame_id)
          point_id_list.append(point) 


      if n_points > 0: 

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
          raise Exception ("ERROR: The number of transformed points {0} does not equal the number of input points {1}".format(n_transformed_points, n_poinst))
 

        #Write the results back in the dictionary 
        for original_frame_id, point_id, transformed_indexes in zip(original_frame_id_list, point_id_list, output_points):
          self.add_transformed_point_to_dictionary(float_id, original_frame_id, point_id, transformed_indexes)


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
      raise Exception ('ERROR: The transformix call "{0}" returned with code {1}'.format(elastix_call, status))
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
        raise Exception('ERROR: The file "{0}" does not exist.'.format(file_name))
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

    print len(work_list)
    from multiprocessing.pool import ThreadPool
    nthreads = multiprocessing.cpu_count()
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
    """Tracks object and plots movement
        Parameters:
          np_array: generated from transformix_results_to_array
          memory: maximum number of frames an object can be absent before it is declared a separate object
          distance: search range between frames for the movement of one object
    """
    t = tp.link_df(np_array,search_range = distance, memory = memory)
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
    """Creates a directory even if it exists 
      Prameters:
        dir_name: The directory path

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
      regex = re.compile(".*OutputIndexFixed\s*=\s*\[\s*(?P<x>[0-9]*)\s*(?P<y>[0-9]*)\s*\]")
    else:
      regex = re.compile(".*OutputPoint\s*=\s*\[\s*(?P<x>\d+\.\d+)\s*(?P<y>\d+\.\d+)\s*\]") 

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


  def add_original_point_to_dictionary(self, original_frame_id, indexes):
    """Adds an original point to the frames_dictionary.
      Parameters:
        original_frame_id: The original frame where the point belongs. 
        indexes: the indexes of the point in the destination frame id.
    """

    #Make sure indexes contains two number for X and Y.
    if len(indexes) != 2:
      print indexes
      raise Exception ("ERROR: the lenght of indexes does not equal 2.")
    
    if not original_frame_id in self.frames_dictionary:
      self.frames_dictionary[original_frame_id] = {}

    if not 'original_points' in self.frames_dictionary[original_frame_id]:
      self.frames_dictionary[original_frame_id]['original_points']={}
      self.frames_dictionary[original_frame_id]['n_points']=0
  
    #Increament 
    self.frames_dictionary[original_frame_id]['n_points'] += 1
    point_id = self.frames_dictionary[original_frame_id]['n_points']
    self.frames_dictionary[original_frame_id]['original_points'][point_id] = indexes  


  def frame_points_to_dataframe(self, frame_id):
    """Collect the original and transformed points from a frame and returns a data frame
      Parameters:
        frame_id: The ID of the frame where the points will be collected 

      Returns:
        data_frame: The DataFrame object, None if there are points in the frame.
    """

    if frame_id == None: 
      raise Exception ("Invalid frame_id")

    #frame_id should be converted to int for trackpy
    int_frame_id = int(frame_id)

    if not frame_id in self.frames_dictionary:
      print "WARNING: the frame_id {0} is not part of the frames dictionary"
      return None

    temp_np = np.empty((0,4))
    #Add original points
    if 'original_points' in  self.frames_dictionary[frame_id]:
      for point_id , indexes in self.frames_dictionary[frame_id]['original_points'].items():
        temp_np = np.append(temp_np,[[int_frame_id,indexes[0], indexes[1], point_id]], axis=0)
          
    #Add transformed points
    if 'transformed_points' in  self.frames_dictionary[frame_id]:
      for original_frame_id , points_ids in self.frames_dictionary[frame_id]['transformed_points'].items():
        for point_id , indexes in self.frames_dictionary[frame_id]['transformed_points'][original_frame_id].items():
          temp_np = np.append(temp_np,[[int(original_frame_id),indexes[0], indexes[1], point_id]], axis=0)

  
    if temp_np.shape[0] == 0:
      return None
    else:    
      return DataFrame(data = temp_np, columns = ['frame', 'x', 'y', 'particle_id'])

    return temp_np

