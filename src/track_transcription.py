
# coding: utf-8

# In[1]:

#Python files for taking in point files and returning tracking values

#from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

#import matplotlib as mpl
#import matplotlib.pyplot as plt

# Optionally, tweak styles.
#mpl.rc('figure',  figsize=(10, 6))
#mpl.rc('image', cmap='gray')

import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import pims
import trackpy as tp
import re
import os

class track_transcription:
    

    
    #####PUBLIC CLASS FUNCTIONS
    
    def preprocess_track_file(self, track, index): 
        """Create a numpy array that contains the frame numer and location of the object in only the frames where the object was physically present


        Parameters: 
            track: file that indicates the spacial location of the object in each frame
            index: file that indicates which objects are physically present in each frame
            
            
        Returns:
            [length(file),3] np array
        """
        
        
        location_file_1 = open(track)
        location_file_2 = open(track)
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
        x = [float(x[5:11]) for x in location_file_1]
        y = [float(y[20:26]) for y in location_file_2]
        x.pop(0)
        y.pop(0)
        locations[:,1] = x
        locations[:,2] = y
        locations[:,0] = list(range(1,length+1))
    
        positions = locations[indices]
        return positions
    
    

    def frames_id_to_elastix(self, list_np_array, output_frames_file_location):
        """Creates a list of frames where there is at least one particle present and writes it to an output file
        
        
        Parameters:
            list_np_array: list of numpy arrays generated from the preprocess_track_file function
            output_frames_file_location: file location for frames to be written to
            
            
        Returns:
            list of frames
        """
        
        
        frames_id = []
        if len(list_np_array) == 1:
            for i in range(len(list_np_array[0])):
                frames_id.append(int(list_np_array[0][i,0]))
        else:
            for i in range(len(list_np_array)):
                for j in range(len(list_np_array[i][:,0])):
                    frames_id.append(int(list_np_array[i][j,0]))
                
        t = list(set(frames_id))
        file = open(output_frames_file_location + 'elastix_frames_file.txt','w')
        for i in t:
            file.write(str(i) + '\n')
        return t
    
    
    def run_elastix_2(script_location, output_directory, format_, image_directory, parameter_file, slice_prefix, frames_input):
        os.system('mkdir ' + output_directory + '/all_images')
        os.system('mkdir ' + output_directory + '/elastix')

        for i in frames_input:
            moving = ' -m ' + image_directory + '/' + slice_prefix + '2.tif '
            fixed = '-f ' + image_directory + '/' + slice_prefix + str(i) + '.tif '
            parameter = '-p ' + parameter_file + ' '
            output = '-out ' + output_directory + '/frame_' + str(i)
        
            line = ' "' + moving + fixed+ parameter + output + '"'
            print(script_location + line)
            os.system('mkdir '+ output_directory +'/frame_' + str(i))
            os.system(script_location + line)
            os.system('cp ' + output_directory + '/frame_' + str(i) + '/result.0.tif ' +  output_directory + '/all_images/frame.' + str(i) + '.tif')
  
    
    def run_transformix_2(self, run_transformix_location, transformix_input_points_directory, output_directory, frames_input, object_name):
        os.system('mkdir '+ output_directory + '/transformix/')
        os.system('mkdir '+ output_directory + '/transformix/' + object_name)
        os.system('mkdir ' + output_directory + '/transformix/' + object_name + '/all_output_points')

    
        for i in frames_input:
            i = int(i)
            os.system('mkdir ' + output_directory + '/transformix/' + object_name + '/frame_' + str(i))
        
        
            input_points = ' -def ' + transformix_input_points_directory + '/inputPointsFrame' + str(i) + '.0.txt'
            output = ' -out ' + output_directory + '/transformix/' + object_name + '/frame_' + str(i)
            parameter = ' -tp ' + output_directory + '/elastix/frame_' + str(i) + '/TransformParameters.0.txt'
        
            transformix_call = input_points + output + parameter
            os.system(run_transformix_location + ' "'+ transformix_call + '"')
        
       

#    def run_elastix(self, run_elastix_location, output_directory, format_, image_directory, parameter_file, slice_prefix, frames_input):
#        """Runs elastix for entire cell
#        
#        
#        Parameters:
#            run_elastix_location: where the bash script is located
#            output_directory: where output will be placed
#            format_: image file type returned - preferrably .tif
#            image_directory: folder where the input images are located
#            parameter_file: parameter file elastix uses to register each image
#            slice_prefix: prefix for each image from the image directory
#            frames_input: file where the frames to run elastix are indicated
#        """
#        
#        
#        os.system(run_elastix_location + ' ' + output_directory + ' '+  format_ + ' ' + image_directory + ' ' + 
#                  parameter_file + ' ' + slice_prefix + ' ' + frames_input)
#    
#    
#
#    def run_transformix(self, run_transformix_location, transformix_input_points_directory, trans_output_folder, trans_parameter_output_folder, frames_input):
#        """Runs transformix for a given particle
#        
#        Parameters:
#            run_transformix_location: where the bash script is located
#            input_points_folder: where the input points to run transformix are located
#            trans_output_folder: where the output points files will be written
#            frames_input: txt file that contains which frames to run transformix
#        """
#        
#        
#        os.system(run_transformix_location + ' ' + transformix_input_points_directory + ' ' + trans_output_folder + ' ' + 
#                  trans_parameter_output_folder + ' ' + frames_input)               
#    
 

    def generate_input_point_files(self, np_array, transformix_input_points_directory):
        """Generates input files for transformix parameters
        
        Parameters:
            np_array: numpy array genertated from preprocess_track_file
            transformix_input_points_directory: directory to copy transformix input points files into
        """
        
        
        for i in range(len(np_array[:,0].tolist())):
            input_points_file = open(transformix_input_points_directory + '/inputPointsFrame.' + str(np_array[i,0]) + '.txt', 'w')
            input_points_file.write('index\n1\n')
            input_points_file.write(str(np_array[i,1]) + ' ' + str(np_array[i,2]))
    
   

    def transformix_results_to_array(self, np_array, transformix_points_directory):
        """Reads Transformix results and copies the frames and output points to a np array
        
        
        Arguments:
            np_array: generated from preprocess_track_file
            transformix_points_directory: location where transformix points were written into
            
            
        Returns:
            [length(file),3] np array
        """
        
        
        path = transformix_points_directory + '/Frame_'
        particle = np.zeros((len(np_array[:,0].tolist()),3))
        for i in range(len(np_array[:,0].tolist())):
            points = self.find_points_from_transformix_output(path + str(int(np_array[i,0])) + '/outputpoints.txt')
            particle[i,1] = float(points[0])
            particle[i,2] = float(points[1])
            particle[:,0] = np_array[:,0]
        p = DataFrame(data = particle, columns = ['frame','x','y'])
        return p
            

        
    def track(self, np_array, memory, distance):
        """Tracks object and plots movement
        
        
        Parameters:
            np_array: generated from transformix_results_to_array
            memory: maximum number of frames an object can be absent before it is declared a separate object
            distance: search range between frames for the movement of one object
        """
        
        t = tp.link_df(np_array,search_range = distance, memory = memory)
        return t # tp.plot_traj(t)
        
        
    
    ###PRIVATE CLASS FUNCTIONS   
    def find_points_from_transformix_output(self,file):
        """Finds output points in a transformix outputpoints.txt file
        
        
        Parameters:
            Transformix outputpoints.txt file
        
        
        Returns:
            size 2 list
        """
        points = open(file)
        a = re.findall("OutputIndexFixed\s*=\s*\[\s*[0-9]*\s*[0-9]*\s*\]", points.read())
        b = []
        for j in range(len(a)):
            b.append([int(s) for s in a[j].split() if s.isdigit()])
        c = []
        for i in range(len(b[0])):
            c.append(b[0][i])
        return c


