from track_transcription  import Track
import numpy as np
import ConfigParser
import argparse
import glob
import os
import sys
import skimage.io
from PIL import Image
from util_functions import error_exit, read_config





parser = argparse.ArgumentParser(description="A script to run the iterative registration and tracking algorithm")
parser.add_argument('input_master', help="The master configuration file")
args = parser.parse_args()
config_file = args.input_master


options_dic =  read_config(config_file)

images_stack = options_dic['stack-dir']
images_prefix = options_dic['images-prefix'] 
images_suffix= options_dic['images-suffix']
processing_dir = options_dic['processing-dir']
override_registration = options_dic['override-registration'] 
segmented_rnas = options_dic['segmented-rnas'] 
point_index = options_dic['point-or-index']
registration_file = options_dic['registration-file'] 
destination = int(options_dic['destination'])
stride = int(options_dic['stride'])

try: 
  resolution = eval(options_dic['resolution'])
except Exception as e:
  error_exit("Can not convert {0} to a tuple".format(options_dic['resolution']))


#Get the dimensions of the images
images_regex= images_prefix + "*" + images_suffix
full_regex = os.path.join(images_stack,images_regex)
image_files = glob.glob(full_regex)

if image_files == None or len(image_files) == 0 :
  error_exit("There are no files that matches the regex:{0}".format(full_regex))

#sample_image = skimage.io.imread(image_files[0])
#dimensions = sample_image.shape

sample_image = Image.open(image_files[0])
dimensions =  sample_image.size

print "The image dimensions are:{0}".format(dimensions)

track_obj =  Track(processing_dir, override_results=override_registration)

track_obj.set_image_information(images_stack, images_prefix, images_suffix)
track_obj.set_point_index(point_index)
#added_points =  track_obj.preprocess_rna_segmentation(segmented_rnas, resolution, dimensions)
#added_points =  track_obj.preprocess_rna_segmentation(segmented_rnas, resolution, dimensions)
frames = track_obj.populate_rna_segmentation(segmented_rnas, resolution, dimensions)


#trak_obj.particles_df.to_csv("csv_particles.csv")
#populate the object with the segmented spots and get the frames
#frames = [point[0] for point in added_points ]



track_obj.iterative_tracking(destination, stride, frames, registration_file)
print "The image dimensions are:{0}".format(dimensions)


dest_str = track_obj.id2str(destination)
df = track_obj.frame_points_to_dataframe(dest_str)


#Track and merge the tracked data_frame with the original dataframe 
try:
  memory = int(options_dic['memory'])
except:
  memory = len(image_files)
  print "No memory for the tracking found. Using the maximum value of:{0}".format(memory)

try:
  distance = int(options_dic['distance'])
except:
  print "No maximum distance the partcile can move in is found. Using distance=20"
  distance = 20

track_obj.track(df,memory , distance)
print "Tracked the particles using  memory={0}, distance={1}".format(memory, distance)

df.to_csv("tracker-results.csv")

#Overlay the results


try:
  adjusted_dir = options_dic['adjusted-dir'] 
except:
  adjusted_dir=None
  pass

track_obj.overlay_tracking_results(os.path.join(processing_dir,"overlayed"), df, source_dir=adjusted_dir)

#stride_str = str(stride)
#track_id = 'to_frame_' + dest_str + '_stride_' + stride_str
#track_obj.track_results_to_dictionary(df, track_id)
#track_obj.assess_tracking_result(0,track_id)
#track_obj.assess_tracking_result(1, track_id)
#track_obj.assess_tracking_result(2, track_id)

