import ConfigParser
import glob
import argparse
import os


parser = argparse.ArgumentParser(description="A script to generate configuration files for cells, and to generate a swarm that can run on Biowulf.")

parser.add_argument('cells_dir', help="The top level directory that contains the cells images directories, and the corresponding segmented rna.csv. The rna file should have the same prefix as the cell images directory.")
parser.add_argument('input_master', help="The generic configuration file")
parser.add_argument('-o', '--overwrite', action="store_true", help="Overwirte the configuration files if they exist.")


args = parser.parse_args()

general_config = args.input_master
cells_dir = args.cells_dir

if not os.path.isfile(general_config):
  raise Exception ('ERROR: The configuration file "{0}" does not exist.'.format(general_config))


if not os.path.isdir(cells_dir):
  raise Exception ('ERROR: The cells directory "{0}" does not exist.'.format(cells_dir))

cells_dir = os.path.abspath(cells_dir)


if not "GENE_TRAP_ROOT" in os.environ:
  raise Exception ('ERROR: The environment variable "GENE_TRAP_ROOT" is not defined.')

biowulf_runme = os.path.join(os.environ["GENE_TRAP_ROOT"],"utils/biowulf_runme")
if not os.path.isfile(biowulf_runme):
  raise Exception("ERROR: Can not locate the runme script {0}".format(biowulf_runme))
    


#Parse the general configuration file.
config = ConfigParser.ConfigParser()
try:
  config.read(general_config)
except ConfigParser.Error as error:
  error_exit('ERROR: unable to parse the configuration file {0}.\nERROR details: {1}'.format(general_config,error))

#Get the list of directories
cells = [x for x in os.listdir(cells_dir) if os.path.isdir(x)]


#This is the top level tracking direcotries for all the cells.
base_tracking_dir = config.get('arguments', 'processing-dir')

print "found {0} cell directories".format(len(cells))
#Modify the config file to include the cell directory
processing_dirs=[]
for cell_name in cells: 

    #I should find a corresponding <cell_name>.csv file    
    #These names will come from the conventions
    stack_dir = os.path.abspath(os.path.join(cells_dir, cell_name))
    segmented_rnas = os.path.abspath(os.path.join(cells_dir, cell_name + '.csv'))
    tracking_dir = os.path.abspath(os.path.join(base_tracking_dir, cell_name + '-processing'))


    if not os.path.exists(segmented_rnas):
      print "WARNING: Can not find the segmented rnas files {0}".format(segmented_rnas)
      continue


    if not os.path.exists(tracking_dir):
      os.mkdir(tracking_dir) 

    #Get the number of frames for the cell:
    image_prefix = config.get('arguments', 'images-prefix')
    search_prefix = os.path.join(stack_dir, image_prefix) + '*' 
    images = glob.glob(search_prefix)
    n_frames = len(images)
    if n_frames == 0:
      raise Exception('ERROR: Can not find any image that match the prefix "{0}"'.format(search_prefix))

    destination_frame = n_frames / 2
    config.set('arguments', 'stack-dir', stack_dir)
    config.set('arguments', 'segmented-RNAs', segmented_rnas)
    config.set('arguments', 'processing-dir', tracking_dir)
    config.set('arguments', 'destination', destination_frame)
    config.set('arguments', 'memory', n_frames)

    config_path = os.path.join(tracking_dir, "config.cfg")
    if not args.overwrite:
      if os.path.isfile(config_path):
        raise Exception ('Error, the configuration file "{0}" already exists'.format(config_path))

    with open(config_path, 'wb') as configfile:
        config.write(configfile)

    processing_dirs.append(tracking_dir)
    print "Added  {0} to the joblist".format(tracking_dir)
   

#Generate a swarm file:
commands = [] 
for job in processing_dirs:
   commands.append(biowulf_runme + " " + job)


swarm_file_name = "tracking-jobs.swarm"
with open(swarm_file_name, 'wb') as swarm_file:
    swarm_file.write("\n".join(commands))


print "\n\nGenerated the swarm file {0}".format(swarm_file_name)
print "To execute the swarm, run the command:"
print "swarm -f {0} --time [requested-time] ".format(swarm_file_name)
