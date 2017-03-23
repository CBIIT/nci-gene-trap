import sys
import os
import ConfigParser


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


