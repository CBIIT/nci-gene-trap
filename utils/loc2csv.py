
import pandas
import glob
import argparse
import os


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="A script to convert .loc files generated using localize to .csv files. The columns of the .loc file should be sapce separated and in the order of x y intensity frame")
  parser.add_argument('input', help='A direcotry that contains .loc file, or a single .loc filename. The columns of the .loc file should be sapce separated and in the order of "x" "y" "intensity" "frame"')
  args = parser.parse_args()

  if not os.path.exists(args.input):
    raise Exception ("Can not find the input file/directory {0}".format(args.input))

  input_files_list = []

  #Check if  input is a directory or a file 
  if os.path.isfile(args.input):
    input_files_list.append(args.input)


  #If the input is a directory, get all the .loc file within that directory 
  if os.path.isdir(args.input):
    input_files_list = glob.glob(args.input + "*.loc")
    print "found {0} loc files in {1}".format(len(input_files_list), args.input)


  # A file to store the prepended loc file
  temp_file = "loc2csv_temp_file.csv"

  for loc_file in input_files_list: 
    csv_file_name = loc_file.replace(".loc", ".csv")
  
    #Append a line of x y intensity frame to the loc file
    with open(loc_file, 'r') as file_stream: 
      lines = file_stream.readlines()
      lines.insert(0, " x y  intensity frame \n")

    with open(temp_file, 'w') as temp_stream:
      temp_stream.write("".join(lines))

    #Read the temp file using pandas
    particles = pandas.read_csv(temp_file, sep= ' *')

    #Get rid of the zero values
    true_particles = particles[(particles.x != 0) & (particles.y != 0)]
    true_particles.to_csv(csv_file_name,index=False)
    print "convert {0} to {1}".format(loc_file, csv_file_name)

  #Cleanup the temp file if it exists
  if os.path.exists(temp_file):
    os.remove(temp_file)


