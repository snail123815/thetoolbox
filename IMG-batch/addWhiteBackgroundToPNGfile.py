import os
import subprocess
import sys

# Uses ImageMagick's convert command. Make sure it is installed in the system.
# Will set background color to white and remove alpha layer

# Input and output folder and  can be specified in command line or here in the script

input = "/mnt/d/Desktop/New folder"
output = input

if len([arg for arg in sys.argv if not arg.endswith('.py')]) !=0:
	input = sys.argv[1]
	try:
		output = sys.argv[2]
	except:
		pass
else:
	pass

removeAlpha = False

try:
	os.mkdir(output)
except:
	pass

for file in os.listdir(input):
	if file.endswith('png'):
		args = ['convert', os.path.join(input, file), '-background', 'white', '-alpha', 'remove', os.path.join(output, f'{file.split(".png")[0]}_whitebackground.png')]
		if not removeAlpha:
			args.pop[4,5]
		print(args)