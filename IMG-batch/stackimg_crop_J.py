
# For opening a set of images as a stack, crop, save cropped image.

import os
from ij.io import DirectoryChooser
from ij import IJ, ImagePlus, VirtualStack

def OPENFolder():
	try:
		srcDir = DirectoryChooser("Choose!").getDirectory()
	except:
		return	
	return(srcDir)

def OPENFolder_as_stack(ext, NoSub, srcDir, name):
	# ext should be like '.jpg'
	# Not including sub directories
	NoSub = 'NoSub'
	#srcDir = OPENFolder()
	vs = None
	for root, directories, filenames in os.walk(srcDir):
		for filename in filenames:
			if not filename.endswith(ext):
				continue
			path = os.path.join(root, filename)
			print(path)
			if vs is None:
				imp = IJ.openImage(path)
				vs = VirtualStack(imp.width, imp.height, None, srcDir)
			vs.addSlice(path[len(srcDir):])
		if NoSub == 'NoSub':
			break
	return ImagePlus(name, vs)

def crop():
	for i in range(2,80):
		path = 'E:\\BIOBIO\\%s\\solidgrowth_%s_row_3.tif' % (setdir,str(i).zfill(2))
		print(path)
		im = IJ.openImage(path)
		im.show()
		process('solidgrowth_%s_row_3' % (str(i).zfill(2)))
		im.close()
def process(filename):
	# just to pass pramaters down
	#im = OPENFolder_as_stack()
	#im.show()
	cropping([0, 0, 2160, 2208],'%s_JA'%filename)
	cropping([2160, 0, 2160, 2208],'%s_MM'%filename)
	"""
	cropping([1341, 441, 266, 266],'%s_5_1' %output,targetdir)
	cropping([867, 474, 266, 266],'%s_4_1' %output,targetdir)
	cropping([1392, 936, 266, 266],'%s_4_2' %output,targetdir)
	cropping([405, 486, 266, 266],'%s_3_1' %output,targetdir)
	cropping([902, 988, 266, 266],'%s_3_2' %output
	cropping([1409,1487, 266, 266],'%s_3_3' %output,targetdir)
	cropping([422, 1001, 266, 266],'%s_2_1' %output,targetdir)
	cropping([935, 1490, 266, 266],'%s_2_2' %output,targetdir)
	cropping([467, 1526, 266, 266],'%s_1_1' %output,targetdir)
	"""
	#im.close()
def cropping(rec,filename_final):
	IJ.makeRectangle(rec[0],rec[1],rec[2],rec[3])
	IJ.run("Duplicate...", "duplicate")
	im_sub = IJ.getImage()
	IJ.run("Save", "save=[E:\\BIOBIO\\JASvs_QL109_2.slices\\plates\\%s.tif]" %(filename_final))
	im_sub.close()


def savestacks():
	for root, dirs, files in os.walk(setdir):
		for dir in dirs:
			path = os.path.join(rootpath, root, dir)
			im = OPENFolder_as_stack('.tif', 'NoSub', path, path[len(path)-4:len(path)])
			im.show()
			IJ.saveAs(im, 'tif', path)
			im.close()
		break

def openstacks():
	for targetdir in targetdirs[3:4]:
		for target in targets:
			IJ.open('E:\\BIOBIO\\JASvs.QL109.slices\\plates\\Singlecolonystacks\\%s_%s.tif'%(targetdir,target))
		
targets = ['1_1','2_1','2_2','3_1','3_2','3_3','4_1','4_2','5_1']
targetdirs = ['row_1_JA','row_1_MM','row_2_JA','row_2_MM','row_3_JA','row_3_MM']
setdir = 'JASvs.QL109.slices\plates'
rootpath = os.getcwd()


#savestacks()
#crop()
#openstacks()
im = OPENFolder_as_stack('.tif', 'NoSub', OPENFolder(), '5-1')
im.show()
IJ.makeOval(90, 94, 90, 90);
