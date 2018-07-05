import os
import sys

def main(pts_dirname):
	pts_fname = []
	for file in os.listdir(pts_dirname):
		parts = file.split('.')
		if len(parts) == 2:
			fname, fext = parts
			if fext == 'pts':
				pts_fname.append(int(fname))
	pts_fname.sort()
	for i in xrange(len(pts_fname)):
		filename = os.path.join(pts_dirname, '.'.join([str(pts_fname[i]), 'pts']))
		newname = os.path.join(pts_dirname, '.'.join([str(i), 'pts']))
		os.rename(filename, newname)




if __name__ == '__main__':
	if len(sys.argv) != 2:
			print("Usage: rename_pts.py [pts_dirname]")
			exit()
	main(sys.argv[1])