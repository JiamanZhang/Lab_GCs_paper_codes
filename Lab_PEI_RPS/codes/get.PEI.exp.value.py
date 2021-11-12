import sys
import glob
import os
import numpy as np
import math
def main():
	inp = sys.argv[1]
	oud2 = sys.argv[2]
	cutFDR = float(sys.argv[3])
	distancenum = int(sys.argv[4])
	
	out2 = open(oud2,'w')
	out2.write('geneid\tEnhancer_num\tRPS\n')

	f1 = open(inp,'r')
	geneid_list = []
	mystrength_dict = {}
	for line1 in f1:
		info1 = line1.strip().split('\t')
		chr = info1[0]
		geneid = info1[-3].split(':')[0]
		mystrength = float(info1[-1]) - float(info1[-3].split(':')[-2])
		mystrength=math.log(mystrength+1)/math.log(10)
		distance = float(info1[-3].split(':')[1][:-2])

		if abs(distance) < distancenum:
			# print "this is wrong!!!"
			# sys.exit(0)
			continue
		FDR = float(info1[-3].split(':')[2])
		if FDR < cutFDR:
			try:
				mystrength_dict[geneid].append(mystrength)
			except KeyError as reason:
				mystrength_dict[geneid] = []
				geneid_list.append(geneid)
				mystrength_dict[geneid].append(mystrength)
			
	f1.close()

	for geneid in geneid_list:
		info_list = [geneid,len(mystrength_dict[geneid]),np.sum(mystrength_dict[geneid])]
		out2.write('{0}\n'.format('\t'.join(map(str,info_list))))
	out2.close()


if __name__ == '__main__':
	main()