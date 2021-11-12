import sys
import os 
import numpy as np
import gzip
def main():
	tadfile = sys.argv[1]
	chr = sys.argv[2]
	contactfile = sys.argv[3]
	faifile = sys.argv[4]
	oup = sys.argv[5]

	f1 = open(faifile,'r')
	fai_dict = {}
	for line1 in f1:
		info1 = line1.strip().split('\t')
		fai_dict[info1[0]] = int(info1[1])
	f1.close()

	chr_len = fai_dict[chr]

	if chr_len%20000 == 0:
		bin_num = chr_len / 20000
	else:
		bin_num = chr_len / 20000 + 1

	bin1 = 0
	f2 = gzip.open(contactfile,'r')
	contact_dict ={}
	for line2 in f2:
		info2 = line2.strip().split('\t')
		if len(info2) != bin_num:
			print('this is wrong!!!~')
			sys.exit(0)

		for bin2 in range(0,len(info2),1):
			if bin1 > bin2:
				continue
			contact_dict[(bin1,bin2)] = float(info2[bin2])
		bin1+=1
	f2.close()

	out = open(oup,'w')
	f3 = open(tadfile,'r')
	for line3 in f3:
		info3 = line3.strip().split('\t')
		if info3[0] != chr:
			continue
		start = int(info3[1])
		end = int(info3[2])

		if end - start <= 20000:
			info_list = info3 + ['NA']
			out.write('{0}\n'.format('\t'.join([str(i) for i in info_list])))
			continue

		bin_s = start / 20000
		bin_e = end / 20000

		lines = bin_e - bin_s

		Dscore_dict = {}
		for i in range(1,lines,1):
			Dscore_dict[i] = [[],[],[]]

		for bin1 in range(bin_s,bin_e,1):
			for bin2 in range(bin_s-lines,bin_s,1):
				if bin1-lines >= bin2:
					continue
				i = abs(bin2 - bin1)
				try:
					Dscore_dict[i][0].append(contact_dict[(bin2,bin1)])
				except KeyError as reason:
					continue

		for bin1 in range(bin_s,bin_e,1):
			for bin2 in range(bin_e,bin_e+lines,1):
				if bin1 <= bin2-lines:
					continue
				i = abs(bin2 - bin1)
				try:
					Dscore_dict[i][1].append(contact_dict[(bin1,bin2)])
				except KeyError as reason:
					continue

		for bin1 in range(bin_s,bin_e,1):
			for bin2 in range(bin_s,bin_e,1):
				if bin1 >= bin2:
					continue
				i = abs(bin2 - bin1)
				Dscore_dict[i][2].append(contact_dict[(bin1,bin2)])

		a_list = []
		for i in range(1,lines,1):
			# print(len(Dscore_dict[i][0]),len(Dscore_dict[i][1]),len(Dscore_dict[i][2]))
			if len(Dscore_dict[i][0]) == 0:
				mean_val = np.mean(Dscore_dict[i][1])
			elif len(Dscore_dict[i][1]) ==0:
				mean_val = np.mean(Dscore_dict[i][0])
			else:
				mean_val = np.mean([np.mean(Dscore_dict[i][0]),np.mean(Dscore_dict[i][1])])
			
			a_list.append(np.median(Dscore_dict[i][2])/mean_val)

		value = np.median(a_list)

		info_list = info3 + [value]
		out.write('{0}\n'.format('\t'.join([str(i) for i in info_list])))
	out.close()
	# write.close()
	f3.close()


if __name__ == '__main__':
	main()
