import sys
import os
def main():
	borderfile = sys.argv[1]
	ind = sys.argv[2]
	samplefile = sys.argv[3]
	oup = sys.argv[4]
	chr = sys.argv[5]

	f1 = open(borderfile,'r')
	border_list = []
	for line1 in f1:
		info1 = line1.strip().split('\t')
		# chr = info1[0]
		border = int(info1[1])
		if info1[0] == chr:
			border_list.append(border)
		# border_dict[chr].append(border)
	f1.close()

	f1 = open('chr1.size.txt','r')
	chr_l_dict={}
	for line1 in f1:
		info1 = line1.strip().split('\t')
		chr_l_dict[info1[0]] = int(info1[1])
	f1.close()

	f2 = open(samplefile,'r')
	sample_list = []
	alltad_dict = {}
	for line2 in f2:
		sample = line2.strip()
		alltad_dict[sample] = {}
		sample_list.append(sample)
		alltad_dict[sample][chr] = {}
		inp = ind+'/'+sample+'.chr1.tad.bed'
		f2 = open(inp,'r')
		for line2 in f2:
			info2 = line2.strip().split('\t')
			start = int(info2[1])
			end = int(info2[2])
			if info2[0] == chr:
				alltad_dict[sample][chr][start] = end
		f2.close()

	out = open(oup,'w')
	# f1 = open(borderfile,'r')
	# for line1 in f1:
	# 	info1 = line1.strip().split('\t')
	for left in sorted(border_list):
		# if info1[0] != chr:
		# 	continue
		# left = int(info1[1])
		a_list = []
		a_dict = {}
		for sample in sample_list:
			for start in sorted(alltad_dict[sample][chr].keys()):
				end = alltad_dict[sample][chr][start]
				# if start <= left and left < end:
				if left >= start and left < end:
					a_list.append(end)
					a_dict[end] = 0
					break

		if len(a_list) <= len(sample_list)/2:
			continue

		for i in sorted(a_dict.keys()):
			for j in sorted(a_list):
				if j >= i:
					a_dict[i] += 1

		b_list = []
		for i in sorted(a_dict.keys()):
			if a_dict[i] >= len(sample_list)/2:
				b_list.append(i)

		if len(b_list) == 0:
			continue
		else:
			if (b_list[-1] - left < 100000) or (chr_l_dict[chr] - left < 100000):
				continue
			if b_list[-1] != chr_l_dict[chr]:
				info_list = [chr,left,b_list[-1]-1]
			else:
				info_list =[chr,left,chr_l_dict[chr]]
			out.write('{0}\n'.format('\t'.join([str(i) for i in info_list])))
	f1.close()
	out.close()

if __name__ == '__main__':
	main()