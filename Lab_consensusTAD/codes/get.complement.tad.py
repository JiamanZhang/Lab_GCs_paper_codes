import sys
import os
def main():
	IStad = sys.argv[1]
	fai = sys.argv[2]
	oup = sys.argv[3]

	f2 = open(fai,'r')
	chrl_dict = {}
	for line2 in f2:
		info2 = line2.strip().split('\t')
		chrl_dict[info2[0]] = int(info2[1])
	f2.close()

	f1 = open(IStad,'r')
	out = open(oup,'w')
	for line1 in f1:
		info1 = line1.strip().split('\t')
		chr = info1[0]
		start = int(info1[1])
		end = int(info1[2])

		if start%20000 != 0:
			start += 1

		if end == chrl_dict[chr]:
			end = chrl_dict[chr]
		elif end%20000 != 0:
			end +=1
		out.write('{0}\t{1}\t{2}\n'.format(chr,start,end))
	f1.close()
	out.close()

if __name__ == '__main__':
	main()
