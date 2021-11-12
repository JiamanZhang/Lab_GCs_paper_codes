import sys
import os
def main():

	inp = sys.argv[1]
	f1 = open(inp,'r')
	CTCF_list = []
	for line1 in f1:
		info1 = line1.strip().split('\t')
		if len(info1) != 10 or info1[0] == 'motif_id':
			continue
		chrname = info1[2].replace('chr','')
		start = int(info1[3])
		end = int(info1[4])
		strand = info1[5]

		CTCF_list.append([chrname,start,end,strand])
	f1.close()

	oup = sys.argv[3]
	out = open(oup,'w')

	inp = sys.argv[2]
	f1 = open(inp,'r')
	for line1 in f1:
		info1 = line1.strip().split('\t')

		first_list = []
		start1 = int(info1[1])-2500
		end1 = int(info1[1])+2500

		second_list = []
		start2 = int(info1[2])-2500
		end2 = int(info1[2])+2500

		a_list = []
		b_list = []
		for pos_list in CTCF_list:
			CTCFstart = pos_list[1]
			CTCFend = pos_list[2]
			CTCFstrand = pos_list[3]

			if CTCFstart >= start1 and CTCFend <= end1:
				first_list.append(CTCFstrand)
				if CTCFstrand == '+':
					a_list.append(pos_list)

			if CTCFstart >= start2 and CTCFend <= end2:
				second_list.append(CTCFstrand)
				if CTCFstrand == '-':
					b_list.append(pos_list)

		if '+' in first_list and '-' in second_list:
			info_list = info1[:4]
			out.write('{0}\n'.format('\t'.join(map(str,info_list))))
	f1.close()
	out.close()






if __name__ == '__main__':
	main()
