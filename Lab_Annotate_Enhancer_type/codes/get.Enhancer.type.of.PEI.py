import sys
import os
def main():

	inp = sys.argv[1]
	f1 = open(inp,'r')
	Ehancer_dict = {}
	for line1 in f1:
		info1 = line1.strip().split('\t')
		if info1[0] == 'REGION_ID':
			continue

		chrname = info1[1].replace('chr','')
		start = int(info1[2])
		end = int(info1[3])

		tname = info1[-1]
		try:
			Ehancer_dict[chrname].append([chrname,start,end,tname])
		except KeyError as reason:
			Ehancer_dict[chrname] = []
			Ehancer_dict[chrname].append([chrname,start,end,tname])
	f1.close()

	oup = sys.argv[3]
	out = open(oup,'w')

	inp = sys.argv[2]
	f1 = open(inp,'r')
	for line1 in f1:
		info1 = line1.strip().split('\t')
		if info1[0] == 'chrname':
			continue

		chrname = info1[0]
		Estart = int(info1[1])
		Eend = Estart + 5000

		RE='NA'
		SE='NA'
		for a_list in Ehancer_dict[chrname]:
			start = a_list[1]
			end = a_list[2]
			tname = a_list[3]

			ovlap = min(Eend,end) - max(start,Estart)
			if ovlap <= 0:
				continue

			if tname == '0':
				if ovlap >= 1:
					RE='RE'
					REinfo = '_'.join([str(i) for i in a_list])
			elif tname == '1':
				if ovlap >=2500:
					SE='SE'
					SEinfo = '_'.join([str(i) for i in a_list])
			else:
				print 'wrong'
				sys.exit(0)

		if SE =='SE':
			Ehancer_name = 'SE'
			posinfo = SEinfo
		elif RE == 'RE':
			Ehancer_name = 'RE'
			posinfo = REinfo
		else:
			Ehancer_name = 'inactive'
			posinfo = 'NA'

		info_list = info1+[posinfo,Ehancer_name]
		out.write('{0}\n'.format('\t'.join(map(str,info_list))))
	f1.close()
	out.close()


if __name__ == '__main__':
	main()
