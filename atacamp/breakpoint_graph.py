#import networks as nx
import sys
from interval import Interval

#打开bkamplicon文件
f_input = open(sys.argv[1],'r')
f_out = open(sys.argv[2],'w')
f_out1 = open(sys.argv[3],'w')
input_file = f_input.readlines()
i=0
chr_dir = {}
dif_bk = []
ecdna_line = {}
del_key = []


def combine_segments(c_dict,bkline,ii):
	# line_list【10】是左端染色体名，判断这个染色体有没有在字典中，这个字典是以染色体作为key，其上的扩增区域信息构建的字典作为key的value
	if bkline[10] in c_dict.keys():
		index1 = 1
		#遍历字典
		for index,interv in c_dict[bkline[10]].items():
			if Interval(bkline[2],bkline[1]).overlaps(interv):
				c_dict[bkline[10]][index] = Interval(bkline[2],bkline[1]).join(interv)
				#如果这个区域和已有扩增区域重合，则合并，id和已有扩增区域相同
				lbk_id = index
				#为跳出这个循环，在找到已存在的扩增时，不执行下面的append语句
				index1 = 0
		if index1 == 1:
			#如果在该染色体上这个部分还没有已知扩增，则新建一个二级字典，扩增id为对应行号
			lbk_id = ii
			c_dict[bkline[10]][str(i)]=Interval(bkline[2],bkline[1])

	#不在字典中，说明是第一次遇到这个染色体的扩增片段，下面是把这个区域加进字典
	else:
		lbk_id = ii
		#构建一个字典中的字典
		c_dict[bkline[10]]={}
		#二级字典的key是首次出现的行号（一行加两次，左右两个端点），作为这个扩增区域的id，value是扩增区域位置区间
		c_dict[bkline[10]][str(ii)]=Interval(bkline[2],bkline[1])

	#这里加一次是为了让一行有两个号，左右两个端点
	ii=ii+1

	# 检测右端点
	if bkline[12] in c_dict.keys():
		index2 = 1
		for index, interv in c_dict[bkline[12]].items():
			if Interval(bkline[6], bkline[7]).overlaps(interv):
				rbk_id = index
				c_dict[bkline[12]][index] = Interval(bkline[6], bkline[7]).join(interv)
				index2 = 0
				
		if index2 == 1:
			rbk_id = ii
			c_dict[bkline[12]][str(ii)] = Interval(bkline[6], bkline[7])
	else:

		rbk_id = ii
		c_dict[bkline[12]] = {}
		c_dict[bkline[12]][str(ii)] = Interval(bkline[6], bkline[7])


	return (c_dict,ii,lbk_id,rbk_id)

#	line_out = str(lbk_id) + '\t' + str(rbk_id) + '\t' + line

#遍历断点区域对文件
for line in input_file:
	i=i+1
	line_list = line.split('\t')
	tuple_r = combine_segments(chr_dir,line_list,i)

	chr_dir = tuple_r[0]
	i = tuple_r[1]
	lbk_id = tuple_r[2]
	rbk_id = tuple_r[3]


	line_out = str(lbk_id) + '\t' +str(rbk_id) + '\t' + line

	#筛选断点两个端点不在一个扩增区域的
	if int(lbk_id) != int(rbk_id):
	#	print(line_list[10], chr_dir[line_list[10]][str(lbk_id)],line_list[12],str(rbk_id))
		print(line_list[10],chr_dir[line_list[10]][str(lbk_id)],line_list[12],chr_dir[line_list[12]][str(rbk_id)],line_list[16])
		line_num = int((i)/2)
		print(line_num)

		f_out1.write(line_out)
	f_out.write(line_out)










