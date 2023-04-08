import pysam
from interval import Interval
import time

#samfile = pysam.AlignmentFile('../SRR8236758.sorted.rmdup.mapq20.bam','rb')
#f_bkline = open(sys.argv[1],'r')
#interval=int(sys.argv[2])
#bklines = f_bkline.readlines()

#s_chr = sys.argv[1]
#give_site = int(sys.argv[2])
#interval = int(sys.argv[3])
class fina_amp:
	def __init__(self,bam,discbk,interval):
		self.discbk = discbk
		self.interval = interval
		self.bam = bam
		self.samfile = pysam.AlignmentFile(bam,'rb')
		self.chrom_dir = {}
		for ref in self.samfile.references:
			self.chrom_dir[ref] = []
			self.chrom_dir[ref].append(Interval(0, 0))

#计算位点左侧指定长度区域测序深度，默认长度是interval
	def interval_cov(self,samfile_name,chr_n,site,interval_size=None):
		reads=0
		interval_size=self.interval
		llimit = site-interval_size
		if llimit<=0:
			return 0
	#	print(chr_n,samfile_name.get_reference_length(chr_n),llimit,site)
		for i,k in enumerate(samfile_name.count_coverage(chr_n,llimit,site)):

			reads=reads+sum(k)
		return(reads/interval_size)

	#同上，计算右侧区域测序深度
	def interval_cov1(self,samfile_name,chr_n,site,interval_size=None):
		reads=0
		interval_size = self.interval

		end_site = site+interval_size
		#print(chr_n,samfile_name.get_reference_length(chr_n),site)
		if end_site > samfile_name.get_reference_length(chr_n):
			return 0
	#	print(chr_n, samfile_name.get_reference_length(chr_n),site,end_site)
		for i,k in enumerate(samfile_name.count_coverage(chr_n,site,end_site)):
			reads=reads+sum(k)
		return(reads/interval_size)

	#找连续有测序深度的区域
	def find_amplicon(self,sam,start_c,start_s,cov_t=5):
		interval=self.interval
		l_site = start_s
		if l_site >= sam.get_reference_length(start_c):
			l_site = sam.get_reference_length(start_c) - 1
		while self.interval_cov(sam,start_c,l_site)>cov_t or self.interval_cov(sam,start_c,l_site,interval_size=10*interval)>cov_t:

			l_site=l_site-interval
			if l_site <= 0:
				break
		r_site = start_s
		while self.interval_cov1(sam,start_c,r_site)>cov_t or self.interval_cov1(sam,start_c,r_site,interval_size=10*interval)>cov_t:
			r_site = r_site+interval
			if r_site >= sam.get_reference_length(start_c):
				break
	#返回左端点和右端点
		return(l_site,r_site)

	#判断断点在不在已检测扩增中
	def isinchrom_dir(self,chr,pos):
		chrom_dir = self.chrom_dir
		for line in chrom_dir[chr]:
			if pos in line:
	#在已知扩增区域的话，返回布尔值True和扩增区域信息
				return (True,line)

		return (False,line)

	#主要功能，遍历断点文件进行处理
	def process_coverage(self,bam=None,support_reads=5):
		bk_file = self.discbk
		chrom_dir = self.chrom_dir
		bam = self.bam
		samfile = self.samfile
		s_t = time.time()
		f_bkline = open(bk_file,'r')
		bklines = f_bkline.readlines()
		#创建一个结果文件
		outfile_n = bk_file+'.amplicon'
		f_out = open(outfile_n,'w')
		#count3计数读取了多少行
		count3 = 0
		for line in bklines:
			count3 = count3 + 1
			if count3%1000 == 0:
				print(count3)
				print(time.time()-s_t)
			#去除比对到不确定染色体的断点对
			if 'chrUn' in line or 'random' in line:
				continue

			bkline_list = line.split('\t')

			#小于3个支持reads的断点对不要，也就是过滤多数只有一对，2条reads支持的断点
			if int(bkline_list[4]) <= support_reads:
				continue

			#读取断点对的两端位于的染色体
			s_chr = bkline_list[0]
			e_chr = bkline_list[2]

			#读取断点对两端区间，暂时用区间中点表示端点位置
			start_site = int((int(bkline_list[1].strip('[').strip(']').split('.')[-1])+int(bkline_list[1].strip('[').strip(']').split('.')[0])+1)/2)
			end_site = int((int(bkline_list[3].strip('[').strip(']').split('.')[-1])+int(bkline_list[3].strip('[').strip(']').split('.')[0])+1)/2)

			#端点坐标小于0的不要
			if start_site <= 0 or end_site <= 0:
				continue

			#调用isinchrom_dir函数，判断左端是否在已收集扩增区域
			if self.isinchrom_dir(s_chr,start_site)[0]:
				#如果在已知扩增区域，把区域位置拿出来
				first_l,first_r = self.isinchrom_dir(s_chr,start_site)[1].lower_bound,self.isinchrom_dir(s_chr,start_site)[1].upper_bound

			else:
				#如果不在，调用find_amplicon以端点为起始搜索扩增区域，并把这个区域加入扩增区域字典

				first_l,first_r = self.find_amplicon(samfile,s_chr,start_site)
				chrom_dir[s_chr].append(Interval(first_r,first_l))

			#同上，判断右端点
			if self.isinchrom_dir(e_chr,end_site)[0]:
				second_l,second_r = self.isinchrom_dir(e_chr,end_site)[1].lower_bound,self.isinchrom_dir(e_chr,end_site)[1].upper_bound
			else:
				second_l,second_r = self.find_amplicon(samfile,e_chr,end_site)
				chrom_dir[e_chr].append(Interval(second_l,second_r))

			#左右端点都位于大于某个大小的区间才保留
			if (first_r - first_l) > 2000 and (second_r - second_l)>2000:
				out_line = 'first_bk_interval'+'\t'+str(first_l)+'\t'+str(first_r)+'\t'+'first_len'+'\t'+str(first_r - first_l)+'\t''second_bk_interval'+'\t'+str(second_l)+'\t'+str(second_r)+'\t'+'second_len'+'\t'+str(second_r - second_l)+'\t'+line
				f_out.write(out_line)
			
		return(outfile_n)

#fina_amp('COLO320HSR_rep2_atac_possorted_bam.bam',interval=1000).process_coverage()
##print(find_amplicon(samfile,s_chr,start_site))
	#	print(find_amplicon(samfile,s_chr,start_site)[1]-find_amplicon(samfile,s_chr,start_site)[0])
	#	print(find_amplicon(samfile,e_chr,end_site))
	#	print(find_amplicon(samfile,e_chr,end_site)[1]-find_amplicon(samfile,e_chr,end_site)[0])


