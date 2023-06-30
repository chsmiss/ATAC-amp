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
class find_amp:
	def __init__(self,bam,discbk,interval,cov):
		self.cov = cov
		print(cov)
		self.discbk = discbk
		self.interval = interval
		self.bam = bam
		self.samfile = pysam.AlignmentFile(bam,'rb')
		self.chrom_dir = {}
		for ref in self.samfile.references:
			self.chrom_dir[ref] = []
			self.chrom_dir[ref].append(Interval(0, 0))

#Calculates the sequencing depth of the specified length region to the left of the locus, the default length is interval
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

	#As above, calculate the sequencing depth of the right-hand region
	def interval_cov1(self,samfile_name,chr_n,site,interval_size=None):
		reads=0
		interval_size = self.interval

		rsite = site+interval_size
		#print(chr_n,samfile_name.get_reference_length(chr_n),site)
		if rsite > samfile_name.get_reference_length(chr_n):
			return 0
	#	print(chr_n, samfile_name.get_reference_length(chr_n),site,end_site)
		for i,k in enumerate(samfile_name.count_coverage(chr_n,site,rsite)):
			reads=reads+sum(k)
		return(reads/interval_size)

	#Finding contiguous regions with sequencing depth
	def find_amplicon(self,sam,start_c,start_s,cov_t=5):
		interval=self.interval
		cov_t = self.cov
		l_site = start_s
		if l_site >= sam.get_reference_length(start_c):
			l_site = sam.get_reference_length(start_c) - 1

			#Continuous area on the left
		while self.interval_cov(sam,start_c,l_site)>cov_t or self.interval_cov(sam,start_c,l_site,interval_size=10*interval)>cov_t:
			l_site=l_site-interval
			if l_site <= 0:
				break
		r_site = start_s
		#Continuous area on the right
		while self.interval_cov1(sam,start_c,r_site)>cov_t or self.interval_cov1(sam,start_c,r_site,interval_size=10*interval)>cov_t:
			r_site = r_site+interval
			if r_site >= sam.get_reference_length(start_c):
				break
	#Return to left endpoint and right endpoint
		return(l_site,r_site)

	#Determine if the breakpoint is in the detected amplification
	def isinchrom_dir(self,chr,pos):
		chrom_dir = self.chrom_dir
		for line in chrom_dir[chr]:
			if pos in line:
	#If the amplification region is known, the Boolean value True and the amplification region information are returned.
				return (True,line)

		return (False,line)

	#Main functions, traversing breakpoint files for processing
	def process_coverage(self,bam=None,support_reads=5):
		bk_file = self.discbk
		chrom_dir = self.chrom_dir
		bam = self.bam
		samfile = self.samfile
		s_t = time.time()
		f_bkline = open(bk_file,'r')
		bklines = f_bkline.readlines()
		#Create a results file
		outfile_n = bk_file+'.amplicon'
		f_out = open(outfile_n,'w')
		#count3 counts how many rows have been read
		count3 = 0
		for line in bklines:
			count3 = count3 + 1
			if count3%1000 == 0:
				print(count3)
				print(time.time()-s_t)
			#Removal of breakpoint pairs that match to uncertain chromosomes
			if 'chrUn' in line or 'random' in line:
				continue

			bkline_list = line.split('\t')

			#Less than 5 breakpoint pairs with reads support do not
			if int(bkline_list[4]) <= support_reads:
				continue

			#Read the chromosomes located at both ends of the breakpoint pair
			s_chr = bkline_list[0]
			e_chr = bkline_list[2]

			#Read the interval between the ends of the breakpoint pair, temporarily using the midpoint of the interval to indicate the endpoint position
			start_site = int((int(bkline_list[1].strip('[').strip(']').split('.')[-1])+int(bkline_list[1].strip('[').strip(']').split('.')[0])+1)/2)
			end_site = int((int(bkline_list[3].strip('[').strip(']').split('.')[-1])+int(bkline_list[3].strip('[').strip(']').split('.')[0])+1)/2)

			#Do not have endpoints with coordinates less than 0
			if start_site <= 0 or end_site <= 0:
				continue

			#Call the isinchrom_dir function to determine if the left end is in the collected amplification region
			if self.isinchrom_dir(s_chr,start_site)[0]:
				#If in a known amplification region, take out the region location
				first_l,first_r = self.isinchrom_dir(s_chr,start_site)[1].lower_bound,self.isinchrom_dir(s_chr,start_site)[1].upper_bound

			else:
				#If not, call find_amplicon to search for the augmented region starting with the endpoint and add this region to the augmented region dictionary

				first_l,first_r = self.find_amplicon(samfile,s_chr,start_site)
				chrom_dir[s_chr].append(Interval(first_r,first_l))

			#As above, determine the right endpoint
			if self.isinchrom_dir(e_chr,end_site)[0]:
				second_l,second_r = self.isinchrom_dir(e_chr,end_site)[1].lower_bound,self.isinchrom_dir(e_chr,end_site)[1].upper_bound
			else:
				second_l,second_r = self.find_amplicon(samfile,e_chr,end_site)
				chrom_dir[e_chr].append(Interval(second_l,second_r))

			#Left and right endpoints are only retained if they lie in an interval greater than a certain size
			if (first_r - first_l) > 2000 and (second_r - second_l)>2000:
				out_line = 'first_bk_interval'+'\t'+str(first_l)+'\t'+str(first_r)+'\t'+'first_len'+'\t'+str(first_r - first_l)+'\t''second_bk_interval'+'\t'+str(second_l)+'\t'+str(second_r)+'\t'+'second_len'+'\t'+str(second_r - second_l)+'\t'+line
				f_out.write(out_line)


		return(outfile_n)

#fina_amp('COLO320HSR_rep2_atac_possorted_bam.bam',interval=1000).process_coverage()
##print(find_amplicon(samfile,s_chr,start_site))
	#	print(find_amplicon(samfile,s_chr,start_site)[1]-find_amplicon(samfile,s_chr,start_site)[0])
	#	print(find_amplicon(samfile,e_chr,end_site))
	#	print(find_amplicon(samfile,e_chr,end_site)[1]-find_amplicon(samfile,e_chr,end_site)[0])


