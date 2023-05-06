import sys
import networkx as nx
import matplotlib.pyplot as plt
import pysam
import interval
import time
from multiprocessing import Pool
import concurrent.futures


s_t = time.time()
#打开不同区域的断点对
f_dif_bk = open(sys.argv[1],'r')
f_out = open(sys.argv[2],'w')
mode = sys.argv[3]
print('This run is in '+mode+'mode')
gtf_file = sys.argv[4]
threads = int(sys.argv[5])
def annotate():
    gtf = open(gtf_file, 'r')
    chromosomes = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    gtf_dict = {chromosome: [] for chromosome in chromosomes}
    for gtf_line in gtf.readlines():
        gtfs = gtf_line.split('\t')
        if gtfs[2] == 'transcript':
            if gtfs[0] in gtf_dict:
                gtf_dict[gtfs[0]].append({interval.Interval(int(gtfs[3]),
                                                            int(gtfs[4])):
                                              gtfs[8].split(';')[-2].split(
                                                  '\"')[1]})
    return (gtf_dict)


gene_annotation = annotate()


#print(gene_annotation)


def getbarcode(args):
    srr_l,bam_f = args
    bam = pysam.AlignmentFile(bam_f, 'rb')
    cblist = []

    for read in bam:
        if read.qname in srr_l:
            if read.has_tag('CB'):
                cblist.append(read.get_tag('CB'))
    bam.close()
    return list(set(cblist))

def split_list(lst, num_parts):
    # 计算每个子列表的长度
    n = len(lst)
    len_part = n // num_parts
    # 将列表分割成指定数量的子列表
    parts = [lst[i * len_part:(i + 1) * len_part] for i in
             range(num_parts)]
    # 如果有余数，则将余数部分添加到最后一个子列表中
    if n % num_parts != 0:
        parts[-1].extend(lst[num_parts * len_part:])
    return parts

def getbarcode_multiprocess(srr_name,bam_name,f_out,threads):
    srrs = split_list(srr_name, threads)
    pool = Pool(threads)
    #         for srr in srrs:
    #          results = pool.map(getbarcode,[(srr,bam_name)])
    results = pool.map(getbarcode, [(srr, bam_name) for srr in srrs])
    for r_list in results:
        f_out.write('\t'.join(r_list))


def dfs(start, curr, path, visited,G_dfs):
    """
    用于进行深度优先遍历的函数
    """
    visited.add(curr)
    path.append(curr)
    for neighbor in G_dfs.neighbors(curr):
        if neighbor == start and len(path) > 2:
            path.append(start)
            yield path.copy()  # 找到了一个环
            path.pop()
        elif neighbor not in visited:
            yield from dfs(start, neighbor, path, visited,G_dfs)
    path.pop()
    visited.remove(curr)

bam_list = []

if sys.argv[3] == 'sc':
        bam_name = sys.argv[1].strip('.bkline_dif_interval')+'.discordant.bam'
        #读取bam文件，要有索引，获取细胞barcode
        f_bam = pysam.AlignmentFile(bam_name,'rb')
        for line_bam in f_bam:
                bam_list.append(line_bam)
        f_bam.close()

print('fuction finished')
print(time.time()-s_t)
#所有区域id列表
seg_list = []
#断点对左右两个区域id的列表转成元组，再存在列表
edge_list = []
#node字典，key为区域id，值为染色体号，起始终止位置和长度
node_dict = {}
#edge字典，key为断点对的左右两个区域id，值为支持这个断点对的reads id
edge_dict = {}
#遍历断点对文件
for line in f_dif_bk.readlines():
        line_list = line.split('\t')
        seg_list = seg_list+line_list[:2]
        edge_list.append(tuple(line_list[:2]))
        node_dict[line_list[0]] = line_list[12] + '\t'+line_list[3] +'\t'+ line_list[4]+'\t'+line_list[6]
        node_dict[line_list[1]] = line_list[14] + '\t'+line_list[8] +'\t'+ line_list[9]+'\t'+line_list[11]
        edge_key = str(min(int(line_list[0]),int(line_list[1])))+'\t'+str(max(int(line_list[0]),int(line_list[1])))
        if edge_key in edge_dict:
              #  edge_dict[edge_key] = edge_dict[edge_key]  + '\t'.join(line_list[17::2])
                edge_dict[edge_key] = edge_dict[edge_key] + '\t'+'\t'.join(line_list[17::2])
        else:
                edge_dict[edge_key] = '\t'.join(line_list[17::2])
#print(node_dict)
#nodes = list(set(seg_list))
#print(nodes)

print('edge_dict established')
print(time.time()-s_t)


amp_gene = {}
#注释扩增区域
for key in node_dict.keys():
        amp_gene[key]=[]
        amp_line = node_dict[key].split('\t')
        if amp_line[0].startswith('chr'):
            chr_name = amp_line[0]
        else:
            chr_name = 'chr'+amp_line[0]
        if chr_name in gene_annotation:
                for gene_line in gene_annotation[chr_name]:
                        if interval.Interval(int(amp_line[1]),int(amp_line[
                            2])).overlaps(list(gene_line)[0]):
                                amp_gene[key].append(str(list(gene_line.values())[0]))
                amp_list = list(set(amp_gene[key]))
                amp_gene[key] = ','.join(amp_list)

print('annotate finished')
print(time.time()-s_t)

#print(amp_gene)

#创建空的网络图
G=nx.Graph()
#H = nx.path_graph(len(nodes))
#G.add_nodes_from(H)
#从元组列表中获得边和节点
G.add_edges_from(edge_list)


#判断是否是连通图
#print(nx.is_connected(G))




#输入一个图，输出每条边对应的支持reads id
def out_edge(subgraph):
        srr_num = 0
        srr_id = ''
        for lien in subgraph.edges():
                temp = str(min(int(lien[0]), int(lien[1]))) + '\t' + str(max(int(lien[0]), int(lien[1])))
                srr_id = srr_id+'\t'+edge_dict[temp]
                srr_num = srr_num + len(edge_dict[temp].split('\t'))

        return (srr_num,srr_id)

#输出连通分量的节点集合
def node_set(G_e):
        sum_len = 0
        for ii in G_e:

                sum_len = sum_len + int(node_dict[ii].split('\t')[-1])

        return (sum_len)


#输出图的连通分量数量
#print(nx.number_connected_components(G))

#输出图中的环
#print(nx.cycle_basis(G))
#print(nx.find_cycle(G))

#
largest = list(nx.connected_components(G))
sorted_subgraph = []
for line1 in largest:
        sort_i = out_edge(G.subgraph(line1))[0]
        sorted_subgraph.append([line1,sort_i])

#按支持reads数排序扩增区域
sorted_subgraph = sorted(sorted_subgraph, key=lambda x:x[1],reverse=True)


index_i = 0
count4 = 0
for line in sorted_subgraph:
        count4 = count4 +1
        print('processing subgraph:'+str(count4))
        index_c = 0
        index_i = index_i+1
        print(line[0],out_edge(G.subgraph(line[0]))[0])
        info = sys.argv[2]+'_'+'amplicon interval sets'+str(index_i)+': '
        amp_len = node_set(line[0])
        f_out.write('\n'+info+','.join(line[0])+'\t'+str(line[1])+'\t'+str(amp_len)+'\n')
        for line3 in line[0]:
                f_out.write(line3+'\t'+node_dict[line3]+'\t'+str(amp_gene[line3])+'\n')


        if mode == 'sc':
                # count5 = 0
                srr_name = out_edge(G.subgraph(line[0]))[1].split('\t')
                srrs = split_list(srr_name,threads)
                pool = Pool(threads)
                results = pool.map(getbarcode,[(srr, bam_name) for srr in
                                                srrs])
                for r_list in results:
                    f_out.write('\t'.join(r_list))

                f_out.write('\n')
                print(time.time()-s_t)

        circle = nx.cycle_basis(G.subgraph(line[0]))

        if circle != []:

                max_cycle = []

                for start_node in G.subgraph(line[0]).nodes():
                        for cycle in dfs(start_node, start_node, [], set(),G.subgraph(line[0])):
                                if len(cycle) > len(max_cycle):
                                        max_cycle = cycle

                if len(max_cycle) > 0:
                        print("max cycle is：", max_cycle)
                        f_out.write("max cycle："+','.join(max_cycle)+'\n')


                plt.clf ()
                plt.rcParams['figure.figsize'] = (500, 500)
                pos = nx.spring_layout(G.subgraph(line[0]))

                #设置出来的图的节点文字
                label = {}
                for line10 in (pos.keys()):
                        label[line10]=' '.join(node_dict[line10].split('\t')[:-1])


                nx.draw(G.subgraph(line[0]),pos,node_size=300)
                nx.draw_networkx_labels(G.subgraph(line[0]),pos,label, font_size=5, font_color='black')
         #       plt.gcf().subplots_adjust(left=0.3,top=0.91,bottom = 0.1)

                picture_name = info.strip(': ').replace(' ','_')+'.png'
                plt.savefig(picture_name,dpi=300,bbox_inches = 'tight')
                for line4 in circle:
                        index_c = index_c + 1
                        f_out.write('cycle '+ str(index_c)+': '+','.join(line4)+'\n')
                     #   print(circle)
                        for line5 in G.subgraph(line4).edges():
                                f_out.write(','.join(line5)+'\n')

# for lien in largest_connected_subgraph.edges():
#         temp = str(min(int(lien[0]),int(lien[1])))+' '+str(max(int(lien[0]),int(lien[1])))
#         print(temp)
#         print(node_dict[lien[0]],node_dict[lien[1]])
#         #print(edge_dict[temp])


#plt.rcParams['figure.figsize']= (6, 4)
#
# #nx.draw(G, with_labels=True)
# #nx.draw(largest_connected_subgraph,pos = nx.random_layout(largest_connected_subgraph),node_color = 'b',edge_color = 'r',with_labels = True,font_size =18,node_size =20)
#nx.draw(largest_connected_subgraph, with_labels=True)
# plt.show()
# plt.savefig('./generated_image.png')
