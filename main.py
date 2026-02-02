#!/usr/bin/env python3
import sys
import math
import json


def parser_snp_info(data_file):
    snp_info_dict = {}
    with open(data_file) as F:
        F.readline()
        for line in F:
            rsid, effect_allele, maf, beta = line.strip().split(",")
            beta = float(beta)
            maf = float(maf)
            snp_info_dict[rsid] = [effect_allele, maf, beta]
    return snp_info_dict

# 从stdin读取输入数据
body = sys.stdin.read()
try:
    inputs = json.loads(body)['inputs']
    # 读取自己收集好的snp 位点信息
    snp_info_dict = parser_snp_info("snp_info.txt")
    
    # 开始计算输入用户的prs值和人群平均prs值
    prs_value = 0
    population_prs_value = 0
    #prs_value = 1
    #population_prs_value = 1

    for rsid, (effect_allele, maf, beta) in snp_info_dict.items():
        genotype = inputs.get(rsid.upper(), "--")
        # 不包含微基因芯片或者微基因芯片检出结果为no call的位点的值用人群频率来折合代替:
        if genotype in ("--", "__"):
            prs_value += 2 * maf * beta
            #prs_value *= math.pow(_or, 2*maf) 
        else:
            prs_value += genotype.count(effect_allele) * beta
            #prs_value *= math.pow(_or, genotype.count(effect_allele)) 

        population_prs_value += 2 * maf * beta
        #population_prs_value *= math.pow(_or, 2*maf)
    
     
    # 将累加的beta值转为OR值,如果本身就是or值的累乘得到的结果，则不需再转换。
    prs_value = math.exp(prs_value)
    population_prs_value = math.exp(population_prs_value)
    
    #输出两者的比值
    result = f'您患XX疾病的风险是人群平均值的{prs_value/population_prs_value:.2f}倍数'

    # 输出给用户的结果只需要通过 print 输出即可，print只可调用一次
    print(result)
except Exception as e:
    # 错误信息需要被从 stderr 中输出，否则会作为正常结果输出
    sys.stderr.write('这是一段错误信息，可能对用户可见，建议友好输出')
    print(e)