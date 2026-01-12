---
title: 构建 SAW 分析参考索引
date: 2026-01-12
updated: 2026-01-12
description: SAW-ST-V8 make reference.
top_img: https://www.vangoghgallery.com/img/starry_night_full.jpg
cover: https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/29/1/10.1093_bioinformatics_bts635/3/m_bioinformatics_29_1_15_f2.jpeg?Expires=1771021321&Signature=fqt36-Ujsft5yhbtrBM0IVIOgNYIJdnH2nPfHkwb~JZD1RMAQ-Khdm6Qkf6X2m~MbpxQ4D~otvV37upgiKEy7dCA5r5wky4G9PL84IZxKK08t6oBzY1w-LB-gjzFS1gg4RMojlD~WYIgMjkBFkcAWdrrY3ZcoL1NwCkURaodOQbVmyn80P4SD89j3F9ylCvcKp6Dgn3Nj1yLcqgvxu9s-Rgs4cvbsJRWHSxHqpo0JmjtDDyFkzvDXHpwxsMaF60BvwMpvw5hMGBH0D2tBWm7Z~5jOkdpv0QySk~WoALpsNsgJndzXFF3jftYEyqN4h7fyDBN6DJREleM5lyhTwM2ow__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA
categories:
  - Workflow
tags: 
  - BGI
  - Genomics
  - Spatial Transcriptomics
comments: true
---

在运行 SAW (Stereo-seq Analysis Workflow) count 分析前，需要先生成参考基因组索引 (index).这一步可以由 `SAW makeRef` 流程分析实现.我们需要提供的输入文件为参考基因组序列 `reference.fna` 和对应的注释文件 `reference.gff/gtf`.

## miniprot 基因注释

由于原先的参考基因组 `Larus_argentatus.fna` 并未注释出关键的三个视蛋白 (Opsins) 基因 *OPN1SW*、*OPN2SW* 与 *OPN1LW*，因此需要重新对参考基因组进行注释，选用的软件为 **miniprot**，具体注释方法在此不展开论述.最终只注释出位于 Chr1 的 *OPN2SW* 基因，需要添加的另外两个基因序列所在染色体分别为 [CAYQHM010000238.1.fna](https://www.ncbi.nlm.nih.gov/nuccore/CAYQHM010000238.1) 和 [OZ207386.1.fna](https://www.ncbi.nlm.nih.gov/nuccore/OZ207386.1/)，**miniprot** 注释输出的 `GFF3` 格式文件分别为 `CAYQHM010000238.1.gff` 和 `OZ207386.1.gff`.

## 注释文件格式转换

原先的基因组注释文件为 `.gtf` 格式，因此需将 `CAYQHM010000238.1.gff`和 `OZ207386.1.gff` 转换为 GTF 格式.所用脚本如下：

``` python
#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys

def parse_attributes(attr_str):
    """解析 GFF 属性列为字典"""
    attrs = {}
    for item in attr_str.strip(';').split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            attrs[k] = v
        elif ' ' in item: # 处理某些类 GTF 格式
            parts = item.strip().split(' ', 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip('"')
    return attrs

def convert_gff_to_gtf(input_file):
    genes = {} # 存储基因信息
    transcripts = {} # 存储转录本信息 {rna_id: [feature_lines]}

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
                
            feature_type = cols[2]
            attributes = parse_attributes(cols[8])
            
            if feature_type == 'mRNA':
                rna_id = attributes.get('ID')
                # 如果没有 Parent，则用 rna_id 作为 gene_id
                gene_id = attributes.get('Parent', rna_id)
                transcripts[rna_id] = {
                    'info': cols,
                    'gene_id': gene_id,
                    'features': []
                }
            elif feature_type in ['CDS', 'exon', 'five_prime_utr', 'three_prime_utr']:
                rna_id = attributes.get('Parent')
                if rna_id in transcripts:
                    transcripts[rna_id]['features'].append(cols)

    # 开始输出
    for rna_id, data in transcripts.items():
        t_info = data['info']
        gene_id = data['gene_id']
        chrom, source, _, start, end, score, strand, frame = t_info[:8]

        # 1. 输出 gene 行
        gene_attr = f'gene_id "{gene_id}"; gene_name "{gene_id}"; gene_biotype "protein-coding";'
        print(f"{chrom}\t{source}\tgene\t{start}\t{end}\t{score}\t{strand}\t.\t{gene_attr}")

        # 2. 输出 transcript 行
        tran_attr = f'gene_id "{gene_id}"; transcript_id "{rna_id}"; gene_name "{gene_id}"; gene_biotype "protein-coding";'
        print(f"{chrom}\t{source}\ttranscript\t{start}\t{end}\t{score}\t{strand}\t.\t{tran_attr}")

        # 3. 输出子特征 (CDS/Exon)
        # 如果 GFF 里只有 CDS 没有 exon，GTF 必须补全 exon
        has_exon = any(f[2] == 'exon' for f in data['features'])
        
        exon_count = 0
        cds_count = 0
        
        # 排序特征坐标 (正链升序，负链降序)
        sorted_features = sorted(data['features'], key=lambda x: int(x[3]))
        
        for feat in sorted_features:
            f_type = feat[2]
            f_start, f_end, f_score, f_strand, f_frame = feat[3], feat[4], feat[5], feat[6], feat[7]
            
            # 生成通用属性前缀
            base_attr = f'gene_id "{gene_id}"; transcript_id "{rna_id}"; gene_name "{gene_id}";'

            # 处理 CDS 和 自动生成 exon
            if f_type == 'CDS':
                cds_count += 1
                # 如果原始数据没写 exon，把 CDS 当作 exon 输出一次
                if not has_exon:
                    exon_count += 1
                    ex_attr = base_attr + f' exon_id "{rna_id}.exon.{exon_count}";'
                    print(f"{chrom}\t{source}\texon\t{f_start}\t{f_end}\t{f_score}\t{f_strand}\t.\t{ex_attr}")
                
                # 输出 CDS
                c_attr = base_attr + f' exon_id "{rna_id}.cds.{cds_count}";'
                print(f"{chrom}\t{source}\tCDS\t{f_start}\t{f_end}\t{f_score}\t{f_strand}\t{f_frame}\t{c_attr}")
            
            elif f_type == 'exon':
                exon_count += 1
                ex_attr = base_attr + f' exon_id "{rna_id}.exon.{exon_count}";'
                print(f"{chrom}\t{source}\texon\t{f_start}\t{f_end}\t{f_score}\t{f_strand}\t.\t{ex_attr}")
            
            else:
                # 处理 UTR 等其他特征
                print(f"{chrom}\t{source}\t{f_type}\t{f_start}\t{f_end}\t{f_score}\t{f_strand}\t{f_frame}\t{base_attr}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python gff2gtf.py input.gff")
    else:
        convert_gff_to_gtf(sys.argv[1])
```

```bash
# Example
python gff2gtf.py CAYQHM010000238.1.gff > OPN1LW.gtf
```

格式转换成功后，建议将 `gene_id`、`gene_name`、`transcript_id` 等改为基因名（**miniprot** 自动输出的标识为 `MPXXXXXX`），其余命名格式与原先 GTF 文件相同，如一个转录本（transcript）有多个 exon/CDS 时，用 exon.1/.2/... 和 cds.1/.2/... 加以区分.

## 添加序列与注释文件

在 **Linux** 终端，用如下命令将基因序列与注释信息添加进原有的基因组序列与注释文件中：

``` bash
# Example
# Chr1 放最前面
cat OPN2SW.gtf Larus_argentatus.gtf OPN1LW.gtf OPN1SW.gtf > Larus_argentatus_v1.gtf

# 删除 '>' 后多余信息
cat Larus_argentatus.fna CAYQHM010000238.1.fna OZ207386.1.fna | sed -e "s/>OZ207386.1 Larus argentatus genome assembly, chromosome: 1/>OZ207386.1/g" -e "s/>CAYQHM010000238.1 Larus argentatus genome assembly, contig: HAP1_SCAFFOLD_368, whole genome shotgun sequence/>CAYQHM010000238.1/g" > Larus_argentatus_v1.fna
```

## 运行 SAW makeRef

由于需进行转录组比对，所选用的 Mode 为 [**STAR**](https://academic.oup.com/bioinformatics/article/29/1/15/272537).若在终端运行，则命令如下：

``` bash
# 'transcript' 为存放结果文件创建的文件夹
saw makeRef \
    --mode=STAR \
    --fasta=Larus_argentatus_v1.fna \
    --gtf=Larus_argentatus_v1.gtf \
    --genome=./transcript
```

关于 `SAW makeRef` 更多用法，可见参考文献 1.

## ⚠️ 注意事项

1. `fna` 文件中的染色体/Scaffold 等标识符必须唯一，不能重复，且需和 `gtf` 文件中的对应.
2. 在拼接完基因组序列后，建议检查文件大小/行数，防止有碱基序列丢失.可用 `wc -l` 打印行数，或者 `grep -v '>' filename | tr -d '\n' | wc -w` 精确统计碱基个数（去除 Header 信息）.
3. 注释文件中的起止坐标不能超过所在染色体/Scaffold 的总长度，否则构建索引会报错.

## 参考文献

1. [构建索引文件](https://www.stomics.tech/service/saw_8_2/docs/shi-yong-jiao-cheng/bi-dui-wen-jian/gou-jian-suo-yin-wen-jian.html)
2. [STAR: ultrafast universal RNA-seq aligner](https://academic.oup.com/bioinformatics/article/29/1/15/272537)