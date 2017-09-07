#!/usr/bin/env python
# pip install pillow
# pip install selenium
import sys
import numpy as np
import random
import pysam
import os
import argparse
from bokeh.plotting import figure, output_file, show, gridplot, save
from bokeh.layouts import column
from bokeh.resources import CDN
from bokeh.embed import autoload_static

parser=argparse.ArgumentParser()

parser.add_argument("-b",
                  dest="bams",
                  help="Bam file names",
                  nargs="+",
                  required=True)

parser.add_argument("-o",
                  dest="output_file",
                  help="Output file name",
                  #type=lambda e:file_choices(("csv","tab"),e)
                  required=True)

parser.add_argument("--embed",
                  dest="embedded_path",
                  help="Embedded path",
                  #type=lambda e:file_choices(("csv","tab"),e)
                  required=False)

parser.add_argument("-d","--dynamo-config",
                  dest="dynamo_config",
                  help="Config file to store embed info in DynamoDB table. If not present, no data sent to DynamoDB",
                  required=False)

parser.add_argument("-s",
                  dest="start",
                  help="Start range")

parser.add_argument("-e",
                  dest="end",
                  help="End range")

parser.add_argument("-c",
                  dest="chrom",
                  help="Chromosome range")

parser.add_argument("-w",
                  dest="window",
                  type=int,
                  default=1000,
                  help="Window, default(1000)")
                  
args = parser.parse_args()

vals = []
all_plot_pairs = []
all_plot_reads = []

for bam_file_name in args.bams:
    bam_file = pysam.AlignmentFile(bam_file_name, "rb")
    pairs = {}

    plot_reads = []

    for read in bam_file.fetch(args.chrom,
                               int(args.start) - args.window,
                               int(args.end) + args.window):
        if (read.is_paired):
            if (read.reference_end) :
                if read.query_name not in pairs:
                    pairs[read.query_name] = []
                pairs[read.query_name].append(read)
        else:
            if (read.reference_end) :
                plot_reads.append([read.reference_start, read.reference_end])


    plot_pairs = []
    
    for pair in pairs:
        if len(pairs[pair]) == 2:

            plot_pairs.append([[pairs[pair][0].reference_start,
                                pairs[pair][0].reference_end,
                                not(pairs[pair][0].is_reverse)],
                               [pairs[pair][1].reference_start,
                                pairs[pair][1].reference_end,
                                not(pairs[pair][1].is_reverse)]])

            vals.append(pairs[pair][0].reference_start)
            if pairs[pair][0].reference_end:
                vals.append(pairs[pair][0].reference_end)

            vals.append(pairs[pair][1].reference_start)
            if pairs[pair][1].reference_end:
                vals.append(pairs[pair][1].reference_end)


    plot_pairs.sort(key=lambda x: x[1][1] - x[0][0])
    #plot_pairs.sort(key=lambda x: x[0][0], reverse=True)
    all_plot_pairs.append(plot_pairs)
    all_plot_reads.append(plot_reads)



x_range=[int(args.start) - args.window, int(args.end) + args.window]
tools = "pan,box_zoom,reset"

plots = []


p = figure(plot_width=1000, plot_height=50, x_range=x_range, title=args.chrom + ':' + args.start + '-' + args.end, tools=[])
p.title.text_font_style = "normal"
p.multi_line(xs=[[int(args.start),int(args.end)]], ys=[[0,0]], line_width=10, color=['black'])
p.xaxis.visible = False
p.yaxis.visible = False
p.xgrid.grid_line_color = None
p.ygrid.grid_line_color = None
plots.append(p)


p_i = 0
for plot_pairs in all_plot_pairs:

    pair_ends_x = []
    pair_ends_y = []
    links_x = []
    links_y = []
    colors_ends = []
    colors_links = []

    color_map = {
        (True, False): 'black',
        (False, True): 'red',
        (False, False): 'blue',
        (True, True): 'green',
    }

    for p in plot_pairs:
        insert_size = p[1][1] - p[0][0]
        new_y = [insert_size,insert_size]
        new_end = p[0][0:2]
        pair_ends_x.append(new_end)
        pair_ends_y.append(new_y)
        colors_ends.append(color_map[(p[0][2], p[1][2])])
        new_end = p[1][0:2]
        pair_ends_x.append(new_end)
        pair_ends_y.append(new_y)
        colors_ends.append(color_map[(p[0][2], p[1][2])])
        
        new_link=[p[0][1],p[1][0]]
        links_x.append(new_link)
        links_y.append(new_y)
        colors_links.append(color_map[(p[0][2], p[1][2])])


    p = figure(plot_width=1000, plot_height=200, x_range=x_range, title=os.path.basename(args.bams[p_i]), tools=tools)
    p.title.text_font_style = "normal"
    p.multi_line(xs=pair_ends_x, ys=pair_ends_y, line_width=4, color=colors_ends, alpha=0.25)
    p.multi_line(xs=links_x, ys=links_y, line_width=1, color=colors_links, alpha=0.25)
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.yaxis.axis_label = "Insert size"
    p.yaxis.axis_label_text_font_style= 'normal'
    if p_i + 1 < len(all_plot_pairs):
        p.xaxis.visible = False
    else:
        p.xaxis.axis_label = "Position on chromosome " + args.chrom
        p.xaxis.axis_label_text_font_style= 'normal'

    plots.append(p)
    p_i += 1

p = column(*plots)

if args.embedded_path:
    js, tag = autoload_static(p, CDN, args.embedded_path + '/' + args.output_file)
    f = open( args.embedded_path + '/' + args.output_file , 'w')
    f.write(js)
    f.close()
    print (tag)
        

    if args.dynamo_config:
        import boto3
        from boto3.s3.transfer import S3Transfer
        from boto3.dynamodb.conditions import Key, Attr
        from botocore.exceptions import ClientError
        import json
        import ntpath
        with open(args.dynamo_config,'r') as config_file:
            config_data = json.load(config_file)
        key = config_data['folderName'] + '/' + args.output_file
        client = boto3.client('s3')
        transfer = S3Transfer(client)
        transfer.upload_file(
                f.name,
                config_data['bucketName'],
                key,
                extra_args={'ACL': 'public-read'})
        file_url = '%s/%s/%s' % (client.meta.endpoint_url, config_data['bucketName'], key)
        
        script_fields = tag.strip().split('\n')
        script_fields[1] = 'src="' + file_url + '"'
        script = " ".join(script_fields)

        dynamodb = boto3.resource('dynamodb', 
                region_name=config_data['region'], 
                endpoint_url=config_data['dynamoEndpoint'])
        js_info_table = dynamodb.Table(config_data['dynamoTable'])
        try:
            response = js_info_table.put_item(
                Item = {
                    'id' : ("__").join([ntpath.basename(x) for x in args.bams]) + \
                            "__" + args.chrom + ':' + args.start + '-' + args.end,
                    'chr' : args.chrom,
                    'start' : args.start,
                    'end' : args.end,
                    'bams' : [ntpath.basename(x) for x in args.bams],
                    'script' : script
                })
        except ClientError as e:
            print (e.response['Error']['Message'])

else:
    if args.output_file.split('.')[-1] == "html":
        output_file(args.output_file)
        save(p)
    elif args.output_file.split('.')[-1] == "png":
        from bokeh.io import export_png
        export_png(p, filename=args.output_file)
