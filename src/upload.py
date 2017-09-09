# pip install pillow
# pip install selenium
import sys
import os
import argparse
import boto3
from boto3.s3.transfer import S3Transfer
from boto3.dynamodb.conditions import Key, Attr
from botocore.exceptions import ClientError
import json

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--directory",
                  dest="directory",
                  help="Directory to upload",
                  required=True)

parser.add_argument("-c","--config",
                  dest="config",
                  help="Config file to store embed info in DynamoDB table and S3 bucket",
                  required=True)
args = parser.parse_args()

with open(args.config,'r') as config_file:
    config_data = json.load(config_file)
client = boto3.client('s3')
transfer = S3Transfer(client)

sv_args = {}
dir_files = os.listdir(args.directory)
for filename in dir_files:
    basename, ext = os.path.splitext(filename)
    if basename + ".js" not in dir_files or basename + ".args" not in dir_files:
        print ("Warning: mismatched file '" + filename + "' found in '" + args.directory + "'")
        continue
    if basename not in sv_args.keys(): sv_args[basename] = {}
    if ext  == '.js':
        to_store = args.directory + '/' + filename
        key = config_data['folderName'] + '/' + filename
        transfer.upload_file(
                to_store,
                config_data['bucketName'],
                key,
                extra_args={'ACL': 'public-read'})
        sv_args[basename]['file_url'] = str('%s/%s/%s' % (client.meta.endpoint_url, config_data['bucketName'], key))
    elif ext == ".args":
        with open(args.directory + '/' +filename, 'r') as arg_file:
            keys = arg_file.readline().strip().replace("#", '').split('\t')
            values = arg_file.readline().strip().split('\t')
            sv_args[basename] = dict(zip(keys, values))
            sv_args[basename]['bams'] = sv_args[basename]['bams'].split(',')
            script_fields = sv_args[basename]['script'].split()
            script_info = {}
            for i in range(len(script_fields)):
                if "=" in script_fields[i]:
                    pair = script_fields[i].split("=")
                    pair[1] = pair[1].strip('"')
                    script_info[pair[0]] = pair[1]
            sv_args[basename]['script'] = script_info

dynamodb = boto3.resource('dynamodb', 
        region_name=config_data['region'], 
        endpoint_url=config_data['dynamoEndpoint'])
js_info_table = dynamodb.Table(config_data['dynamoTable'])

with js_info_table.batch_writer() as batch:
    for key in sv_args:
        sv_args[key]['script']['src'] = sv_args[key]['file_url']
        batch.put_item(
            Item = {
                'id': sv_args[key]['file_url'],
                'chr' : sv_args[key]['chrom'],
                'start' : sv_args[key]['start'],
                'end' : sv_args[key]['end'],
                'bams' : sv_args[key]['bams'],
                'script': sv_args[key]['script'],
            }
        )
