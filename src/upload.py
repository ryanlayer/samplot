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
client = boto3.client('s3',
    aws_access_key_id=config_data['accessKey'], 
    aws_secret_access_key=config_data['secretAccessKey'])
transfer = S3Transfer(client)

sv_args = {}
dir_files = {}
for filename in os.listdir(args.directory):
    raw_filename, ext = os.path.splitext(filename)
    if raw_filename not in sv_args: sv_args[raw_filename] = {}
    if raw_filename not in dir_files: dir_files[raw_filename] = {}
    if ext == ".args":
        dir_files[raw_filename]['args'] = filename
    else:
        dir_files[raw_filename]['img'] = filename
for raw_filename in dir_files:
    # filter for unmatched files
    if not dir_files[raw_filename]['args'] and dir_files[raw_filename]['img']:
        print ("Warning: mismatched file with prefix'" + raw_filename + "' found in '" + args.directory + "'")
        continue

    # upload images or js to S3
    to_store = args.directory + '/' + dir_files[raw_filename]['img']
    key = config_data['folderName'] + '/' + dir_files[raw_filename]['img']
    transfer.upload_file(
            to_store,
            config_data['bucketName'],
            key,
            extra_args={'ACL': 'public-read'})
    sv_args[raw_filename]['file_url'] = str('%s/%s/%s' % (client.meta.endpoint_url, config_data['bucketName'], key))

    # upload entries for image use to DynamoDB
    with open(args.directory + '/' + dir_files[raw_filename]['args'], 'r') as arg_file:
        keys = arg_file.readline().strip().replace("#", '').split('\t')
        values = arg_file.readline().strip().split('\t')
        if "file_url" in sv_args[raw_filename]:
            keys.append("file_url")
            values.append(sv_args[raw_filename]['file_url'])
        sv_args[raw_filename] = dict(zip(keys, values))
        sv_args[raw_filename]['bams'] = sv_args[raw_filename]['bams'].split(',')
        if (len(sv_args[raw_filename]['script']) > 1):
            script_fields = sv_args[raw_filename]['script'].split()
            script_info = {}
            for i in range(len(script_fields)):
                if "=" in script_fields[i]:
                    pair = script_fields[i].split("=")
                    pair[1] = pair[1].strip('"')
                    script_info[pair[0]] = pair[1]
                sv_args[raw_filename]['script'] = script_info

dynamodb = boto3.resource('dynamodb', 
        region_name=config_data['region'], 
        endpoint_url=config_data['dynamoEndpoint'],
        aws_access_key_id=config_data['accessKey'], 
        aws_secret_access_key=config_data['secretAccessKey']
        )
js_info_table = dynamodb.Table(config_data['dynamoTable'])

with js_info_table.batch_writer() as batch:
    for key in sv_args:
        if type(sv_args[key]['script']) == str:
            sv_args[key]['script'] = sv_args[key]['file_url']
        else:
            sv_args[key]['script']['src'] = sv_args[key]['file_url']
        batch.put_item(
            Item = {
                'id': sv_args[key]['file_url'],
                'chr' : sv_args[key]['chrom'],
                'start' : sv_args[key]['start'],
                'end' : sv_args[key]['end'],
                'bams' : sv_args[key]['bams'],
                'inc_info': sv_args[key]['script'],
            }
        )
