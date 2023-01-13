"""
gefsv12ww3uploadaws.py

VERSION AND LAST UPDATE:
 v1.0  09/16/2022

PURPOSE:
 Upload GEFv12 reforecast to AWS using boto3
 https://noaa-nws-gefswaves-reforecast-pds.s3.amazonaws.com/index.html

USAGE:
 See the 3 input arguments below and enter the path strings

DEPENDENCIES:
 See dependencies.py and the imports below.

AUTHOR and DATE:
 09/16/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import numpy as np
from pylab import *
from boto3.s3.transfer import S3Transfer
import boto3
import warnings
warnings.filterwarnings("ignore")

# https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/s3.html#client
# https://noaa-nws-gefswaves-reforecast-pds.s3.amazonaws.com/index.html

bucket_name = ""
awsricakid=''
awsricsak=''

# start
client = boto3.client('s3', aws_access_key_id=awsricakid,aws_secret_access_key=awsricsak)
transfer = S3Transfer(client)

bdir="" #path

ftime=str(sys.argv[1]) # ftime="20000108" 
fst=str(sys.argv[2]) # fst="1"
fem=str(sys.argv[3]) # fem="c00"

year=str(ftime[0:4])
sfx=np.array([".global.0p25.grib2",".global.0p25.idx",".spec.nc",".tab.nc"]).astype('str')
smin=np.array([1000000000,50000,300000000,5000000]).astype('int')  # minimum file size

# File info
print(" ")
print("-----------------------")
print(" ")
for i in range(0,size(sfx)):
	filepath = bdir+fst+"/"+fem+"/gefs.wave."+ftime+"."+fem+sfx[i]
	if i<2:
		folder_name = "GEFSv12/reforecast/"+str(year)+"/"+ftime+"/gridded"
	else:
		folder_name = "GEFSv12/reforecast/"+str(year)+"/"+ftime+"/station"

	filename = "gefs.wave."+ftime+"."+fem+sfx[i]

	try:
		fsize = boto3.resource('s3', aws_access_key_id=awsricakid,aws_secret_access_key=awsricsak).Bucket(bucket_name).Object(folder_name+"/"+filename).content_length
	except:
		print("  no file "+filename+" on AWS yet. Uploading ...")
		# transfer file
		transfer.upload_file(filepath, bucket_name, folder_name+"/"+filename)
	else:
		# Check size
		fsize = boto3.resource('s3', aws_access_key_id=awsricakid,aws_secret_access_key=awsricsak).Bucket(bucket_name).Object(folder_name+"/"+filename).content_length
		if fsize>smin[i]:
			print("  ok, file "+filename+" on AWS.")
		else:
			# delete file
			client.delete_object(Bucket=bucket_name, Key=folder_name+"/"+filename)
			# transfer file
			transfer.upload_file(filepath, bucket_name, folder_name+"/"+filename)

	# Check size
	fsize = boto3.resource('s3', aws_access_key_id=awsricakid,aws_secret_access_key=awsricsak).Bucket(bucket_name).Object(folder_name+"/"+filename).content_length
	print(" Done "+filename+" Size: "+np.str(fsize)); del fsize
	print(" ")


