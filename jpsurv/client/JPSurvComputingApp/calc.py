import rpy2.robjects as robjects
import os
import json
import urllib
import requests

# Input files
UPLOAD_DIR = os.path.join(os.getcwd(), 'tmp')
file_control_filename = "Breast_RelativeSurvival.dic"
file_data_filename = "Breast_RelativeSurvival.txt"
tokenId='636307'

response = requests.get('http://localhost:8080/queue/restapi/jpsurvqueue')
 
if response.status_code == 200:
	#Init the R Source
	rSource = robjects.r('source')
	rSource('./JPSurvWrapper.R')
	jpsurvDataString = json.dumps(response.json())
	getFittedResults = robjects.r['getFittedResults']
	StrVector = getFittedResults(UPLOAD_DIR, jpsurvDataString)
	print StrVector[0]
else:
	print 'queue is empty' 








