#################################################################
#################################################################
############### Dubois RNA-seq Analysis #########################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#######################################################
#######################################################
########## 1. App Configuration
#######################################################
#######################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. Flask modules #####
from flask import Flask, render_template, url_for
from flask_basicauth import BasicAuth
import os
import pandas as pd
entry_point = '/dubois'
app = Flask(__name__, static_url_path=os.path.join('/app/static'))

##### 2. Authentication #####
if os.environ.get('AUTHENTICATION'):
	for key in ['BASIC_AUTH_USERNAME', 'BASIC_AUTH_PASSWORD']:
		app.config[key] = os.environ.get(key)
	app.config['BASIC_AUTH_FORCE'] = True
	basic_auth = BasicAuth(app)

##### 3. Prefix middleware #####
class PrefixMiddleware(object):

	def __init__(self, app, prefix=''):
		self.app = app
		self.prefix = prefix

	def __call__(self, environ, start_response):
		if environ['PATH_INFO'].startswith(self.prefix):
			environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
			environ['SCRIPT_NAME'] = self.prefix
			return self.app(environ, start_response)
		else:
			start_response('404', [('Content-Type', 'text/plain')])
			return ["This url does not belong to the app.".encode()]
app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix=entry_point)

#######################################################
#######################################################
########## 2. Routes
#######################################################
#######################################################
##### Handles routes used to generate notebooks.

##################################################
########## 2.1 Webpages
##################################################

#############################################
########## 1. Home
#############################################
### Landing page for the website. Links to analyze() endpoint.
### Links to: analyze().

@app.route('/')
def index():
	# Read samples
	print(os.getcwd())
	sample_dataframe = pd.read_table('app/static/dubois-metadata.txt')
	conditions = sample_dataframe['Condition'].unique()
	sample_conditions = ["Control", " + delta", " + delta, 1d", " + delta, 2d", " + delta, 3d", " + delta, 7d"]
	conditions_dict = {
		'show': sample_conditions,
		'hide': list(set(conditions)-set(sample_conditions))
	}
	print(sample_dataframe)
	return render_template('index.html', sample_dataframe=sample_dataframe, conditions_dict=conditions_dict)

if __name__ == "__main__":
	app.run(debug=True, host='0.0.0.0')
