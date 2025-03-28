#!/usr/bin/env python
"""
Calculate Pdet values from classifier.
"""
__author__ = "Alex Drlica-Wagner, Kabelo Tsiane"

import time
import os
import pickle
import itertools
import yaml
import warnings
import sys
sys.path.append(os.path.expandvars('$HOME/software/'))

import numpy as np
import healpy as hp

import sklearn
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
import xgboost as xgb
import pandas as pd

import dc2_satellite_census.code.utils as utils
from IPython.core.debugger import set_trace

class SurveySelectionFunction(object):

    def __init__(self, config_file, lsst_version):
        """
        Parameters
        ----------
        config_file : classifier configuration file

        Returns
        -------
        ssf : instance of survey selection function object
        """
        self.config = yaml.safe_load(open(config_file))
        self.survey = self.config['survey']
        self.load(lsst_version)

    def load_mask(self):
        """ Load the sky area mask

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        try:
            self.mask = utils.load_mask(self.survey)
        except:
            print("WARNING: Failed to load mask.")
            self.mask = None

    def load_sims(self, lsst_version):
        """ Load the input simulations for training

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # try:
        #     self.sims = utils.load_sims(self.survey)
        # except:
        #     print("WARNING: Failed to load sims.")
        #     self.sims = None
        if lsst_version is None:
            print("Version not specified.")
        elif lsst_version == 6: 
            self.sims = pd.read_csv('/home/kb/software/simple_adl/notebooks/results_dir/v6_results.csv')
        elif lsst_version == 7: 
            self.sims = pd.read_csv('/home/kb/software/simple_adl/notebooks/results_dir/v7_results.csv')
        return
            
    def load_density(self):
        """ Load stellar density maps """
        try:
            self.stellar_density = utils.load_density(self.survey)
        except:
            print("WARNING: Failed to load stellar density.")
            self.stellar_density = None
            return

        # Update stellar density
        stellar_density = self.get_stellar_density(self.sims['RA'],self.sims['DEC'])
        self.sims['STELLAR_DENSITY'] = stellar_density

    def load_classifier(self,filename=None):
        """ Load classifier from a pickle file.

        Parameters
        ----------
        filename : file containing the pickled classifier

        Returns
        -------
        classifier : the classifier
        """
        if filename is None: 
            filename = os.path.join(utils.get_datadir(),self.config['filename'])
        if os.path.exists(filename + '.gz') and not os.path.exists(filename):
            print("Unzipping classifier...")
            os.system('gunzip -c %s.gz > %s'%(filename,filename))

        if not os.path.exists(filename):
            print("WARNING: Classifier %s not found."%filename)
            self.classifier = None
            return

        print('Loading classifier from %s ...'%(filename))
        try: 
            #self.classifier = xgb.XGBClassifier(**self.config['classifier']['params'])
            self.classifier = xgb.XGBClassifier()
            self.classifier.load_model(filename)  # load data
        except AttributeError:
            # Ugh, pickles are for delis...
            self.classifier = pickle.loads(open(filename,'r').read())
        return self.classifier

    def load(self, lsst_version):
        """ Load all components """
        self.load_sims(lsst_version)
        self.load_mask()
        #self.load_density()
        #self.load_classifier()
        self.data_sim = None

    def create_train_test_sample(self):
        """ Create the training and testing samples """
        if self.sims is None: raise Exception("No simulations provided.")
        if self.mask is None: raise Exception("No mask provided.")

        sel = utils.select(self.sims,self.mask)
        pix = hp.ang2pix(hp.get_nside(self.mask), self.sims['RA'], 
                         self.sims['DEC'], lonlat=True)
        #sel &= (self.mask[pix] == 0)     #This doesnt seem necessary
        sel &= (self.sims['DIFFICULTY'] == 0)
        self.data_sim = self.sims[sel]
        self.sel_detect = utils.detect(self.data_sim)
        mc_source_id_detect = self.data_sim['MC_SOURCE_ID'][self.sel_detect]

        #Construct dataset
        # x = []
        # for key, operation in self.config['classifier']['params_intrinsic']:
        #     key = key.upper()
        #     assert operation.lower() in ['linear', 'log'], 'ERROR'
        #     if operation.lower() == 'linear':
        #         x.append(self.data_sim[key])
        #     else:
        #         x.append(np.log10(self.data_sim[key]))
        x = []
        # hack for intrinsic satellite params
        x.append(np.log10(self.data_sim['DISTANCE']))
        x.append(self.data_sim['ABS_MAG'])
        x.append(np.log10(self.data_sim['R_PHYSICAL']))
        
        X = np.vstack(x).T
        Y = self.sel_detect

        # Use the same random state for reproducibility
        index = np.arange(len(X))
        samples = train_test_split(X,Y,index,test_size=0.1,random_state=0)
        self.X_train, self.X_test = samples[0:2]
        self.Y_train, self.Y_test = samples[2:4]
        self.sel_train, self.sel_test = samples[4:6]

        print("Ntrain =  %i; Ntest = %i "%(len(self.X_train),len(self.X_test)))

    def train_classifier(self,n_jobs=4):
        """ Train the classification algorithm

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.create_train_test_sample()
        #Train MLP classifier
        t_start = time.time()

        if False:
            print("Training MLPClassifier...")
            model = MLPClassifier(random_state=0)
            parameter_space = {'alpha': [0.001, 0.005, 0.01, 0.05], 
                               'hidden_layer_sizes': list(itertools.product((50,45,40),repeat=3))}
        else:
            print("Training XGBClassifier...")
            #model = xgb.XGBClassifier(self.config['classifier']['params'])
            model = xgb.XGBClassifier()
            parameter_space = {'learning_rate': [0.01,0.05,0.1],
                               'max_depth': [6,7,8],
                               'n_estimators': [100,250,500],
                               'max_delta_step': [1],
                               'seed': [1337]}

        clf = GridSearchCV(model, parameter_space, cv=3, verbose=1, n_jobs=n_jobs)
        # Some DepricationWarnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",DeprecationWarning)
            clf.fit(self.X_train, self.Y_train)
        self.classifier = clf.best_estimator_

        # Print the best hyperparameters:
        print(clf.best_params_)
        t_end = time.time()
        print('  ... training took %.2f seconds'%(t_end - t_start))

    def predict(self, mask=True, **kwargs):
        """ Get the classifier prediction value.

        Note that the BADBITS mask value is applied.

        Parameters
        ----------
        distance       : heliocentric distance (kpc)
        r_physical     : physical half-light radius (kpc)
        abs_mag        : absolute V-band magnitude
        ra [optional]  : right ascension (deg)
        dec [optional] : declination (deg)

        Returns
        -------
        Pdet           : detection probability
        """
        if self.classifier is None: raise Exception("No classifier loaded.")

        # Get the ra,dec values
        ra  = kwargs.pop('ra',None)
        dec = kwargs.pop('dec',None)

        # if ra is not None and dec is not None:
        #     if ('stellar_density' not in kwargs) and \
        #        (self.stellar_density is not None):
        #         kwargs['stellar_density'] = self.get_stellar_density(ra,dec)

        pred = self.predict_proba(**kwargs)

        if mask and (ra is not None) and (dec is not None):
            val = self.get_mask_value(ra,dec)
            pred[(val & utils.BADBITS) != 0] = 0

        return pred

    def predict_proba(self, **kwargs):
        """ Call underlying predictor given intrinsic parameters.
        """
        # x_eval = []
        # for key, operation in self.config['classifier']['params_intrinsic']:
        #     assert operation.lower() in ['linear', 'log'], 'ERROR'
        #     if operation.lower() == 'linear':
        #         x_eval.append(kwargs[key])
        #     else:
        #         x_eval.append(np.log10(kwargs[key]))
        
        x_eval = []
        # hack for intrinsic satellite params
        x_eval.append(np.log10(self.data_sim['DISTANCE']))
        x_eval.append(self.data_sim['ABS_MAG'])
        x_eval.append(np.log10(self.data_sim['R_PHYSICAL']))
        
        x_eval = np.vstack(x_eval).T
        pred = self.classifier.predict_proba(x_eval)[:,1]

        return pred

    def get_stellar_density(self, ra, dec):
        nside = hp.get_nside(self.stellar_density)
        pix = hp.ang2pix(nside, ra, dec, lonlat=True)
        return self.stellar_density[pix]

    def get_mask_value(self, ra, dec):
        nside = hp.get_nside(self.mask)
        pix = hp.ang2pix(nside,ra,dec,lonlat=True)
        return self.mask[pix]

    def write_classifier(self,filename=None, force=False):
        """ Save the classifier to disk.

        Parameters
        ----------
        filename : the name of the output file
        force    : overwrite output file

        Returns
        -------
        None
        """

        if filename is None: filename = self.config['filename']
        if os.path.exists(filename):
            if force:
                print("WARNING: Overwriting %s..."%filename)
            else:
                print("File already exists; use `force` to overwrite")
                return
        try:
            print("Saving model...")
            self.classifier.save_model(filename)
        except AttributeError:
            print("Pickling classifier...")
            # Ugh, pickles are for delis...
            classifier_pkl = pickle.dumps(self.classifier)
            writer = open(filename, 'w')
            writer.write(classifier_pkl)
            writer.close()
