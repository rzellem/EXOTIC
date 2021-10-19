import os
import json
import urllib
import numpy as np

base_uri = "https://exoplanets.nasa.gov/watch_results/"

class ExoplanetWatch():
    url = 'https://exoplanets.nasa.gov/api/v1/exoplanet_watch_results/?per_page=-1&order=planet_name+asc'
    def __init__(self) -> None:
        '''An interface to the Exoplanet Watch Results Database

        Parameters
        ----------

        Returns
        -------
        None
        '''        

        # download the results JSON
        r = urllib.request.urlopen(self.url)
        jdata = json.loads(r.read().decode(r.info().get_param('charset') or 'utf-8'))
        self.alldata = jdata['items']

        # list targets
        self.target_list = [list(item.keys())[0] for item in self.alldata]
        print(f"Loaded {len(self.target_list)} targets")

    def get(self, target):
        """ Returns a dictionary of the result for the given target

        Parameters
        ----------
        target : str
            The target name (case + space sensitive)
        
        Returns
        -------
        result : dict
            The result dictionary with keys:
            ['host', 'name', 'priors', 'ephemeris', 'timestamp', 
            'identifier', 'observations', 'update_action', 
            'reduction_count', 'observation_count']
        """
        try:
            ti = self.target_list.index(target)
            return ExoplanetWatchTarget(self.alldata[ti][target])
        except ValueError:
            print(f"{target} not in list")
            return {}

class ExoplanetWatchTarget():
    def __init__(self,result) -> None:
        self.raw_result = result
        self.host = result['host']
        self.name = result['name']
        self.priors = result['priors']
        self.timestamp = result['timestamp']
        self.identifier = result['identifier']
        self.reduction_count = result['reduction_count']
        self.observation_count = result['observation_count']
        self.ephemeris = check4floats(result['ephemeris'])

        self.observations = []
        for i in range(len(result['observations'])):
            self.observations.append(
                ExoplanetWatchObservation(self.name, result['observations'][i]))

        self.ephemeris['ephemeris_url'] = os.path.join(base_uri,
            self.ephemeris['files']['file_oc_png'][2:])

    def __str__(self) -> str:
        return f"{self.name}: {self.reduction_count} light curves"

    def __repr__(self) -> str:
        return f"{self.name}: {self.reduction_count} light curves"

    def __rstr__(self) -> str:
        return f"{self.name}: {self.reduction_count} light curves"


def check4floats(dictionary):
    for key in dictionary:
        try:
            dictionary[key] = float(dictionary[key])
        except:
            pass
    return dictionary

class ExoplanetWatchObservation():
    def __init__(self, name, observation) -> None:
        self.raw_observation = observation
        self.name = name
        self.files = observation['files']
        self.filter = observation['filter']
        self.obscode = observation['obscode']
        self.identifier = observation['identifier']
        self.secondary_obscodes = observation['secondary_obscodes']
        self.errors = self.translate(check4floats(observation['errors']))
        self.parameters = self.translate(check4floats(observation['parameters']))
        self.lightcurve_url = os.path.join(base_uri,self.files['file_lc_png'][2:])
        self.posterior_url = os.path.join(base_uri,self.files['file_po_png'][2:])

    def __str__(self) -> str:
        return f"{self.name} light curve @ {self.parameters['Tc']:.4f}"

    def __repr__(self) -> str:
        return f"{self.name} light curve @ {self.parameters['Tc']:.4f}"

    def __rstr__(self) -> str:
        return f"{self.name} light curve @ {self.parameters['Tc']:.4f}"

    def get_data(self):
        """ Download the light curve data

        Parameters
        ----------
        
        Returns
        -------
        (time, flux, fluxerr, airmass, airmass_detrend) : tuple of numpy arrays
        """
        r = urllib.request.urlopen(
            os.path.join(base_uri,self.files['file_data_json'][2:]))
        jdata = json.loads(r.read().decode(r.info().get_param('charset') or 'utf-8'))
        return np.array(jdata[1:],dtype=np.float).T

    def translate(self,rdict):
        """ Translates the keys to a compatible format for EXOTIC/ELCA

        Parameters
        ----------
        rdict : dict
            The dictionary of parameters
        
        Returns
        -------
        lc_pars : dict
            A dictionary with keys corresponding to exotic's lc format
        """
        lc_pars = dict(rdict)
        translate_keys = {
            'Tc':'tmid',
            'Am1': 'a1',
            'Am2': 'a2',
            'a/R*': 'ars',
            'Rp/R*': 'rprs',
            'Period': 'per'
        }
        for k in translate_keys:
            if k in lc_pars:
                lc_pars[translate_keys[k]] = lc_pars[k]
                del lc_pars[k]
        for k in lc_pars:
            lc_pars[k] = float(lc_pars[k])
        return lc_pars

if __name__ == "__main__":

    EW = ExoplanetWatch()
    print(EW.target_list)
    result = EW.get('WASP-33 b')
    print(result)
    time, flux, fluxerr, airmass, airmasscorr = result.observations[0].get_data()
    print(result.observations[0].parameters)