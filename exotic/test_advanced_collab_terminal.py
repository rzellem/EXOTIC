import os
import re
import json
from api.nea import NASAExoplanetArchive, ExoplanetNotFoundError
from api.exofop import ExoFOP

def get_target():
    while True:
        target = input('Please enter the name of your exoplanet or exoplanet candidate target (i.e. "HAT-P-32 b", "TOI-4479 b", or "TIC 51094586") and press return: ')
        if target.strip() == "":
            print('Exoplanet target may not be blank.')
        else:
            return target

def confirm_target(target):
    if re.match(r'^TIC \d{8}$', target):
        print(f'You have entered {target}. This appears to be a TESS Input Catalog ID.')
    elif re.match(r'^TOI-\d{4} \w$', target):
        print(f'You have entered {target}. This appears to be a TESS Target of Interest.')
    else:
        print(f'You have entered {target}. This appears to be a standard exoplanet target.')
    confirmation = input('Is this correct (1 for yes or 2 for no)? ')
    return confirmation == '1'

def fix_planetary_params(p_param_dict):
    # Implement this function according to your needs
    return p_param_dict

def main():
    planetary_params = ""
    while not planetary_params:
        target = get_target()
        while not confirm_target(target):
            target = get_target()

        print(f'Searching NASA Exoplanet Archive for "{target}"...')

        targ = None
        try:
            targ = NASAExoplanetArchive(planet=target)
            target = targ.planet_info()[0]
            print(f'Found target "{target}" in the NASA Exoplanet Archive')
        except ExoplanetNotFoundError as e:
            print(f"Caught an exception: {e}")
            print('Sorry, we can\'t find your target in the NASA Exoplanet Archive. Trying ExoFOP...')
            try:
                targ = ExoFOP(tic_code=target)
                # targ.query_exofop()
                target = targ.planet_info()[0]
                print(f'Found target "{target}" in the ExoFOP')
            except Exception as e:
                print(f"Caught an exception: {e}")
                print('''Sorry, we can't find your target in the ExoFOP either. Unfortunately, this
                isn't going to work until we can find it. Please try
                different formats for your target name, until the target is located.
                Looking it up in the NASA Exoplanet Archive at https://exoplanetarchive.ipac.caltech.edu/
                or the ExoFOP at https://exofop.ipac.caltech.edu/
                might help you know where to put the spaces and hyphens and such.
                ''')
                print('''If your target is a candidate, you may need to create your own inits.json file and
                add it to the folder with your FITS images.
                ''')
                target = None
        if target:
            p_param_string = targ.planet_info(fancy=True)
            planetary_params = '"planetary_parameters": ' + p_param_string
            p_param_dict = json.loads(p_param_string)
            planetary_params = fix_planetary_params(p_param_dict)
            print(f'Loading {"NASA Exoplanet Archive" if isinstance(targ, NASAExoplanetArchive) else "ExoFOP"} planetary parameters for {target}')
            print(f'{planetary_params}')
        else:
            print('Unable to find target in either NASA Exoplanet Archive or ExoFOP. Please try again with a different target name.')

    # Prompt for AAVSO code
    aavso_obs_code = input("Enter your AAVSO Observer code or press return to skip: ")
    if aavso_obs_code:
        sec_obs_code = input("Enter a secondary AAVSO Observer code or press return to skip: ")
    else:
        sec_obs_code = ""

if __name__ == "__main__":
    main()