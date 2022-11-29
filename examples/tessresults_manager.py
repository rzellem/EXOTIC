from time import sleep
import subprocess
import argparse
import atexit
import signal
import json
import os
from bs4 import BeautifulSoup as bs
import pandas as pd
import time
import requests
import random

from exotic4tess import tap_query

def kill_all(ps):
    for p in ps:
        if isinstance(ps,dict):
            os.kill(ps[p].pid, signal.SIGTERM)
        else:
            os.kill(p.pid, signal.SIGTERM)

def spawn_bot(target):
    return subprocess.Popen([
        "python3.10","exotic4tess.py", 
        "-t", target ])

def kill_bot(ps):
    os.kill(ps.pid, signal.SIGTERM)

def steal_epwResults():
    planet_list = []
    num = 0
    page = requests.get("https://exoplanets.nasa.gov/exoplanet-watch/results/")
    soup = bs(page.content)
    epw_targets = [i for i in soup.find_all(class_='module exoplanet_watch_results')]
    targ_str = str(epw_targets[0])
    pl_loc = targ_str.find("planet_name&")
    pl_substr = targ_str[pl_loc::]
    while(pl_substr.find("planet_name&")!=-1):
        pname = pl_substr[pl_substr.find(":"):pl_substr.find(",")]
        pname = pname[pname.find(";")+1:pname.rfind("&")]
        planet_list.append(pname)
        pl_substr = pl_substr[pl_substr.find(",")::]
        pl_substr=pl_substr[pl_substr.find("planet_name&")::]
        num+=1
        print(pname)
        print("num of planets: "+str(num))
    return planet_list

if __name__ == "__main__":
    epw_targets = steal_epwResults()
    import pdb; pdb.set_trace()
    for i in range(len(epw_targets)):
        target = epw_targets[i]
        print(target)
        try:
            ps = subprocess.call([
                "python3.10","exotic4tess.py",
                "-t", target])
        except:
            print("failed")