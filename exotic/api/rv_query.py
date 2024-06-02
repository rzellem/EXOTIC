import os
import json
from dace_query.spectroscopy import Spectroscopy

class RadialVelocityQuery:

    def __init__(self, target=None):
        self.target = target
        self.data = {
            "target": target,
            "time": [],
            "instrument": [],
            "rv": [],
            "rv_err": []
        }

        if target is not None:
            self.query_dace()

            # more queries can be added here

    def query_dace(self):

        try:
            # query dace for number of RV observations
            observedTargets = Spectroscopy.query_database(
                limit=1000,
                filters={"public": {"is": True}, "obj_id_catname": {"contains": [self.target]}},
                output_format="pandas"
            )
        except:
            print(f"error for {self.target}")

        if len(observedTargets) == 0:
            print(f"No RV observations for {self.target}")
        else:
            print(f"Found {len(observedTargets)} RV observations for {self.target} from DACE")

            # extract time, instrument name, rv, error
            self.data["time"].extend(observedTargets["obj_date_bjd"] + 2450000) # TODO check this conversion
            self.data["instrument"].extend(observedTargets["ins_name"])
            self.data["rv"].extend(observedTargets["spectro_ccf_rv"])
            self.data["rv_err"].extend(observedTargets["spectro_ccf_rv_err"])

            # TODO add activity index? other metrics?

    ###############################
    # add more functions here
    ###############################

    def save_to_disk(self, path):
        # save data as json
        with open(path, 'w') as f:
            json.dump(self.data, f, ident=4)

    @staticmethod
    def load_from_disk(path):
        # load data from json
        with open(path, 'r') as f:
            data = json.load(f)
        rv = RadialVelocityQuery(target=None)
        rv.data = data
        rv.target = data["target"]
        return rv
