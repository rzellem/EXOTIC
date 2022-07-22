from json import dump, dumps
from numpy import mean, median, std
from pathlib import Path

try:
    from utils import round_to_2
except ImportError:
    from .utils import round_to_2
try:
    from version import __version__
except ImportError:
    from .version import __version__


class VSPOutputFiles:
    def __init__(self, fit, p_dict, i_dict, durs, vsp_params):
        self.fit = fit
        self.p_dict = p_dict
        self.i_dict = i_dict
        self.durs = durs
        self.dir = Path(self.i_dict['save'])
        self.vsp_params = vsp_params

    def aavso(self, comp_star, airmasses, ld0, ld1, ld2, ld3):

        params_file = self.dir / f"vspAAVSO_{self.p_dict['sName']}_{self.i_dict['date']}.txt"
        with params_file.open('w') as f:
            f.write("#TYPE=EXTENDED\n"  # fixed
                    f"#OBSCODE={self.i_dict['aavso_num']}\n"  # UI
                    f"#SECONDARY_OBSCODES={self.i_dict['second_obs']}\n"  # UI
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE_TYPE=BJD_TDB\n"  # fixed
                    f"#OBSTYPE={self.i_dict['camera']}\n")
            f.write(
                "# EXOTIC is developed by Exoplanet Watch (exoplanets.nasa.gov/exoplanet-watch/), a citizen science "
                "project managed by NASA's Jet Propulsion Laboratory on behalf of NASA's Universe of Learning. "
                "This work is supported by NASA under award number NNX16AC65A to the "
                "Space Telescope Science Institute.\n"
                "# Use of this data is governed by the AAVSO Data Usage Guidelines: "
                "aavso.org/data-usage-guidelines\n")

            f.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")
            for aavsoC in range(0, len(self.vsp_params['time'][self.vsp_params['OOT']])):
                f.write(f"{self.p_dict['sName']},"f"{round(self.vsp_params['time'][self.vsp_params['OOT']][aavsoC], 8)}," f"{round(self.vsp_params['mag'][aavsoC], 6)}," f"{round(self.vsp_params['mag_err'][aavsoC],5)},"
                        "V,NO,STD," f"{self.vsp_params['cname']}," f"{self.vsp_params['cmag']}," "na,na," f"{round(airmasses[aavsoC], 7)}," "na," 
                        f"{self.vsp_params['chart_id']}," "na\n")


