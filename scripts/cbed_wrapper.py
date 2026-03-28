
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter

def bla():

    def _pars2r(pars_dict):
        """Convert a Python dict to an R named list."""
        r_list = ro.ListVector(
            {k: ro.FloatVector([v]) if isinstance(v, (int, float))
                else ro.FloatVector(v)
                for k, v in pars_dict.items()}
        )
        return r_list



    def r2dict(varbs):
        result = {}
        for name, val in zip(varbs.names, varbs):
            # nested ListVector — recurse
            if isinstance(val, ro.vectors.ListVector):
                result[name] = r2dict(val)
            else:
                try:
                    arr = np.array(val, dtype=float)
                    result[name] = arr[0] if arr.size == 1 else arr
                except (ValueError, TypeError):
                    result[name] = val
        return result 

    POC_flux = 0
    btm_O2   = 0
    btm_no3  = 0
    btm_nh4  = 0
    btm_temp = 0


    script_path = 'src/cbed_R/cbed_v1_func.R'


    source = lambda x: f"source(\"{x}\")"
    r = ro.r
    r(source(script_path))
    r(source('scripts/default_CBED_pars.R'))  # contains get_default_pars function

    # set up parameters
    r_pars = r('get_default_pars()')
    _default_pars = dict(zip(r_pars.names, [list(v) if len(v) > 1 else v[0] for v in r_pars]))
    # edit _default_pars as needed

    _default_pars["J.OM"] = 11.05*365/10
    # pars["J.OM"]  = np.nanmean(df["POC_flux"]) * 365 / 10
    # pars["O2.w"]  = np.nanmean(df["btm_O2"])   / 1000
    # pars["NO3.w"] = np.nanmean(df["btm_no3"])  / 1000
    # pars["NH4.w"] = np.nanmean(df["btm_nh4"])  / 1000
    # pars["temp"]  = np.nanmean(df["btm_temp"])

    # por0 = pars["por.0"]
    # pars["w"] = (
    #     np.nanmean(df["j_caco3_arag"]) * 100 / 2.71 +
    #     np.nanmean(df["j_caco3_calc"]) * 100 / 2.94 +
    #     np.nanmean(df["j_sidet"])      *  60 / 2.65 +
    #     np.nanmean(df["j_fedet"])      * 160 / 5.24 +
    #     np.nanmean(df["j_pdet"])       * 120 / 2.3  +
    #     np.nanmean(df["j_lithdet"])         / 2.65  +
    #     pars["J.OM"] * (1e4 / 1e6 / SEC_TO_YEAR) * 22.4 / 0.9
    # ) / 10000 * SEC_TO_YEAR / (1 - por0)


    # _default_pars
    r_pars_updated = _pars2r(_default_pars)

    r("""
    manual_control_bioturbation <- F
    manual_control_bioirrigation <- F
    manual_OM_decay_rate <- F
    """)

    ro.globalenv["py_pars"] = r_pars_updated
    r_out_ss = r("cbed_model(py_pars)")



    test = r2dict(r_out_ss)

%time bla()