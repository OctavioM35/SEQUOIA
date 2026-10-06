import os
import numpy as np
import json


def diference(outdir,label,dicc):
        '''
        Generates a JSON file containing the posterior medians and the lower and upper bounds of their 68% credible intervals, 
        together with the injection parameters of the SXS event.

        INPUTS:
        ------

        -Outdir: The outdir directory where the PE results are.
        -dicc: The created dicctionay in config.py 
        '''
        
        title = label + '_result.json'
        inyection_diccionary = {'mass_ratio' : [0,0,0] , 'chirp_mass': [0,0,0] , 'luminosity_distance': [0,0,0] ,'chi_1' : [0,0,0], 'chi_2': [0,0,0]}

        
        with open(os.path.join(outdir, title), "r") as f:
            data = json.load(f)
            content = data["posterior"]["content"]
            for key in inyection_diccionary:
                if key in content:
                    values = np.asarray(content[key], dtype=float)

                    median = np.median(values)
                    lower = np.percentile(values, 16)
                    upper = np.percentile(values, 84)

                    inyection_diccionary[key] = [
                        median,
                        upper - median,   # error superior
                        median - lower    # error inferior
                    ]


        inyection_diccionary = {
                    "obtained_values": inyection_diccionary,
                    "injected_values": dicc["inyection"]
                    } 

        
        output_path = os.path.join(
            outdir,
            "Obtained vs inyected parameters for surrogate " + label + '.json'
        )


        with open(output_path, "w") as f:
            json.dump(inyection_diccionary, f, indent=4)