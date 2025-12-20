from pyteomics import mgf
import argparse
import sys
import os
import glob
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Generate a summary of a library')
    parser.add_argument('library', help='library file')
    parser.add_argument('output_summary')

    args = parser.parse_args()

    output_list = []


    # Pyteomics will fail silently if there is no SCANS= line, so we will warn the user
    with open(args.library, "r") as f:
        for line in f:
            line = line.strip()

            # If we hit SCANS=, we're done
            if line.startswith("SCANS="):
                break

            # If we hit END IONS before SCANS=
            if line == "END IONS":
                raise ValueError('No Scan Numbers')

    with mgf.read(args.library, index_by_scans=True) as reader:
        for spectrum in reader:
            # We will use the spectrum ID as the key for the dictionary
            spectrum_id = spectrum['params'].get('spectrumid', spectrum['params'].get('title', "No ID"))
            compound_name = spectrum['params'].get('compound_name', spectrum['params'].get("name", "No Compound"))
            smiles = spectrum['params'].get('smiles', "")
            collision_energy = spectrum['params'].get('collision_energy', 0)
            instrument = spectrum['params'].get('instrument', "")
            ion_source = spectrum['params'].get('ion_source', "")
            charge = spectrum['params'].get('charge', 0)
            adduct = spectrum['params'].get('adduct', "NA")
            try:
                precursormz = spectrum['params'].get('pepmass', [0])[0]
            except:
                precursormz = 0

            output_dictionary = {}
            output_dictionary["spectrum_id"] = spectrum_id
            output_dictionary["compound_name"] = compound_name
            output_dictionary["smiles"] = smiles
            output_dictionary["collision_energy"] = collision_energy
            output_dictionary["instrument"] = instrument
            output_dictionary["ion_source"] = ion_source
            output_dictionary["charge"] = charge
            output_dictionary["adduct"] = adduct
            output_dictionary["precursormz"] = precursormz
            
            output_list.append(output_dictionary)

    # creating an output df
    output_df = pd.DataFrame(output_list)
    if output_df.empty:
        # Create an empty DataFrame with the correct columns
        output_df = pd.DataFrame(columns=["spectrum_id", "compound_name", "smiles", "collision_energy", "instrument", "ion_source", "charge", "adduct", "precursormz"])
    output_df.to_csv(args.output_summary, sep="\t", index=False)

if __name__ == "__main__":
    main()