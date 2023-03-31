import loompy as lp
import pickle
import os
import argparse

def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Write loom",)
    parser.add_argument('--input_cto', '-i', type=str, required=True,
                        help='Path to cisTopic object pickle file.')

    parser.add_argument('--output_loom', '-o', type=str, required=True,
                    help='Path to out loom.')

    return parser


def main():
    """
    The main executable function
    """

    parser = make_argument_parser()
    args = parser.parse_args()


    cto_path = args.input_cto
    print('input file:', cto_path)

    loom_path = args.output_loom
    print('output file:', loom_path)
    
    if not os.path.exists(loom_path):
        with open(cto_path, 'rb') as f:
            cto = pickle.load(f)

        print(f"Loaded filtered cistopic object {cto_path}")
        lp.create(
            filename = loom_path,
            layers=cto.fragment_matrix,
            row_attrs={ 'Gene': cto.region_names },
            col_attrs={ 'CellID': [ x.split('__')[0]  for x in cto.cell_names ] }
        )

        print(f"\tFinished {loom_path} loom writing")
        
    else:
        print(f'{loom_path} already exists, skipping...')


if __name__ == "__main__":
    main()