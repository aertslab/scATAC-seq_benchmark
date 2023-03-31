import loompy as lp
import pickle
import argparse
import os
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling, get_consensus_peaks
import pyranges as pr
import requests
import ray


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
    print('input file:', loom_path)

    ####
    
    chromsizes = pd.read_csv(chromsizes_path)
  
    bed_paths = {x.split('/')[-1].split('__')[0].split('.bed.gz')[0]:x for x in glob.glob(bed_path_dict[sample] + "/*")}
    peak_path = os.path.join('SCREEN_peaks', f'{sample}__SCREEN_consensus_peaks')
    peak_path_pkl = os.path.join('SCREEN_peaks', f'{sample}__SCREEN_consensus_peaks.pkl')

    if not os.path.exists(peak_path):
        os.makedirs(peak_path)
        
    if not len(os.listdir(peak_path)) > 0:
        # Run peak calling
        narrow_peaks_dict = peak_calling('macs2',
                                 bed_paths,
                                 peak_path,
                                 genome_size='hs',
                                 n_cpu=20,
                                 input_format='BEDPE',
                                 shift=73, 
                                 ext_size=146,
                                 keep_dup = 'all',
                                 q_value = 0.05,
                                 )
        
        with open(peak_path_pkl, "wb") as f:
            pickle.dump(narrow_peaks_dict, f, protocol=4)
            
    else:
        print(f"already ran {sample}, delete {peak_path}")
        
    if not os.path.exists(consensus_peaks_path_bed):
        with open(peak_dict_path, 'rb') as f:
            peak_dict = pickle.load(f)

        consensus_peaks = get_consensus_peaks(
            peak_dict,
            peak_half_width,
            chromsizes=chromsizes,
            path_to_blacklist=path_to_blacklist
        )   
        
        consensus_peaks.to_bed(
            path=consensus_peaks_path_bed,
            keep=True,
            compression='infer',
            chain=False
        )
    
    
if __name__ == "__main__":
    main()
