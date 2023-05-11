import glob
from IPython.display import HTML
import logging
import os
import pandas as pd
from pybiomart import Dataset
import pyranges as pr 
import ray
import subprocess
import shutil
import sys
from typing import Union, Dict, Sequence, Optional

from .utils import *

# Set stderr to null when using ray.init to avoid ray printing Broken pipe million times
_stderr = sys.stderr                                                         
null = open(os.devnull,'wb') 

class Homer(): 
    """
    Homer class.
    :class:`Homer` contains Homer for motif enrichment analysis on sets of regions. 
    
    Attributes
    ---------
    homer_path: str
        Path to Homer bin folder.
    bed_path: str
        Path to bed file containing region set to be analyzed with Homer.
    name: str
        Analysis name.
    outdir: str
        Path to folder to output Homer results.
    genome: str
        Homer genome label to use.
    size: str, optional
        Fragment size to use for motif finding. Default: 'given' [uses the exact regions you give it]
    mask: bool, optional
        Whether to mask repeats or not. Default: True
    denovo : bool, optional
        Whether to infer overrepresented motifs de novo. Default: False
    length: str, optional
        Motif length values. Default: 8,10,12
    meme_path: str, optional
        Path to meme bin folder. Meme will be used if given for motif annotation. Default: None
    meme_collection_path: str, optional
        Path to motif collection (in .cb format) to compare homer motifs with. Default: None
    cistrome_annotation: List, optional
        Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
        linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
        (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
        Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']
    motif_similarity_fdr: float, optional
        Minimal motif similarity value to consider two motifs similar. Default: 0.001
    orthologous_identity_threshold: float, optional
        Minimal orthology value for considering two TFs orthologous. Default: 0.0
    known_motifs: pd.DataFrame
        A dataframe containing known motif enrichment results.
    denovo_motifs: pd.DataFrame
        A dataframe containing de novo motif enrichment results.
    known_motif_hits: Dict
        A dictionary containing regions with motif hits for each known motif.
    denovo_motif_hits: Dict
        A dictionary containing regions with motif hits for each de novo motif.
    known_cistromes: Dict
        A dictionary containing regions with motif hits for each TF found with known motifs.
    denovo_motif_hits: Dict
        A dictionary containing regions with motif hits for each TF found de novo.

    References
    ---------
    Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining 
    Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities.
    Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
    """
    def __init__(self,
                 homer_path: str,
                 bed_path: str,
                 name: str,
                 outdir: str,
                 genome: str,
                 size: str = 'given',
                 mask: bool = True,
                 denovo: bool = False,
                 length: str = '8,10,12',
                 meme_path: str = None,
                 meme_collection_path: str = None,
                 cistrome_annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
                 motif_similarity_fdr: float = 0.001,
                 orthologous_identity_threshold: float = 0.0):
        """
        Initialize Homer class. 

        Parameters
        ---------
        homer_path: str
            Path to Homer bin folder.
        bed_path: str
            Path to bed file containing region set to be analyzed with Homer.
        name: str
            Analysis name.
        outdir: str
            Path to folder to output Homer results.
        genome: str
            Homer genome label to use.
        size: str, optional
            Fragment size to use for motif finding. Default: 'given' [uses the exact regions you give it]
        mask: bool, optional
            Whether to mask repeats or not. Default: True
        denovo : bool, optional
            Whether to infer overrepresented motifs de novo. Default: False
        length: str, optional
            Motif length values. Default: 8,10,12
        meme_path: str, optional
            Path to meme bin folder. Meme will be used if given for motif annotation. Default: None
        meme_collection_path: str, optional
            Path to motif collection (in .cb format) to compare homer motifs with. Default: None
        cistrome_annotation: List, optional
            Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
            linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
            (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
            Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']
        motif_similarity_fdr: float, optional
            Minimal motif similarity value to consider two motifs similar. Default: 0.001
        orthologous_identity_threshold: float, optional
            Minimal orthology value for considering two TFs orthologous. Default: 0.0

        References
        ---------
        Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining 
        Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities.
        Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
        """
        
        self.homer_path = homer_path
        self.bed_path = bed_path
        self.genome = genome
        self.outdir = outdir
        self.size = size
        self.len = length
        self.mask = mask
        self.denovo = denovo
        self.name = name
        self.meme_path = meme_path
        self.meme_collection_path = meme_collection_path
        self.cistrome_annotation = cistrome_annotation
        self.motif_similarity_fdr = motif_similarity_fdr
        self.orthologous_identity_threshold = orthologous_identity_threshold
        self.known_motifs = None
        self.denovo_motifs = None
        self.known_motif_hits = None
        self.denovo_motif_hits = None
        self.known_cistromes = None
        self.denovo_cistromes = None

    def run(self):
        """
        Run Homer
        """
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('Homer')
        
        if self.mask == True and self.denovo == False:
            cmd = os.path.join(self.homer_path, 'findMotifsGenome.pl') + ' %s %s %s -preparsedDir %s -size %s -len %s -mask -nomotif -keepFiles'
        if self.mask == True and self.denovo == True:
            cmd = os.path.join(self.homer_path, 'findMotifsGenome.pl') + ' %s %s %s -preparsedDir %s -size %s -len %s -mask -keepFiles'
        if self.mask == False and self.denovo == False:
            cmd = os.path.join(self.homer_path, 'findMotifsGenome.pl') + ' %s %s %s -preparsedDir %s -size %s -len %s -nomotif -keepFiles'
        if self.mask == False and self.denovo == True:
            cmd = os.path.join(self.homer_path, 'findMotifsGenome.pl') + ' %s %s %s -preparsedDir %s -size %s -len %s -keepFiles'
            
        cmd = cmd % (self.bed_path, self.genome, self.outdir, self.outdir, self.size, self.len)
        log.info("Running Homer for " + self.name + " with %s", cmd)
        try:
            subprocess.check_output(args=cmd, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
        
        try:
            self.known_motifs = self.load_known()
        except:
            log.info('No known results found')
        if self.denovo == True:
            try:
                self.denovo_motifs = self.load_denovo()
            except:
                log.info('No de novo results found')
           
        log.info("Annotating motifs for " + self.name)     
        self.add_motif_annotation_homer()
        log.info("Finding motif hits for " + self.name)   
        self.find_motif_hits(n_cpu=1)
        log.info("Getting cistromes for " + self.name) 
        self.get_cistromes(self.cistrome_annotation)
        
    def load_known(self):
        """
        Load known motif enrichment results from file.
        """
        known = pd.read_csv(os.path.join(self.outdir, 'knownResults.txt'), sep='\t')
        return known
    
    def load_denovo(self):
        """
        Load de novo motif enrichment results from file.
        """
        denovo = pd.read_html(os.path.join(self.outdir, 'homerResults.html'), header=0)[0].iloc[:,[7,2,3,4,5,6]]
        denovo.iloc[:,0] = [x.split('More')[0] for x in denovo.iloc[:,0]]
        denovo.to_csv(os.path.join(self.outdir, 'homerResults.txt'), sep='\t', index=False)
        return denovo
    
    def add_motif_annotation_homer(self):
        """
        Add motif annotations (based on Homer, cisTarget and meme if specified)
        """
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('Homer')
        
        if self.known_motifs is not None:
            if self.known_motifs.shape[0] != 0:
                log.info('Annotating known motifs')
                # Prepare cistarget annotation
                if 'mm' in self.genome:
                    species = 'mus_musculus'
                if 'dm' in self.genome:
                    species = 'drosophila_melanogaster'
                if 'hg' in self.genome:
                    species = 'homo_sapiens'
                ctx_motif_annotation = load_motif_annotations(species,
                                         motif_similarity_fdr= self.motif_similarity_fdr,
                                          orthologous_identity_threshold=self.orthologous_identity_threshold)
                motifs = self.known_motifs
                homer_motifs = 'homer__' + motifs['Consensus'] + '_' + [x.split('(')[0] for x in motifs['Motif Name']]

                motifs['MotifID'] = homer_motifs
                homer_motifs = [x for x in homer_motifs if x in ctx_motif_annotation.index.tolist()]
                ctx_motif_annotation = ctx_motif_annotation.loc[list(set(homer_motifs))].reset_index()

                # Prepare homer annotation
                try:
                    homer_motif_annotation = pd.read_csv(os.path.abspath(os.path.join(os.path.dirname(self.homer_path), '..', 'motifs/extras/table.txt')),
                                                    sep='\t', error_bad_lines=False).iloc[:,[1,11]].dropna()
                except:
                    homer_motif_annotation = pd.read_csv(os.path.abspath(os.path.join(os.path.dirname(self.homer_path), '..', 'motifs/extras/motifTable.txt')),
                                                    sep='\t', error_bad_lines=False).iloc[:,[1,11]].dropna()
                    homer_motif_annotation.columns = ['Name', 'Symbol'] 
                # If not human, convert by homology
                if species != 'homo_sapiens':
                    dataset = Dataset(name='hsapiens_gene_ensembl',
                                  host='http://www.ensembl.org')
                    if species == 'mus_musculus':
                        biomart_query = 'mmusculus_homolog_associated_gene_name'
                    if species == 'drosophila_melanogaster':
                        biomart_query = 'dmelanogaster_homolog_associated_gene_name'

                    human2specie = dataset.query(attributes=['external_gene_name', biomart_query])
                    human2specie.index = human2specie['Gene name']
                    # Check that the TF has homolog
                    TF_names = [x for x in homer_motif_annotation.iloc[:,1].tolist() if x in human2specie.index.tolist()]
                    human2specie = human2specie.loc[TF_names,:]
                    human2specie.columns = ['Symbol', 'Homolog']
                    df = pd.merge(homer_motif_annotation, human2specie, on='Symbol', how='left')
                    homer_motif_annotation = df.iloc[:,[0,2]]
                    homer_motif_annotation.columns = ['Name', 'Symbol']

                # We first bind the cisTarget annotation
                motifs = pd.merge(motifs, ctx_motif_annotation, on='MotifID', how='left')
                # We now bind the Homer annotation
                homer_motif_annotation.columns = ['Motif Name', 'Homer_annot']
                motifs = pd.merge(motifs, homer_motif_annotation, on='Motif Name', how='left')
                # If Homer_annot is not in Direct_annot we will add it
                # Concatenate
                motifs.Direct_annot = [str(motifs.Direct_annot.tolist()[x]) + ', ' + str(motifs.Homer_annot.tolist()[x]) 
                                   if (str(motifs.Homer_annot.tolist()[x]) not in str(motifs.Direct_annot.tolist()[x]))
                                   else motifs.Direct_annot.tolist()[x] for x in range(motifs.shape[0])]
                motifs.Direct_annot = motifs.Direct_annot.replace('nan, ', '', regex=True)
                motifs.Direct_annot = motifs.Direct_annot.replace(', nan', '', regex=True)
                motifs = motifs.drop(['MotifID', 'Homer_annot'], axis=1)
                self.known_motifs = motifs
        
        if self.denovo_motifs is not None:
            if self.denovo_motifs.shape[0] != 0:
                if self.meme_path is None:
                    log.info('Parameter meme_path is not provided. Skipping annotation of de novo motifs')
                elif self.meme_collection_path is None:
                    log.info('Parameter meme_collection_path is not provided. Skipping annotation of de novo motifs')
                else:
                    # Find closest match for denovo motifs in cistarget database (as meme)
                    log.info('Comparing de novo motifs with given motif collection with tomtom')
                    homer_motif_paths = glob.glob(os.path.join(self.outdir, 'homerResults', 'motif*[0-9$].motif'))
                    homer_motif_paths = [x for x in homer_motif_paths if 'similar' not in x] 
                    tomtom_pd = pd.concat([tomtom(x, self.meme_path, self.meme_collection_path) for x in homer_motif_paths])
                    ctx_motif_annotation = load_motif_annotations(species,
                                         motif_similarity_fdr= self.motif_similarity_fdr,
                                          orthologous_identity_threshold=self.orthologous_identity_threshold)
                    homer_motifs = [x for x in tomtom_pd.iloc[:,1].tolist() if x in ctx_motif_annotation.index.tolist()]
                    ctx_motif_annotation = ctx_motif_annotation.loc[list(set(homer_motifs))].reset_index()
                    ctx_motif_annotation = ctx_motif_annotation.rename(columns={'MotifID': 'Best Match/Tomtom'})
                    # Bind cisTarget annotation
                    tomtom_pd = pd.merge(tomtom_pd, ctx_motif_annotation, on='Best Match/Tomtom', how='left')
                    motifs = pd.merge(self.denovo_motifs, tomtom_pd, on='Best Match/Details', how='left')
                    self.denovo_motifs = motifs

    def find_motif_hits(self, n_cpu=1):
        """
        Find motif hits with `homer2 find`
        
        Parameters
        ---------
        n_cpu: int
            Number of cores to use.
        """
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('Homer')
        
        if self.known_motifs is not None:
            if self.known_motifs.shape[0] != 0:
                # Merge all motifs to file
                log.info('Retrieving enriched regions per known motif')
                if os.path.exists(os.path.join(self.outdir, 'knownResults', 'all_motifs.motif')):
                    os.remove(os.path.join(self.outdir, 'knownResults', 'all_motifs.motif'))
                for f in glob.glob(os.path.join(self.outdir, 'knownResults', '*.motif')):
                    os.system("cat "+f+" >> "+os.path.join(self.outdir, 'knownResults', 'all_motifs.motif'))
                cmd = os.path.join(self.homer_path, 'homer2 find') + ' -s %s -m %s -o %s -p %s'
                cmd = cmd % (os.path.join(self.outdir, 'targetgiven.seq'), os.path.join(self.outdir, 'knownResults', 'all_motifs.motif'), os.path.join(self.outdir, 'knownResults_motif_hits.bed'), n_cpu)
                try:
                    subprocess.check_output(args=cmd, shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

                known_motif_hits = pd.read_csv(os.path.join(self.outdir, 'knownResults_motif_hits.bed'), sep='\t', header=None)
                self.known_motif_hits = known_motif_hits.groupby(3)[0].apply(lambda g: list(set(g.values.tolist()))).to_dict()
        if self.denovo_motifs is not None:
            if self.denovo_motifs.shape[0] != 0:
                # Merge all motifs to file
                log.info('Retrieving enriched regions per de novo motif')
                if os.path.exists(os.path.join(self.outdir, 'homerResults', 'all_motifs.motif')):
                    os.remove(os.path.join(self.outdir, 'homerResults', 'all_motifs.motif'))
                for f in glob.glob(os.path.join(self.outdir, 'homerResults', '*.motif')):
                    os.system("cat "+f+" >> "+os.path.join(self.outdir, 'homerResults', 'all_motifs.motif'))
                os.system("sed -i 's/\t.*BestGuess:/\t/g' "+os.path.join(self.outdir, 'homerResults', 'all_motifs.motif'))
                cmd = os.path.join(self.homer_path, 'homer2 find') + ' -s %s -m %s -o %s -p %s'
                cmd = cmd % (os.path.join(self.outdir, 'targetgiven.seq'), os.path.join(self.outdir, 'homerResults', 'all_motifs.motif'), os.path.join(self.outdir, 'homerResults_motif_hits.bed'), n_cpu)
                try:
                    subprocess.check_output(args=cmd, shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

                denovo_motif_hits = pd.read_csv(os.path.join(self.outdir, 'homerResults_motif_hits.bed'), sep='\t', header=None)
                denovo_motif_hits  = denovo_motif_hits.groupby(3)[0].apply(lambda g: list(set(g.values.tolist()))).to_dict()
                self.denovo_motif_hits = {k:denovo_motif_hits[k] for k in denovo_motif_hits.keys() if not k[0].isdigit()}
                
    def get_cistromes(self, annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']):
        """
        Format cistromes per TF
        
        Parameters
        ---------
        cistrome_annotation: List, optional
            Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
            linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
            (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
            Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']

        """
        if self.known_motif_hits is not None:
            if 'Direct_annot' in annotation:
                tfs = get_TF_list(self.known_motifs, annotation=['Direct_annot'])
                cistrome_dict_direct = {tf: get_cistrome_per_TF(self.known_motif_hits,  get_motifs_per_TF(self.known_motifs, tf, motif_column = 'Motif Name', annotation=['Direct_annot'])) for tf in tfs}
            else:
                cistrome_dict_direct = {}
                
            if not 'Direct_annot' in annotation or len(annotation) > 1:
                tfs = get_TF_list(self.known_motifs, annotation=annotation)
                cistrome_dict_extended = {tf+'_extended': get_cistrome_per_TF(self.known_motif_hits,  get_motifs_per_TF(self.known_motifs, tf, motif_column = 'Motif Name', annotation=annotation)) for tf in tfs}
            else:
                cistrome_dict_extended = {}
                
            cistrome_dict = {**cistrome_dict_direct, **cistrome_dict_extended}
            cistrome_dict = {x + ' (' + str(len(cistrome_dict[x])) + 'r)': cistrome_dict[x] for x in cistrome_dict.keys()}
            self.known_cistromes = cistrome_dict 
        if self.denovo_motif_hits is not None:
            if 'Direct_annot' in annotation:
                tfs = get_TF_list(self.denovo_motifs, annotation=['Direct_annot'])
                cistrome_dict_direct = {tf: get_cistrome_per_TF(self.denovo_motif_hits,  get_motifs_per_TF(self.denovo_motifs, tf, motif_column = 'Motif Name', annotation=['Direct_annot'])) for tf in tfs}
            else:
                cistrome_dict_direct = {}
                
            if not 'Direct_annot' in annotation or len(annotation) > 1:
                tfs = get_TF_list(self.denovo_motifs, annotation=annotation)
                cistrome_dict_extended = {tf+'_extended': get_cistrome_per_TF(self.denovo_motif_hits,  get_motifs_per_TF(self.denovo_motifs, tf, motif_column = 'Motif Name', annotation=annotation)) for tf in tfs}
            else:
                cistrome_dict_extended = {}
                
            cistrome_dict = {**cistrome_dict_direct, **cistrome_dict_extended}
            cistrome_dict = {x + '_(' + str(len(cistrome_dict[x])) + 'r)': cistrome_dict[x] for x in cistrome_dict.keys()}
            self.denovo_cistromes = cistrome_dict
 
# Run Homer           
def run_homer(homer_path: str,
                             region_sets: Dict[str, pr.PyRanges],
                             outdir: str,
                             genome: str,
                             size: str = 'given',
                             mask: bool = True,
                             denovo: bool = False,
                             length: str = '8,10,12',
                             n_cpu: int = 1,
                             meme_path: str = None,
                             meme_collection_path: str = None,
                             cistrome_annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
                             motif_similarity_fdr: float = 0.001,
                             orthologous_identity_threshold: float = 0.0,
                             **kwargs):
    """
    Run Homer 

    Parameters
    ---------
    homer_path: str
        Path to Homer bin folder.
    region_sets: Dict
        A dictionary of PyRanges containing region coordinates for the region sets to be analyzed.
    outdir: str
        Path to folder to output Homer results.
    genome: str
        Homer genome label to use.
    size: str, optional
        Fragment size to use for motif finding. Default: 'given' [uses the exact regions you give it]
    mask: bool, optional
        Whether to mask repeats or not. Default: True
    denovo : bool, optional
        Whether to infer overrepresented motifs de novo. Default: False
    length: str, optional
        Motif length values. Default: 8,10,12
    n_cpu: int
        Number of cores to use.
    meme_path: str, optional
        Path to meme bin folder. Meme will be used if given for motif annotation. Default: None
    meme_collection_path: str, optional
        Path to motif collection (in .cb format) to compare homer motifs with. Default: None
    cistrome_annotation: List, optional
        Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
        linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
        (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
        Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']
    motif_similarity_fdr: float, optional
        Minimal motif similarity value to consider two motifs similar. Default: 0.001
    orthologous_identity_threshold: float, optional
        Minimal orthology value for considering two TFs orthologous. Default: 0.0
    **kwargs
        Extra parameters to pass to `ray.init()`.

    References
    ---------
    Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining 
    Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities.
    Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
    """
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('Homer')
    # Save regions in dict to the output dir
    bed_paths={}
    bed_dir = os.path.join(outdir, 'regions_bed')
    # Create bed directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(bed_dir):
        os.mkdir(bed_dir)
    # Create bed files for Homer
    for key in region_sets.keys():
        bed_path = os.path.join(bed_dir, key+'.bed')
        region_sets[key].Name =  coord_to_region_names(region_sets[key])
        region_sets[key].to_bed(path=bed_path, keep=False, compression='infer', chain=False)
        bed_paths[key] = bed_path
    # Run Homer
    ray.init(num_cpus=n_cpu, **kwargs)
    sys.stderr = null
    homer_dict = ray.get([homer_ray.remote(homer_path,
                                bed_paths[name],
                                name,
                                outdir + name, 
                                genome,
                                size,
                                mask,
                                denovo,
                                length, 
                                meme_path,
                                meme_collection_path,
                                cistrome_annotation,
                                motif_similarity_fdr,
                                orthologous_identity_threshold) for name in list(bed_paths.keys())])
    ray.shutdown()
    sys.stderr = sys.__stderr__
    homer_dict={list(bed_paths.keys())[i]: homer_dict[i] for i in range(len(homer_dict))}
    return homer_dict

@ray.remote
def homer_ray(homer_path: str,
              bed_path: str,
              name: str,
              outdir: str,
              genome: str,
              size: str = 'given',
              mask: bool = True,
              denovo: bool = False,
              length: str = '8,10,12',
              meme_path: str = None,
              meme_collection_path: str = None,
              cistrome_annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
              motif_similarity_fdr: float = 0.001,
              orthologous_identity_threshold: float = 0.0):
    """
    Ray method to run Homer.
    """
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('Homer')
    
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)
    
    log.info('Running '+ name)
    Homer_res = Homer(homer_path,
                bed_path,
                name,
                outdir, 
                genome, 
                size, 
                mask, 
                denovo, 
                length, 
                meme_path, 
                meme_collection_path, 
                cistrome_annotation, 
                motif_similarity_fdr)
    Homer_res.run()
    log.info(name + ' done!')
    return Homer_res
            
# Utils
## Show results
def homer_results(homer_dict, name, results='known'):
    """
    A function to show Homer results in jupyter notebooks.
    
    Parameters
    ---------
    Homer_dict: Dict
        A dictionary with one :class:`Homer` object per slot.
    name: str
        Dictionary key of the analysis result to show. Default: None (All)
    results: str
        Whether to show know or de novo results. Default: 'known'
    """
    if results == 'known':
        file = os.path.join(homer_dict[name].outdir, 'knownResults.html')
    if results == 'denovo':
        file = os.path.join(homer_dict[name].outdir, 'homerResults.html')
    inplace_change(file, 'width="505" height="50"', 'width="1010" height="200"')
    return HTML(file)
                