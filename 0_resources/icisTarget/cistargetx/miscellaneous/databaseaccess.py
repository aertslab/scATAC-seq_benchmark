import sqlite3, operator, numpy, os
import cistargetx.common.utils as utils


CREATE_TABLE_STATEMENTS = r"""
DROP TABLE IF EXISTS rankings;
DROP TABLE IF EXISTS motifs;
CREATE TABLE rankings (geneID VARCHAR(255), ranking BLOB);
CREATE TABLE motifs (motifName VARCHAR(255), idx INTEGER); """
CREATE_INDEX_STATEMENT   = """CREATE UNIQUE INDEX id ON rankings (geneID);"""
INSERT_MOTIF_STATEMENT   = "INSERT INTO motifs (idx, motifName) VALUES (?, ?)"
INSERT_RANKING_STATEMENT = "INSERT INTO rankings (geneID, ranking) VALUES (?, ?)"
GENECOUNT_QUERY          = r"SELECT COUNT(*) FROM rankings"
GENESET_QUERY            = r"SELECT geneID FROM rankings WHERE geneID IN ({0:s}) ORDER BY geneID"
ALLGENESET_QUERY         = r"SELECT geneID FROM rankings ORDER BY geneID"
RANKINGS_QUERY           = r"SELECT geneID, ranking FROM rankings WHERE geneID IN ({0:s}) ORDER BY geneID"
ALLRANKINGS_QUERY        = r"SELECT geneID, ranking FROM rankings ORDER BY geneID"
MOTIFS_QUERY             = r"SELECT motifName FROM motifs ORDER BY idx"


def fetch_motifs(database_filename):
    '''Fetch motifs from SQLite3 database with supplied filename.'''
    if not os.path.exists(database_filename):
        raise ValueError("Database {0:s} doesn't exist.".format(database_filename))
    with sqlite3.connect(database_filename) as db:
        cursor = db.cursor()
        totalmotifs = map(operator.itemgetter(0), cursor.execute(MOTIFS_QUERY).fetchall())
        cursor.close()
    return totalmotifs


def _quoted_csv(values):
    # Use double quotes, because sometimes ID's contain a single quote ...
    def quote(value): return '"' + value + '"'
    return ','.join(map(quote, values))


def fetch_rankings(database_filename, geneids=None, motifs=None):
    '''  Fetch rankings from SQLite3 database with supplied filename.
    0-based rankings are returned.
    '''
    if not os.path.exists(database_filename):
        raise ValueError("Database {0:s} doesn't exist.".format(database_filename))
    
    with sqlite3.connect(database_filename) as db:
        cursor = db.cursor()
        
        # Get motif information ...
        totalmotifs = map(operator.itemgetter(0), cursor.execute(MOTIFS_QUERY).fetchall())
        totalmotifcount = len(totalmotifs)
        if motifs:
            selectedmotifs = [ (totalmotifs.index(motif), motif) for motif in motifs if motif in totalmotifs ]
            if selectedmotifs: motifidxs, motifs = zip(*selectedmotifs)
            else: motifidxs, motifs = None, []
        else: motifidxs, motifs = None, totalmotifs
        motifs = list(motifs)
        motifcount = len(motifs)
        
        # Get total geneset information ...
        totalgenecount = int(cursor.execute(GENECOUNT_QUERY).fetchone()[0])
        dtype = utils.derive_dtype(totalgenecount)
        if geneids:
            cursor.execute(GENESET_QUERY.format(_quoted_csv(geneids)))
            geneids = set(id for id, in cursor)
            genecount = len(geneids)
        else:
            geneids = []
            genecount = totalgenecount  

        # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
        rankings = numpy.full(shape=(genecount, motifcount), fill_value=0, dtype=dtype)
        
        if motifcount == 0:
            if genecount == totalgenecount and len(geneids) == 0:
                cursor.execute(ALLGENESET_QUERY)
                geneids = [id for id, in cursor]
            return motifs, numpy.array(geneids, dtype='|S255'), rankings, totalgenecount, dtype
        
        # Load all requested rankings ...
        if geneids:
            cursor.execute(RANKINGS_QUERY.format(_quoted_csv(geneids)))
        else:
            cursor.execute(ALLRANKINGS_QUERY)
        rowidx = 0; geneids = []
        if motifidxs:
            for ID, ranking in cursor:
                rankings[rowidx, :] = numpy.frombuffer(ranking, dtype=dtype)[motifidxs, :]
                geneids.append(ID)
                rowidx += 1
        else:
            for ID, ranking in cursor:
                rankings[rowidx, :] = numpy.frombuffer(ranking, dtype=dtype)
                geneids.append(ID)
                rowidx += 1
        
        cursor.close()
    return motifs, numpy.array(geneids, dtype='|S255'), rankings, totalgenecount, dtype


def fetch_all_rankings(database_filename):
    '''  Fetch rankings from SQLite3 database with supplied filename.
    0-based rankings are returned.
    '''
    return fetch_rankings(database_filename)


def write_rankings(motifs, geneids, rankings, database_filename):
    ''' Write 0-based rankings to SQLite3 database. '''
    with sqlite3.connect(database_filename) as db:
        db.text_factory = str
        cursor = db.cursor()
        cursor.executescript(CREATE_TABLE_STATEMENTS)
        for idx, motifname in enumerate(motifs):
            cursor.execute(INSERT_MOTIF_STATEMENT, (idx, motifname) )
        for rowidx, geneid in enumerate(geneids):
            cursor.execute(INSERT_RANKING_STATEMENT, (geneid, buffer(rankings[rowidx, :])) )
        cursor.execute(CREATE_INDEX_STATEMENT)
        cursor.close()


#def fetch_rankings(database_filename, geneids):
#    # Returns rankings 0-based ...
#    if not os.path.exists(database_filename):
#        raise ValueError("Database {0:s} doesn't exist.".format(database_filename))
#    with sqlite3.connect(database_filename) as db:
#        cursor = db.cursor()
#        # Get motif information ...
#        motifs = map(operator.itemgetter(0), cursor.execute(MOTIFS_QUERY).fetchall())
#        motifcount = len(motifs)
#        # Get total geneset information ...
#        totalgenecount = int(cursor.execute(GENECOUNT_QUERY).fetchone()[0])
#        dtype = utils.derive_dtype(totalgenecount)
#        # Query geneset information ...
#        cursor.execute(GENESET_QUERY.format(_quoted_csv(geneids)))
#        geneids = set(id for id, in cursor)
#        # Pre-allocation ...
#        rankings = numpy.empty(shape=(len(geneids), motifcount), dtype=dtype)
#        # Load all requested rankings ...
#        cursor.execute(RANKINGS_QUERY.format(_quoted_csv(geneids)))
#        idx = 0; geneids = []
#        for ID, ranking in cursor:
#            rankings[idx, :] = numpy.frombuffer(ranking, dtype=dtype)
#            geneids.append(ID)
#            idx += 1
#        cursor.close()
#    return motifs, numpy.array(geneids, dtype='|S255'), rankings, totalgenecount, dtype
