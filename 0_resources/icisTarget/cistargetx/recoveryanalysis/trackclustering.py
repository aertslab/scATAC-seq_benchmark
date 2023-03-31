import logging
import operator


def get_track_clusters(track_ids, track_cluster_filename):
    found_tracks = set()
    track_clusters = []
    clust_number = None
    track_clusters_tmp = []

    with open(track_cluster_filename, 'r') as input_fh:
        for line in input_fh:
            line = line.split()
            clust_number_tmp = line[1]
            cur_track = line[0]

            if cur_track not in track_ids:
                continue
            else:
                found_tracks.add(cur_track)

            if clust_number == clust_number_tmp:
                track_clusters_tmp.append(line[0])
            else:
                if track_clusters_tmp != []:
                    track_clusters.append(track_clusters_tmp)

                track_clusters_tmp = []
                track_clusters_tmp.append(line[0])
                clust_number = clust_number_tmp

        if track_clusters_tmp != []:
            track_clusters.append(track_clusters_tmp)

    if len(found_tracks) != len(track_ids):
        logging.warning("CLUSTERING: only {0:d} out of {1:d} tracks were found.".format(len(found_tracks),
                                                                                        len(track_ids)))

        return track_clusters

    return track_clusters


def track_cluster(track_cluster_filename, track_names):
    # Create subset of track database with only the enriched tracks.
    track_clusters = get_track_clusters(track_names, track_cluster_filename)

    if len(track_clusters) < 1:
        logging.warning("Clustering: no tracks where found so skipping clustering.")

        return [], None

    return track_clusters, None


def track_cluster_table_generator(ranked_tracks, track_clusters):
    for idx, cluster in enumerate(track_clusters):
        ranks = (ranked_tracks.index(track) + 1 for track in cluster)

        for rank, track in sorted(zip(ranks, cluster), key=operator.itemgetter(0)):
            yield idx + 1, rank, track
