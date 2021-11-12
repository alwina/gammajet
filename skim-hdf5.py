import h5py
import numpy as np


def skimHDF5Centrality(filename, centrange, eventgroupsize=100000, chunksize=2000):
    centrality_index = 3
    fullfile = h5py.File(filename, 'r')

    eventshape = fullfile['event'].shape
    clustershape = fullfile['cluster'][:1].shape
    trackshape = fullfile['track'][:1].shape
    jetshape = fullfile['jet'][:1].shape
    neventsfull = eventshape[0]

    skimevents = []
    skimclusters = []
    skimtracks = []
    skimjets = []

    for startevent in range(0, neventsfull, eventgroupsize):
        event = fullfile['event'][startevent:min(startevent + eventgroupsize, neventsfull)]
        cluster = fullfile['cluster'][startevent:min(startevent + eventgroupsize, neventsfull)]
        track = fullfile['track'][startevent:min(startevent + eventgroupsize, neventsfull)]
        jet = fullfile['jet'][startevent:min(startevent + eventgroupsize, neventsfull)]

        event_pass_indices = np.logical_and(centrange[0] < event[:, centrality_index], event[:, centrality_index] < centrange[1])

        skimevents.append(event[event_pass_indices])
        skimclusters.append(cluster[event_pass_indices])
        skimtracks.append(track[event_pass_indices])
        skimjets.append(jet[event_pass_indices])

        print('Keeping {0}/{1} events'.format(skimevents[-1].shape[0], event.shape[0]))

    skimfile = h5py.File(filename.replace('.hdf5', '_skimcent{0}{1}.hdf5'.format(*centrange)), 'w')
    skimfile.create_dataset('event', data=np.concatenate(skimevents), chunks=(chunksize, eventshape[1]), compression='gzip')
    skimfile.create_dataset('cluster', data=np.concatenate(skimclusters), chunks=(chunksize, clustershape[1], clustershape[2]), compression='gzip')
    skimfile.create_dataset('track', data=np.concatenate(skimtracks), chunks=(chunksize, trackshape[1], trackshape[2]), compression='gzip')
    skimfile.create_dataset('jet', data=np.concatenate(skimjets), chunks=(chunksize, jetshape[1], jetshape[2]), compression='gzip')

    neventsskim = np.concatenate(skimevents).shape[0]
    print('Kept {0}/{1} events'.format(neventsskim, neventsfull))

    skimfile.close()
    fullfile.close()
