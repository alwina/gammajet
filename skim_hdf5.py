import h5py
import numpy as np


def skimHDF5Centrality(filename, centrange, eventgroupsize=100000, chunksize=2000):
    centrality_index = 3
    fullfile = h5py.File(filename, 'r')
    skimfile = h5py.File(filename.replace('.hdf5', '_skimcent{0}{1}.hdf5'.format(*centrange)), 'w')

    eventshape = fullfile['event'].shape
    clustershape = fullfile['cluster'][:1].shape
    trackshape = fullfile['track'][:1].shape
    jetshape = fullfile['jet'][:1].shape
    neventsfull = eventshape[0]

    print('Skimming from {0} events'.format(neventsfull))

    skimevent = None
    skimcluster = None
    skimtrack = None
    skimjet = None
    neventsskimtotal = 0

    for startevent in range(0, neventsfull, eventgroupsize):
        event = fullfile['event'][startevent:min(startevent + eventgroupsize, neventsfull)]
        cluster = fullfile['cluster'][startevent:min(startevent + eventgroupsize, neventsfull)]
        track = fullfile['track'][startevent:min(startevent + eventgroupsize, neventsfull)]
        jet = fullfile['jet'][startevent:min(startevent + eventgroupsize, neventsfull)]

        event_pass_indices = np.logical_and(centrange[0] < event[:, centrality_index], event[:, centrality_index] < centrange[1])
        neventsskim = np.count_nonzero(event_pass_indices)

        if skimevent is None:
            skimevent = skimfile.create_dataset('event', data=event[event_pass_indices], chunks=(chunksize, eventshape[1]), compression='gzip', maxshape=(neventsfull, eventshape[1]))
            skimcluster = skimfile.create_dataset('cluster', data=cluster[event_pass_indices], chunks=(chunksize, clustershape[1], clustershape[2]), compression='gzip', maxshape=(neventsfull, clustershape[1], clustershape[2]))
            skimtrack = skimfile.create_dataset('track', data=track[event_pass_indices], chunks=(chunksize, trackshape[1], trackshape[2]), compression='gzip', maxshape=(neventsfull, trackshape[1], trackshape[2]))
            skimjet = skimfile.create_dataset('jet', data=jet[event_pass_indices], chunks=(chunksize, jetshape[1], jetshape[2]), compression='gzip', maxshape=(neventsfull, jetshape[1], jetshape[2]))
        else:
            skimevent.resize(skimevent.shape[0] + neventsskim, axis=0)
            skimcluster.resize(skimcluster.shape[0] + neventsskim, axis=0)
            skimtrack.resize(skimtrack.shape[0] + neventsskim, axis=0)
            skimjet.resize(skimjet.shape[0] + neventsskim, axis=0)

            skimevent[-neventsskim:] = event[event_pass_indices]
            skimcluster[-neventsskim:] = cluster[event_pass_indices]
            skimtrack[-neventsskim:] = track[event_pass_indices]
            skimjet[-neventsskim:] = jet[event_pass_indices]

        print('Keeping {0}/{1} events'.format(neventsskim, event.shape[0]))
        neventsskimtotal = neventsskimtotal + neventsskim

    print('Kept {0}/{1} events'.format(neventsskimtotal, neventsfull))

    skimfile.close()
    fullfile.close()
