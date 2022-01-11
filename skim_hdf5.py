import h5py
import numpy as np


def skimHDF5Centrality(filename, centrange, eventgroupsize=100000, chunksize=2000):
    centrality_index = 3
    fullfile = h5py.File(filename, 'r')
    skimfile = h5py.File(filename.replace('.hdf5', '_skimcent{0}{1}.hdf5'.format(*centrange)), 'w')

    eventshape = fullfile['event'].shape
    clustershape = fullfile['cluster'][:1].shape
    trackshape = fullfile['track'][:1].shape
    jet04tpcshape = fullfile['jet_ak04tpc'][:1].shape
    jet02tpcshape = fullfile['jet_ak02tpc'][:1].shape
    neventsfull = eventshape[0]

    print('Skimming from {0} events'.format(neventsfull))

    skimevent = None
    skimcluster = None
    skimtrack = None
    skimjet04tpc = None
    skimjet02tpc = None
    neventsskimtotal = 0

    for startevent in range(0, neventsfull, eventgroupsize):
        event = fullfile['event'][startevent:min(startevent + eventgroupsize, neventsfull)]
        cluster = fullfile['cluster'][startevent:min(startevent + eventgroupsize, neventsfull)]
        track = fullfile['track'][startevent:min(startevent + eventgroupsize, neventsfull)]
        jet04tpc = fullfile['jet_ak04tpc'][startevent:min(startevent + eventgroupsize, neventsfull)]
        jet02tpc = fullfile['jet_ak02tpc'][startevent:min(startevent + eventgroupsize, neventsfull)]

        event_pass_indices = np.logical_and(centrange[0] <= event[:, centrality_index], event[:, centrality_index] < centrange[1])
        neventsskim = np.count_nonzero(event_pass_indices)

        if skimevent is None:
            skimevent = skimfile.create_dataset('event', data=event[event_pass_indices], chunks=(chunksize, eventshape[1]), compression='gzip', maxshape=(neventsfull, eventshape[1]))
            skimcluster = skimfile.create_dataset('cluster', data=cluster[event_pass_indices], chunks=(chunksize, clustershape[1], clustershape[2]), compression='gzip', maxshape=(neventsfull, clustershape[1], clustershape[2]))
            skimtrack = skimfile.create_dataset('track', data=track[event_pass_indices], chunks=(chunksize, trackshape[1], trackshape[2]), compression='gzip', maxshape=(neventsfull, trackshape[1], trackshape[2]))
            skimjet04tpc = skimfile.create_dataset('jet_ak04tpc', data=jet04tpc[event_pass_indices], chunks=(chunksize, jet04tpcshape[1], jet04tpcshape[2]), compression='gzip', maxshape=(neventsfull, jet04tpcshape[1], jet04tpcshape[2]))
            skimjet02tpc = skimfile.create_dataset('jet_ak02tpc', data=jet02tpc[event_pass_indices], chunks=(chunksize, jet02tpcshape[1], jet02tpcshape[2]), compression='gzip', maxshape=(neventsfull, jet02tpcshape[1], jet02tpcshape[2]))
        else:
            skimevent.resize(skimevent.shape[0] + neventsskim, axis=0)
            skimcluster.resize(skimcluster.shape[0] + neventsskim, axis=0)
            skimtrack.resize(skimtrack.shape[0] + neventsskim, axis=0)
            skimjet04tpc.resize(skimjet04tpc.shape[0] + neventsskim, axis=0)
            skimjet02tpc.resize(skimjet02tpc.shape[0] + neventsskim, axis=0)

            skimevent[-neventsskim:] = event[event_pass_indices]
            skimcluster[-neventsskim:] = cluster[event_pass_indices]
            skimtrack[-neventsskim:] = track[event_pass_indices]
            skimjet04tpc[-neventsskim:] = jet04tpc[event_pass_indices]
            skimjet02tpc[-neventsskim:] = jet02tpc[event_pass_indices]

        print('Keeping {0}/{1} events'.format(neventsskim, event.shape[0]))
        neventsskimtotal = neventsskimtotal + neventsskim

    print('Kept {0}/{1} events'.format(neventsskimtotal, neventsfull))

    skimfile.close()
    fullfile.close()
