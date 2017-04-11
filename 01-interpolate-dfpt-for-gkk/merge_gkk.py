from merge_gkk_nc import merge_gkk_nc

nqpt = 29
npert = 6

for iqpt in range(nqpt):

    out_fname = 'GkkInterp-g2-g4/GKK/qpt-{:0=4}/GKK/out_data/odat_GKK.nc'.format(iqpt+1)

    fnames = list()
    for ipert in range(npert):
        fnames.append('GkkInterp-g2-g4/GKK/qpt-{:0=4}/GKK/out_data/odat_GKK{}.nc'.format(iqpt+1, ipert+1))

    merge_gkk_nc(out_fname, fnames)
