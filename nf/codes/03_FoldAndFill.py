# +
# :Description:
#
# :AUTHOR: Katarzyna Ewa Lewinska
# :DATE: 17 March 2022
# :VERSION: 1.0     based on Fold2Days
# :UPDATES:         2023-03-17  Monthly composites
#                   2023-03-20  Filling in from RBF output 16 days
#                   2023-03-22  Filling in from RBF output 24 days for wider gaps signa 16, 48, 96
#                   2023-04-21  pass parameters from the command line; iterate through all tiles in the cube
# -


# :ENVIRONMENT: #

from osgeo import gdal
import os
import pandas as pd
import numpy as np
import glob
import rasterio
from datetime import datetime
from datetime import timedelta
from dateutil.relativedelta import relativedelta
import warnings
import sys

startTime = datetime.now()
# :INPUTS: #
if sys.argv[1] is None:
#    cube = r'J:\fonda\grassdata\force\unm'  # parent directory
  print('No input dir')
else:
    cube = (sys.argv[1])
    print(cube)

cubedir = list(os.listdir(os.path.join(cube, 'gv')))
tiles = [x for x in cubedir if x.startswith('X')]

os.chdir(cube)  # parent directory
cdir = os.getcwd()

gvdir = 'gv'
npvdir = 'npv'
soildir = 'soil'
shadedir = 'shade'

gvdirW = 'gv_wide'
npvdirW = 'npv_wide'
soildirW = 'soil_wide'
shadedirW = 'shade_wide'

foldMonths = 1      # number of months to base the composite on

# :CODE: #

# with warnings.catch_warnings():
    # warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')
    # warnings.filterwarnings(action='ignore', message='invalid value encountered in cast')
warnings.simplefilter("ignore", category=RuntimeWarning)

nos = len(tiles)

print(tiles)
print(nos)

# iterate through tiles
for i in range(0, nos):
    # itile = tiles.iat[i, 0]
    itile = tiles[i]
    print(itile)

    # get paths to endmembers
    gvtile = os.path.join(cdir, gvdir, itile)
    npvtile = os.path.join(cdir, npvdir, itile)
    soiltile = os.path.join(cdir, soildir, itile)
    shtile = os.path.join(cdir, shadedir, itile)

    # get paths to wide TSIs
    gvtileW = os.path.join(cdir, gvdirW, itile)
    npvtileW = os.path.join(cdir, npvdirW, itile)
    soiltileW = os.path.join(cdir, soildirW, itile)
    shtileW = os.path.join(cdir, shadedirW, itile)

    # unmixing data
    rmse = glob.glob(gvtile+'/*SMA_RMS.tif')
    gvtss = glob.glob(gvtile+'/*TSS.tif')
    npvtss = glob.glob(npvtile+'/*TSS.tif')
    soiltss = glob.glob(soiltile + '/*TSS.tif')
    shtss = glob.glob(shtile + '/*TSS.tif')

    # 16-day RBF interpolation
    gvrbfp = glob.glob(gvtile+'/*TSI.tif')
    npvrbfp = glob.glob(npvtile + '/*TSI.tif')
    soilrbfp = glob.glob(soiltile + '/*TSI.tif')
    shrbfp = glob.glob(shtile + '/*TSI.tif')

    # 32-day RBF interpolation
    gvrbfpW = glob.glob(gvtileW+'/*TSI.tif')
    npvrbfpW = glob.glob(npvtileW + '/*TSI.tif')
    soilrbfpW = glob.glob(soiltileW + '/*TSI.tif')
    shrbfpW = glob.glob(shtileW + '/*TSI.tif')
    # sosStak = rasterio.open(rmse[0], "r", tiled=True, BIGTIFF='YES')
    # sos = sosStak.read()
    # sosStak.descriptions    # band names

    rms = gdal.Open(rmse[0])
    gv = gdal.Open(gvtss[0])
    npv = gdal.Open(npvtss[0])
    soil = gdal.Open(soiltss[0])
    sh = gdal.Open(shtss[0])

    gvrbf = gdal.Open(gvrbfp[0])
    npvrbf = gdal.Open(npvrbfp[0])
    soilrbf = gdal.Open(soilrbfp[0])
    shrbf = gdal.Open(shrbfp[0])

    gvrbfW = gdal.Open(gvrbfpW[0])
    npvrbfW = gdal.Open(npvrbfpW[0])
    soilrbfW = gdal.Open(soilrbfpW[0])
    shrbfW = gdal.Open(shrbfpW[0])
    # print(gdal.Info(gv))

    # get time info from Band information in the stack
    datDes = []
    for b in range(1, gv.RasterCount+1):
        bt = gv.GetRasterBand(b)
        dest = bt.GetDescription()
        datDes.append(dest)

    datDesRBF = []
    for b2 in range(1, gvrbf.RasterCount+1):
        bt2 = gvrbf.GetRasterBand(b2)
        dest2 = bt2.GetDescription()
        datDesRBF.append(dest2)

    datDesRBFW = []
    for b2W in range(1, gvrbfW.RasterCount+1):
        bt2W = gvrbfW.GetRasterBand(b2W)
        dest2W = bt2W.GetDescription()
        datDesRBFW.append(dest2W)

    # dates = []
    datesYmd = []
    for n in range(0, len(datDes)):
        # dates.append(datDes[n][0:8])
        datesYmd.append(datetime.strptime(datDes[n][0:8], "%Y%m%d"))

    # get dates from RBF layerstack
    datesYmdRBF = []
    for n2 in range(0, len(datDesRBF)):
        datesYmdRBF.append(datetime.strptime(datDesRBF[n2][0:8], "%Y%m%d"))

    # get dates from wide RBF layerstack
    datesYmdRBFW = []
    for n3 in range(0, len(datDesRBFW)):
        datesYmdRBFW.append(datetime.strptime(datDesRBFW[n3][0:8], "%Y%m%d"))

    firstYear = int(datesYmd[0].strftime('%Y'))
    lastYear = int(datesYmd[-1].strftime('%Y'))

    gvLS = np.empty((rms.GetRasterBand(1).ReadAsArray().shape[0],
                    rms.GetRasterBand(1).ReadAsArray().shape[1],1),
                    like = rms.GetRasterBand(1).ReadAsArray())
    npvLS = gvLS.copy()
    soilLS = gvLS.copy()
    shLS = gvLS.copy()
    bandNames = []
    for y in range(firstYear, lastYear+1):

        # iterate through composites in a single year
        for mm in range(0, 12):                # m is not used in the for loop
            fday = datetime(year=y, month=mm+1, day=1)  # the first day in the time series
            lday = (fday + relativedelta(months=foldMonths)) - timedelta(days=1)

            bname = fday.strftime('%Y%m%d')
            bandNames.append(bname)
            # ind = [index for index, adate in enumerate(datesYmd) if fday <= adate <= lday]
            ind = [index for index in range(len(datesYmd)) if fday <= datesYmd[index] <= lday]

            # prepare fill in values
            indRBF = [index for index in range(len(datesYmdRBF)) if fday <= datesYmdRBF[index] <= lday]
            indRBFW = [index for index in range(len(datesYmdRBFW)) if fday <= datesYmdRBFW[index] <= lday]

            empty = np.ones((rms.RasterXSize, rms.RasterYSize, len(indRBF))) * (-9999)
            gvRBFarr = empty.copy()
            npvRBFarr = empty.copy()
            soilRBFarr = empty.copy()
            shRBFarr = empty.copy()

            gvRBFWarr = empty.copy()
            npvRBFWarr = empty.copy()
            soilRBFWarr = empty.copy()
            shRBFWarr = empty.copy()

            for bR, bandR in enumerate(indRBF):  # will use the value passed by ind
                # print('b is:', bR, 'band is: ', bandR)
                gvRBFarr[:, :, bR] = gvrbf.GetRasterBand(bandR + 1).ReadAsArray()
                npvRBFarr[:, :, bR] = npvrbf.GetRasterBand(bandR + 1).ReadAsArray()
                soilRBFarr[:, :, bR] = soilrbf.GetRasterBand(bandR + 1).ReadAsArray()
                shRBFarr[:, :, bR] = shrbf.GetRasterBand(bandR + 1).ReadAsArray()

            for bRW, bandRW in enumerate(indRBFW):  # will use the value passed by ind
                gvRBFWarr[:, :, bRW] = gvrbfW.GetRasterBand(bandRW + 1).ReadAsArray()
                npvRBFWarr[:, :, bRW] = npvrbfW.GetRasterBand(bandRW + 1).ReadAsArray()
                soilRBFWarr[:, :, bRW] = soilrbfW.GetRasterBand(bandRW + 1).ReadAsArray()
                shRBFWarr[:, :, bRW] = shrbfW.GetRasterBand(bandRW + 1).ReadAsArray()

            gvRBFarr[gvRBFarr == -9999] = np.nan
            npvRBFarr[npvRBFarr == -9999] = np.nan
            soilRBFarr[soilRBFarr == -9999] = np.nan
            shRBFarr[shRBFarr == -9999] = np.nan

            gvRBFWarr[gvRBFWarr == -9999] = np.nan
            npvRBFWarr[npvRBFWarr == -9999] = np.nan
            soilRBFWarr[soilRBFWarr == -9999] = np.nan
            shRBFWarr[shRBFWarr == -9999] = np.nan

            gvRBFComp = np.nanmean(gvRBFarr[:, :, :], axis=2)
            npvRBFComp = np.nanmean(npvRBFarr[:, :, :], axis=2)
            soilRBFComp = np.nanmean(soilRBFarr[:, :, :], axis=2)
            shRBFComp = np.nanmean(shRBFarr[:, :, :], axis=2)

            gvRBFWmean = np.nanmean(gvRBFWarr[:, :, :], axis=2)
            npvRBFWmean = np.nanmean(npvRBFWarr[:, :, :], axis=2)
            soilRBFWmean = np.nanmean(soilRBFWarr[:, :, :], axis=2)
            shRBFWmean = np.nanmean(shRBFWarr[:, :, :], axis=2)

            # fill in the 0s and nan in 8-day RBF with 16-day RBF
            gvMa = np.squeeze(np.any([[gvRBFComp == 0], [np.isnan(gvRBFComp)]], axis=0), axis=(0))
            npvMa = np.squeeze(np.any([[npvRBFComp == 0], [np.isnan(npvRBFComp)]], axis=0), axis=(0))
            soilMa = np.squeeze(np.any([[soilRBFComp == 0], [np.isnan(soilRBFComp)]], axis=0), axis=(0))
            shMa = np.squeeze(np.any([[shRBFComp == 0], [np.isnan(shRBFComp)]], axis=0), axis=(0))

            gvRBFComp[gvMa] = gvRBFWmean[gvMa]
            # gvRBFComp[gvRBFComp == 0 | np.isnan(gvRBFComp)] = gvRBFWmean[gvRBFComp == 0 | np.isnan(gvRBFComp)]
            npvRBFComp[npvMa] = npvRBFWmean[npvMa]
            soilRBFComp[soilMa] = soilRBFWmean[soilMa]
            shRBFComp[shMa] = shRBFWmean[shMa]

            if len(ind) > 0:

                if len(ind) > 1:    # exclusive False for len([1])
                    print(y, mm, ': multiple images: look for lowest RMSE')
                    empty = np.ones((rms.RasterXSize, rms.RasterYSize, len(ind))) * (-9999)
                    rmsearr = empty.copy()
                    rmsearrM = empty.copy()
                    gvarr = empty.copy()
                    npvarr = empty.copy()
                    soilarr = empty.copy()
                    sharr = empty.copy()

                    for b, band in enumerate(ind):          # will use the value passed by ind
                        # print('b is:', b, 'band is: ', band)
                        rmsearr[:, :, b] = rms.GetRasterBand(band + 1).ReadAsArray()
                        gvarr[:, :, b] = gv.GetRasterBand(band + 1).ReadAsArray()
                        npvarr[:, :, b] = npv.GetRasterBand(band + 1).ReadAsArray()
                        soilarr[:, :, b] = soil.GetRasterBand(band + 1).ReadAsArray()
                        sharr[:, :, b] = sh.GetRasterBand(band + 1).ReadAsArray()

                    rmsearr[gvarr == -9999] = np.nan
                    rmsearr[rmsearr == -9999] = np.nan
                    minrmse = np.nanmin(rmsearr[:, :, :], axis=2)
                    # minrmse = rmsearr[:, :, :].min(axis=2)

                    for m in range(0, len(ind)):  # 1 change to b
                        rmsearrM[:, :, m] = np.where(rmsearr[:, :, m] == minrmse, 1, 0)

                    gvM = empty.copy()
                    npvM = empty.copy()
                    soilM = empty.copy()
                    shM = empty.copy()
                    for m1 in range(0, len(ind)):  # 1 change to b
                        gvM[:, :, m1] = np.where(rmsearrM[:, :, m1] == 1, gvarr[:, :, m1], -9999)
                        npvM[:, :, m1] = np.where(rmsearrM[:, :, m1] == 1, npvarr[:, :, m1], -9999)
                        soilM[:, :, m1] = np.where(rmsearrM[:, :, m1] == 1, soilarr[:, :, m1], -9999)
                        shM[:, :, m1] = np.where(rmsearrM[:, :, m1] == 1, sharr[:, :, m1], -9999)

                    gvM[gvM == -9999] = np.nan
                    npvM[npvM == -9999] = np.nan
                    soilM[soilM == -9999] = np.nan
                    shM[shM == -9999] = np.nan

                    gvComp = np.nanmean(gvM, axis=2)
                    npvComp = np.nanmean(npvM, axis=2)
                    soilComp = np.nanmean(soilM, axis=2)
                    shComp = np.nanmean(shM, axis=2)

                else:
                    print(y, mm, ': single image in the time window: copy as is')

                    gvComp = gv.GetRasterBand(ind[0]+1).ReadAsArray()
                    npvComp = npv.GetRasterBand(ind[0]+1).ReadAsArray()
                    soilComp = soil.GetRasterBand(ind[0]+1).ReadAsArray()
                    shComp = sh.GetRasterBand(ind[0]+1).ReadAsArray()

                # fill in all nan with RBD and RBFW
                # first identify the masks
                gvMs = np.squeeze(np.any([[gvComp == 0], [np.isnan(gvComp)], [gvComp == -9999]], axis = 0),axis=(0))
                npvMs = np.squeeze(np.any([[npvComp == 0], [np.isnan(npvComp)], [npvComp == -9999]], axis=0), axis=(0))
                soilMs = np.squeeze(np.any([[soilComp == 0], [np.isnan(soilComp)], [soilComp == -9999]], axis=0), axis=(0))
                shMs = np.squeeze(np.any([[shComp == 0], [np.isnan(shComp)], [shComp == -9999]], axis=0), axis=(0))
                # and apply them
                gvComp[gvMs] = gvRBFComp[gvMs]
                npvComp[npvMs] = npvRBFComp[npvMs]
                soilComp[soilMs] = soilRBFComp[soilMs]
                shComp[shMs] = shRBFComp[shMs]

            if len(ind) == 0:
            # else:
                print(y, mm, ': empty range, take RBF')

                gvComp = gvRBFComp
                npvComp = npvRBFComp
                soilComp = soilRBFComp
                shComp = shRBFComp

            # add band to the out LS
            gvLS = np.append(gvLS, np.atleast_3d(gvComp), axis=2)
            npvLS = np.append(npvLS, np.atleast_3d(npvComp), axis=2)
            soilLS = np.append(soilLS, np.atleast_3d(soilComp), axis=2)
            shLS = np.append(shLS, np.atleast_3d(shComp), axis=2)

    # create the output filename
    nOut = os.path.basename(os.path.splitext((npvtss[0]))[0])

    # output paths
    gvOutP = os.path.join(gvtile, nOut + '_FnF.tif')
    npvOutP = os.path.join(npvtile, nOut + '_FnF.tif')
    soilOutP = os.path.join(soiltile, nOut + '_FnF.tif')
    shOutP = os.path.join(shtile, nOut + '_FnF.tif')

    # get rid of the first empty array at axis 2
    gvLS = gvLS[:,:,1:gvLS.shape[2]]
    npvLS = npvLS[:, :, 1:npvLS.shape[2]]
    soilLS = soilLS[:, :, 1:soilLS.shape[2]]
    shLS = shLS[:, :, 1:shLS.shape[2]]

    gvStak = rasterio.open(gvrbfp[0], "r", tiled=True, BIGTIFF='YES')

    gvProfile = gvStak.profile
    gvProfile.update(count=gvLS.shape[2])
    with rasterio.open(gvOutP, 'w', **gvProfile) as dst:
        for ith, outbn in enumerate(bandNames):
            dst.set_band_description(ith+1, outbn)
            # dst.update_tags(ith, ns='TIMESERIES', dates=outbn[0:4]+'-'+outbn[4:6]+'-'+outbn[6:8]+'T00:00:00')
            # dst.update_tags(ith, ns='TIMESERIES',dates=outbnDate[ith])
            # dst.update_tags(ith, ns='FORCE', Date=datetime.strptime(outbn, '%Y%m%d' ) )
            # dst.update_tags(ith+1, ns='FORCE', Date=outbnDate[ith])
            dst.update_tags(ith+1, ns='FORCE', Date=outbn[0:4]+'-'+outbn[4:6]+'-'+outbn[6:8])
            dst.update_tags(ith+1, ns='FORCE', Wavelength_unit='decimal year')
            dst.update_tags(ith+1, ns='FORCE', Domain='SMA_TSA')
            # dst.update_tags(ith+1, start_time=outbnDate[ith]+'T00:00:00')
            # dst.update_tags(ith, DATE_TIME=outbn[0:4] + '-' + outbn[4:6] + '-' + outbn[6:8]+'T12:00')
            # dst.update_tags(ith, dates=outbn[0:4] + '-' + outbn[4:6] + '-' + outbn[6:8])
        # rasterio needs the band axis first.
        dst.write(np.moveaxis(gvLS,-1,0), list(range(1,gvLS.shape[2]+1)))
        dst.close()

    with rasterio.open(npvOutP, 'w', **gvProfile) as dst:
        for ith, outbn in enumerate(bandNames):
            dst.set_band_description(ith+1, outbn)
            dst.update_tags(ith + 1, ns='FORCE', Date=outbn[0:4] + '-' + outbn[4:6] + '-' + outbn[6:8])
            dst.update_tags(ith + 1, ns='FORCE', Wavelength_unit='decimal year')
            dst.update_tags(ith + 1, ns='FORCE', Domain='SMA_TSA')
        # rasterio needs the band axis first.
        dst.write(np.moveaxis(npvLS,-1,0), list(range(1,npvLS.shape[2]+1)))
        dst.close()

    with rasterio.open(soilOutP, 'w', **gvProfile) as dst:
        for ith, outbn in enumerate(bandNames):
            dst.set_band_description(ith+1, outbn)
            dst.update_tags(ith + 1, ns='FORCE', Date=outbn[0:4] + '-' + outbn[4:6] + '-' + outbn[6:8])
            dst.update_tags(ith + 1, ns='FORCE', Wavelength_unit='decimal year')
            dst.update_tags(ith + 1, ns='FORCE', Domain='SMA_TSA')
        # rasterio needs the band axis first.
        dst.write(np.moveaxis(soilLS,-1,0), list(range(1,soilLS.shape[2]+1)))
        dst.close()

    with rasterio.open(shOutP, 'w', **gvProfile) as dst:
        for ith, outbn in enumerate(bandNames):
            dst.set_band_description(ith+1, outbn)
            dst.update_tags(ith + 1, ns='FORCE', Date=outbn[0:4] + '-' + outbn[4:6] + '-' + outbn[6:8])
            dst.update_tags(ith + 1, ns='FORCE', Wavelength_unit='decimal year')
            dst.update_tags(ith + 1, ns='FORCE', Domain='SMA_TSA')
        # rasterio needs the band axis first.
        dst.write(np.moveaxis(shLS,-1,0), list(range(1,shLS.shape[2]+1)))
        dst.close()

print('execution took ', datetime.now() - startTime)

# :END OF THE CODE: #


# debug
#
# gvOutP = gvtile + '\\' + nOut + '_gv2020_03_stack_gvM.tif'
# gvProfile = gvStak.profile
# gvProfile.update(count=4)
# with rasterio.open(gvOutP, 'w', **gvProfile) as dst:
#     dst.set_band_description(1, xnames[0])
#     dst.set_band_description(2, xnames[1])
#     dst.set_band_description(3, xnames[2])
#     dst.set_band_description(4, xnames[3])
#     # dst.update_tags(1, ns='FORCE', Date=bandNames[0][0:4] + '-' + bandNames[0][4:6] + '-' + bandNames[0][6:8])
#     dst.update_tags(1, ns='FORCE', Date=xnames[0][0:4] + '-' + xnames[0][4:6] + '-' + xnames[0][6:8])
#     dst.update_tags(2, ns='FORCE', Date=xnames[1][0:4] + '-' + xnames[1][4:6] + '-' + xnames[1][6:8])
#     dst.update_tags(3, ns='FORCE', Date=xnames[2][0:4] + '-' + xnames[2][4:6] + '-' + xnames[2][6:8])
#     dst.update_tags(4, ns='FORCE', Date=xnames[3][0:4] + '-' + xnames[3][4:6] + '-' + xnames[3][6:8])
#     dst.update_tags(1, ns='FORCE', Wavelength_unit='decimal year')
#     dst.update_tags(1, ns='FORCE', Domain='SMA_TSA')
#     dst.update_tags(2, ns='FORCE', Domain='SMA_TSA')
#     dst.update_tags(3, ns='FORCE', Domain='SMA_TSA')
#     dst.update_tags(4, ns='FORCE', Domain='SMA_TSA')
#     # rasterio needs the band axis first.
#     dst.write(np.moveaxis(gvM,-1,0), [1,2,3,4])
#     # dst.write(np.moveaxis(gvComp, -1, 0), 1)
#     dst.close()
