import ee
import sys
import time
import math
import datetime
import pandas as pd
import traceback

ee.Initialize()

# SPECIFY SEVERAL VARIABLES
BUFFER_DIST = 100  # set it to 0 if buffer is not needed, else it will make script timed-out
CRS = 'EPSG:4326'
SCALE = 30
NULL_VALUE = 0
MASK_VALUE = 0  # used for ADDON image
startJulian = 152  # start april
endJulian = 243  # early oct

# define region of interest
sites = ee.FeatureCollection('users/loganberner/tundra_shrub_cluster_sites_wgs84')
OUTFILE = 'C:/tmp/lsat_tundra_shrub_cluster_sites_buf100m_20190507.csv'

# parts to divide into if 'too many values' error
n = 15

# wait time
wait = 10  # seconds


# values to copy from fusion table/feature collection
feat_properties = ['site']
KEYS = ee.List(feat_properties)

# properties to retrieve from scene
scene_properties = ['CLOUD_COVER',
                    'GEOMETRIC_RMSE_MODEL',
                    'LANDSAT_ID',
                    'SOLAR_ZENITH_ANGLE']
PROPERTY_LIST = ee.List(scene_properties)

# Bands to retrieve
bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'pixel_qa', 'radsat_qa']
BAND_LIST = ee.List(bands)

# addon asset and bands
ADDON = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').float().unmask(MASK_VALUE)
ADDON_BANDLIST = ee.List(['max_extent'])

# Landsat Surface Reflectance collections
ls5_1 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR")
ls5_2 = ee.ImageCollection("LANDSAT/LT05/C01/T2_SR")
ls7_1 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
ls7_2 = ee.ImageCollection("LANDSAT/LE07/C01/T2_SR")
ls8_1 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
ls8_2 = ee.ImageCollection("LANDSAT/LC08/C01/T2_SR")

# ------------------------------------------------------------------------------------
# BE CAREFUL MODIFYING THE STUFF BELOW
# ------------------------------------------------------------------------------------


# Functions ----------------------------------------

# function to concatenate two GEE lists
def add_lists(list1, list2):
    first = ee.List(list1).add(list2.get(0))
    n = list2.size()

    def _combine(this, prev):
        return ee.List(prev).add(this)

    return ee.List(list2).slice(1, n).iterate(_combine, first)


# print and flush function
def cprint(text):
    sys.stdout.write(str(text) + '\n')
    sys.stdout.flush()


# print exception with traceback without breaking the code
def eprint(excep):
    print(excep)
    print('---------------------------------------------------------------------')
    print('Error Traceback:')
    print(traceback.format_exc())
    print('---------------------------------------------------------------------')


# Function calls and usage ----------------------------
ALL_BANDS = add_lists(BAND_LIST, ADDON_BANDLIST)

# boundary of the region
boundary = sites.geometry().bounds()

try:
    n_sites = sites.size().getInfo()
except Exception as e:
    print(e)
    exit(1)

# convert sites collection to list
sites_list = sites.toList(n_sites)


# merge all collections in one
def ls_coll(start_day=startJulian,
            end_day=endJulian):

    return ee.ImageCollection(ls5_1
                                 .merge(ls7_1
                                        .merge(ls8_1
                                               .merge(ls5_2
                                                      .merge(ls7_2
                                                             .merge(ls8_2)))))) \
        .filter(ee.Filter.calendarRange(start_day, end_day, "day_of_year")) \
        .filterBounds(boundary) \
        .map(lambda image: image.addBands(ADDON, ADDON_BANDLIST)) \
        .select(ALL_BANDS) \
        .map(lambda image: image.float())


# function to extract properties from collection as list of dictionaries
def get_coll_dict(image_collection):
    n_images = ee.ImageCollection(image_collection).size()
    coll_list = ee.ImageCollection(image_collection).toList(n_images)

    # extract properties in PROPERTY_LIST from one image
    def get_prop(image):
        image_dict = ee.Image(image).toDictionary()
        val_list = PROPERTY_LIST.map(lambda prop: image_dict.get(prop)).add(ee.Image(image).id())
        return ee.Dictionary.fromLists(PROPERTY_LIST.add('id'), val_list)

    return coll_list.map(get_prop)


LS_COLL = ls_coll()


# function to extract a feature from a non-null collection
def get_region(feature,
               collection=LS_COLL):

    coll = collection.filterBounds(ee.Feature(feature).geometry())
    ncoll = coll.size().getInfo()

    # feature to be extracted
    feat = ee.Algorithms.If(ee.Number(BUFFER_DIST).gt(0),
                            ee.Feature(feature).buffer(BUFFER_DIST).bounds(),
                            ee.Feature(feature))

    feat_dict = feat.getInfo()
    feat_props = feat_dict['properties']

    pt_list = list()

    if ncoll > 0:

        # list of properties of each image in the collection as dictionary
        coll_dicts = get_coll_dict(coll).getInfo()

        # extract pixel values from the collection
        # the output is a list or lists with the first row as column names
        temp_list = coll.getRegion(ee.Feature(feat).geometry(), SCALE).getInfo()
        n = len(temp_list)

        elem_names = temp_list[0]
        id_column = elem_names.index('id')

        for j in range(1, n):
            pt_dict = dict()
            for k in range(0, len(elem_names)):
                pt_dict[elem_names[k]] = temp_list[j][k]
                id = temp_list[j][id_column]

                scene_prop_list = list(coll_dict for coll_dict in coll_dicts if coll_dict['id'] == id)
                scene_prop = scene_prop_list[0]

                for prop in scene_properties:
                    pt_dict[prop] = scene_prop[prop]

                for prop in feat_properties:
                    pt_dict[prop] = feat_props[prop]

            pt_list.append(pt_dict)

    else:
        pt_dict = dict()
        pt_dict['id'] = NULL_VALUE
        pt_dict['time'] = NULL_VALUE

        for prop in feat_properties:
            pt_dict[prop] = feat_props[prop]

        for prop in scene_properties:
            pt_dict[prop] = NULL_VALUE

        for band in bands:
            pt_dict[band] = NULL_VALUE

        pt_list.append(pt_dict)

    return pt_list


time0 = datetime.datetime.now()

# iterate thru site list
i=0
while i < n_sites:

    time1 = datetime.datetime.now()
    site = sites_list.get(i)

    try:
        temp_dicts = get_region(site)

    # if any error occurs, print it and break the loop
    except Exception as e:
        print(e)

        if 'Earth Engine memory capacity exceeded' in e.args[0]:
            cprint('Waiting 30 secs...'.format(str(wait)))
            time.sleep(wait)
            continue

        elif 'Too many values' in e.args[0]:
            cprint('Dividing the julian date range in {} parts and retrying...'.format(str(n)))

            # divide into n parts

            dates = list(range(startJulian, endJulian, int(math.floor((endJulian-startJulian)/n))))
            dates[-1] = endJulian
            date_tuples = [(dates[i], dates[i+1]-1)
                           if i != len(dates)-2
                           else (dates[i], dates[i+1])
                           for i in range(0, len(dates)-1)]

            part_dicts = list()

            for date_tuple in date_tuples:
                coll = ls_coll(start_day=date_tuple[0],
                               end_day=date_tuple[1])

                part_dict = get_region(site,
                                       collection=coll)

                for part in part_dict:
                    part_dicts.append(part)

            temp_dicts = part_dicts

        elif 'ServerNotFoundError' in e.args[0] or \
                'Unable to find the server' in e.args[0] or \
                'getaddrinfo failed' in e.args[0] or \
                'connection attempt failed' in e.args[0]:

            cprint('Waiting 30 secs...'.format(str(wait)))
            time.sleep(wait)
            continue

        else:
            cprint('Exiting...')
            break

    # convert to pandas data frame
    df = pd.DataFrame(temp_dicts)

    # all extracted values to file
    if i == 0:
        with open(OUTFILE, 'w') as f:
            df.to_csv(f, index=False, header=True)
    else:
        with open(OUTFILE, 'a') as f:
            df.to_csv(f, index=False, header=False)

    time2 = datetime.datetime.now()

    cprint('Time taken for sample {s}: {t} seconds'.format(s=str(i + 1),
                                                           t=str(round((time2 - time1).total_seconds(), 1))))

    i = i + 1

time0_ = datetime.datetime.now()
cprint('Total Time taken: {t} seconds'.format(t=str(round((time0_ - time0).total_seconds(), 1))))
