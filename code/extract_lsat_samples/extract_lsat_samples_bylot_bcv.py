import ee
import sys
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
endJulian = 243
#startJulian = 201  # start april
#endJulian = 243  # early oct


# values to copy from fusion table/feature collection
feat_properties = ['Site']
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

OUTFILE = 'C:/tmp/lsat_bylot_bcv_subsite_samples_buf100m.csv'

# Landsat Surface Reflectance collections
ls5_1 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR")
ls5_2 = ee.ImageCollection("LANDSAT/LT05/C01/T2_SR")
ls7_1 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
ls7_2 = ee.ImageCollection("LANDSAT/LE07/C01/T2_SR")
ls8_1 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
ls8_2 = ee.ImageCollection("LANDSAT/LC08/C01/T2_SR")

# define region of interest
sites = ee.FeatureCollection('users/loganberner/bylot_bcv_subsites')

boundary = sites.geometry().bounds()

# ------------------------------------------------------------------------------------
# BE CAREFUL MODIFYING THE STUFF BELOW
# ------------------------------------------------------------------------------------


# function to concatenate two GEE lists
def add_lists(list1, list2):
    first = ee.List(list1).add(list2.get(0))
    n = list2.size()

    def _combine(this, prev):
        return ee.List(prev).add(this)

    return ee.List(list2).slice(1, n).iterate(_combine, first)


ALL_BANDS = add_lists(BAND_LIST, ADDON_BANDLIST)


# merge all collections in one
LS_COLL = ee.ImageCollection(ls5_1
                             .merge(ls7_1
                                    .merge(ls8_1
                                           .merge(ls5_2
                                                  .merge(ls7_2
                                                         .merge(ls8_2)))))) \
    .filter(ee.Filter.calendarRange(startJulian,endJulian, "day_of_year")) \
    .filterBounds(boundary) \
    .map(lambda image: image.addBands(ADDON, ADDON_BANDLIST)) \
    .select(ALL_BANDS) \
    .map(lambda image: image.float())


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


# function to extract a feature from a non-null collection
def get_region(feature):

    coll = LS_COLL.filterBounds(ee.Feature(feature).geometry())
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

# convert sites collection to list
n_sites = sites.size().getInfo()
sites_list = sites.toList(n_sites)

# iterate thru site list
for i in range(0, n_sites):
    time1 = datetime.datetime.now()
    site = sites_list.get(i)

    try:
        temp_dicts = get_region(site)

    # if any error occurs, print it and break the loop
    except Exception as e:
        eprint(e)
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

time0_ = datetime.datetime.now()
cprint('Total Time taken: {t} seconds'.format(t=str(round((time0_ - time0).total_seconds(), 1))))