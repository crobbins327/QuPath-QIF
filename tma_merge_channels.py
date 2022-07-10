import numpy as np
# import tifffile
import pyvips
import os
import sys
import argparse
import re
import multiprocessing as mp

def main(argv=sys.argv):


    parser = argparse.ArgumentParser(
        description='Merge IF channels for TMA spot images (from HistoRX microscopes)',
        formatter_class=HelpFormatter,
    )

    parser.add_argument(
        'filepaths', metavar='FILE', nargs='+',
        help='Image file(s) or directory containing channel images to be processed',
    )
    parser.add_argument(
        '-o', '--output', dest='output', default='QIF_images/Merged_{spotID:05}.tif',
        metavar='PATH',
        help=("Output files. Uses regex pattern of spotID to name images, store in QIF_images folder. "
              "Can be renamed but filenames MUST contain {spotID:#}."),
    )

    parser.add_argument(
        '-f', '--filename_format', dest='filename_format', default='{channel}_{spotID:05}.tif',
        metavar='str',
        help=("Input file name format/pattern. Channel images require consistent naming pattern to be combined together."
              "MUST contain {channel} and {spotID:05}.")
    )

    parser.add_argument(
        '-c', '--channels', dest='channels', default='all',
        metavar='str',
        help=("'all' or comma separated string of channel names to include (depends on input filename-format). Eg. 'DAPI, CK, HER2'.")
    )

    parser.add_argument(
        '-k', '--max_cores', dest="max_cores", default = 6,
        metavar = int,
        help=("maximum cpu cores to use for multiprocessing.")
    )

    args = parser.parse_args(argv[1:])

    #verify that filename_format contains {channel} and {spotID:05}
    sc = re.search("{channel}", args.filename_format)
    ss = re.search("{spotID:05}", args.filename_format)
    if not sc:
        raise ValueError("filename_format MUST contain {channel}!")
    if not ss:
        raise ValueError("filename_format MUST contain {spotID:05}!")

    # check if filepaths are .tif/.tiff files or a directory
    filepaths = args.filepaths
    print(filepaths)
    if len(filepaths) == 1 and os.path.isdir(filepaths[0]):
        print("reading .tif/.tiff channel images from directory {}".format(filepaths[0]))
        main_path = filepaths[0]
        c_imgfiles = [t for t in os.listdir(filepaths[0]) if t.endswith('.tif') or t.endswith('.tiff')]
    else:
        main_path = os.path.split(filepaths[0])[0]
        c_imgfiles = [t for t in filepaths if t.endswith('.tif') or t.endswith('.tiff')]
    if len(c_imgfiles) <= 0:
        raise FileNotFoundError("no .tif/.tiff files were found in filepaths provided...")

    # extract channels and spotID from filename_format pattern
    output_pattern = args.output
    input_pattern = args.filename_format
    input_regex = re.sub(r'{([^:}]+)(?:[^}]*)}', r'(?P<\1>.*?)',
                   input_pattern.replace('.', '\.'))

    gds = []
    channels = []
    spotIDs = []
    for f in c_imgfiles:
        m = re.match(input_regex, f)
        if m:
            gd = m.groupdict()
            gd['name'] = f
            channels.append(gd['channel'])
            spotIDs.append(gd['spotID'])
            gds.append(gd)
        else:
            print("file {} does not match input filename-format {}... please check that this is correct.".format(f, input_pattern))

    if args.channels == 'all':
        channels = np.unique(channels)
    elif isinstance(args.channels, str):
        input_channels = [c.strip() for c in args.channels.split(',')]
        channels = list(set(channels).intersection(set(input_channels)))
        if len(channels) < 1:
            raise ValueError("None of the input channels were found. input channels: {}".format(input_channels))
    else:
        raise ValueError("Channels need to be a comma separated string or 'all'. {} not acceptable".format(args.channels))
    print("channels to be used: {}".format(channels))

    spotIDs = np.unique(spotIDs)
    print("found {} unique spots...".format(len(spotIDs)))

    out_dir = os.path.split(output_pattern)[0]
    if out_dir and not os.path.exists(out_dir):
        print("making {} directory..".format(out_dir))
        os.makedirs(out_dir)
    elif os.path.exists(out_dir):
        print("using {} as output directory...".format(out_dir))

    #use multiprocessing to merge all the images fast
    # Multiprocessing pool to load and write images
    num_workers = mp.cpu_count()
    max_cores = args.max_cores
    if num_workers > max_cores:
        num_workers = max_cores
    print('Using CPUs:', num_workers)

    with mp.Pool(num_workers) as pool:
        # pool = mp.Pool(num_workers)
        iterable = [(spotID, channels, input_pattern, main_path, output_pattern) for spotID in spotIDs]
        res = pool.starmap(mergeChannelImgs, iterable)
        pool.close()
        pool.join()

    # trials = 0
    # while trials <= 10:
    #     try:
    #         res = pool.starmap(mergeChannelImgs, iterable)
    #         if res.count(None) == len(res):
    #             pool.close()
    #             break
    #         trials += 1
    #     except Exception as e:
    #         print(e)
    #         trials += 1

    # mergeChannelImgs(spotIDs[0], channels, input_pattern, main_path, output_pattern)



def mergeChannelImgs(spotID, channels, input_pattern, main_path, output_pattern):
    filenames = [os.path.join(main_path, input_pattern.format(channel=c, spotID=int(spotID))) for c in channels]

    tiff_files = [pyvips.Image.new_from_file(file, access="sequential")[0]
                  for file in filenames]
    final_stack_filename = os.path.abspath(output_pattern.format(spotID=int(spotID)))
    # join all list elements
    if len(tiff_files) > 1:
        final_image = tiff_files[0].bandjoin(tiff_files[1:])
    else:
        raise Warning("merge resulted in a single channel image, is this intentional?")
        final_image = tiff_files[0]
    final_image.set_type(pyvips.GValue.gint_type, "page-height", final_image.height)
    final_image.set_type(pyvips.GValue.gstr_type, "image-description",
                f"""<?xml version="1.0" encoding="UTF-8"?>
      <OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">
        <Image ID="Image:0">
          <!-- Minimum required fields about image dimensions -->
          <Pixels DimensionOrder="XYCZT"
                  ID="Pixels:0"
                  SizeC="{final_image.image_bands}"
                  SizeT="1"
                  SizeX="{final_image.width}"
                  SizeY="{final_image.height}"
                  SizeZ="1"
                  Type="uint8">
          </Pixels>
        </Image>
      </OME>""")
    print("saving merge {} of spotID: {}, channels: {}".format(final_stack_filename, spotID, channels))
    final_image.cast("uchar").tiffsave(final_stack_filename,
                                       compression="lzw", Q=100, tile=True,
                                       tile_width=512, tile_height=512,
                                       pyramid=True, bigtiff=True, rgbjpeg=False)

class HelpFormatter(argparse.HelpFormatter):
    """Help message formatter which adds default values to argument help.
    Based on argparse.ArgumentDefaultsHelpFormatter.
    """

    def _get_help_string(self, action):
        help = action.help
        if isinstance(action, (argparse._HelpAction, argparse._VersionAction)):
            help = help.capitalize()
        elif (
            not isinstance(action, argparse._StoreTrueAction)
            and "%(default)" not in help
            and "(default:" not in help
            and action.default is not argparse.SUPPRESS
        ):
            defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
            if action.option_strings or action.nargs in defaulting_nargs:
                help += " (default: %(default)s)"
        return help

if __name__ == '__main__':
    sys.exit(main())
