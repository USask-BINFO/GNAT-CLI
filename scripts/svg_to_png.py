#!/usr/bin/env python

if __name__ == "__main__":
    # try:
    from argparse import ArgumentParser
    import cairosvg
    import os.path

    parser = ArgumentParser("Convert an SVG image of a plot to PNG to display in the dashboard more responsively. Saves a copy of `svg_file` as PNG.")
    parser.add_argument("svg_file", help="Path to SVG file to convert.")
    args = parser.parse_args()

    if os.path.exists(args.svg_file):
        if args.svg_file.endswith(".svg"):
            out_file = args.svg_file[:-4] + '.png'

            cairosvg.svg2png(url=args.svg_file, write_to=out_file)

            # with open(args.svg_file, 'rb') as svg_file:
            #     with wand.image.Image() as image:
            #         with wand.color.Color('transparent') as background_color:
            #             library.MagickSetBackgroundColor(image.wand, background_color.resource)
                    
            #         image.read(blob=svg_file.read(), format='svg')
            #         png_image = image.make_blob('png32')

            #     with open(out_file, 'wb') as f:
            #         f.write(png_image)

        else:
            raise TypeError(f"Invalid file type. Expecting a 'svg' instead got '{args.svg_file.split('.')[-1]}'.")
    else:
        raise FileNotFoundError(f"Provided file not found at {args.svg_file}")
    # except ImportError as IE:
    #     print(f"Warning: ImageMagick not installed on system, cannot convert SVGs to PNGs.")