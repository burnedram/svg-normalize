#! /usr/bin/env python3
import svgpathtools
import svgwrite
import sys
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description='normalize svg images')
    parser.add_argument('inFile',
        help='input svg')
    parser.add_argument('-o', '--outFile', dest='outFile',
        action='store',
        help='output svg. default same as input')
    parser.add_argument('-f', '--fit-view-box', dest='fit',
        action='store_true',
        help='change viewBox such that the image perfectly fits inside it')
    parser.add_argument('-b', '--bounding-box', dest='bb',
        action='store_true',
        help='move the bounding box to the center of the viewBox')
    parser.add_argument('-s', '--scale', dest='scale',
        action='store_true',
        help='scale paths to fit inside the viewBox (radial)')
    parser.add_argument('-v', '--view-box', dest='vb',
        action='store',
        help='resize viewBox ("x y w h")')
    parser.add_argument('-p', '--keep-precision', dest='precision',
        action='store_true',
        help='use same precision in output as in input svg')
    args = parser.parse_args()
    if args.outFile is None:
        args.outFile = args.inFile
    return args

def main(args):
    paths, attributes, svg_attributes = svgpathtools.svg2paths2(args.inFile)
    if not paths:
        raise Exception('no <path>:s in svg')
    if not 'viewBox' in svg_attributes:
        raise Exception("no 'viewBox' attribute on svg root element")

    if args.precision:
        original_precision = max(map(get_precision, iter_points(paths)))

    if args.fit:
        (xmin, xmax, ymin, ymax) = svgpathtools.paths2svg.big_bounding_box(paths)
        width = xmax - xmin
        height = ymax - ymin
        svg_attributes['viewBox'] = ' '.join(map(format_float, (xmin, ymin, width, height)))
    else:
        (xmin, ymin, width, height) = list(map(float, svg_attributes['viewBox'].split()))
        xmax = xmin + width
        ymax = ymin + height
    cx = xmin + width/2
    cy = ymin + height/2
    center = complex(cx, cy)
    if args.bb:
        (bb_xmin, bb_xmax, bb_ymin, bb_ymax) = svgpathtools.paths2svg.big_bounding_box(paths);
        bb_width = bb_xmax - bb_xmin
        bb_height = bb_ymax - bb_ymin
        bb_cx = bb_xmin + bb_width/2
        bb_cy = bb_ymin + bb_height/2
        dx = bb_cx - cx
        dy = bb_cy - cy
        # move center of paths to center of viewBox
        translate(paths, complex(-dx, -dy))

    if args.scale:
        # r_max is the farthest distance of any point from the center of the viewBox/bounding box
        r_max = 0
        for path in paths:
            for segment in path:
                if isinstance(segment, svgpathtools.Arc):
                    print("WARNING: Arc is not supported by radialrange, skipping")
                    continue
                ((d_min, t_min), (d_max, t_max)) = segment.radialrange(center)
                r_max = max(r_max, d_max)

        if r_max > 0:
            # rescale all paths
            diameter = r_max*2
            ratio = min(width/diameter, height/diameter)
            rescale(paths, ratio, center)

    if not args.vb is None:
        (new_xmin, new_ymin, new_width, new_height) = list(map(float, args.vb.split()))
        new_cx = new_xmin + new_width/2
        new_cy = new_ymin + new_height/2
        dx = cx - new_cx
        dy = cy - new_cy
        translate(paths, complex(-dx, -dy))
        ratio = min(new_width/width, new_height/height)
        rescale(paths, ratio, complex(new_cx, new_cy))
        svg_attributes['viewBox'] = args.vb

        tlate = complex(new_xmin + new_width / (xmin + width),
                        new_ymin + new_height / (ymin + height))
        for i in range(len(paths)):
            if 'transform' not in attributes[i]:
                continue
            transform = attributes[i]['transform']
            new_transform = ' '.join(fix_transform(transform, tlate))
            attributes[i]['transform'] = new_transform
            
    if args.precision:
        round_paths(paths, original_precision)
        svg_attributes['viewBox'] = ' '.join(map(format_float, map(lambda x: round(x, original_precision), map(float, svg_attributes['viewBox'].split()))))

    svgpathtools.wsvg(paths, attributes=attributes, svg_attributes=svg_attributes, filename=args.outFile)

def format_float(f):
    return '{:f}'.format(f).rstrip('0').rstrip('.')

def fix_transform(transform, ratio):
    p = re.compile('(\w+)\(((?:-?(?:\d|\.)+(?:\s|,)*)+)\)')
    for m in p.finditer(transform):
        name = m.group(1)
        if name == 'translate' or name == 'matrix':
            t_args = list(map(float, re.split('\s|,', m.group(2))))
        else:
            yield m.group(0)
            continue
        t_args[-2] = t_args[-2] * ratio.real
        t_args[-1] = t_args[-1] * ratio.imag
        yield '{}({})'.format(name, ' '.join(map(str, t_args)))

def translate(paths, xy):
    for i in range(len(paths)):
        paths[i] = paths[i].translated(xy)

def rescale(paths, ratio, origin):
    def _r(x, independant=False):
        if independant:
            return complex(x.real*ratio, x.imag*ratio)
        else:
            return complex(
                (x.real - origin.real)*ratio + origin.real, 
                (x.imag - origin.imag)*ratio + origin.imag)
    for path in paths:
        for segment in path:
            segment.start = _r(segment.start)
            segment.end = _r(segment.end)
            if isinstance(segment, svgpathtools.Line):
                pass
            elif isinstance(segment, svgpathtools.QuadraticBezier):
                segment.control = _r(segment.control)
            elif isinstance(segment, svgpathtools.CubicBezier):
                segment.control1 = _r(segment.control1)
                segment.control2 = _r(segment.control2)
            elif isinstance(segment, svgpathtools.Arc):
                segment.radius = _r(segment.radius, True)
            else:
                raise Exception('unknown segment {}'.format(segment))

def round_complex(c, precision):
    return complex(
        round(c.real, precision),
        round(c.imag, precision))

def round_paths(paths, precision):
    for path in paths:
        for segment in path:
            segment.start = round_complex(segment.start, precision)
            segment.end = round_complex(segment.end, precision)
            if isinstance(segment, svgpathtools.Line):
                pass
            elif isinstance(segment, svgpathtools.QuadraticBezier):
                segment.control = round_complex(segment.control, precision)
            elif isinstance(segment, svgpathtools.CubicBezier):
                segment.control1 = round_complex(segment.control1, precision)
                segment.control2 = round_complex(segment.control2, precision)
            elif isinstance(segment, svgpathtools.Arc):
                segment.radius = round_complex(segment.radius, precision)
            else:
                raise Exception('unknown segment {}'.format(segment))

def get_precision(f):
    if isinstance(f, complex):
        return max(get_precision(f.real), get_precision(f.imag))
    f = format_float(f)
    if not '.' in f:
        return 0
    return len(f) - f.index('.') - 1

def iter_points(paths):
    for path in paths:
        for segment in path:
            if isinstance(segment, svgpathtools.Arc):
                yield segment.start
                yield segment.radius
                yield segment.end
            else:
                for point in segment:
                    yield point

if __name__ == '__main__':
    main(parse_args())
