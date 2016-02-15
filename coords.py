#!/usr/bin/env python

import argparse
import zipfile
import xml.etree.ElementTree
import sys

parser = argparse.ArgumentParser(description="Convert KML file into SVG")
parser.add_argument('input', help="Filename to process (KML or KMZ supported)", nargs='?', type=argparse.FileType("rb"), default=sys.stdin)
parser.add_argument('output', help="Filename to write", nargs='?', type=argparse.FileType("w"), default=sys.stdout)
parser.add_argument('-m', '--massstab', help="Specify 200 for 1:200 etc", type=int, default=1000)
args = parser.parse_args()

import numpy as np

class MapReference:
    def __init__(self, A, B):
        self.A = A
        self.B = B
        self.boundingbox = [(0,0), (0,0), (0,0), (0,0)]
        self.R = np.asarray([[1,0], [0,1]])

    def orientate(self, points_plain):
        # deprecated
        points = [self.get_xy(p) for p in points_plain]
        min_x = min([x for x, y in points])
        max_x = max([x for x, y in points])
        min_y = min([y for x, y in points])
        max_y = max([y for x, y in points])

        self.boundingbox = [(max_x, min_y), (max_x, max_y), (min_x, max_y), (min_x, min_y)]

        # see http://gis.stackexchange.com/questions/22895/how-to-find-the-minimum-area-rectangle-for-given-points
        from scipy.spatial import ConvexHull

        points = np.asarray([self.get_xy(p) for p in points_plain])

        # get the convex hull for the points
        hull = points[ConvexHull(points).vertices]

        # calculate edge angles
        edges = np.zeros((len(hull)-1, 2))
        edges = hull[1:] - hull[:-1]

        angles = np.zeros((len(edges)))
        angles = np.arctan2(edges[:, 1], edges[:, 0])

        angles = np.abs(np.mod(angles, np.pi / 2))
        angles = np.unique(angles)

        # find rotation matrices (given transposed for later use)
        rotations = np.vstack([
          np.cos(angles),
          np.sin(angles),
          -np.sin(angles),
          np.cos(angles)]).T

        rotations = rotations.reshape((-1, 2, 2))

        # apply rotations to the hull
        rot_points = np.dot(rotations, hull.T)

        # find the bounding points
        min_x = np.nanmin(rot_points[:, 0], axis=1)
        max_x = np.nanmax(rot_points[:, 0], axis=1)
        min_y = np.nanmin(rot_points[:, 1], axis=1)
        max_y = np.nanmax(rot_points[:, 1], axis=1)

        # find the box with the best area
        areas = (max_x - min_x) * (max_y - min_y)
        best_idx = np.argmin(areas)

        # return the best box
        x1 = max_x[best_idx]
        y1 = max_y[best_idx]
        x2 = min_x[best_idx]
        y2 = min_y[best_idx]
        r = rotations[best_idx]

        self.rotate(r)
        self.boundingbox = [np.dot([x1, y2], r), np.dot([x1, y1], r), np.dot([x2, y1], r), np.dot([x2, y2], r)]

        if self.max_height() > self.max_width():
            self.rotate(np.asarray([[0, 1], [1, 0]]))
            self.boundingbox = [self.boundingbox[1], self.boundingbox[2], self.boundingbox[3], self.boundingbox[0]]

        if self.get_north() > 90 and self.get_north() < 270:
            self.rotate(np.asarray([[-1, 0], [0, -1]]))
            self.boundingbox = [self.boundingbox[2], self.boundingbox[3], self.boundingbox[0], self.boundingbox[1]]

    def max_width(self):
        return np.linalg.norm(self.boundingbox[1]-self.boundingbox[2])

    def max_height(self):
        return np.linalg.norm(self.boundingbox[0]-self.boundingbox[1])

    # TODO only depend on self.boundingbox

    def get_north(self):
        # TODO
        assert(self.A.lat == self.B.lat)
#        print(np.degrees(np.arcsin(np.linalg.norm((self.boundingbox[2][0], self.boundingbox[1][1]) - self.boundingbox[1]) / self.max_width())) - 90 - -1 * np.degrees(np.arccos(self.R[0,0])))
#        return np.degrees(np.arcsin(np.linalg.norm((self.boundingbox[2][0], self.boundingbox[1][1]) - self.boundingbox[1]) / self.max_width())) - 90
        return -1 * np.degrees(np.arccos(self.R[0,0]))

    def rotate(self, R):
#        self.boundingbox = np.dot(self.boundingbox, R.T)
        self.R = np.dot(self.R, R.T)

    def get_xy(self, C):
        a = self.B.distance(C)
        b = self.A.distance(C)
        c = self.A.distance(self.B)

        # TODO kosinussatz für kugeloberflächen
        alpha = np.arccos((b*b + c*c - a*a) / (2*b*c))

        x = b * np.cos(alpha)
        y = b * np.sin(alpha)

        refpoint = np.dot(self.boundingbox[3], self.R)
        return tuple(np.dot([x, y], self.R) - refpoint)

class MapPoint:
    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon

    def distance(self, B):
        lat1, lon1, lat2, lon2 = map(np.radians, [self.lat, self.lon, B.lat, B.lon])
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))
        m = 6367000 * c
        return m

class KmlFigure:
    pass

class KmlPoint(KmlFigure):
    def __init__(self, coordinate, label):
        self.coordinate = coordinate
        self.label = label

    def get_coordinates(self):
        return [self.coordinate]

class KmlLineString(KmlFigure):
    def __init__(self, coordinates):
        self.coordinates = coordinates

    def get_coordinates(self):
        return self.coordinates

class KmlPolygon(KmlFigure):
    def __init__(self, coordinates):
        self.coordinates = coordinates

    def get_coordinates(self):
        return self.coordinates

class KmlFile:
    def __init__(self, file):
        # Detect KMZ-files
        if zipfile.is_zipfile(file):
            kmz_file = zipfile.ZipFile(file)
            kml_file = kmz_file.open("doc.kml")
        else:
            # input was already seeked by is_zipfile
            file.seek(0)
            kml_file = file
        self.kml = xml.etree.ElementTree.parse(kml_file).getroot()

    def parse_coordinate(self, coord):
        ordinates = coord.split(",")
        # lat and long is switched in kml
        return MapPoint(float(ordinates[1]), float(ordinates[0]))

    def get_name(self):
        return self.kml.find("{http://www.opengis.net/kml/2.2}Document").find("{http://www.opengis.net/kml/2.2}name").text.strip()

    def parse_point(self, placemark, point):
        label = placemark.find("{http://www.opengis.net/kml/2.2}name").text
        coordinate = self.parse_coordinate(point.find("{http://www.opengis.net/kml/2.2}coordinates").text)
        return KmlPoint(coordinate, label)

    def parse_linestring(self, placemark, linestring):
        coordinates = list(map(self.parse_coordinate, linestring.find("{http://www.opengis.net/kml/2.2}coordinates").text.strip().split(" ")))
        return KmlLineString(coordinates)

    def parse_polygon(self, placemark, polygon):
        coordinates = list(map(self.parse_coordinate, polygon.find("{http://www.opengis.net/kml/2.2}outerBoundaryIs").find("{http://www.opengis.net/kml/2.2}LinearRing").find("{http://www.opengis.net/kml/2.2}coordinates").text.strip().split(" ")))
        return KmlPolygon(coordinates)

    def parse_placemark(self, placemark):
        point = placemark.find("{http://www.opengis.net/kml/2.2}Point")
        if point:
            return self.parse_point(placemark, point)

        linestring = placemark.find("{http://www.opengis.net/kml/2.2}LineString")
        if linestring:
            return self.parse_linestring(placemark, linestring)

        polygon = placemark.find("{http://www.opengis.net/kml/2.2}Polygon")
        if polygon:
            return self.parse_polygon(placemark, polygon)

    def get_figures(self):
        return [self.parse_placemark(placemark) for placemark in self.kml.iter("{http://www.opengis.net/kml/2.2}Placemark")]

class Paper:
    def __init__(self, ref, massstab, width, height, drawingoffsetx, drawingoffsety, drawingwidth, drawingheight, legendoffsetx, legendoffsety, foldingx=[], foldingy=[]):
        self.ref = ref
        self.massstab = massstab
        self.width = width
        self.height = height
        self.drawingoffsetx = drawingoffsetx
        self.drawingoffsety = drawingoffsety
        self.drawingwidth = drawingwidth
        self.drawingheight = drawingheight
        self.legendoffsetx = legendoffsetx
        self.legendoffsety = legendoffsety
        self.foldingx = foldingx
        self.foldingy = foldingy

    def get_drawing_offset(self):
        return (self.drawingoffsetx + (self.drawingwidth - ref.max_width() * 1000 / self.massstab) / 2, self.drawingoffsety + (self.drawingheight - ref.max_height() * 1000 / self.massstab) / 2)

    def get_drawing_point(self, point):
        x, y = ref.get_xy(point)
        offsetx, offsety = self.get_drawing_offset()
        return (offsetx + x * 1000 / self.massstab, offsety + y * 1000 / self.massstab)

    def match(self, ref):
        return self.drawingwidth > ref.max_width() * 1000 / self.massstab and self.drawingheight > ref.max_height() * 1000 / self.massstab

    def svg_points(self, coords):
        return " ".join(["{},{}".format(x, y) for x, y in coords])

    def move_points(self, delta_x, delta_y, points):
        return [(x + delta_x, y + delta_y) for x, y in points]

    def scale_points(self, scale, points):
        return [(x * scale, y * scale) for x, y in points]

    def draw_folding(self, output):
        def draw_line(x1, y1, x2, y2):
            output.write('<line x1="{}" y1="{}" x2="{}" y2="{}" style="stroke:black;stroke-width:0.5" />'.format(x1, y1, x2, y2))

        for x in self.foldingx:
            draw_line(x, 0, x, 5)
            draw_line(x, self.height - 5, x, self.height)

        for y in self.foldingy:
            draw_line(0, y, 5, y)
            draw_line(self.width - 5, y, self.width, y)

    def draw_north(self, output, pos_x, pos_y, size=10, a=0.4):
        points = [(0,1), (0, -1 + a * np.sin(np.pi/3)), (a * np.cos(np.pi/3), -1 + a * np.sin(np.pi/3)), (0,-1), (-a * np.cos(np.pi/3),-1 + a * np.sin(np.pi/3)), (0,-1 + a * np.sin(np.pi/3))]
        output.write('<g transform="rotate({}, {}, {})">'.format(self.ref.get_north(), pos_x + size, pos_y + size))
        output.write('<text x="{}" y="{}" style="font-family:sans-serif; font-size:{}; text-anchor:middle">N</text>'.format(pos_x + size, pos_y + size * 1.5, size))
        output.write('<polygon points="{}" style="fill:black;stroke:black;stroke-width:{}" />'.format(self.svg_points(self.move_points(pos_x + size, pos_y + size, self.scale_points(size, points))), size / 10))
        output.write('</g>')

    def draw_scale(self, output, x, y, max_width=100, scale_count=5, scale_height=5, font_family="HUN-din 1451", font_size=5):
        round_digits = lambda x, digits: round(x, digits) if digits > 0 else int(x)
        max_scale_width = max_width / scale_count
        digits = np.floor(max_scale_width / self.massstab)
        scale_unit = round_digits(max_scale_width * self.massstab / 1000, digits)
        scale_width = scale_unit * 1000 / self.massstab

        if scale_unit >= 1000:
            display_divisor, display_unit = (1000, "km")
        else:
            display_divisor, display_unit = (1, "m")

        for i in range(0,scale_count):
            output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:middle">{}</text>'.format(x + scale_width * i, y + font_size + font_size + scale_height + font_size, font_family, font_size, str(round_digits(scale_unit * i / display_divisor, digits))))
            output.write('<rect x="{}" y="{}" width="{}" height="{}" style="fill:{};stroke:black;stroke-width:0.5" />'.format(x + scale_width * i, y + font_size + font_size, scale_width, scale_height, ["white", "black"][i % 2]))
        output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:middle">{}</text>'.format(x + scale_width * scale_count, y + font_size + font_size + scale_height + font_size, font_family, font_size, str(round_digits(scale_unit * scale_count / display_divisor, digits)) + " " + display_unit))
        output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:middle" xml:space="preserve">{}</text>'.format(x + scale_width * scale_count / 2, y + font_size, font_family, font_size, "Maßstab 1 : {}     1 cm ≙ {} {}".format(self.massstab, self.massstab / 100 / display_divisor, display_unit)))

    def draw_titleblock(self, output, x, y, lines, width=180, font_family="HUN-din 1451", font_size=10):
        height = (font_size + 2) * len(lines)

        output.write('<rect x="{}" y="{}" width="{}" height="{}" style="fill:none;stroke:black;stroke-width:2" />'.format(x, y, width, height))

        for i in range(0, len(lines)):
            output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}">{}</text>'.format(x + 2, y - 2 + (font_size + 2) * (i + 1), font_family, font_size, lines[i]))
            output.write('<line x1="{}" y1="{}" x2="{}" y2="{}" style="stroke:black;stroke-width:1" />'.format(x, y + (font_size + 2) * (i + 1), x + width, y + (font_size + 2) * (i + 1)))

    def draw_point(self, output, point, radius=1, text_position="topleft", font_family="HUN-din 1451", font_size=7):
        x, y = self.get_drawing_point(point.coordinate)
        output.write('<circle cx="{}" cy="{}" r="{}" style="fill:black;stroke:black;stroke-width:1" />'.format(x, y, radius))
        if text_position == "topleft":
            output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:end">{}</text>'.format(x - radius, y - radius, font_family, font_size, point.label))
        elif text_position == "topright":
            output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:start">{}</text>'.format(x + radius, y - radius, font_family, font_size, point.label))
        elif text_position == "bottomleft":
            output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:end">{}</text>'.format(x - radius, y + font_size + radius, font_family, font_size, point.label))
        elif text_position == "bottomright":
            output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:start">{}</text>'.format(x + radius, y + font_size + radius, font_family, font_size, point.label))

    def draw_linestring(self, output, linestring, font_family="HUN-din 1451", font_size=7):
        output.write('<polyline points="{}" style="fill:none;stroke:black;stroke-width:1" />'.format(self.svg_points([self.get_drawing_point(coord) for coord in linestring.coordinates])))
        perimeter = 0
        for i in range(1, len(linestring.coordinates)):
            a = linestring.coordinates[(i-1) % len(linestring.coordinates)]
            b = linestring.coordinates[i % len(linestring.coordinates)]
            a_x, a_y = ref.get_xy(a)
            b_x, b_y = ref.get_xy(b)
            perimeter += np.sqrt((a_x - b_x)**2 + (a_y - b_y)**2)
        pos_x, pos_y = self.get_drawing_point(linestring.coordinates[0])
        output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:middle">{} m</text>'.format(pos_x, pos_y, font_family, font_size, round(perimeter)))

    def draw_polygon(self, output, polygon, font_family="HUN-din 1451", font_size=7):
        output.write('<polygon points="{}" style="fill:none;stroke:black;stroke-width:1" />'.format(self.svg_points([self.get_drawing_point(coord) for coord in polygon.coordinates])))
        area = 0
        perimeter = 0
        for i in range(0, len(polygon.coordinates)):
            a = polygon.coordinates[(i-1) % len(polygon.coordinates)]
            b = polygon.coordinates[i % len(polygon.coordinates)]
            a_x, a_y = ref.get_xy(a)
            b_x, b_y = ref.get_xy(b)
            area += a_x * b_y - a_y * b_x
            perimeter += np.sqrt((a_x - b_x)**2 + (a_y - b_y)**2)
        area = np.abs(area / 2) / 10000
        pos_x, pos_y = self.get_drawing_point(MapPoint(sum([c.lat for c in polygon.coordinates]) / len(polygon.coordinates), sum([c.lon for c in polygon.coordinates]) / len(polygon.coordinates)))
        output.write('<text x="{}" y="{}" style="font-family:\'{}\'; font-size:{}; text-anchor:middle">{} ha / {} m</text>'.format(pos_x, pos_y, font_family, font_size, round(area, 4), round(perimeter)))

    def render(self, figures, output):
        # viewBox for polygons etc
        output.write('<svg width="{}mm" height="{}mm" viewBox="0 0 {} {}">'.format(paper.width, paper.height, paper.width, paper.height))

        if False:
            output.write('<rect x="{}" y="{}" width="{}" height="{}" style="fill:none;stroke:black;stroke-width:1" />'.format(self.drawingoffsetx, self.drawingoffsety, self.drawingwidth, self.drawingheight))
            offsetx, offsety = self.get_drawing_offset()
            output.write('<rect x="{}" y="{}" width="{}" height="{}" style="fill:none;stroke:black;stroke-width:1;stroke-dasharray:10,10" />'.format(offsetx, offsety, ref.max_width() * 1000 / self.massstab, ref.max_height() * 1000 / self.massstab))

        for figure in figures:
            if isinstance(figure, KmlPoint):
                self.draw_point(output, figure)
            elif isinstance(figure, KmlLineString):
                self.draw_linestring(output, figure)
            elif isinstance(figure, KmlPolygon):
                self.draw_polygon(output, figure)
            else:
                raise Exception("Unsupported figure")

        self.draw_folding(output)
        self.draw_scale(output, self.legendoffsetx, self.legendoffsety)
        self.draw_north(output, self.legendoffsetx + 110, self.legendoffsety)

        import datetime
        now = datetime.datetime.now()
        lines = [
          name,
          "GPS-Koordinaten: {}, {}".format(round(sum([c.lat for c in all_coords]) / len(all_coords), 5), round(sum([c.lon for c in all_coords]) / len(all_coords), 5)),
          "Projekt: ",
          "Datum: {}".format(now.strftime("%d.%m.%Y"))
        ]
        self.draw_titleblock(output, self.width - 180 - 10, self.height - 50 - 10, lines)

        output.write('</svg>')

def select_paper(ref, massstab):
    # width, height, drawingoffsetx, drawingoffsety, drawingwidth, drawingheight, legendoffsetx, legendoffsety, foldingx, foldingy
    # MUST be sorted (smallest to largest)
    papersizes = [
      Paper(ref, massstab,  210, 297, 20,10,  180, 187, 20, 207, [], []), # a4
      Paper(ref, massstab,  420, 297, 20,10,  390, 217, 20, 247, [125, 125+105], []), # a3
      Paper(ref, massstab,  594, 420, 25,10,  559, 340, 25, 370, [210, 210+192], [420-297]), # a2
      Paper(ref, massstab,  841, 594, 25,10,  806, 514, 25, 544, [210, 210+190, 210+190+125, 841-190], [297]), # a1
      Paper(ref, massstab, 1189, 841, 25,10, 1154, 761, 25, 791, [210, 210+190, 210+190+190, 210+190+190+190, 210+190+190+190+110, 1189-190], [841-297-297, 841-297]), # a0
      Paper(ref, massstab, 1682,1189, 25,10, 1652,1109, 25,1139, [186, 1682-7*190, 1682-6*190, 1682-5*190, 1682-4*190, 1682-3*190, 1682-2*190, 1682-190], [1189-3*297, 1189-2*297, 1189-297]), # 2a0
    ]

    for paper in papersizes:
        if paper.match(ref):
            return paper

    # if no generic papersize could be found we just generate our own
    drawingwidth, drawingheight = np.ceil(ref.max_width() * 1000 / massstab), np.ceil(ref.max_height() * 1000 / massstab)
    return Paper(ref, massstab, 10 + drawingwidth + 10, 10 + drawingheight + 10 + 50 + 10, 10, 10, drawingwidth, drawingheight, 10, 10 + drawingheight + 10)

kml_file = KmlFile(args.input)
name = kml_file.get_name()
figures = kml_file.get_figures()

all_coords = sum([figure.get_coordinates() for figure in figures], [])
lats = [coord.lat for coord in all_coords]
lons = [coord.lon for coord in all_coords]
ref = MapReference(MapPoint(max(lats), min(lons)), MapPoint(max(lats), max(lons)))
ref.orientate(all_coords)

paper = select_paper(ref, args.massstab)
paper.render(figures, args.output)