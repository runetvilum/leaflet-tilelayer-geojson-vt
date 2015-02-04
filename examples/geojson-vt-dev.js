!function(e){if("object"==typeof exports&&"undefined"!=typeof module)module.exports=e();else if("function"==typeof define&&define.amd)define([],e);else{var o;"undefined"!=typeof window?o=window:"undefined"!=typeof global?o=global:"undefined"!=typeof self&&(o=self),o.geojsonvt=e()}}(function(){var define,module,exports;return (function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
'use strict';

module.exports = clip;

/* clip features between two axis-parallel lines:
 *     |        |
 *  ___|___     |     /
 * /   |   \____|____/
 *     |        |
 */

function clip(features, scale, k1, k2, axis, intersect) {

    k1 /= scale;
    k2 /= scale;

    var clipped = [];

    for (var i = 0; i < features.length; i++) {

        var feature = features[i],
            geometry = feature.geometry,
            type = feature.type,
            min, max;

        if (feature.min) {
            min = feature.min[axis];
            max = feature.max[axis];

            if (min >= k1 && max <= k2) { // trivial accept
                clipped.push(feature);
                continue;
            } else if (min > k2 || max < k1) continue; // trivial reject
        }

        var slices = type === 1 ?
                clipPoints(geometry, k1, k2, axis) :
                clipGeometry(geometry, k1, k2, axis, intersect, type === 3);

        if (slices.length) {
            // if a feature got clipped, it will likely get clipped on the next zoom level as well,
            // so there's no need to recalculate bboxes
            clipped.push({
                geometry: slices,
                type: type,
                tags: features[i].tags || null
            });
        }
    }

    return clipped.length ? clipped : null;
}

function clipPoints(geometry, k1, k2, axis) {
    var slice = [];

    for (var i = 0; i < geometry.length; i++) {
        var a = geometry[i],
            ak = a[axis];

        if (ak >= k1 && ak <= k2) slice.push(a);
    }
    return slice;
}

function clipGeometry(geometry, k1, k2, axis, intersect, closed) {

    var slices = [];

    for (var i = 0; i < geometry.length; i++) {

        var ak = 0,
            bk = 0,
            b = null,
            points = geometry[i],
            area = points.area,
            dist = points.dist,
            len = points.length,
            a, j;

        var slice = [];

        for (j = 0; j < len - 1; j++) {
            a = b || points[j];
            b = points[j + 1];
            ak = bk || a[axis];
            bk = b[axis];

            if (ak < k1) {

                if ((bk > k2)) { // ---|-----|-->
                    slice.push(intersect(a, b, k1), intersect(a, b, k2));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk >= k1) slice.push(intersect(a, b, k1)); // ---|-->  |

            } else if (ak > k2) {

                if ((bk < k1)) { // <--|-----|---
                    slice.push(intersect(a, b, k2), intersect(a, b, k1));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk <= k2) slice.push(intersect(a, b, k2)); // |  <--|---

            } else {

                slice.push(a);

                if (bk < k1) { // <--|---  |
                    slice.push(intersect(a, b, k1));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk > k2) { // |  ---|-->
                    slice.push(intersect(a, b, k2));
                    if (!closed) slice = newSlice(slices, slice, area, dist);
                }
                // | --> |
            }
        }

        // add the last point
        a = points[len - 1];
        ak = a[axis];
        if (ak >= k1 && ak <= k2) slice.push(a);

        // close the polygon if its endpoints are not the same after clipping
        if (closed && slice[0] !== slice[slice.length - 1]) slice.push(slice[0]);

        // add the final slice
        newSlice(slices, slice, area, dist);
    }

    return slices;
}

function newSlice(slices, slice, area, dist) {
    if (slice.length) {
        // we don't recalculate the area/length of the unclipped geometry because the case where it goes
        // below the visibility threshold as a result of clipping is rare, so we avoid doing unnecessary work
        slice.area = area;
        slice.dist = dist;

        slices.push(slice);
    }
    return [];
}

},{}],2:[function(require,module,exports){
'use strict';

module.exports = convert;

var simplify = require('./simplify');

// converts GeoJSON feature into an intermediate projected JSON vector format with simplification data

function convert(data, tolerance) {
    var features = [];

    if (data.type === 'FeatureCollection') {
        for (var i = 0; i < data.features.length; i++) {
            convertFeature(features, data.features[i], tolerance);
        }
    } else if (data.type === 'Feature') {
        convertFeature(features, data, tolerance);

    } else {
        // single geometry or a geometry collection
        convertFeature(features, {geometry: data}, tolerance);
    }
    return features;
}

function convertFeature(features, feature, tolerance) {
    var geom = feature.geometry,
        type = geom.type,
        coords = geom.coordinates,
        tags = feature.properties,
        i, j, rings;

    if (type === 'Point') {
        features.push(create(tags, 1, [projectPoint(coords)]));

    } else if (type === 'MultiPoint') {
        features.push(create(tags, 1, project(coords)));

    } else if (type === 'LineString') {
        features.push(create(tags, 2, [project(coords, tolerance)]));

    } else if (type === 'MultiLineString' || type === 'Polygon') {
        rings = [];
        for (i = 0; i < coords.length; i++) {
            rings.push(project(coords[i], tolerance));
        }
        features.push(create(tags, type === 'Polygon' ? 3 : 2, rings));

    } else if (type === 'MultiPolygon') {
        rings = [];
        for (i = 0; i < coords.length; i++) {
            for (j = 0; j < coords[i].length; j++) {
                rings.push(project(coords[i][j], tolerance));
            }
        }
        features.push(create(tags, 3, rings));

    } else if (type === 'GeometryCollection') {
        for (i = 0; i < geom.geometries.length; i++) {
            convertFeature(features, {
                geometry: geom.geometries[i],
                properties: tags
            }, tolerance);
        }

    } else {
        console.warn('Unsupported GeoJSON type: ' + geom.type);
    }
}

function create(tags, type, geometry) {
    var feature = {
        geometry: geometry,
        type: type,
        tags: tags || null,
        min: [1, 1], // initial bbox values;
        max: [0, 0]  // note that all coords are in [0..1] range
    };
    calcBBox(feature);
    return feature;
}

function project(lonlats, tolerance) {
    var projected = [];
    for (var i = 0; i < lonlats.length; i++) {
        projected.push(projectPoint(lonlats[i]));
    }
    if (tolerance) {
        simplify(projected, tolerance);
        calcSize(projected);
    }
    return projected;
}

function projectPoint(p) {
    var sin = Math.sin(p[1] * Math.PI / 180),
        x = (p[0] / 360 + 0.5),
        y = (0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI);
    return [x, y, 0];
}

// calculate area and length of the poly
function calcSize(points) {
    var area = 0,
        dist = 0;

    for (var i = 0, a, b; i < points.length - 1; i++) {
        a = b || points[i];
        b = points[i + 1];

        area += a[0] * b[1] - b[0] * a[1];

        // use Manhattan distance instead of Euclidian one to avoid expensive square root computation
        dist += Math.abs(b[0] - a[0]) + Math.abs(b[1] - a[1]);
    }
    points.area = Math.abs(area / 2);
    points.dist = dist;
}

// calculate the feature bounding box for faster clipping later
function calcBBox(feature) {
    var geometry = feature.geometry,
        min = feature.min,
        max = feature.max;

    if (feature.type === 1) {
        calcRingBBox(min, max, geometry);
    } else {
        for (var i = 0; i < geometry.length; i++) {
            calcRingBBox(min, max, geometry[i]);
        }
    }
    return feature;
}

function calcRingBBox(min, max, points) {
    for (var i = 0, p; i < points.length; i++) {
        p = points[i];
        min[0] = Math.min(p[0], min[0]);
        max[0] = Math.max(p[0], max[0]);
        min[1] = Math.min(p[1], min[1]);
        max[1] = Math.max(p[1], max[1]);
    }
}

},{"./simplify":4}],3:[function(require,module,exports){
'use strict';

module.exports = geojsonvt;

var convert = require('./convert'), // GeoJSON conversion and preprocessing
    clip = require('./clip'),       // stripe clipping algorithm
    createTile = require('./tile'); // final simplified tile generation


function geojsonvt(data, options) {
    return new GeoJSONVT(data, options);
}

function GeoJSONVT(data, options) {
    options = this.options = extend(Object.create(this.options), options);

    var debug = options.debug;

    if (debug) console.time('preprocess data');

    var z2 = 1 << options.baseZoom, // 2^z
        features = convert(data, options.tolerance / (z2 * options.extent));

    this.tiles = {};

    if (debug) {
        console.timeEnd('preprocess data');
        console.time('generate tiles up to z' + options.maxZoom);
        this.stats = {};
        this.total = 0;
    }

    // start slicing from the top tile down
    this.splitTile(features, 0, 0, 0);

    if (debug) {
        console.log('features: %d, points: %d', this.tiles[0].numFeatures, this.tiles[0].numPoints);
        console.timeEnd('generate tiles up to z' + options.maxZoom);
        console.log('tiles generated:', this.total, JSON.stringify(this.stats));
    }
}

GeoJSONVT.prototype.options = {
    baseZoom: 14,   // max zoom to preserve detail on
    maxZoom: 4,     // zoom to slice down to on first pass
    maxPoints: 100, // stop slicing a tile below this number of points
    tolerance: 3,   // simplification tolerance (higher means simpler)
    extent: 4096,   // tile extent
    buffer: 64,     // tile buffer on each side
    debug: 0        // logging level (0, 1 or 2)
};

GeoJSONVT.prototype.splitTile = function (features, z, x, y, cz, cx, cy) {

    var stack = [features, z, x, y],
        options = this.options,
        debug = options.debug,
        extent = options.extent,
        buffer = options.buffer;

    // avoid recursion by using a processing queue
    while (stack.length) {
        features = stack.shift();
        z = stack.shift();
        x = stack.shift();
        y = stack.shift();

        var z2 = 1 << z,
            id = toID(z, x, y),
            tile = this.tiles[id],
            tileTolerance = z === options.baseZoom ? 0 : options.tolerance / (z2 * extent);

        if (!tile) {
            if (debug > 1) console.time('creation');

            tile = this.tiles[id] = createTile(features, z2, x, y, tileTolerance, extent, z === options.baseZoom);

            if (debug) {
                if (debug > 1) {
                    console.log('tile z%d-%d-%d (features: %d, points: %d, simplified: %d)',
                        z, x, y, tile.numFeatures, tile.numPoints, tile.numSimplified);
                    console.timeEnd('creation');
                }
                var key = 'z' + z + ':';
                this.stats[key] = (this.stats[key] || 0) + 1;
                this.total++;
            }
        }

        if (!cz && (z === options.maxZoom || tile.numPoints <= options.maxPoints ||
                isClippedSquare(tile.features, extent, buffer)) || z === options.baseZoom || z === cz) {
            tile.source = features;
            continue; // stop tiling
        }

        if (cz) tile.source = features;
        else tile.source = null;

        if (debug > 1) console.time('clipping');

        // values we'll use for clipping
        var k1 = 0.5 * buffer / extent,
            k2 = 0.5 - k1,
            k3 = 0.5 + k1,
            k4 = 1 + k1,

            tl, bl, tr, br, left, right,
            m, goLeft, goTop;

        if (cz) { // if we have a specific tile to drill down to, calculate where to go
            m = 1 << (cz - z);
            goLeft = cx / m - x < 0.5;
            goTop = cy / m - y < 0.5;
        }

        tl = bl = tr = br = left = right = null;

        if (!cz ||  goLeft) left  = clip(features, z2, x - k1, x + k3, 0, intersectX);
        if (!cz || !goLeft) right = clip(features, z2, x + k2, x + k4, 0, intersectX);

        if (left) {
            if (!cz ||  goTop) tl = clip(left, z2, y - k1, y + k3, 1, intersectY);
            if (!cz || !goTop) bl = clip(left, z2, y + k2, y + k4, 1, intersectY);
        }

        if (right) {
            if (!cz ||  goTop) tr = clip(right, z2, y - k1, y + k3, 1, intersectY);
            if (!cz || !goTop) br = clip(right, z2, y + k2, y + k4, 1, intersectY);
        }

        if (debug > 1) console.timeEnd('clipping');

        if (tl) stack.push(tl, z + 1, x * 2,     y * 2);
        if (bl) stack.push(bl, z + 1, x * 2,     y * 2 + 1);
        if (tr) stack.push(tr, z + 1, x * 2 + 1, y * 2);
        if (br) stack.push(br, z + 1, x * 2 + 1, y * 2 + 1);
    }
};

GeoJSONVT.prototype.getTile = function (z, x, y) {
    var id = toID(z, x, y);
    if (this.tiles[id]) return this.tiles[id];

    var options = this.options,
        debug = options.debug;

    if (debug > 1) console.log('drilling down to z%d-%d-%d', z, x, y);

    var z0 = z,
        x0 = x,
        y0 = y,
        parent;

    while (!parent && z0 > 0) {
        z0--;
        x0 = Math.floor(x0 / 2);
        y0 = Math.floor(y0 / 2);
        parent = this.tiles[toID(z0, x0, y0)];
    }

    if (debug > 1) console.log('found parent tile z%d-%d-%d', z0, x0, y0);

    // if we found a parent tile containing the original geometry, we can drill down from it
    if (parent.source) {
        if (isClippedSquare(parent.features, options.extent, options.buffer)) return parent;

        if (debug) console.time('drilling down');
        this.splitTile(parent.source, z0, x0, y0, z, x, y);
        if (debug) console.timeEnd('drilling down');
    }

    return this.tiles[id];
};

// checks whether a tile is a whole-area fill after clipping; if it is, there's no sense slicing it further
function isClippedSquare(features, extent, buffer) {
    if (features.length !== 1) return false;

    var feature = features[0];
    if (feature.type !== 3 || feature.geometry.length > 1) return false;

    for (var i = 0; i < feature.geometry[0].length; i++) {
        var p = feature.geometry[0][i];
        if ((p[0] !== -buffer && p[0] !== extent + buffer) ||
            (p[1] !== -buffer && p[1] !== extent + buffer)) return false;
    }
    return true;
}

function toID(z, x, y) {
    return (((1 << z) * y + x) * 32) + z;
}

function intersectX(a, b, x) {
    return [x, (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0]) + a[1], 1];
}
function intersectY(a, b, y) {
    return [(y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) + a[0], y, 1];
}

function extend(dest, src) {
    for (var i in src) {
        dest[i] = src[i];
    }
    return dest;
}

},{"./clip":1,"./convert":2,"./tile":5}],4:[function(require,module,exports){
'use strict';

module.exports = simplify;

// calculate simplification data using optimized Douglas-Peucker algorithm

function simplify(points, tolerance) {

    var sqTolerance = tolerance * tolerance,
        len = points.length,
        first = 0,
        last = len - 1,
        stack = [],
        i, maxSqDist, sqDist, index;

    // always retain the endpoints (1 is the max value)
    points[first][2] = 1;
    points[last][2] = 1;

    // avoid recursion by using a stack
    while (last) {

        maxSqDist = 0;

        for (i = first + 1; i < last; i++) {
            sqDist = getSqSegDist(points[i], points[first], points[last]);

            if (sqDist > maxSqDist) {
                index = i;
                maxSqDist = sqDist;
            }
        }

        if (maxSqDist > sqTolerance) {
            points[index][2] = maxSqDist; // save the point importance in squared pixels as a z coordinate
            stack.push(first, index, index, last);
        }

        last = stack.pop();
        first = stack.pop();
    }
}

// square distance from a point to a segment
function getSqSegDist(p, a, b) {

    var x = a[0], y = a[1],
        bx = b[0], by = b[1],
        px = p[0], py = p[1],
        dx = bx - x,
        dy = by - y;

    if (dx !== 0 || dy !== 0) {

        var t = ((px - x) * dx + (py - y) * dy) / (dx * dx + dy * dy);

        if (t > 1) {
            x = bx;
            y = by;

        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }

    dx = px - x;
    dy = py - y;

    return dx * dx + dy * dy;
}

},{}],5:[function(require,module,exports){
'use strict';

module.exports = createTile;

function createTile(features, z2, tx, ty, tolerance, extent, noSimplify) {
    var tile = {
        features: [],
        numPoints: 0,
        numSimplified: 0,
        numFeatures: 0,
        source: null
    };
    for (var i = 0; i < features.length; i++) {
        tile.numFeatures++;
        addFeature(tile, features[i], z2, tx, ty, tolerance, extent, noSimplify);
    }
    return tile;
}

function addFeature(tile, feature, z2, tx, ty, tolerance, extent, noSimplify) {

    var geom = feature.geometry,
        type = feature.type,
        transformed = [],
        sqTolerance = tolerance * tolerance,
        i, j, ring, p;

    if (type === 1) {
        for (i = 0; i < geom.length; i++) {
            transformed.push(transformPoint(geom[i], z2, tx, ty, extent));
            tile.numPoints++;
            tile.numSimplified++;
        }

    } else {

        // simplify and transform projected coordinates for tile geometry
        for (i = 0; i < geom.length; i++) {
            ring = geom[i];

            // filter out tiny polylines & polygons
            if (!noSimplify && ((type === 2 && ring.dist < tolerance) ||
                                (type === 3 && ring.area < sqTolerance))) {
                tile.numPoints += ring.length;
                continue;
            }

            var transformedRing = [];

            for (j = 0; j < ring.length; j++) {
                p = ring[j];
                // keep points with importance > tolerance
                if (noSimplify || p[2] > sqTolerance) {
                    transformedRing.push(transformPoint(p, z2, tx, ty, extent));
                    tile.numSimplified++;
                }
                tile.numPoints++;
            }

            transformed.push(transformedRing);
        }
    }

    if (transformed.length) {
        tile.features.push({
            geometry: transformed,
            type: type,
            tags: feature.tags || null
        });
    }
}

function transformPoint(p, z2, tx, ty, extent) {
    var x = Math.round(extent * (p[0] * z2 - tx)),
        y = Math.round(extent * (p[1] * z2 - ty));
    return [x, y];
}

},{}]},{},[3])(3)
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJzcmMvY2xpcC5qcyIsInNyYy9jb252ZXJ0LmpzIiwic3JjL2luZGV4LmpzIiwic3JjL3NpbXBsaWZ5LmpzIiwic3JjL3RpbGUuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbEpBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOU1BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2RUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3ZhciBmPW5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIik7dGhyb3cgZi5jb2RlPVwiTU9EVUxFX05PVF9GT1VORFwiLGZ9dmFyIGw9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGwuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sbCxsLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIid1c2Ugc3RyaWN0JztcblxubW9kdWxlLmV4cG9ydHMgPSBjbGlwO1xuXG4vKiBjbGlwIGZlYXR1cmVzIGJldHdlZW4gdHdvIGF4aXMtcGFyYWxsZWwgbGluZXM6XG4gKiAgICAgfCAgICAgICAgfFxuICogIF9fX3xfX18gICAgIHwgICAgIC9cbiAqIC8gICB8ICAgXFxfX19ffF9fX18vXG4gKiAgICAgfCAgICAgICAgfFxuICovXG5cbmZ1bmN0aW9uIGNsaXAoZmVhdHVyZXMsIHNjYWxlLCBrMSwgazIsIGF4aXMsIGludGVyc2VjdCkge1xuXG4gICAgazEgLz0gc2NhbGU7XG4gICAgazIgLz0gc2NhbGU7XG5cbiAgICB2YXIgY2xpcHBlZCA9IFtdO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBmZWF0dXJlcy5sZW5ndGg7IGkrKykge1xuXG4gICAgICAgIHZhciBmZWF0dXJlID0gZmVhdHVyZXNbaV0sXG4gICAgICAgICAgICBnZW9tZXRyeSA9IGZlYXR1cmUuZ2VvbWV0cnksXG4gICAgICAgICAgICB0eXBlID0gZmVhdHVyZS50eXBlLFxuICAgICAgICAgICAgbWluLCBtYXg7XG5cbiAgICAgICAgaWYgKGZlYXR1cmUubWluKSB7XG4gICAgICAgICAgICBtaW4gPSBmZWF0dXJlLm1pbltheGlzXTtcbiAgICAgICAgICAgIG1heCA9IGZlYXR1cmUubWF4W2F4aXNdO1xuXG4gICAgICAgICAgICBpZiAobWluID49IGsxICYmIG1heCA8PSBrMikgeyAvLyB0cml2aWFsIGFjY2VwdFxuICAgICAgICAgICAgICAgIGNsaXBwZWQucHVzaChmZWF0dXJlKTtcbiAgICAgICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgICAgIH0gZWxzZSBpZiAobWluID4gazIgfHwgbWF4IDwgazEpIGNvbnRpbnVlOyAvLyB0cml2aWFsIHJlamVjdFxuICAgICAgICB9XG5cbiAgICAgICAgdmFyIHNsaWNlcyA9IHR5cGUgPT09IDEgP1xuICAgICAgICAgICAgICAgIGNsaXBQb2ludHMoZ2VvbWV0cnksIGsxLCBrMiwgYXhpcykgOlxuICAgICAgICAgICAgICAgIGNsaXBHZW9tZXRyeShnZW9tZXRyeSwgazEsIGsyLCBheGlzLCBpbnRlcnNlY3QsIHR5cGUgPT09IDMpO1xuXG4gICAgICAgIGlmIChzbGljZXMubGVuZ3RoKSB7XG4gICAgICAgICAgICAvLyBpZiBhIGZlYXR1cmUgZ290IGNsaXBwZWQsIGl0IHdpbGwgbGlrZWx5IGdldCBjbGlwcGVkIG9uIHRoZSBuZXh0IHpvb20gbGV2ZWwgYXMgd2VsbCxcbiAgICAgICAgICAgIC8vIHNvIHRoZXJlJ3Mgbm8gbmVlZCB0byByZWNhbGN1bGF0ZSBiYm94ZXNcbiAgICAgICAgICAgIGNsaXBwZWQucHVzaCh7XG4gICAgICAgICAgICAgICAgZ2VvbWV0cnk6IHNsaWNlcyxcbiAgICAgICAgICAgICAgICB0eXBlOiB0eXBlLFxuICAgICAgICAgICAgICAgIHRhZ3M6IGZlYXR1cmVzW2ldLnRhZ3MgfHwgbnVsbFxuICAgICAgICAgICAgfSk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICByZXR1cm4gY2xpcHBlZC5sZW5ndGggPyBjbGlwcGVkIDogbnVsbDtcbn1cblxuZnVuY3Rpb24gY2xpcFBvaW50cyhnZW9tZXRyeSwgazEsIGsyLCBheGlzKSB7XG4gICAgdmFyIHNsaWNlID0gW107XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGdlb21ldHJ5Lmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIHZhciBhID0gZ2VvbWV0cnlbaV0sXG4gICAgICAgICAgICBhayA9IGFbYXhpc107XG5cbiAgICAgICAgaWYgKGFrID49IGsxICYmIGFrIDw9IGsyKSBzbGljZS5wdXNoKGEpO1xuICAgIH1cbiAgICByZXR1cm4gc2xpY2U7XG59XG5cbmZ1bmN0aW9uIGNsaXBHZW9tZXRyeShnZW9tZXRyeSwgazEsIGsyLCBheGlzLCBpbnRlcnNlY3QsIGNsb3NlZCkge1xuXG4gICAgdmFyIHNsaWNlcyA9IFtdO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBnZW9tZXRyeS5sZW5ndGg7IGkrKykge1xuXG4gICAgICAgIHZhciBhayA9IDAsXG4gICAgICAgICAgICBiayA9IDAsXG4gICAgICAgICAgICBiID0gbnVsbCxcbiAgICAgICAgICAgIHBvaW50cyA9IGdlb21ldHJ5W2ldLFxuICAgICAgICAgICAgYXJlYSA9IHBvaW50cy5hcmVhLFxuICAgICAgICAgICAgZGlzdCA9IHBvaW50cy5kaXN0LFxuICAgICAgICAgICAgbGVuID0gcG9pbnRzLmxlbmd0aCxcbiAgICAgICAgICAgIGEsIGo7XG5cbiAgICAgICAgdmFyIHNsaWNlID0gW107XG5cbiAgICAgICAgZm9yIChqID0gMDsgaiA8IGxlbiAtIDE7IGorKykge1xuICAgICAgICAgICAgYSA9IGIgfHwgcG9pbnRzW2pdO1xuICAgICAgICAgICAgYiA9IHBvaW50c1tqICsgMV07XG4gICAgICAgICAgICBhayA9IGJrIHx8IGFbYXhpc107XG4gICAgICAgICAgICBiayA9IGJbYXhpc107XG5cbiAgICAgICAgICAgIGlmIChhayA8IGsxKSB7XG5cbiAgICAgICAgICAgICAgICBpZiAoKGJrID4gazIpKSB7IC8vIC0tLXwtLS0tLXwtLT5cbiAgICAgICAgICAgICAgICAgICAgc2xpY2UucHVzaChpbnRlcnNlY3QoYSwgYiwgazEpLCBpbnRlcnNlY3QoYSwgYiwgazIpKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKCFjbG9zZWQpIHNsaWNlID0gbmV3U2xpY2Uoc2xpY2VzLCBzbGljZSwgYXJlYSwgZGlzdCk7XG5cbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKGJrID49IGsxKSBzbGljZS5wdXNoKGludGVyc2VjdChhLCBiLCBrMSkpOyAvLyAtLS18LS0+ICB8XG5cbiAgICAgICAgICAgIH0gZWxzZSBpZiAoYWsgPiBrMikge1xuXG4gICAgICAgICAgICAgICAgaWYgKChiayA8IGsxKSkgeyAvLyA8LS18LS0tLS18LS0tXG4gICAgICAgICAgICAgICAgICAgIHNsaWNlLnB1c2goaW50ZXJzZWN0KGEsIGIsIGsyKSwgaW50ZXJzZWN0KGEsIGIsIGsxKSk7XG4gICAgICAgICAgICAgICAgICAgIGlmICghY2xvc2VkKSBzbGljZSA9IG5ld1NsaWNlKHNsaWNlcywgc2xpY2UsIGFyZWEsIGRpc3QpO1xuXG4gICAgICAgICAgICAgICAgfSBlbHNlIGlmIChiayA8PSBrMikgc2xpY2UucHVzaChpbnRlcnNlY3QoYSwgYiwgazIpKTsgLy8gfCAgPC0tfC0tLVxuXG4gICAgICAgICAgICB9IGVsc2Uge1xuXG4gICAgICAgICAgICAgICAgc2xpY2UucHVzaChhKTtcblxuICAgICAgICAgICAgICAgIGlmIChiayA8IGsxKSB7IC8vIDwtLXwtLS0gIHxcbiAgICAgICAgICAgICAgICAgICAgc2xpY2UucHVzaChpbnRlcnNlY3QoYSwgYiwgazEpKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKCFjbG9zZWQpIHNsaWNlID0gbmV3U2xpY2Uoc2xpY2VzLCBzbGljZSwgYXJlYSwgZGlzdCk7XG5cbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKGJrID4gazIpIHsgLy8gfCAgLS0tfC0tPlxuICAgICAgICAgICAgICAgICAgICBzbGljZS5wdXNoKGludGVyc2VjdChhLCBiLCBrMikpO1xuICAgICAgICAgICAgICAgICAgICBpZiAoIWNsb3NlZCkgc2xpY2UgPSBuZXdTbGljZShzbGljZXMsIHNsaWNlLCBhcmVhLCBkaXN0KTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgLy8gfCAtLT4gfFxuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgLy8gYWRkIHRoZSBsYXN0IHBvaW50XG4gICAgICAgIGEgPSBwb2ludHNbbGVuIC0gMV07XG4gICAgICAgIGFrID0gYVtheGlzXTtcbiAgICAgICAgaWYgKGFrID49IGsxICYmIGFrIDw9IGsyKSBzbGljZS5wdXNoKGEpO1xuXG4gICAgICAgIC8vIGNsb3NlIHRoZSBwb2x5Z29uIGlmIGl0cyBlbmRwb2ludHMgYXJlIG5vdCB0aGUgc2FtZSBhZnRlciBjbGlwcGluZ1xuICAgICAgICBpZiAoY2xvc2VkICYmIHNsaWNlWzBdICE9PSBzbGljZVtzbGljZS5sZW5ndGggLSAxXSkgc2xpY2UucHVzaChzbGljZVswXSk7XG5cbiAgICAgICAgLy8gYWRkIHRoZSBmaW5hbCBzbGljZVxuICAgICAgICBuZXdTbGljZShzbGljZXMsIHNsaWNlLCBhcmVhLCBkaXN0KTtcbiAgICB9XG5cbiAgICByZXR1cm4gc2xpY2VzO1xufVxuXG5mdW5jdGlvbiBuZXdTbGljZShzbGljZXMsIHNsaWNlLCBhcmVhLCBkaXN0KSB7XG4gICAgaWYgKHNsaWNlLmxlbmd0aCkge1xuICAgICAgICAvLyB3ZSBkb24ndCByZWNhbGN1bGF0ZSB0aGUgYXJlYS9sZW5ndGggb2YgdGhlIHVuY2xpcHBlZCBnZW9tZXRyeSBiZWNhdXNlIHRoZSBjYXNlIHdoZXJlIGl0IGdvZXNcbiAgICAgICAgLy8gYmVsb3cgdGhlIHZpc2liaWxpdHkgdGhyZXNob2xkIGFzIGEgcmVzdWx0IG9mIGNsaXBwaW5nIGlzIHJhcmUsIHNvIHdlIGF2b2lkIGRvaW5nIHVubmVjZXNzYXJ5IHdvcmtcbiAgICAgICAgc2xpY2UuYXJlYSA9IGFyZWE7XG4gICAgICAgIHNsaWNlLmRpc3QgPSBkaXN0O1xuXG4gICAgICAgIHNsaWNlcy5wdXNoKHNsaWNlKTtcbiAgICB9XG4gICAgcmV0dXJuIFtdO1xufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGNvbnZlcnQ7XG5cbnZhciBzaW1wbGlmeSA9IHJlcXVpcmUoJy4vc2ltcGxpZnknKTtcblxuLy8gY29udmVydHMgR2VvSlNPTiBmZWF0dXJlIGludG8gYW4gaW50ZXJtZWRpYXRlIHByb2plY3RlZCBKU09OIHZlY3RvciBmb3JtYXQgd2l0aCBzaW1wbGlmaWNhdGlvbiBkYXRhXG5cbmZ1bmN0aW9uIGNvbnZlcnQoZGF0YSwgdG9sZXJhbmNlKSB7XG4gICAgdmFyIGZlYXR1cmVzID0gW107XG5cbiAgICBpZiAoZGF0YS50eXBlID09PSAnRmVhdHVyZUNvbGxlY3Rpb24nKSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZGF0YS5mZWF0dXJlcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY29udmVydEZlYXR1cmUoZmVhdHVyZXMsIGRhdGEuZmVhdHVyZXNbaV0sIHRvbGVyYW5jZSk7XG4gICAgICAgIH1cbiAgICB9IGVsc2UgaWYgKGRhdGEudHlwZSA9PT0gJ0ZlYXR1cmUnKSB7XG4gICAgICAgIGNvbnZlcnRGZWF0dXJlKGZlYXR1cmVzLCBkYXRhLCB0b2xlcmFuY2UpO1xuXG4gICAgfSBlbHNlIHtcbiAgICAgICAgLy8gc2luZ2xlIGdlb21ldHJ5IG9yIGEgZ2VvbWV0cnkgY29sbGVjdGlvblxuICAgICAgICBjb252ZXJ0RmVhdHVyZShmZWF0dXJlcywge2dlb21ldHJ5OiBkYXRhfSwgdG9sZXJhbmNlKTtcbiAgICB9XG4gICAgcmV0dXJuIGZlYXR1cmVzO1xufVxuXG5mdW5jdGlvbiBjb252ZXJ0RmVhdHVyZShmZWF0dXJlcywgZmVhdHVyZSwgdG9sZXJhbmNlKSB7XG4gICAgdmFyIGdlb20gPSBmZWF0dXJlLmdlb21ldHJ5LFxuICAgICAgICB0eXBlID0gZ2VvbS50eXBlLFxuICAgICAgICBjb29yZHMgPSBnZW9tLmNvb3JkaW5hdGVzLFxuICAgICAgICB0YWdzID0gZmVhdHVyZS5wcm9wZXJ0aWVzLFxuICAgICAgICBpLCBqLCByaW5ncztcblxuICAgIGlmICh0eXBlID09PSAnUG9pbnQnKSB7XG4gICAgICAgIGZlYXR1cmVzLnB1c2goY3JlYXRlKHRhZ3MsIDEsIFtwcm9qZWN0UG9pbnQoY29vcmRzKV0pKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ011bHRpUG9pbnQnKSB7XG4gICAgICAgIGZlYXR1cmVzLnB1c2goY3JlYXRlKHRhZ3MsIDEsIHByb2plY3QoY29vcmRzKSkpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnTGluZVN0cmluZycpIHtcbiAgICAgICAgZmVhdHVyZXMucHVzaChjcmVhdGUodGFncywgMiwgW3Byb2plY3QoY29vcmRzLCB0b2xlcmFuY2UpXSkpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnTXVsdGlMaW5lU3RyaW5nJyB8fCB0eXBlID09PSAnUG9seWdvbicpIHtcbiAgICAgICAgcmluZ3MgPSBbXTtcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IGNvb3Jkcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgcmluZ3MucHVzaChwcm9qZWN0KGNvb3Jkc1tpXSwgdG9sZXJhbmNlKSk7XG4gICAgICAgIH1cbiAgICAgICAgZmVhdHVyZXMucHVzaChjcmVhdGUodGFncywgdHlwZSA9PT0gJ1BvbHlnb24nID8gMyA6IDIsIHJpbmdzKSk7XG5cbiAgICB9IGVsc2UgaWYgKHR5cGUgPT09ICdNdWx0aVBvbHlnb24nKSB7XG4gICAgICAgIHJpbmdzID0gW107XG4gICAgICAgIGZvciAoaSA9IDA7IGkgPCBjb29yZHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGZvciAoaiA9IDA7IGogPCBjb29yZHNbaV0ubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICByaW5ncy5wdXNoKHByb2plY3QoY29vcmRzW2ldW2pdLCB0b2xlcmFuY2UpKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBmZWF0dXJlcy5wdXNoKGNyZWF0ZSh0YWdzLCAzLCByaW5ncykpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnR2VvbWV0cnlDb2xsZWN0aW9uJykge1xuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZ2VvbS5nZW9tZXRyaWVzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgICAgICBjb252ZXJ0RmVhdHVyZShmZWF0dXJlcywge1xuICAgICAgICAgICAgICAgIGdlb21ldHJ5OiBnZW9tLmdlb21ldHJpZXNbaV0sXG4gICAgICAgICAgICAgICAgcHJvcGVydGllczogdGFnc1xuICAgICAgICAgICAgfSwgdG9sZXJhbmNlKTtcbiAgICAgICAgfVxuXG4gICAgfSBlbHNlIHtcbiAgICAgICAgY29uc29sZS53YXJuKCdVbnN1cHBvcnRlZCBHZW9KU09OIHR5cGU6ICcgKyBnZW9tLnR5cGUpO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gY3JlYXRlKHRhZ3MsIHR5cGUsIGdlb21ldHJ5KSB7XG4gICAgdmFyIGZlYXR1cmUgPSB7XG4gICAgICAgIGdlb21ldHJ5OiBnZW9tZXRyeSxcbiAgICAgICAgdHlwZTogdHlwZSxcbiAgICAgICAgdGFnczogdGFncyB8fCBudWxsLFxuICAgICAgICBtaW46IFsxLCAxXSwgLy8gaW5pdGlhbCBiYm94IHZhbHVlcztcbiAgICAgICAgbWF4OiBbMCwgMF0gIC8vIG5vdGUgdGhhdCBhbGwgY29vcmRzIGFyZSBpbiBbMC4uMV0gcmFuZ2VcbiAgICB9O1xuICAgIGNhbGNCQm94KGZlYXR1cmUpO1xuICAgIHJldHVybiBmZWF0dXJlO1xufVxuXG5mdW5jdGlvbiBwcm9qZWN0KGxvbmxhdHMsIHRvbGVyYW5jZSkge1xuICAgIHZhciBwcm9qZWN0ZWQgPSBbXTtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGxvbmxhdHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgcHJvamVjdGVkLnB1c2gocHJvamVjdFBvaW50KGxvbmxhdHNbaV0pKTtcbiAgICB9XG4gICAgaWYgKHRvbGVyYW5jZSkge1xuICAgICAgICBzaW1wbGlmeShwcm9qZWN0ZWQsIHRvbGVyYW5jZSk7XG4gICAgICAgIGNhbGNTaXplKHByb2plY3RlZCk7XG4gICAgfVxuICAgIHJldHVybiBwcm9qZWN0ZWQ7XG59XG5cbmZ1bmN0aW9uIHByb2plY3RQb2ludChwKSB7XG4gICAgdmFyIHNpbiA9IE1hdGguc2luKHBbMV0gKiBNYXRoLlBJIC8gMTgwKSxcbiAgICAgICAgeCA9IChwWzBdIC8gMzYwICsgMC41KSxcbiAgICAgICAgeSA9ICgwLjUgLSAwLjI1ICogTWF0aC5sb2coKDEgKyBzaW4pIC8gKDEgLSBzaW4pKSAvIE1hdGguUEkpO1xuICAgIHJldHVybiBbeCwgeSwgMF07XG59XG5cbi8vIGNhbGN1bGF0ZSBhcmVhIGFuZCBsZW5ndGggb2YgdGhlIHBvbHlcbmZ1bmN0aW9uIGNhbGNTaXplKHBvaW50cykge1xuICAgIHZhciBhcmVhID0gMCxcbiAgICAgICAgZGlzdCA9IDA7XG5cbiAgICBmb3IgKHZhciBpID0gMCwgYSwgYjsgaSA8IHBvaW50cy5sZW5ndGggLSAxOyBpKyspIHtcbiAgICAgICAgYSA9IGIgfHwgcG9pbnRzW2ldO1xuICAgICAgICBiID0gcG9pbnRzW2kgKyAxXTtcblxuICAgICAgICBhcmVhICs9IGFbMF0gKiBiWzFdIC0gYlswXSAqIGFbMV07XG5cbiAgICAgICAgLy8gdXNlIE1hbmhhdHRhbiBkaXN0YW5jZSBpbnN0ZWFkIG9mIEV1Y2xpZGlhbiBvbmUgdG8gYXZvaWQgZXhwZW5zaXZlIHNxdWFyZSByb290IGNvbXB1dGF0aW9uXG4gICAgICAgIGRpc3QgKz0gTWF0aC5hYnMoYlswXSAtIGFbMF0pICsgTWF0aC5hYnMoYlsxXSAtIGFbMV0pO1xuICAgIH1cbiAgICBwb2ludHMuYXJlYSA9IE1hdGguYWJzKGFyZWEgLyAyKTtcbiAgICBwb2ludHMuZGlzdCA9IGRpc3Q7XG59XG5cbi8vIGNhbGN1bGF0ZSB0aGUgZmVhdHVyZSBib3VuZGluZyBib3ggZm9yIGZhc3RlciBjbGlwcGluZyBsYXRlclxuZnVuY3Rpb24gY2FsY0JCb3goZmVhdHVyZSkge1xuICAgIHZhciBnZW9tZXRyeSA9IGZlYXR1cmUuZ2VvbWV0cnksXG4gICAgICAgIG1pbiA9IGZlYXR1cmUubWluLFxuICAgICAgICBtYXggPSBmZWF0dXJlLm1heDtcblxuICAgIGlmIChmZWF0dXJlLnR5cGUgPT09IDEpIHtcbiAgICAgICAgY2FsY1JpbmdCQm94KG1pbiwgbWF4LCBnZW9tZXRyeSk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBnZW9tZXRyeS5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY2FsY1JpbmdCQm94KG1pbiwgbWF4LCBnZW9tZXRyeVtpXSk7XG4gICAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuIGZlYXR1cmU7XG59XG5cbmZ1bmN0aW9uIGNhbGNSaW5nQkJveChtaW4sIG1heCwgcG9pbnRzKSB7XG4gICAgZm9yICh2YXIgaSA9IDAsIHA7IGkgPCBwb2ludHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgcCA9IHBvaW50c1tpXTtcbiAgICAgICAgbWluWzBdID0gTWF0aC5taW4ocFswXSwgbWluWzBdKTtcbiAgICAgICAgbWF4WzBdID0gTWF0aC5tYXgocFswXSwgbWF4WzBdKTtcbiAgICAgICAgbWluWzFdID0gTWF0aC5taW4ocFsxXSwgbWluWzFdKTtcbiAgICAgICAgbWF4WzFdID0gTWF0aC5tYXgocFsxXSwgbWF4WzFdKTtcbiAgICB9XG59XG4iLCIndXNlIHN0cmljdCc7XG5cbm1vZHVsZS5leHBvcnRzID0gZ2VvanNvbnZ0O1xuXG52YXIgY29udmVydCA9IHJlcXVpcmUoJy4vY29udmVydCcpLCAvLyBHZW9KU09OIGNvbnZlcnNpb24gYW5kIHByZXByb2Nlc3NpbmdcbiAgICBjbGlwID0gcmVxdWlyZSgnLi9jbGlwJyksICAgICAgIC8vIHN0cmlwZSBjbGlwcGluZyBhbGdvcml0aG1cbiAgICBjcmVhdGVUaWxlID0gcmVxdWlyZSgnLi90aWxlJyk7IC8vIGZpbmFsIHNpbXBsaWZpZWQgdGlsZSBnZW5lcmF0aW9uXG5cblxuZnVuY3Rpb24gZ2VvanNvbnZ0KGRhdGEsIG9wdGlvbnMpIHtcbiAgICByZXR1cm4gbmV3IEdlb0pTT05WVChkYXRhLCBvcHRpb25zKTtcbn1cblxuZnVuY3Rpb24gR2VvSlNPTlZUKGRhdGEsIG9wdGlvbnMpIHtcbiAgICBvcHRpb25zID0gdGhpcy5vcHRpb25zID0gZXh0ZW5kKE9iamVjdC5jcmVhdGUodGhpcy5vcHRpb25zKSwgb3B0aW9ucyk7XG5cbiAgICB2YXIgZGVidWcgPSBvcHRpb25zLmRlYnVnO1xuXG4gICAgaWYgKGRlYnVnKSBjb25zb2xlLnRpbWUoJ3ByZXByb2Nlc3MgZGF0YScpO1xuXG4gICAgdmFyIHoyID0gMSA8PCBvcHRpb25zLmJhc2Vab29tLCAvLyAyXnpcbiAgICAgICAgZmVhdHVyZXMgPSBjb252ZXJ0KGRhdGEsIG9wdGlvbnMudG9sZXJhbmNlIC8gKHoyICogb3B0aW9ucy5leHRlbnQpKTtcblxuICAgIHRoaXMudGlsZXMgPSB7fTtcblxuICAgIGlmIChkZWJ1Zykge1xuICAgICAgICBjb25zb2xlLnRpbWVFbmQoJ3ByZXByb2Nlc3MgZGF0YScpO1xuICAgICAgICBjb25zb2xlLnRpbWUoJ2dlbmVyYXRlIHRpbGVzIHVwIHRvIHonICsgb3B0aW9ucy5tYXhab29tKTtcbiAgICAgICAgdGhpcy5zdGF0cyA9IHt9O1xuICAgICAgICB0aGlzLnRvdGFsID0gMDtcbiAgICB9XG5cbiAgICAvLyBzdGFydCBzbGljaW5nIGZyb20gdGhlIHRvcCB0aWxlIGRvd25cbiAgICB0aGlzLnNwbGl0VGlsZShmZWF0dXJlcywgMCwgMCwgMCk7XG5cbiAgICBpZiAoZGVidWcpIHtcbiAgICAgICAgY29uc29sZS5sb2coJ2ZlYXR1cmVzOiAlZCwgcG9pbnRzOiAlZCcsIHRoaXMudGlsZXNbMF0ubnVtRmVhdHVyZXMsIHRoaXMudGlsZXNbMF0ubnVtUG9pbnRzKTtcbiAgICAgICAgY29uc29sZS50aW1lRW5kKCdnZW5lcmF0ZSB0aWxlcyB1cCB0byB6JyArIG9wdGlvbnMubWF4Wm9vbSk7XG4gICAgICAgIGNvbnNvbGUubG9nKCd0aWxlcyBnZW5lcmF0ZWQ6JywgdGhpcy50b3RhbCwgSlNPTi5zdHJpbmdpZnkodGhpcy5zdGF0cykpO1xuICAgIH1cbn1cblxuR2VvSlNPTlZULnByb3RvdHlwZS5vcHRpb25zID0ge1xuICAgIGJhc2Vab29tOiAxNCwgICAvLyBtYXggem9vbSB0byBwcmVzZXJ2ZSBkZXRhaWwgb25cbiAgICBtYXhab29tOiA0LCAgICAgLy8gem9vbSB0byBzbGljZSBkb3duIHRvIG9uIGZpcnN0IHBhc3NcbiAgICBtYXhQb2ludHM6IDEwMCwgLy8gc3RvcCBzbGljaW5nIGEgdGlsZSBiZWxvdyB0aGlzIG51bWJlciBvZiBwb2ludHNcbiAgICB0b2xlcmFuY2U6IDMsICAgLy8gc2ltcGxpZmljYXRpb24gdG9sZXJhbmNlIChoaWdoZXIgbWVhbnMgc2ltcGxlcilcbiAgICBleHRlbnQ6IDQwOTYsICAgLy8gdGlsZSBleHRlbnRcbiAgICBidWZmZXI6IDY0LCAgICAgLy8gdGlsZSBidWZmZXIgb24gZWFjaCBzaWRlXG4gICAgZGVidWc6IDAgICAgICAgIC8vIGxvZ2dpbmcgbGV2ZWwgKDAsIDEgb3IgMilcbn07XG5cbkdlb0pTT05WVC5wcm90b3R5cGUuc3BsaXRUaWxlID0gZnVuY3Rpb24gKGZlYXR1cmVzLCB6LCB4LCB5LCBjeiwgY3gsIGN5KSB7XG5cbiAgICB2YXIgc3RhY2sgPSBbZmVhdHVyZXMsIHosIHgsIHldLFxuICAgICAgICBvcHRpb25zID0gdGhpcy5vcHRpb25zLFxuICAgICAgICBkZWJ1ZyA9IG9wdGlvbnMuZGVidWcsXG4gICAgICAgIGV4dGVudCA9IG9wdGlvbnMuZXh0ZW50LFxuICAgICAgICBidWZmZXIgPSBvcHRpb25zLmJ1ZmZlcjtcblxuICAgIC8vIGF2b2lkIHJlY3Vyc2lvbiBieSB1c2luZyBhIHByb2Nlc3NpbmcgcXVldWVcbiAgICB3aGlsZSAoc3RhY2subGVuZ3RoKSB7XG4gICAgICAgIGZlYXR1cmVzID0gc3RhY2suc2hpZnQoKTtcbiAgICAgICAgeiA9IHN0YWNrLnNoaWZ0KCk7XG4gICAgICAgIHggPSBzdGFjay5zaGlmdCgpO1xuICAgICAgICB5ID0gc3RhY2suc2hpZnQoKTtcblxuICAgICAgICB2YXIgejIgPSAxIDw8IHosXG4gICAgICAgICAgICBpZCA9IHRvSUQoeiwgeCwgeSksXG4gICAgICAgICAgICB0aWxlID0gdGhpcy50aWxlc1tpZF0sXG4gICAgICAgICAgICB0aWxlVG9sZXJhbmNlID0geiA9PT0gb3B0aW9ucy5iYXNlWm9vbSA/IDAgOiBvcHRpb25zLnRvbGVyYW5jZSAvICh6MiAqIGV4dGVudCk7XG5cbiAgICAgICAgaWYgKCF0aWxlKSB7XG4gICAgICAgICAgICBpZiAoZGVidWcgPiAxKSBjb25zb2xlLnRpbWUoJ2NyZWF0aW9uJyk7XG5cbiAgICAgICAgICAgIHRpbGUgPSB0aGlzLnRpbGVzW2lkXSA9IGNyZWF0ZVRpbGUoZmVhdHVyZXMsIHoyLCB4LCB5LCB0aWxlVG9sZXJhbmNlLCBleHRlbnQsIHogPT09IG9wdGlvbnMuYmFzZVpvb20pO1xuXG4gICAgICAgICAgICBpZiAoZGVidWcpIHtcbiAgICAgICAgICAgICAgICBpZiAoZGVidWcgPiAxKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnNvbGUubG9nKCd0aWxlIHolZC0lZC0lZCAoZmVhdHVyZXM6ICVkLCBwb2ludHM6ICVkLCBzaW1wbGlmaWVkOiAlZCknLFxuICAgICAgICAgICAgICAgICAgICAgICAgeiwgeCwgeSwgdGlsZS5udW1GZWF0dXJlcywgdGlsZS5udW1Qb2ludHMsIHRpbGUubnVtU2ltcGxpZmllZCk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnNvbGUudGltZUVuZCgnY3JlYXRpb24nKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgdmFyIGtleSA9ICd6JyArIHogKyAnOic7XG4gICAgICAgICAgICAgICAgdGhpcy5zdGF0c1trZXldID0gKHRoaXMuc3RhdHNba2V5XSB8fCAwKSArIDE7XG4gICAgICAgICAgICAgICAgdGhpcy50b3RhbCsrO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgaWYgKCFjeiAmJiAoeiA9PT0gb3B0aW9ucy5tYXhab29tIHx8IHRpbGUubnVtUG9pbnRzIDw9IG9wdGlvbnMubWF4UG9pbnRzIHx8XG4gICAgICAgICAgICAgICAgaXNDbGlwcGVkU3F1YXJlKHRpbGUuZmVhdHVyZXMsIGV4dGVudCwgYnVmZmVyKSkgfHwgeiA9PT0gb3B0aW9ucy5iYXNlWm9vbSB8fCB6ID09PSBjeikge1xuICAgICAgICAgICAgdGlsZS5zb3VyY2UgPSBmZWF0dXJlcztcbiAgICAgICAgICAgIGNvbnRpbnVlOyAvLyBzdG9wIHRpbGluZ1xuICAgICAgICB9XG5cbiAgICAgICAgaWYgKGN6KSB0aWxlLnNvdXJjZSA9IGZlYXR1cmVzO1xuICAgICAgICBlbHNlIHRpbGUuc291cmNlID0gbnVsbDtcblxuICAgICAgICBpZiAoZGVidWcgPiAxKSBjb25zb2xlLnRpbWUoJ2NsaXBwaW5nJyk7XG5cbiAgICAgICAgLy8gdmFsdWVzIHdlJ2xsIHVzZSBmb3IgY2xpcHBpbmdcbiAgICAgICAgdmFyIGsxID0gMC41ICogYnVmZmVyIC8gZXh0ZW50LFxuICAgICAgICAgICAgazIgPSAwLjUgLSBrMSxcbiAgICAgICAgICAgIGszID0gMC41ICsgazEsXG4gICAgICAgICAgICBrNCA9IDEgKyBrMSxcblxuICAgICAgICAgICAgdGwsIGJsLCB0ciwgYnIsIGxlZnQsIHJpZ2h0LFxuICAgICAgICAgICAgbSwgZ29MZWZ0LCBnb1RvcDtcblxuICAgICAgICBpZiAoY3opIHsgLy8gaWYgd2UgaGF2ZSBhIHNwZWNpZmljIHRpbGUgdG8gZHJpbGwgZG93biB0bywgY2FsY3VsYXRlIHdoZXJlIHRvIGdvXG4gICAgICAgICAgICBtID0gMSA8PCAoY3ogLSB6KTtcbiAgICAgICAgICAgIGdvTGVmdCA9IGN4IC8gbSAtIHggPCAwLjU7XG4gICAgICAgICAgICBnb1RvcCA9IGN5IC8gbSAtIHkgPCAwLjU7XG4gICAgICAgIH1cblxuICAgICAgICB0bCA9IGJsID0gdHIgPSBiciA9IGxlZnQgPSByaWdodCA9IG51bGw7XG5cbiAgICAgICAgaWYgKCFjeiB8fCAgZ29MZWZ0KSBsZWZ0ICA9IGNsaXAoZmVhdHVyZXMsIHoyLCB4IC0gazEsIHggKyBrMywgMCwgaW50ZXJzZWN0WCk7XG4gICAgICAgIGlmICghY3ogfHwgIWdvTGVmdCkgcmlnaHQgPSBjbGlwKGZlYXR1cmVzLCB6MiwgeCArIGsyLCB4ICsgazQsIDAsIGludGVyc2VjdFgpO1xuXG4gICAgICAgIGlmIChsZWZ0KSB7XG4gICAgICAgICAgICBpZiAoIWN6IHx8ICBnb1RvcCkgdGwgPSBjbGlwKGxlZnQsIHoyLCB5IC0gazEsIHkgKyBrMywgMSwgaW50ZXJzZWN0WSk7XG4gICAgICAgICAgICBpZiAoIWN6IHx8ICFnb1RvcCkgYmwgPSBjbGlwKGxlZnQsIHoyLCB5ICsgazIsIHkgKyBrNCwgMSwgaW50ZXJzZWN0WSk7XG4gICAgICAgIH1cblxuICAgICAgICBpZiAocmlnaHQpIHtcbiAgICAgICAgICAgIGlmICghY3ogfHwgIGdvVG9wKSB0ciA9IGNsaXAocmlnaHQsIHoyLCB5IC0gazEsIHkgKyBrMywgMSwgaW50ZXJzZWN0WSk7XG4gICAgICAgICAgICBpZiAoIWN6IHx8ICFnb1RvcCkgYnIgPSBjbGlwKHJpZ2h0LCB6MiwgeSArIGsyLCB5ICsgazQsIDEsIGludGVyc2VjdFkpO1xuICAgICAgICB9XG5cbiAgICAgICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS50aW1lRW5kKCdjbGlwcGluZycpO1xuXG4gICAgICAgIGlmICh0bCkgc3RhY2sucHVzaCh0bCwgeiArIDEsIHggKiAyLCAgICAgeSAqIDIpO1xuICAgICAgICBpZiAoYmwpIHN0YWNrLnB1c2goYmwsIHogKyAxLCB4ICogMiwgICAgIHkgKiAyICsgMSk7XG4gICAgICAgIGlmICh0cikgc3RhY2sucHVzaCh0ciwgeiArIDEsIHggKiAyICsgMSwgeSAqIDIpO1xuICAgICAgICBpZiAoYnIpIHN0YWNrLnB1c2goYnIsIHogKyAxLCB4ICogMiArIDEsIHkgKiAyICsgMSk7XG4gICAgfVxufTtcblxuR2VvSlNPTlZULnByb3RvdHlwZS5nZXRUaWxlID0gZnVuY3Rpb24gKHosIHgsIHkpIHtcbiAgICB2YXIgaWQgPSB0b0lEKHosIHgsIHkpO1xuICAgIGlmICh0aGlzLnRpbGVzW2lkXSkgcmV0dXJuIHRoaXMudGlsZXNbaWRdO1xuXG4gICAgdmFyIG9wdGlvbnMgPSB0aGlzLm9wdGlvbnMsXG4gICAgICAgIGRlYnVnID0gb3B0aW9ucy5kZWJ1ZztcblxuICAgIGlmIChkZWJ1ZyA+IDEpIGNvbnNvbGUubG9nKCdkcmlsbGluZyBkb3duIHRvIHolZC0lZC0lZCcsIHosIHgsIHkpO1xuXG4gICAgdmFyIHowID0geixcbiAgICAgICAgeDAgPSB4LFxuICAgICAgICB5MCA9IHksXG4gICAgICAgIHBhcmVudDtcblxuICAgIHdoaWxlICghcGFyZW50ICYmIHowID4gMCkge1xuICAgICAgICB6MC0tO1xuICAgICAgICB4MCA9IE1hdGguZmxvb3IoeDAgLyAyKTtcbiAgICAgICAgeTAgPSBNYXRoLmZsb29yKHkwIC8gMik7XG4gICAgICAgIHBhcmVudCA9IHRoaXMudGlsZXNbdG9JRCh6MCwgeDAsIHkwKV07XG4gICAgfVxuXG4gICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS5sb2coJ2ZvdW5kIHBhcmVudCB0aWxlIHolZC0lZC0lZCcsIHowLCB4MCwgeTApO1xuXG4gICAgLy8gaWYgd2UgZm91bmQgYSBwYXJlbnQgdGlsZSBjb250YWluaW5nIHRoZSBvcmlnaW5hbCBnZW9tZXRyeSwgd2UgY2FuIGRyaWxsIGRvd24gZnJvbSBpdFxuICAgIGlmIChwYXJlbnQuc291cmNlKSB7XG4gICAgICAgIGlmIChpc0NsaXBwZWRTcXVhcmUocGFyZW50LmZlYXR1cmVzLCBvcHRpb25zLmV4dGVudCwgb3B0aW9ucy5idWZmZXIpKSByZXR1cm4gcGFyZW50O1xuXG4gICAgICAgIGlmIChkZWJ1ZykgY29uc29sZS50aW1lKCdkcmlsbGluZyBkb3duJyk7XG4gICAgICAgIHRoaXMuc3BsaXRUaWxlKHBhcmVudC5zb3VyY2UsIHowLCB4MCwgeTAsIHosIHgsIHkpO1xuICAgICAgICBpZiAoZGVidWcpIGNvbnNvbGUudGltZUVuZCgnZHJpbGxpbmcgZG93bicpO1xuICAgIH1cblxuICAgIHJldHVybiB0aGlzLnRpbGVzW2lkXTtcbn07XG5cbi8vIGNoZWNrcyB3aGV0aGVyIGEgdGlsZSBpcyBhIHdob2xlLWFyZWEgZmlsbCBhZnRlciBjbGlwcGluZzsgaWYgaXQgaXMsIHRoZXJlJ3Mgbm8gc2Vuc2Ugc2xpY2luZyBpdCBmdXJ0aGVyXG5mdW5jdGlvbiBpc0NsaXBwZWRTcXVhcmUoZmVhdHVyZXMsIGV4dGVudCwgYnVmZmVyKSB7XG4gICAgaWYgKGZlYXR1cmVzLmxlbmd0aCAhPT0gMSkgcmV0dXJuIGZhbHNlO1xuXG4gICAgdmFyIGZlYXR1cmUgPSBmZWF0dXJlc1swXTtcbiAgICBpZiAoZmVhdHVyZS50eXBlICE9PSAzIHx8IGZlYXR1cmUuZ2VvbWV0cnkubGVuZ3RoID4gMSkgcmV0dXJuIGZhbHNlO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBmZWF0dXJlLmdlb21ldHJ5WzBdLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIHZhciBwID0gZmVhdHVyZS5nZW9tZXRyeVswXVtpXTtcbiAgICAgICAgaWYgKChwWzBdICE9PSAtYnVmZmVyICYmIHBbMF0gIT09IGV4dGVudCArIGJ1ZmZlcikgfHxcbiAgICAgICAgICAgIChwWzFdICE9PSAtYnVmZmVyICYmIHBbMV0gIT09IGV4dGVudCArIGJ1ZmZlcikpIHJldHVybiBmYWxzZTtcbiAgICB9XG4gICAgcmV0dXJuIHRydWU7XG59XG5cbmZ1bmN0aW9uIHRvSUQoeiwgeCwgeSkge1xuICAgIHJldHVybiAoKCgxIDw8IHopICogeSArIHgpICogMzIpICsgejtcbn1cblxuZnVuY3Rpb24gaW50ZXJzZWN0WChhLCBiLCB4KSB7XG4gICAgcmV0dXJuIFt4LCAoeCAtIGFbMF0pICogKGJbMV0gLSBhWzFdKSAvIChiWzBdIC0gYVswXSkgKyBhWzFdLCAxXTtcbn1cbmZ1bmN0aW9uIGludGVyc2VjdFkoYSwgYiwgeSkge1xuICAgIHJldHVybiBbKHkgLSBhWzFdKSAqIChiWzBdIC0gYVswXSkgLyAoYlsxXSAtIGFbMV0pICsgYVswXSwgeSwgMV07XG59XG5cbmZ1bmN0aW9uIGV4dGVuZChkZXN0LCBzcmMpIHtcbiAgICBmb3IgKHZhciBpIGluIHNyYykge1xuICAgICAgICBkZXN0W2ldID0gc3JjW2ldO1xuICAgIH1cbiAgICByZXR1cm4gZGVzdDtcbn1cbiIsIid1c2Ugc3RyaWN0JztcblxubW9kdWxlLmV4cG9ydHMgPSBzaW1wbGlmeTtcblxuLy8gY2FsY3VsYXRlIHNpbXBsaWZpY2F0aW9uIGRhdGEgdXNpbmcgb3B0aW1pemVkIERvdWdsYXMtUGV1Y2tlciBhbGdvcml0aG1cblxuZnVuY3Rpb24gc2ltcGxpZnkocG9pbnRzLCB0b2xlcmFuY2UpIHtcblxuICAgIHZhciBzcVRvbGVyYW5jZSA9IHRvbGVyYW5jZSAqIHRvbGVyYW5jZSxcbiAgICAgICAgbGVuID0gcG9pbnRzLmxlbmd0aCxcbiAgICAgICAgZmlyc3QgPSAwLFxuICAgICAgICBsYXN0ID0gbGVuIC0gMSxcbiAgICAgICAgc3RhY2sgPSBbXSxcbiAgICAgICAgaSwgbWF4U3FEaXN0LCBzcURpc3QsIGluZGV4O1xuXG4gICAgLy8gYWx3YXlzIHJldGFpbiB0aGUgZW5kcG9pbnRzICgxIGlzIHRoZSBtYXggdmFsdWUpXG4gICAgcG9pbnRzW2ZpcnN0XVsyXSA9IDE7XG4gICAgcG9pbnRzW2xhc3RdWzJdID0gMTtcblxuICAgIC8vIGF2b2lkIHJlY3Vyc2lvbiBieSB1c2luZyBhIHN0YWNrXG4gICAgd2hpbGUgKGxhc3QpIHtcblxuICAgICAgICBtYXhTcURpc3QgPSAwO1xuXG4gICAgICAgIGZvciAoaSA9IGZpcnN0ICsgMTsgaSA8IGxhc3Q7IGkrKykge1xuICAgICAgICAgICAgc3FEaXN0ID0gZ2V0U3FTZWdEaXN0KHBvaW50c1tpXSwgcG9pbnRzW2ZpcnN0XSwgcG9pbnRzW2xhc3RdKTtcblxuICAgICAgICAgICAgaWYgKHNxRGlzdCA+IG1heFNxRGlzdCkge1xuICAgICAgICAgICAgICAgIGluZGV4ID0gaTtcbiAgICAgICAgICAgICAgICBtYXhTcURpc3QgPSBzcURpc3Q7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBpZiAobWF4U3FEaXN0ID4gc3FUb2xlcmFuY2UpIHtcbiAgICAgICAgICAgIHBvaW50c1tpbmRleF1bMl0gPSBtYXhTcURpc3Q7IC8vIHNhdmUgdGhlIHBvaW50IGltcG9ydGFuY2UgaW4gc3F1YXJlZCBwaXhlbHMgYXMgYSB6IGNvb3JkaW5hdGVcbiAgICAgICAgICAgIHN0YWNrLnB1c2goZmlyc3QsIGluZGV4LCBpbmRleCwgbGFzdCk7XG4gICAgICAgIH1cblxuICAgICAgICBsYXN0ID0gc3RhY2sucG9wKCk7XG4gICAgICAgIGZpcnN0ID0gc3RhY2sucG9wKCk7XG4gICAgfVxufVxuXG4vLyBzcXVhcmUgZGlzdGFuY2UgZnJvbSBhIHBvaW50IHRvIGEgc2VnbWVudFxuZnVuY3Rpb24gZ2V0U3FTZWdEaXN0KHAsIGEsIGIpIHtcblxuICAgIHZhciB4ID0gYVswXSwgeSA9IGFbMV0sXG4gICAgICAgIGJ4ID0gYlswXSwgYnkgPSBiWzFdLFxuICAgICAgICBweCA9IHBbMF0sIHB5ID0gcFsxXSxcbiAgICAgICAgZHggPSBieCAtIHgsXG4gICAgICAgIGR5ID0gYnkgLSB5O1xuXG4gICAgaWYgKGR4ICE9PSAwIHx8IGR5ICE9PSAwKSB7XG5cbiAgICAgICAgdmFyIHQgPSAoKHB4IC0geCkgKiBkeCArIChweSAtIHkpICogZHkpIC8gKGR4ICogZHggKyBkeSAqIGR5KTtcblxuICAgICAgICBpZiAodCA+IDEpIHtcbiAgICAgICAgICAgIHggPSBieDtcbiAgICAgICAgICAgIHkgPSBieTtcblxuICAgICAgICB9IGVsc2UgaWYgKHQgPiAwKSB7XG4gICAgICAgICAgICB4ICs9IGR4ICogdDtcbiAgICAgICAgICAgIHkgKz0gZHkgKiB0O1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgZHggPSBweCAtIHg7XG4gICAgZHkgPSBweSAtIHk7XG5cbiAgICByZXR1cm4gZHggKiBkeCArIGR5ICogZHk7XG59XG4iLCIndXNlIHN0cmljdCc7XG5cbm1vZHVsZS5leHBvcnRzID0gY3JlYXRlVGlsZTtcblxuZnVuY3Rpb24gY3JlYXRlVGlsZShmZWF0dXJlcywgejIsIHR4LCB0eSwgdG9sZXJhbmNlLCBleHRlbnQsIG5vU2ltcGxpZnkpIHtcbiAgICB2YXIgdGlsZSA9IHtcbiAgICAgICAgZmVhdHVyZXM6IFtdLFxuICAgICAgICBudW1Qb2ludHM6IDAsXG4gICAgICAgIG51bVNpbXBsaWZpZWQ6IDAsXG4gICAgICAgIG51bUZlYXR1cmVzOiAwLFxuICAgICAgICBzb3VyY2U6IG51bGxcbiAgICB9O1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZmVhdHVyZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgdGlsZS5udW1GZWF0dXJlcysrO1xuICAgICAgICBhZGRGZWF0dXJlKHRpbGUsIGZlYXR1cmVzW2ldLCB6MiwgdHgsIHR5LCB0b2xlcmFuY2UsIGV4dGVudCwgbm9TaW1wbGlmeSk7XG4gICAgfVxuICAgIHJldHVybiB0aWxlO1xufVxuXG5mdW5jdGlvbiBhZGRGZWF0dXJlKHRpbGUsIGZlYXR1cmUsIHoyLCB0eCwgdHksIHRvbGVyYW5jZSwgZXh0ZW50LCBub1NpbXBsaWZ5KSB7XG5cbiAgICB2YXIgZ2VvbSA9IGZlYXR1cmUuZ2VvbWV0cnksXG4gICAgICAgIHR5cGUgPSBmZWF0dXJlLnR5cGUsXG4gICAgICAgIHRyYW5zZm9ybWVkID0gW10sXG4gICAgICAgIHNxVG9sZXJhbmNlID0gdG9sZXJhbmNlICogdG9sZXJhbmNlLFxuICAgICAgICBpLCBqLCByaW5nLCBwO1xuXG4gICAgaWYgKHR5cGUgPT09IDEpIHtcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IGdlb20ubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIHRyYW5zZm9ybWVkLnB1c2godHJhbnNmb3JtUG9pbnQoZ2VvbVtpXSwgejIsIHR4LCB0eSwgZXh0ZW50KSk7XG4gICAgICAgICAgICB0aWxlLm51bVBvaW50cysrO1xuICAgICAgICAgICAgdGlsZS5udW1TaW1wbGlmaWVkKys7XG4gICAgICAgIH1cblxuICAgIH0gZWxzZSB7XG5cbiAgICAgICAgLy8gc2ltcGxpZnkgYW5kIHRyYW5zZm9ybSBwcm9qZWN0ZWQgY29vcmRpbmF0ZXMgZm9yIHRpbGUgZ2VvbWV0cnlcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IGdlb20ubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIHJpbmcgPSBnZW9tW2ldO1xuXG4gICAgICAgICAgICAvLyBmaWx0ZXIgb3V0IHRpbnkgcG9seWxpbmVzICYgcG9seWdvbnNcbiAgICAgICAgICAgIGlmICghbm9TaW1wbGlmeSAmJiAoKHR5cGUgPT09IDIgJiYgcmluZy5kaXN0IDwgdG9sZXJhbmNlKSB8fFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAodHlwZSA9PT0gMyAmJiByaW5nLmFyZWEgPCBzcVRvbGVyYW5jZSkpKSB7XG4gICAgICAgICAgICAgICAgdGlsZS5udW1Qb2ludHMgKz0gcmluZy5sZW5ndGg7XG4gICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHZhciB0cmFuc2Zvcm1lZFJpbmcgPSBbXTtcblxuICAgICAgICAgICAgZm9yIChqID0gMDsgaiA8IHJpbmcubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICBwID0gcmluZ1tqXTtcbiAgICAgICAgICAgICAgICAvLyBrZWVwIHBvaW50cyB3aXRoIGltcG9ydGFuY2UgPiB0b2xlcmFuY2VcbiAgICAgICAgICAgICAgICBpZiAobm9TaW1wbGlmeSB8fCBwWzJdID4gc3FUb2xlcmFuY2UpIHtcbiAgICAgICAgICAgICAgICAgICAgdHJhbnNmb3JtZWRSaW5nLnB1c2godHJhbnNmb3JtUG9pbnQocCwgejIsIHR4LCB0eSwgZXh0ZW50KSk7XG4gICAgICAgICAgICAgICAgICAgIHRpbGUubnVtU2ltcGxpZmllZCsrO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB0aWxlLm51bVBvaW50cysrO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICB0cmFuc2Zvcm1lZC5wdXNoKHRyYW5zZm9ybWVkUmluZyk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBpZiAodHJhbnNmb3JtZWQubGVuZ3RoKSB7XG4gICAgICAgIHRpbGUuZmVhdHVyZXMucHVzaCh7XG4gICAgICAgICAgICBnZW9tZXRyeTogdHJhbnNmb3JtZWQsXG4gICAgICAgICAgICB0eXBlOiB0eXBlLFxuICAgICAgICAgICAgdGFnczogZmVhdHVyZS50YWdzIHx8IG51bGxcbiAgICAgICAgfSk7XG4gICAgfVxufVxuXG5mdW5jdGlvbiB0cmFuc2Zvcm1Qb2ludChwLCB6MiwgdHgsIHR5LCBleHRlbnQpIHtcbiAgICB2YXIgeCA9IE1hdGgucm91bmQoZXh0ZW50ICogKHBbMF0gKiB6MiAtIHR4KSksXG4gICAgICAgIHkgPSBNYXRoLnJvdW5kKGV4dGVudCAqIChwWzFdICogejIgLSB0eSkpO1xuICAgIHJldHVybiBbeCwgeV07XG59XG4iXX0=
