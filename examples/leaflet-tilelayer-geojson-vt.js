// Load data tiles from an AJAX data source
L.TileLayer.Ajax = L.TileLayer.extend({

    initialize: function (url, options) {
        L.Util.extend(this.options, options);
        if (typeof (url) === 'object') {
            this.tileIndex = geojsonvt(url, {
                baseZoom: 18, // max zoom to preserve detail on
                maxZoom: 4, // zoom to slice down to on first pass
                maxPoints: 100, // stop slicing each tile below this number of points
                tolerance: 3, // simplification tolerance (higher means simpler)
                extent: 4096, // tile extent (both width and height)
                buffer: 0, // tile buffer on each side
                debug: 0 // logging level (0 to disable, 1 or 2)
            });
            L.TileLayer.prototype.initialize.call(this, null, options);

        } else {
            L.TileLayer.prototype.initialize.call(this, url, options);
        }
    },
    _requests: [],
    _addTile: function (tilePoint) {
        var tile = {
            datum: null,
            processed: false
        };
        this._tiles[tilePoint.x + ':' + tilePoint.y] = tile;
        this._loadTile(tile, tilePoint);
    },
    // XMLHttpRequest handler; closure over the XHR object, the layer, and the tile
    _xhrHandler: function (req, layer, tile, tilePoint) {
        return function () {
            if (req.readyState !== 4) {
                return;
            }
            var s = req.status;
            if ((s >= 200 && s < 300) || s === 304) {
                tile.datum = JSON.parse(req.responseText);
                layer._tileLoaded(tile, tilePoint);
            } else {
                layer._tileLoaded(tile, tilePoint);
            }
        };
    },
    sinh: function sinh(x) {
        return (Math.exp(x) - Math.exp(-x)) / 2;
    },
    transformPoint: function (p, z2, tx, ty, extent) {
        var x = (p[0] / extent + tx), //Math.round(extent * (p[0] * z2 - tx)),
            y = (p[1] / extent + ty); //Math.round(extent * (p[1] * z2 - ty));
        var lon = 360 * x / z2 - 180;
        var lat = 180 / Math.PI * Math.atan(this.sinh(Math.PI * (1 - 2 * y / z2)));
        return [lon, lat];

        /*y = (0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI);
        var sin = Math.sin(p[1] * Math.PI / 180),
        x = (p[0] / 360 + 0.5),
        

        return [x, y];*/
    },
    addFeature: function (feature, tilePoint) {

        var geom = feature.geometry,
            type = feature.type,
            geojson = {
                coordinates: []
            },
            i, j, ring, p,
            z2 = 1 << tilePoint.z,
            tx = tilePoint.x,
            ty = tilePoint.y,
            extent = 4096;

        if (type === 1) {
            geojson.type = "Point";
            for (i = 0; i < geom.length; i++) {
                geojson.coordinates.push(this.transformPoint(geom[i], z2, tx, ty, extent));

            }

        } else {
            if (type === 2) {
                geojson.type = "LineString";
            } else if (type === 3) {
                geojson.type = "Polygon"
            }
            // simplify and transform projected coordinates for tile geometry
            for (i = 0; i < geom.length; i++) {
                ring = geom[i];

                var transformedRing = [];

                for (j = 0; j < ring.length; j++) {
                    p = ring[j];
                    transformedRing.push(this.transformPoint(p, z2, tx, ty, extent));
                }

                geojson.coordinates.push(transformedRing);
            }
        }


        return geojson;
    },
    // Load the requested tile via AJAX
    _loadTile: function (tile, tilePoint) {
        this._adjustTilePoint(tilePoint);
        if (this.tileIndex) {
            var data = this.tileIndex.getTile(tilePoint.z, tilePoint.x, tilePoint.y);
            if (data) {
                var geojson = [];
                for (var i = 0; i < data.features.length; i++) {
                    var feature = data.features[i];
                    geojson.push({
                        type: "Feature",
                        geometry: this.addFeature(feature, tilePoint),
                        properties: feature.tags
                    });
                }


                tile.datum = geojson;
                this._tileLoaded(tile, tilePoint);
            }
        } else {
            var layer = this;
            var req = new XMLHttpRequest();
            this._requests.push(req);
            req.onreadystatechange = this._xhrHandler(req, layer, tile, tilePoint);
            req.open('GET', this.getTileUrl(tilePoint), true);
            req.send();
        }
    },
    _reset: function () {
        L.TileLayer.prototype._reset.apply(this, arguments);
        for (var i in this._requests) {
            this._requests[i].abort();
        }
        this._requests = [];
    },
    _update: function () {
        if (this._map && this._map._panTransition && this._map._panTransition._inProgress) {
            return;
        }
        if (this._tilesToLoad < 0) {
            this._tilesToLoad = 0;
        }
        L.TileLayer.prototype._update.apply(this, arguments);
    }
});


L.TileLayer.GeoJSON = L.TileLayer.Ajax.extend({
    // Store each GeometryCollection's layer by key, if options.unique function is present
    _keyLayers: {},

    // Used to calculate svg path string for clip path elements
    _clipPathRectangles: {},

    initialize: function (url, options, geojsonOptions) {
        L.TileLayer.Ajax.prototype.initialize.call(this, url, options);
        this.geojsonLayer = new L.GeoJSON(null, geojsonOptions);
    },
    onAdd: function (map) {
        this._map = map;
        L.TileLayer.Ajax.prototype.onAdd.call(this, map);
        map.addLayer(this.geojsonLayer);
    },
    onRemove: function (map) {
        map.removeLayer(this.geojsonLayer);
        L.TileLayer.Ajax.prototype.onRemove.call(this, map);
    },
    _reset: function () {
        this.geojsonLayer.clearLayers();
        this._keyLayers = {};
        this._removeOldClipPaths();
        L.TileLayer.Ajax.prototype._reset.apply(this, arguments);
    },

    // Remove clip path elements from other earlier zoom levels
    _removeOldClipPaths: function () {
        for (var clipPathId in this._clipPathRectangles) {
            var clipPathZXY = clipPathId.split('_').slice(1);
            var zoom = parseInt(clipPathZXY[0], 10);
            if (zoom !== this._map.getZoom()) {
                var rectangle = this._clipPathRectangles[clipPathId];
                this._map.removeLayer(rectangle);
                var clipPath = document.getElementById(clipPathId);
                if (clipPath !== null) {
                    clipPath.parentNode.removeChild(clipPath);
                }
                delete this._clipPathRectangles[clipPathId];
            }
        }
    },

    // Recurse LayerGroups and call func() on L.Path layer instances
    _recurseLayerUntilPath: function (func, layer) {
        if (layer instanceof L.Path) {
            func(layer);
        } else if (layer instanceof L.LayerGroup) {
            // Recurse each child layer
            layer.getLayers().forEach(this._recurseLayerUntilPath.bind(this, func), this);
        }
    },

    _clipLayerToTileBoundary: function (layer, tilePoint) {
        // Only perform SVG clipping if the browser is using SVG
        if (!L.Path.SVG) {
            return;
        }
        if (!this._map) {
            return;
        }

        if (!this._map._pathRoot) {
            this._map._pathRoot = L.Path.prototype._createElement('svg');
            this._map._panes.overlayPane.appendChild(this._map._pathRoot);
        }
        var svg = this._map._pathRoot;

        // create the defs container if it doesn't exist
        var defs = null;
        if (svg.getElementsByTagName('defs').length === 0) {
            defs = document.createElementNS(L.Path.SVG_NS, 'defs');
            svg.insertBefore(defs, svg.firstChild);
        } else {
            defs = svg.getElementsByTagName('defs')[0];
        }

        // Create the clipPath for the tile if it doesn't exist
        var clipPathId = 'tileClipPath_' + tilePoint.z + '_' + tilePoint.x + '_' + tilePoint.y;
        var clipPath = document.getElementById(clipPathId);
        if (clipPath === null) {
            clipPath = document.createElementNS(L.Path.SVG_NS, 'clipPath');
            clipPath.id = clipPathId;

            // Create a hidden L.Rectangle to represent the tile's area
            var tileSize = this.options.tileSize,
                nwPoint = tilePoint.multiplyBy(tileSize),
                sePoint = nwPoint.add([tileSize, tileSize]),
                nw = this._map.unproject(nwPoint),
                se = this._map.unproject(sePoint);
            this._clipPathRectangles[clipPathId] = new L.Rectangle(new L.LatLngBounds([nw, se]), {
                opacity: 0,
                fillOpacity: 0,
                clickable: false,
                noClip: true
            });
            this._map.addLayer(this._clipPathRectangles[clipPathId]);

            // Add a clip path element to the SVG defs element
            // With a path element that has the hidden rectangle's SVG path string  
            var path = document.createElementNS(L.Path.SVG_NS, 'path');
            var pathString = this._clipPathRectangles[clipPathId].getPathString();
            path.setAttribute('d', pathString);
            clipPath.appendChild(path);
            defs.appendChild(clipPath);
        }

        // Add the clip-path attribute to reference the id of the tile clipPath
        this._recurseLayerUntilPath(function (pathLayer) {
            pathLayer._container.setAttribute('clip-path', 'url(#' + clipPathId + ')');
        }, layer);
    },

    // Add a geojson object from a tile to the GeoJSON layer
    // * If the options.unique function is specified, merge geometries into GeometryCollections
    // grouped by the key returned by options.unique(feature) for each GeoJSON feature
    // * If options.clipTiles is set, and the browser is using SVG, perform SVG clipping on each
    // tile's GeometryCollection 
    addTileData: function (geojson, tilePoint) {
        var features = L.Util.isArray(geojson) ? geojson : geojson.features,
            i, len, feature;

        if (features) {
            for (i = 0, len = features.length; i < len; i++) {
                // Only add this if geometry or geometries are set and not null
                feature = features[i];
                if (feature.geometries || feature.geometry || feature.features || feature.coordinates) {
                    this.addTileData(features[i], tilePoint);
                }
            }
            return this;
        }

        var options = this.geojsonLayer.options;

        if (options.filter && !options.filter(geojson)) {
            return;
        }

        var parentLayer = this.geojsonLayer;
        var incomingLayer = null;
        if (this.options.unique && typeof (this.options.unique) === 'function') {
            var key = this.options.unique(geojson);

            // When creating the layer for a unique key,
            // Force the geojson to be a geometry collection
            if (!(key in this._keyLayers && geojson.geometry.type !== 'GeometryCollection')) {
                geojson.geometry = {
                    type: 'GeometryCollection',
                    geometries: [geojson.geometry]
                };
            }

            // Transform the geojson into a new Layer
            try {
                incomingLayer = L.GeoJSON.geometryToLayer(geojson, options.pointToLayer, options.coordsToLatLng);
            }
            // Ignore GeoJSON objects that could not be parsed
            catch (e) {
                return this;
            }

            incomingLayer.feature = L.GeoJSON.asFeature(geojson);
            // Add the incoming Layer to existing key's GeometryCollection
            if (key in this._keyLayers) {
                parentLayer = this._keyLayers[key];
                parentLayer.feature.geometry.geometries.push(geojson.geometry);
            }
            // Convert the incoming GeoJSON feature into a new GeometryCollection layer
            else {
                this._keyLayers[key] = incomingLayer;
            }
        }
        // Add the incoming geojson feature to the L.GeoJSON Layer
        else {
            // Transform the geojson into a new layer
            try {
                incomingLayer = L.GeoJSON.geometryToLayer(geojson, options.pointToLayer, options.coordsToLatLng);
            }
            // Ignore GeoJSON objects that could not be parsed
            catch (e) {
                return this;
            }
            incomingLayer.feature = L.GeoJSON.asFeature(geojson);
        }
        incomingLayer.defaultOptions = incomingLayer.options;

        this.geojsonLayer.resetStyle(incomingLayer);

        if (options.onEachFeature) {
            options.onEachFeature(geojson, incomingLayer);
        }
        parentLayer.addLayer(incomingLayer);

        // If options.clipTiles is set and the browser is using SVG
        // then clip the layer using SVG clipping
        if (this.options.clipTiles) {
            this._clipLayerToTileBoundary(incomingLayer, tilePoint);
        }
        return this;
    },

    _tileLoaded: function (tile, tilePoint) {
        L.TileLayer.Ajax.prototype._tileLoaded.apply(this, arguments);
        if (tile.datum === null) {
            return null;
        }
        this.addTileData(tile.datum, tilePoint);
    }
});