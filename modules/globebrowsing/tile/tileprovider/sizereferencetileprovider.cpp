/*****************************************************************************************
 *                                                                                       *
 * OpenSpace                                                                             *
 *                                                                                       *
 * Copyright (c) 2014-2017                                                               *
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF  *
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE  *
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                         *
 ****************************************************************************************/

#include <modules/globebrowsing/tile/tileprovider/sizereferencetileprovider.h>

#include <modules/globebrowsing/geometry/geodeticpatch.h>

#include <openspace/engine/openspaceengine.h>

#include <ghoul/filesystem/filesystem.h>
#include <ghoul/font/fontmanager.h>
#include <ghoul/font/fontrenderer.h>
#include <ghoul/io/texture/texturereader.h>
#include <ghoul/misc/dictionary.h>

using namespace ghoul::fontrendering;

namespace openspace {
namespace globebrowsing {
namespace tileprovider {
    
namespace {
    const char* KeyRadii = "Radii";
    const char* KeyBackgroundImagePath = "BackgroundImagePath";
}

SizeReferenceTileProvider::SizeReferenceTileProvider(const ghoul::Dictionary& dictionary)
{
    _fontSize = 50;
    _font = OsEng.fontManager().font("Mono", _fontSize);

    glm::dvec3 radii(1,1,1);
    if (!dictionary.getValue(KeyRadii, radii)) {
        throw std::runtime_error("Must define key '" + std::string(KeyRadii) + "'");
    }
    _ellipsoid = Ellipsoid(radii);

    _backgroundTile.status = Tile::Status::Unavailable;
    std::string backgroundImagePath;
    if (dictionary.getValue(KeyBackgroundImagePath, backgroundImagePath)) {
        using namespace ghoul::io;
        std::string imgAbsPath = absPath(backgroundImagePath);
        _backgroundTile.texture = TextureReader::ref().loadTexture(imgAbsPath);
        _backgroundTile.texture->uploadTexture();
        _backgroundTile.texture->setFilter(ghoul::opengl::Texture::FilterMode::Linear);
        _backgroundTile.status = Tile::Status::OK;
    }
}

void SizeReferenceTileProvider::renderText(const ghoul::fontrendering::FontRenderer&
                                           fontRenderer,
                                           const TileIndex& tileIndex) const
{
    GeodeticPatch patch(tileIndex);
    bool aboveEquator = patch.isNorthern();
        
    double tileLongitudalLength = roundedLongitudalLength(tileIndex);

    std::string unit = "m";
    if (tileLongitudalLength > 9999) {
        tileLongitudalLength *= 0.001;
        unit = "km";
    }

    glm::vec2 textPosition;
    textPosition.x = 0;
    textPosition.y = aboveEquator ? _fontSize / 2 : _textureSize.y - 3 * _fontSize / 2;
    glm::vec4 color(1.0, 1.0, 1.0, 1.0);

    fontRenderer.render(
        *_font,
        textPosition,
        color,
        " %.0f %s",
        tileLongitudalLength, unit.c_str()
        );
}

int SizeReferenceTileProvider::roundedLongitudalLength(const TileIndex& tileIndex) const {
    GeodeticPatch patch(tileIndex);
    bool aboveEquator = patch.isNorthern();
    double lat = aboveEquator ? patch.minLat() : patch.maxLat();
    double lon1 = patch.minLon();
    double lon2 = patch.maxLon();
    int l = static_cast<int>(_ellipsoid.longitudalDistance(lat, lon1, lon2));

    bool useKm = l > 9999;
    if (useKm) {
        l /= 1000;
    }
    l = std::round(l);
    if (useKm) {
        l *= 1000;
    }

    return l;
}

TileIndex::TileHashKey SizeReferenceTileProvider::toHash(const TileIndex& tileIndex) const {
    int l = roundedLongitudalLength(tileIndex);
    TileIndex::TileHashKey key = static_cast<TileIndex::TileHashKey>(l);
    return key;
}

Tile SizeReferenceTileProvider::backgroundTile(const TileIndex& tileIndex) const {
    if (_backgroundTile.status == Tile::Status::OK) {
        Tile tile;
        auto t = _backgroundTile.texture;
        void* pixelData = new char[t->expectedPixelDataSize()];
        memcpy(pixelData, t->pixelData(), t->expectedPixelDataSize());
        tile.texture = std::make_shared<ghoul::opengl::Texture>(
            pixelData,
            t->dimensions(),
            t->format(),
            t->internalFormat(),
            t->dataType(),
            t->filter(),
            t->wrapping()
        );
        tile.texture->uploadTexture();
        tile.texture->setDataOwnership(ghoul::opengl::Texture::TakeOwnership::Yes);
        tile.status = Tile::Status::OK;
        return tile;
    }
    else {
        // use default background
        return TextTileProvider::backgroundTile(tileIndex);
    }
}

} // namespace tileprovider
} // namespace globebrowsing
} // namespace openspace
