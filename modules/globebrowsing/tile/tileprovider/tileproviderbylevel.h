/*****************************************************************************************
*                                                                                       *
* OpenSpace                                                                             *
*                                                                                       *
* Copyright (c) 2014-2016                                                               *
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

#ifndef __TILE_PROVIDER_BY_LEVEL__
#define __TILE_PROVIDER_BY_LEVEL__

#include <modules/globebrowsing/geometry/geodetic2.h>
#include <modules/globebrowsing/tile/tileprovider/tileprovider.h>
#include <modules/globebrowsing/tile/asynctilereader.h>

#include <gdal_priv.h>

#include <openspace/engine/downloadmanager.h>

#include <ghoul/logging/logmanager.h>
#include <ghoul/filesystem/filesystem.h> // absPath
#include <ghoul/opengl/texture.h>
#include <ghoul/io/texture/texturereader.h>
#include <ghoul/font/fontrenderer.h>

#include <set>

namespace openspace {
namespace globebrowsing {

    using namespace ghoul::opengl;    

    class TileProviderByLevel : public TileProvider {
    public:
        
        TileProviderByLevel(const ghoul::Dictionary& dictionary);
        TileProviderByLevel(const std::string& imagePath);
        virtual ~TileProviderByLevel() { }

        virtual Tile getTile(const TileIndex& tileIndex);
        virtual Tile getDefaultTile();
        virtual Tile::Status getTileStatus(const TileIndex& index);
        virtual TileDepthTransform depthTransform();
        virtual void update();
        virtual void reset();
        virtual int maxLevel();
    private:
        inline int providerIndex(int level) const;
        inline TileProvider* levelProvider(int level) const;

        std::vector<int> _providerIndices;
        std::vector<std::shared_ptr<TileProvider>> _levelTileProviders;
    };

} // namespace globebrowsing
} // namespace openspace

#endif  // __TILE_PROVIDER_BY_LEVEL__