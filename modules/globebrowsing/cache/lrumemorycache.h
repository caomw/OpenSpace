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

#ifndef __OPENSPACE_MODULE_GLOBEBROWSING___LRU_MEMORY_CACHE___H__
#define __OPENSPACE_MODULE_GLOBEBROWSING___LRU_MEMORY_CACHE___H__

#include <list>
#include <map>

namespace openspace {
namespace globebrowsing {
namespace cache {

// Templated class implementing a Least-Recently-Used Cache
template<typename KeyType, typename ValueType>
class LRUMemoryCache {
public:
    LRUMemoryCache(long maximumSize);

    void put(const KeyType& key, const ValueType& value);
    void clear();
    bool exist(const KeyType& key) const;
    ValueType get(const KeyType& key);
    long size() const;
    long maximumSize() const;

private:
    void clean();

    std::list<std::pair<KeyType, ValueType> > _itemList;
    std::map<KeyType, decltype(_itemList.begin())> _itemMap;
    long _cacheSize;
    long _maximumCacheSize;
};

} // namespace cache
} // namespace globebrowsing
} // namespace openspace

#include <modules/globebrowsing/cache/lrumemorycache.inl>

#endif // __OPENSPACE_MODULE_GLOBEBROWSING___LRU_MEMORY_CACHE___H__
