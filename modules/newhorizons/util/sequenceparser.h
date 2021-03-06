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

#ifndef __OPENSPACE_MODULE_NEWHORIZONS___SEQUENCEPARSER___H__
#define __OPENSPACE_MODULE_NEWHORIZONS___SEQUENCEPARSER___H__

#include <openspace/network/networkengine.h>
#include <modules/newhorizons/util/decoder.h>
#include <openspace/util/timerange.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace openspace {

struct Image {
    TimeRange timeRange;
    std::string path;
    std::vector<std::string> activeInstruments;
    std::string target;
    bool isPlaceholder = false;
    bool projected = false;
};

struct ImageSubset {
    TimeRange _range;
    std::vector<Image> _subset;
};

class SequenceParser {
public:
    virtual ~SequenceParser() = default;
    virtual bool create() = 0;
    virtual std::map<std::string, ImageSubset> getSubsetMap() final;
    virtual std::vector<std::pair<std::string, TimeRange>> getInstrumentTimes() final;
    virtual std::vector<std::pair<double, std::string>> getTargetTimes() final;
    std::map<std::string, std::unique_ptr<Decoder>>& getTranslation();
    virtual std::vector<double> getCaptureProgression() final;

protected:
    void sendPlaybookInformation(const std::string& name);

    std::map<std::string, ImageSubset> _subsetMap;
    std::vector<std::pair<std::string, TimeRange>> _instrumentTimes;
    std::vector<std::pair<double, std::string>> _targetTimes;
    std::vector<double> _captureProgression;

    std::map<std::string, std::unique_ptr<Decoder>> _fileTranslation;

    NetworkEngine::MessageIdentifier _messageIdentifier;
};

} // namespace openspace

#endif // __OPENSPACE_MODULE_NEWHORIZONS___SEQUENCEPARSER___H__
