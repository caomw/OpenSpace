/*****************************************************************************************
 *                                                                                       *
 * OpenSpace                                                                             *
 *                                                                                       *
 * Copyright (c) 2015                                                                    *
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

#ifndef __TFBRICKSELECTOR_H__
#define __TFBRICKSELECTOR_H__

#include <vector>
#include <modules/multiresvolume/rendering/brickselection.h>
#include <modules/multiresvolume/rendering/brickselector.h>
#include <modules/multiresvolume/rendering/brickcover.h>


namespace openspace {

class TSP;
class ErrorHistogramManager;
class TransferFunction;

class TfBrickSelector : public BrickSelector {
public:
    TfBrickSelector(TSP* tsp, ErrorHistogramManager* hm, TransferFunction* tf, int memoryBudget, int streamingBudget);

    ~TfBrickSelector();

    virtual bool initialize();

    void selectBricks(int timestep, std::vector<int>& bricks);
    void setMemoryBudget(int memoryBudget);
    void setStreamingBudget(int streamingBudget);
    bool calculateBrickErrors();
 private:

    TSP* _tsp;
    ErrorHistogramManager* _histogramManager;
    TransferFunction* _transferFunction;
    std::vector<float> _brickErrors;
    float spatialSplitPoints(unsigned int brickIndex);
    float temporalSplitPoints(unsigned int brickIndex);
    float splitPoints(unsigned int brickIndex, BrickSelection::SplitType& splitType);

    int linearCoords(int x, int y, int z);
    void writeSelection(BrickSelection coveredBricks, std::vector<int>& bricks);

    int _memoryBudget;
    int _streamingBudget;
    
};

} // namespace openspace

#endif // __TFBRICKSELECTOR_H__
