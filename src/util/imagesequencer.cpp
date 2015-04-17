/*****************************************************************************************
*                                                                                       *
* OpenSpace                                                                             *
*                                                                                       *
* Copyright (c) 2014                                                                    *
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

// open space includes
#include <openspace/util/imagesequencer.h>
#include <ghoul/logging/logmanager.h>
#include <ghoul/filesystem/filesystem.h>
#include <ghoul/filesystem/directory.h>
#include <openspace/util/time.h>
#include <ghoul/filesystem/cachemanager.h>

#include <openspace/util/spicemanager.h>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <limits>

namespace {
const std::string _loggerCat = "ImageSequencer";
}

namespace openspace {

ImageSequencer* ImageSequencer::_sequencer = nullptr;

struct ImageParams{
	double startTime;
	double stopTime;
	std::string path;
	std::string activeInstrument;
	bool projected;
};

auto cmp = [](const ImageParams &a, const ImageParams &b)->bool{
	return a.startTime < b.startTime;
};

std::vector<ImageParams> _timeStamps;

ImageSequencer::ImageSequencer() 
    : _nextCapture(0.0)
    , _defaultCaptureImage(absPath("${OPENSPACE_DATA}/scene/common/textures/placeholder.png"))
{}


ImageSequencer& ImageSequencer::ref() {
	assert(_sequencer != nullptr);
	return *_sequencer;
}
void ImageSequencer::initialize() {
	assert(_sequencer == nullptr);
	_sequencer = new ImageSequencer;
}

void ImageSequencer::deinitialize() {
	delete _sequencer;
	_sequencer = nullptr;
}

void ImageSequencer::createImage(double t1, double t2, std::string instrument, std::string path) {
	// insert
	ImageParams image;
	image.startTime = t1;
	image.stopTime  = t2;
	image.path      = path;
	image.activeInstrument = instrument;
	image.projected = false;
	
	_timeStamps.push_back(image);
	// sort
}

double ImageSequencer::getNextCaptureTime(){
	return _nextCapture;
}
double ImageSequencer::nextCaptureTime(double _time){
	auto binary_find = [](std::vector<ImageParams>::iterator begin,
		std::vector<ImageParams>::iterator end,
		const ImageParams &val,
		bool(*cmp)(const ImageParams &a, const ImageParams &b))->std::vector<ImageParams>::iterator{
		// Finds the lower bound in at most log(last - first) + 1 comparisons
		std::vector<ImageParams>::iterator it = std::lower_bound(begin, end, val, cmp);
		if (it != begin){
			return it;
		}
		return end;
	};

	auto it = binary_find(_timeStamps.begin(), _timeStamps.end(), { _time, 0, "", "",  false }, cmp);

	if (it == _timeStamps.end() || _time < _nextCapture) return _nextCapture;

	return it->startTime;
}

std::string ImageSequencer::findActiveInstrument(double time){
	auto binary_find = [](std::vector<ImageParams>::iterator begin,
		std::vector<ImageParams>::iterator end,
		const ImageParams &val,
		bool(*cmp)(const ImageParams &a, const ImageParams &b))->std::vector<ImageParams>::iterator{
		// Finds the lower bound in at most log(last - first) + 1 comparisons
		std::vector<ImageParams>::iterator it = std::lower_bound(begin, end, val, cmp);
		if (it != begin){
			return it;
		}
		return end;
	};

	auto it = binary_find(_timeStamps.begin(), _timeStamps.end(), { time, 0, "", "", false }, cmp);
	
	if ((it == _timeStamps.end())){
		_activeInstrument = "Not found, incufficient playbook-data";
	}else{
		_activeInstrument = std::prev(it)->activeInstrument;
	}
	

	return _activeInstrument;
}

bool ImageSequencer::getImagePath(double& currentTime, std::string& path,  bool closedInterval){
	auto binary_find = [](std::vector<ImageParams>::iterator begin,
						std::vector<ImageParams>::iterator end,
						const ImageParams &val,
						bool(*cmp)(const ImageParams &a, const ImageParams &b))->std::vector<ImageParams>::iterator{
		// Finds the lower bound in at most log(last - first) + 1 comparisons
		std::vector<ImageParams>::iterator it = std::lower_bound(begin, end, val, cmp);
		if (it != begin){
			return std::prev(it);
		}
		return end;
	};

	auto it = binary_find(_timeStamps.begin(), _timeStamps.end(), { currentTime, 0, "", "", false }, cmp);
	//check [start, stop] 
	if (closedInterval && (it == _timeStamps.end() || it->stopTime < currentTime || it->projected)){
		return false;
	}else if (!closedInterval && (it == _timeStamps.end() || it->projected)){
		return false;
	}
	double upcoming = nextCaptureTime(currentTime);
	if (_nextCapture != upcoming){
		_nextCapture = upcoming;
		_intervalLength = upcoming - currentTime;
	}

	it->projected = true;
	path = it->path;
	currentTime = it->startTime;
	
	return true;
}

bool ImageSequencer::sequenceReset(){
	for (auto image : _timeStamps){
		image.projected = false;
	}
	return true;
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = str.find(from);
	if (start_pos == std::string::npos)
		return false;
	str.replace(start_pos, from.length(), to);
	return true;
}

bool ImageSequencer::parsePlaybook(const std::string& dir, const std::string& type, std::string year){
	ghoul::filesystem::Directory playbookDir(dir, true);
	std::vector<std::string> dirlist = playbookDir.read(true, false);
	for (auto path : dirlist) {
        bool success = parsePlaybookFile(path, year);
        if (!success)
            return false;
    }
	return true; // add check
}

double ImageSequencer::getMissionElapsedTime(std::string timestr){
	std::string::size_type sz;     // alias of size_t
	double met = std::stod(timestr, &sz);
	double diff;
	double et;
	//met ref time.
	openspace::SpiceManager::ref().getETfromDate("2015-07-14T11:50:00.00", et);

	diff = abs(met - _metRef);
	if (met > _metRef){
		et += diff;
	}
	else if (met < _metRef){
		et -= diff;
	}

	return et;
}

bool ImageSequencer::parsePlaybookFile(const std::string& fileName, std::string year) {
	if (size_t position = fileName.find_last_of(".") + 1){
		if (position != std::string::npos){
			std::string extension = ghoul::filesystem::File(fileName).fileExtension();
#ifdef BACKUP_PLAYBOOK // Comma separated playbook in case we dont recieve updates from kang.
			if (extension == "csv"){ // comma separated playbook 
				std::cout << "USING COMMA SEPARATED TIMELINE V9F" << std::endl;

                std::string cachedFile = "";
                FileSys.cacheManager()->getCachedFile(fileName, cachedFile, true);

                bool hasCachedFile = FileSys.fileExists(cachedFile);
                if (hasCachedFile) {
                    std::ifstream file(cachedFile);
                    if (!file.good())
                        LERROR("Error loading cached playbook '" << cachedFile << "'");
                    else {
                        do {
                            std::string line;
                            std::getline(file, line);

                            std::stringstream s(line);

                            double start, end;
                            std::string path;

                            s >> start;
                            s >> end;

                            std::getline(s, path);
                            createImage(start, end, _defaultCaptureImage);
                        } while (!file.eof());
                    }                
                } else {
                    std::ifstream file(fileName);
                    if (!file.good()) LERROR("Failed to open csv file '" << fileName << "'");

                    std::string timestr = "";
                    double shutter = 0.01;
                    double et;
                    do{
                        std::getline(file, timestr);
                        auto pos = timestr.find("LORRI image started");
                        if (pos != std::string::npos){
                            timestr = timestr.substr(timestr.find_first_of(",") + 1);
                            timestr = timestr.substr(0, timestr.find_first_of(","));

                            replace(timestr, " ", "::");
                            timestr = year + " " + timestr;

                            openspace::SpiceManager::ref().getETfromDate(timestr, et);
                            createImage(et, et + shutter, _defaultCaptureImage);
                        }
                    } while (!file.eof());

                    std::sort(_timeStamps.begin(), _timeStamps.end(), cmp);

                    std::ofstream cachedFileStream(cachedFile);
                    cachedFileStream << std::setprecision(64);
                    if (cachedFileStream.good()) {
                        for (const ImageParams& i : _timeStamps)
                            cachedFileStream << i.startTime << "\t" << i.stopTime << "\t" << i.path << std::endl;
                    }

                }
			} 
#endif
			if (extension == "txt"){// Hong Kang. pre-parsed playbook
                LINFO("Using Preparsed Playbook V9H");
				std::ifstream file(fileName);
				if (!file.good()) LERROR("Failed to open txt file '" << fileName << "'");

				std::string timestr = "";
				double shutter = 0.01;
				double et;
				
				//@TODO: change similar to renderableFOV
				std::string instruments[3] = { "LORRI", "RALPH_LEISA", "MVIC" };
				std::string id[3] = { "NH_LORRI", "NH_RALPH_LEISA", "MVIC" };
				int isize = sizeof(instruments) / sizeof(instruments[0]);
				
				do{
					std::getline(file, timestr);
					for (int i = 0; i < isize; i++){
						auto pos = timestr.find(instruments[i]);
						if (pos != std::string::npos){
							timestr = timestr.substr(24, 9);
							et = getMissionElapsedTime(timestr);
							createImage(et, et + shutter, id[i], _defaultCaptureImage);
						}
					}
				} while (!file.eof());
                std::sort(_timeStamps.begin(), _timeStamps.end(), cmp);
			}
				
		}
	}
    return true;
}

bool ImageSequencer::loadSequence(const std::string& dir) {	
	ghoul::filesystem::Directory sequenceDir(dir, true);
	std::vector<std::string> sequencePaths = sequenceDir.read(true, false); // check inputs 
	for (auto path : sequencePaths){
		if (size_t position = path.find_last_of(".") + 1){
			if (position != std::string::npos){
				ghoul::filesystem::File currentFile(path);
				std::string extension = currentFile.fileExtension();

				if (extension == "lbl"){ // discovered header file 		
					std::ifstream file(currentFile.path());

					if (!file.good()) LERROR("Failed to open label file '" << currentFile.path() << "'");
				
					// open up label files
					std::string line = "";
					std::string specsOfInterest[2] = { "START_TIME", "STOP_TIME" };
					double timestamps[2] = { 0.0, 0.0 };
					bool found = false;
					do {
						std::getline(file, line);
						for (int i = 0; i < 2; i++){
							auto pos = line.find(specsOfInterest[i]);
							if (pos != std::string::npos){
								std::string time = line.substr(line.find("=") + 2);
								time.erase(std::remove(time.begin(), time.end(), ' '), time.end());
								openspace::SpiceManager::ref().getETfromDate(time, timestamps[i]);
							}
						}
						if (timestamps[0] != 0.0 && timestamps[1] != 0.0){
							found = true;
							std::string ext = "jpg";
							path.replace(path.begin() + position, path.end(), ext);
                            bool fileExists = FileSys.fileExists(path);
							if (fileExists) {
								createImage(timestamps[0], timestamps[1], "NH_LORRI", path); /// fix active instrument!
                                std::sort(_timeStamps.begin(), _timeStamps.end(), cmp);
							}
						}
					} while (!file.eof() && found == false);
				}
			}
		}
	}
	return !_timeStamps.empty();
}

}  // namespace openspace