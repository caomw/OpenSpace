/*****************************************************************************************
*                                                                                       *
* OpenSpace                                                                             *
*                                                                                       *
* Copyright (c) 2014-2015                                                               *
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
#include <modules/iswa/util/iswamanager.h>
#include <modules/iswa/rendering/iswacygnet.h>
#include <ghoul/filesystem/filesystem>
#include <modules/kameleon/include/kameleonwrapper.h>
#include <modules/iswa/rendering/dataplane.h>
#include <modules/iswa/rendering/textureplane.h>
#include <openspace/util/time.h>
#include <modules/iswa/rendering/iswacontainer.h>
#include <modules/iswa/rendering/screenspacecygnet.h>
#include <modules/iswa/ext/json/json.hpp>
#include <fstream>

namespace {
    using json = nlohmann::json;
    const std::string _loggerCat = "ISWAManager";
}

namespace openspace{
    ISWAManager::ISWAManager()
        :_container(nullptr)
    {
        _month["JAN"] = "01";
        _month["FEB"] = "02";
        _month["MAR"] = "03";
        _month["APR"] = "04";
        _month["MAY"] = "05";
        _month["JUN"] = "06";
        _month["JUL"] = "07";
        _month["AUG"] = "08";
        _month["SEP"] = "09";
        _month["OCT"] = "10";
        _month["NOV"] = "11";
        _month["DEC"] = "12";
    }

    ISWAManager::~ISWAManager(){}

    std::shared_ptr<ISWACygnet> ISWAManager::createISWACygnet(std::shared_ptr<Metadata> metadata){
        LDEBUG("Creating ISWACygnet with id " << metadata->id);
        if(metadata->path != ""){
            const std::string& extension = ghoul::filesystem::File(absPath(metadata->path)).fileExtension();
            std::shared_ptr<ISWACygnet> cygnet;

            if(extension == "plain"){
                LWARNING("This cygnet image does not exist");
                return nullptr;
            }else if(extension == "cdf"){

                if(!FileSys.fileExists(absPath(metadata->path))) {
                    LERROR("Could not find cdf file:  " << absPath(metadata->path));
                    return nullptr;
                }

                std::shared_ptr<KameleonWrapper> kw = std::make_shared<KameleonWrapper>(absPath(metadata->path));
                auto parentNode = OsEng.renderEngine().scene()->sceneGraphNode(kw->getParent());

                if(parentNode){
                    ghoul::Dictionary metadataDic = 
                    {
                        // {std::string("Name"),        std::string("DataPlane")},
                        {std::string("Type"),       std::string("DataPlane")},
                        {std::string("StartTime"),  std::string("")},
                        {std::string("EndTime"),    std::string("")},
                        {std::string("Id"),         metadata->id},
                        {std::string("Path"),       metadata->path},
                        {std::string("Scale"),      glm::vec4(kw->getModelScaleScaled())},
                        {std::string("Offset"),     glm::vec4(kw->getModelBarycenterOffsetScaled())},
                        {std::string("Parent"),     kw->getParent()},
                        {std::string("Frame"),      kw->getFrame()},
                        {std::string("KW"),         kw}
                    };

                    
                    ghoul::Dictionary nodeDic = 
                    {
                        {std::string("Name"),       std::string("DataPlane")},
                        {std::string("Parent"),     kw->getParent()},
                        {std::string("Renderable"), metadataDic}
                    };
/*                  SceneGraphNode* cygnetNode = SceneGraphNode::createFromDictionary(nodeDic);
                    cygnetNode->setParent(parentNode);
                    parentNode->addChild(cygnetNode);
                    OsEng.renderEngine().scene()->addSceneGraphNode(cygnetNode);
                    cygnetNode->initialize();*/

                }
                // cygnet = std::make_shared<DataPlane>(metadataDic);               
            }else {
                auto parentNode = OsEng.renderEngine().scene()->sceneGraphNode(metadata->parent);
                if(parentNode){
                    std::string script = "openspace.addSceneGraphNode({"
                        "Name = 'TexturePlane"+std::to_string(metadata->id)+"',"
                        "Parent = '" + metadata->parent + "',"
                        "Renderable = {"
                            "Type = 'TexturePlane',"
                            "Id = "+std::to_string(metadata->id)+","
                            "Path = '" + metadata->path + "',"
                            "StartTime = '',"
                            "EndTime = '',"
                            "Frame = 'GALACTIC',"
                            "Parent = '" + metadata->parent + "',"
                            "Scale = { 3, 3, 3, 10},"
                            "Offset = {0, 0, 0, 1}"
                        "}"
                    "});";
                    std::cout << script << std::endl;
                    OsEng.scriptEngine().queueScript(script);
                    //cygnetNode->initialize();

                }else{
                    OsEng.renderEngine().registerScreenSpaceRenderable(std::make_shared<ScreenSpaceCygnet>(metadata));
                    return nullptr;
                } 
            }
            // cygnet->initialize();
            // return cygnet;
        }

        return nullptr;
    }

    void ISWAManager::addISWACygnet(std::string info){
        std::string token;
        std::stringstream ss(info);
        getline(ss,token,',');
        int cygnetId = std::stoi(token);

        getline(ss,token,',');
        std::string data = token;
        addISWACygnet(cygnetId, data);
    }

    void ISWAManager::addISWACygnet(int id, std::string info){
        if(id > 0){
            // std::shared_ptr<ExtensionFuture> extFuture = fileExtension(id);
            // extFuture->parent = info;
            // _extFutures.push_back(extFuture);
            // _container->addISWACygnet(cygnetId, data);
            createScreenSpace(id);
        }else if(id < 0){
            //download metadata to texture plane
            std::shared_ptr<MetadataFuture> extFuture = downloadMetadata(id);
            extFuture->type = "TEXTURE";
            extFuture->id = -id;
            // extFuture->parent = info;
            // _extFutures.push_back(extFuture);
            // std::shared_ptr<ExtensionFuture> extFuture
            _metaFutures.push_back(extFuture);
        }
        else {
            std::shared_ptr<MetadataFuture> extFuture = downloadMetadata(-1);
            extFuture->type = "DATA";
            extFuture->id = -1;
            // std::shared_ptr<Metadata> mdata = std::make_shared<Metadata>();
            // mdata->id = 0;
            // mdata->path = absPath("${OPENSPACE_DATA}/"+info);
            // createISWACygnet(mdata);
            _metaFutures.push_back(extFuture);
            // createDataPlane(absPath("${OPENSPACE_DATA}/"+info));
        }
    }

    void ISWAManager::deleteISWACygnet(std::string name){
        //_container->deleteISWACygnet(name);
        OsEng.scriptEngine().queueScript("openspace.removeSceneGraphNode('" + name + "')");
        
    }

    std::shared_ptr<DownloadManager::FileFuture> ISWAManager::downloadImage(int id, std::string path){
        return  DlManager.downloadFile(
                    iSWAurl(id),
                    path,
                    true,
                    [path](const DownloadManager::FileFuture& f){
                        LDEBUG("Download finished: " << path);
                    }
                );
    }

    std::shared_ptr<DownloadManager::FileFuture> ISWAManager::downloadImageToMemory(int id, std::string& buffer){
        return  DlManager.downloadToMemory(
                    iSWAurl(id),
                    buffer,
                    [](const DownloadManager::FileFuture& f){
                        LDEBUG("Download to memory finished");
                    }
                );
    }

    std::shared_ptr<DownloadManager::FileFuture> ISWAManager::downloadDataToMemory(int id, std::string& buffer){
        return DlManager.downloadToMemory(
                iSWADataUrl(id),
                buffer,
                [](const DownloadManager::FileFuture& f){
                    LDEBUG("Download data finished");
                }
            );
    }

    std::shared_ptr<MetadataFuture> ISWAManager::downloadMetadata(int id){
        std::shared_ptr<MetadataFuture> metaFuture = std::make_shared<MetadataFuture>();

        metaFuture->id = id;

        std::ifstream file(absPath("${OPENSPACE_DATA}/GM_openspace_Y0_info.txt"));
        if(file.is_open()){
            std::string json( (std::istreambuf_iterator<char>(file) ),
                                (std::istreambuf_iterator<char>()));
            metaFuture->isFinished = true;
            metaFuture->json = json;
        }

        return metaFuture;
    }

    std::shared_ptr<ExtensionFuture> ISWAManager::fileExtension(int id){
        std::shared_ptr<ExtensionFuture> extFuture = std::make_shared<ExtensionFuture>();
        extFuture->isFinished = false;
        extFuture->id = id;
        DlManager.getFileExtension(
                iSWAurl(id),
                [extFuture](std::string extension){
                    std::stringstream ss(extension);
                    std::string token;
                    std::getline(ss, token, '/');
                    std::getline(ss, token);

                    std::string ext = "."+token;
                    extFuture->extension = ext;
                    extFuture->isFinished = true;
                }
            );

        return extFuture;
    }

    void ISWAManager::setContainer(ISWAContainer* container){
        _container = container;
    }


    std::shared_ptr<ISWACygnet> ISWAManager::iSWACygnet(std::string name){
        if(_container)
            return _container->iSWACygnet(name);
        return nullptr;
    }

    std::string ISWAManager::iSWAurl(int id){
        std::string url = "http://iswa2.ccmc.gsfc.nasa.gov/IswaSystemWebApp/iSWACygnetStreamer?timestamp=";
        
        std::string t = Time::ref().currentTimeUTC(); 
        std::stringstream ss(t);
        std::string token;

        std::getline(ss, token, ' ');
        url += token + "-"; 
        std::getline(ss, token, ' ');
        url += _month[token] + "-";
        std::getline(ss, token, 'T');
        url += token + "%20";
        std::getline(ss, token, '.');
        url += token;

        url += "&window=-1&cygnetId=";
        url += std::to_string(id);

        //std::cout << url <<  std::endl;

        return url;
    }

    std::string ISWAManager::iSWADataUrl(int id){
        std::string url = "http://128.183.168.116:3000/data/" + std::to_string(id) + "/"; 
        // /2996-01-23%2000:44:00
        std::string t = Time::ref().currentTimeUTC(); 
        std::stringstream ss(t);
        std::string token;

        std::getline(ss, token, ' ');
        url += token + "-"; 
        std::getline(ss, token, ' ');
        url += _month[token] + "-";
        std::getline(ss, token, 'T');
        url += token + "%20";
        std::getline(ss, token, '.');
        url += token;

        return url;
    }

    void ISWAManager::update(){
        for (auto it = _extFutures.begin(); it != _extFutures.end(); )
        {
            if ((*it)->isFinished) {
                std::string path = "${OPENSPACE_DATA}/scene/iswa/" + std::to_string((*it)->id) + (*it)->extension;
                
                std::shared_ptr<Metadata> data = std::make_shared<Metadata>();
                data->id = (*it)->id;
                data->path = path;
                data->parent = (*it)->parent;

                createISWACygnet(data);
                it = _extFutures.erase( it );
            }
            else {
                ++it;
            }
        }

        for (auto it = _metaFutures.begin(); it != _metaFutures.end(); ){
            if((*it)->isFinished) {
                if((*it)->type == "TEXTURE"){
                    createPlane((*it)->id,(*it)->json,std::string("TexturePlane"));
                }else{
                    createPlane((*it)->id,(*it)->json,std::string("DataPlane"));

                }
                it = _metaFutures.erase( it );
            }else{
                ++it;
            }
        }    
    }

    std::string ISWAManager::getDictionaryTable(int id, std::string path){
        json j;
        std::ifstream file(path);
        if(file.is_open()){
            j = json::parse(file);
        }

        std::string parent = j["Central Body"];
        std::string frame = j["Coordinates"];

        int xmax = j["Plot XMAX"];
        int ymax = j["Plot YMAX"];
        int zmax = j["Plot ZMAX"];
        int xmin = j["Plot XMIN"];
        int ymin = j["Plot YMIN"];
        int zmin = j["Plot ZMIN"];

        float spatScale=1, scalew=10;
        std::string spatialScale = j["Spatial Scale (Custom)"];
        if(spatialScale == "R_E"){
            // spatScale = 6.371f;
            // scalew = 6;
        }

        std::string scale = "{" 
                                + std::to_string(spatScale*(xmax-xmin)) + ","
                                + std::to_string(spatScale*(ymax-ymin)) + ","
                                + std::to_string(spatScale*(zmax-zmin)) + ","
                                + std::to_string(scalew) +
                            "}";

        std::string offset ="{"
                                + std::to_string(spatScale*(xmin + (std::abs(xmin)+std::abs(xmax))/2.0f)) + "," 
                                + std::to_string(spatScale*(ymin + (std::abs(ymin)+std::abs(ymax))/2.0f)) + ","
                                + std::to_string(spatScale*(zmin + (std::abs(zmin)+std::abs(zmax))/2.0f)) + ","
                                + std::to_string(scalew) +
                            "}";

        std::string table = "{"
        "Name = 'TexturePlane' , "
        "Parent = '" + parent + "', "
        "Renderable = {"    
            "Type = 'TexturePlane', "
            "Id = " + std::to_string(id) + ", "
            "Frame = '" + frame + "' , "
            "Scale = " + scale + ", "
            "Offset = " + offset + 
            "}"
        "}"
        ;

        std::cout << table << std::endl;
        // ghoul::Dictionary dic;
        return table;
    }

    std::string ISWAManager::parseJSONToLuaTable(int id, std::string jsonString, std::string type){
        if(jsonString != ""){
            json j = json::parse(jsonString);

            std::string parent = j["Central Body"];
            std::string frame = j["Coordinates"];

            int xmax = j["Plot XMAX"];
            int ymax = j["Plot YMAX"];
            int zmax = j["Plot ZMAX"];
            int xmin = j["Plot XMIN"];
            int ymin = j["Plot YMIN"];
            int zmin = j["Plot ZMIN"];

            float spatScale=1, scalew=10;
            std::string spatialScale = j["Spatial Scale (Custom)"];
            if(spatialScale == "R_E"){
                spatScale = 6.371f;
                scalew = 6;
            }

            glm::vec4 scale(
                spatScale*(xmax-xmin),
                spatScale*(ymax-ymin),
                spatScale*(zmax-zmin),
                scalew
            );
            glm::vec4 offset (
                spatScale*(xmin + (std::abs(xmin)+std::abs(xmax))/2.0f),
                spatScale*(ymin + (std::abs(ymin)+std::abs(ymax))/2.0f),
                spatScale*(zmin + (std::abs(zmin)+std::abs(zmax))/2.0f),
                scalew
            );

            std::string table = "{"
            "Name = 'TexturePlane "+ std::to_string(id) +"' , "
            "Parent = '" + parent + "', "
            "Renderable = {"    
                "Type = '" + type + "', "
                "Id = " + std::to_string(id) + ", "
                "Frame = '" + frame + "' , "
                "Scale = " + std::to_string(scale) + ", "
                "Offset = " + std::to_string(offset) + 
                "}"
            "}"
            ;

            std::cout << table << std::endl;
            // ghoul::Dictionary dic;
            return table;
        }
        return "";
    }

    std::string ISWAManager::parseKWToLuaTable(std::string kwPath){

        if(kwPath != ""){
            const std::string& extension = ghoul::filesystem::File(absPath(kwPath)).fileExtension();
            if(extension == "cdf"){
                KameleonWrapper kw = KameleonWrapper(absPath(kwPath));
         
                std::string parent  = kw.getParent();
                std::string frame   = kw.getFrame();
                glm::vec4 scale     = kw.getModelScaleScaled();
                glm::vec4 offset    = kw.getModelBarycenterOffsetScaled();

                std::string table = "{"
                    "Name = 'DataPlane',"
                    "Parent = '" + parent + "', "
                    "Renderable = {"    
                        "Type = 'DataPlane', "
                        "Id = 0 ,"
                        "Frame = '" + frame + "' , "
                        "Scale = " + std::to_string(scale) + ", "
                        "Offset = " + std::to_string(offset) + ", "
                        "kwPath = '" + kwPath + "'" 
                        "}"
                    "}"
                    ;
                std::cout << table << std::endl;    
                return table;
            }
        }
        return "";
    }

    //Create KameleonPlane?
    // void ISWAManager::createDataPlane(std::string kwPath){
    //     std::string luaTable = parseKWToLuaTable(kwPath);
    //     if(luaTable != ""){
    //         std::string script = "openspace.addSceneGraphNode(" + luaTable + ");";
    //         OsEng.scriptEngine().queueScript(script);
    //     }
    // }


    // void ISWAManager::createTexturePlane(int id, std::string json){
    //     std::string luaTable = parseJSONToLuaTable(id, json);
    //     if(luaTable != ""){
    //         std::string script = "openspace.addSceneGraphNode(" + parseJSONToLuaTable(id, json) + ");";
    //         OsEng.scriptEngine().queueScript(script);
    //     }
    // }

    void ISWAManager::createPlane(int id, std::string json, std::string type){
        std::string luaTable = parseJSONToLuaTable(id, json, type);
        if(luaTable != ""){
            std::string script = "openspace.addSceneGraphNode(" + luaTable + ");";
            OsEng.scriptEngine().queueScript(script);
        }
    }

    void ISWAManager::createScreenSpace(int id){
        OsEng.renderEngine().registerScreenSpaceRenderable(std::make_shared<ScreenSpaceCygnet>(id));
    }
}// namsepace openspace