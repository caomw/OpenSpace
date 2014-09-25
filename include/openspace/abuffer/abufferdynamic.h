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

#ifndef __ABUFFERDYNAMIC_H__
#define __ABUFFERDYNAMIC_H__

#include <openspace/abuffer/abuffer.h>

namespace openspace {

class ABufferDynamic: public ABuffer {
public:

	ABufferDynamic();
	virtual ~ABufferDynamic();
	virtual bool initialize();

	virtual void clear();
	virtual void preRender();
	virtual void postRender();

	virtual std::string settings();

private:

	GLuint *_data;
	GLuint _anchorPointerTexture;
	GLuint _anchorPointerTextureInitializer;
	GLuint _atomicCounterBuffer;
	GLuint _fragmentBuffer;
	GLuint _fragmentTexture;




}; 		// ABufferDynamic
} 		// openspace

#endif 	// __ABUFFERDYNAMIC_H__