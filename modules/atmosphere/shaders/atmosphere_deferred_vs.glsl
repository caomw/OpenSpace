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

//#version __CONTEXT__
#version 400

layout(location = 0) in vec4 in_position;

//#include "PowerScaling/powerScaling_vs.hglsl"

uniform mat4 sgctProjectionMatrix;
uniform mat4 inverseSgctProjectionMatrix;
uniform mat4 objToWorldTransform;
uniform mat4 worldToObjectTransform;
uniform mat4 worldToEyeTransform;
uniform mat4 eyeToWorldTransform;
uniform mat4 eyeToViewTranform;
uniform mat4 viewToEyeTranform;

out vec3 interpolatedNDCPos;
out vec4 vertexPosObjVS;
out vec3 interpolatedRayDirection;

void main()
{
    //viewDirectionVS = normalize( (completeInverse * vec4((projInverse * in_position).xyz, 0.0)).xyz - cameraPosObj.xyz);

    //viewDirectionVS = normalize( (completeInverse * vec4(projInverse * in_position) ).xyz );

    //viewDirectionVS = (completeInverse * projInverse * in_position).xyz;
    interpolatedRayDirection = (viewToEyeTranform * vec4((inverseSgctProjectionMatrix * in_position).xyz, 0.0)).xyz;
    interpolatedNDCPos       = in_position.xyz;
    gl_Position              = in_position;
}