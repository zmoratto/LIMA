// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
//
#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <iostream>

/*
#include <vw/Core/Exception.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>

#include <vw/Cartography.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
*/
using namespace std;

void FindAndReplace( std::string& tInput, std::string tFind, std::string tReplace );
vector<std::string> FindAndSplit(std::string& tInput, std::string tFind);
std::string GetFilenameNoExt(std::string const& filename);
std::string GetFilenameNoPath(std::string const& filename);
std::string GetFilenameExt(std::string const& filename);


#endif

