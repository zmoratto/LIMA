// __BEGIN_LICENSE__
// Copyright (C) 2012 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef ALIGNSETTINGS_H
#define ALIGNSETTINGS_H

boost::program_options::options_description getCommonSettings();
std::string CommonSettings_to_string( const boost::program_options::variables_map& vm );
std::string VariablesMap_to_string( const boost::program_options::variables_map& vm );

#endif
