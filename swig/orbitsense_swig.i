/* -*- c++ -*- */

#define ORBITSENSE_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "orbitsense_swig_doc.i"

%{
#include "orbitsense/detection_engine.h"
%}


%include "orbitsense/detection_engine.h"
GR_SWIG_BLOCK_MAGIC2(orbitsense, detection_engine);
