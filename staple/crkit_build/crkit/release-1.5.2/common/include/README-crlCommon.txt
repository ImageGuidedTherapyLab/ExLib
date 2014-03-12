/*
 * Copyright (c) 2008-2009 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/


================================================================================
 HOW TO INCLUDE crlCommon IN A PROJECT
================================================================================

Just insert the three following lines (and modify the value of CRKIT_SOURCE_DIR):

  SET(CRKIT_SOURCE_DIR "XXXX")
  INCLUDE_DIRECTORIES( ${CRKIT_SOURCE_DIR}/common/include )
  SUBDIRS ( ${CRKIT_SOURCE_DIR}/common/include )

It will creates a new target in your build tree to build crlCommon.

