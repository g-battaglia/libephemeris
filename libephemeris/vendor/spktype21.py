# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
# Copyright (c) 2018 The Python Packaging Authority
# Third-party/adapted code — see file docstring and THIRD_PARTY_NOTICES.md
"""A supporting module for jplephem to handle data type 21 (Version 0.1.0)

VENDORED by libephemeris from spktype21==0.1.0 (PyPI).
Original author: Shushi Uetsuki (whiskie14142), MIT License.
Upstream is unmaintained since 2018.

Changes from upstream:
  - Lines 568, 599: Added .item() to daf.map_array() calls for numpy 2.x
    compatibility (numpy 2.x rejects int() on 1-element arrays).
  - Module header: added the explicit vendoring and provenance record; the
    numerical implementation is otherwise the identified 0.1.0 upstream.

Provenance:
  This file remains MIT-licensed third-party code by Shushi Uetsuki, vendored
  from the exact upstream/version named above. Its scientific input format is
  NASA/JPL/NAIF SPK type 21. LibEphemeris does not claim this implementation as
  project-owned and does not relicense it; see ``THIRD_PARTY_NOTICES.md``.

This module computes position and velocity of a celestial small body, from a 
NASA SPICE SPK ephemeris kernel file of data type 21 (Extended Modified 
Difference Arrays).
http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

You can get SPK files for many solar system small bodies from HORIZONS 
system of NASA/JPL.  See https://ssd.jpl.nasa.gov/?horizons

This module reads SPK files of data type 21, one of the types of binary SPK 
file.  

At the point of Oct. 2018, HORIZONS system provides files of type 21 for 
binary SPK files by default.  You can get type 21 binary SPK file for celestial 
small bodies through TELNET interface by answering back 'Binary' for 
'SPK file format'.  Also you can get type 21 binary SPK file from:
https://ssd.jpl.nasa.gov/x/spk.html

Modules required:
    jplephem (version 2.6 or later)
    numpy

Usage:
    from spktype21 import SPKType21
    kernel = SPKType21.open('path')
    position, velocity = kernel.compute_type21(center, target, jd)
    
    where:
        path - path to the SPK file
        center - SPKID of central body (0 for SSB, 10 for Sun, etc.)
        target - SPKID of target body
        jd - time for computation (Julian date)

Exceptions:
    RuntimeError will be raised when:
        invalid data_type of SPK file, or
        SPK file contains too large table in EMDA record(s)
    ValueError will be raised when:
        invalid parameter(s) of compute_type21 function

Author: Shushi Uetsuki (whiskie14142)
This module has been developed based on jplephem and FORTRAN source 
of the SPICE Toolkit of NASA/JPL/NAIF.
jplephem : https://pypi.org/project/jplephem/
SPICE Toolkit : http://naif.jpl.nasa.gov/naif/toolkit.html
"""

import threading

from numpy import array, zeros, reshape
from jplephem.daf import DAF
from jplephem.names import target_names

T0 = 2451545.0
S_PER_DAY = 86400.0

# Included from 'spk21.inc' on the FORTRAN source 'spke21.f'
MAXTRM = 25

def jd(seconds):
    """Convert a number of seconds since J2000 to a Julian Date.
    """
    return T0 + seconds / S_PER_DAY

class SPKType21(object):
    """Class for SPK kernel to handle data type 21 (Extended Modified Difference Arrays)
    """
    def __init__(self, daf):
        self.daf = daf
        self.segments = [Segment(self.daf, *t) for t in self.daf.summaries()]
        ssec = lambda s : s.start_second
        self.segments.sort(key=ssec)
        
        # initialize arrays for spke21
        self.G = zeros(MAXTRM)
        
        self.REFPOS = zeros(3)
        self.REFVEL = zeros(3)
        
        self.KQ = array([0, 0, 0])
        self.FC = zeros(MAXTRM)
        self.FC[0] = 1.0
        self.WC = zeros(MAXTRM - 1)
        self.W = zeros(MAXTRM + 2)
        
        # initialize for compute_type21
        # Per-body cache of the last MDA record: (target, center) -> (record,
        # lb, ub). The mda_lb/mda_ub window is NOT body-specific, so a single
        # slot keyed only on time would return the wrong body's state when
        # multiple type-21 targets share one .bsp/instance. Keying on
        # (target, center) also avoids cache thrashing when callers interleave
        # several bodies in the same time window. The key space is tiny (a
        # handful of small-body targets), so no eviction is needed.
        self._mda_cache = {}
        self.current_segment_exist = False

        # libephemeris addition: the evaluation work arrays (G/FC/WC/W/KQ,
        # REFPOS/REFVEL, mda_record, current_segment) are instance state
        # mutated by spke21()/compute_type21(); kernels are cached and
        # shared across calc threads, so evaluation must be serialized.
        self._compute_lock = threading.Lock()
        
    @classmethod
    def open(cls, path):
        """Open the file at `path` and return an SPK instance.
        """
        return cls(DAF(open(path, 'rb')))

    def close(self):
        """Close this SPK file."""
        self.daf.file.close()

    def __str__(self):
        daf = self.daf
        d = lambda b: b.decode('latin-1')
        lines = (str(segment) for segment in self.segments)
        return 'File type {0} and format {1} with {2} segments:\n{3}'.format(
            d(daf.locidw), d(daf.locfmt), len(self.segments), '\n'.join(lines))
    
    def comments(self):
        return self.daf.comments()

    def compute_type21(self, center, target, jd1, jd2=0.0):
        """Compute position and velocity of target from SPK data (data type 21).
        Inputs:
            center - SPKID of the coordinate center (0 for Solar System Barycenter, 
                     10 for Sun, etc)
            target - SPKID of the target
            jd1, jd2 - Julian date of epoch for computation.  (jd1 + jd2) will 
                be used for computation.  If you want precise definition of 
                epoch, jd1 should be an integer or a half integer, and jd2 should
                be a relatively small floating point number.
        Returns:
            Position (X, Y, Z) and velocity (XD, YD, ZD) of the target at 
            epoch.  Position and velocity are provided as Numpy arrays 
            respectively.
        """
        eval_sec = (jd1 - T0)
        eval_sec = (eval_sec + jd2) * S_PER_DAY

        with self._compute_lock:
            key = (target, center)
            cached = self._mda_cache.get(key)
            if cached is not None:
                mda_record, mda_lb, mda_ub, mda_segment = cached
                if mda_lb <= eval_sec < mda_ub:
                    # spke21 reads MAXDIM from self.current_segment; restore the
                    # segment this record came from so a cache hit for one body
                    # does not inherit a stale MAXDIM from another body whose
                    # type-21 segments share this same SPKType21 instance.
                    self.current_segment = mda_segment
                    self.current_segment_exist = True
                    result = self.spke21(eval_sec, mda_record)
                    return result[0:3], result[3:]

            mda_record, mda_lb, mda_ub = self.get_MDA_record(
                eval_sec, target, center
            )
            self._mda_cache[key] = (
                mda_record,
                mda_lb,
                mda_ub,
                self.current_segment,
            )

            result = self.spke21(eval_sec, mda_record)
            return result[0:3], result[3:]
                
    def get_MDA_record(self, eval_sec, target, center):
        """Return a EMDA record for defined epoch.
        Inputs:
            eval_sec - epoch for computation, seconds from J2000
            target - body ID of the target
            center - body ID of coordinate center
        Returns:
            EMDA record - a Numpy array of DLSIZE floating point numbers
        Exception:
            ValueError will be raised when:
                eval_sed is outside of SPK data
                target and center are not in SPK data
            RuntimeError will be raised when:
                invalid data type of SPK data
        """
        
        # chech last segment can be used
        if self.current_segment_exist:
            if eval_sec >= self.current_segment.start_second    \
                and eval_sec < self.current_segment.end_second  \
                and target == self.current_segment.target            \
                and center == self.current_segment.center:
                
                return self.current_segment.get_MDA_record(eval_sec)

        # select segments with matched 'target' and 'center'
        matched = []
        for segment in self.segments:
            if segment.target == target and segment.center == center:
                matched.append(segment)
        if len(matched) == 0:
            raise ValueError('Invalid Target and/or Center')
        if eval_sec < matched[0].start_second or eval_sec >= matched[-1].end_second:
            raise ValueError('Invalid Time to evaluate')
        
        # selet a segment based on eval_sec
        found = False
        for segment in matched:
            if eval_sec < segment.end_second:
                found = True
                self.current_segment = segment
                break
        if not found:
            self.current_segment = matched[-1]
        self.current_segment_exist = True
        
        # get the MDA record from selected segment
        if self.current_segment.data_type != 21:
            raise RuntimeError('Invalid data. Data Type must be 21')
        
        return self.current_segment.get_MDA_record(eval_sec)
        


# left this module only 2018/10/12

    def spke21(self, ET, RECORD):
        """Compute position and velocity from a Modified Difference Array record
        
        Inputs:
            ET: Epoch time to evaluate position and velocity (seconds since J2000)
            RECORD: A record of Extended Modified Difference Array
        Returns: STATE
            STATE: A numpy array which contains position and velocity
        """
        
# This method was translated from FORTRAN source code ‘spke21.f’ of SPICE 
# Toolkit and modified by Shushi Uetsuki.
# 
# SPICE Toolkit for FORTRAN : http://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
# SPK Required Reading : http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
# 
# Unfortunately, I found some discrepancies between FORTRAN source code 
# and actual data contained in SPK files.  So, I tried to compose a 
# method that compute positions and velocities correctly by referencing 
# code of spktype01.

# Following comments start with #C were copied from original FORTRAN code.

#C$ Abstract
#C
#C     Evaluate a single SPK data record from a segment of type 21
#C     (Extended Difference Lines).
#C
#C$ Disclaimer
#C
#C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
#C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
#C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
#C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
#C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
#C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
#C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
#C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
#C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
#C
#C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
#C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
#C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
#C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
#C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
#C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
#C
#C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
#C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
#C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
#C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
#C
#C$ Required_Reading
#C
#C     SPK
#C     TIME
#C
#C$ Keywords
#C
#C     EPHEMERIS
#C
#C$ Declarations

        STATE = zeros(6)

#C$ Brief_I/O
#C
#C     Variable  I/O  Description
#C     --------  ---  --------------------------------------------------
#C     ET         I   Evaluation epoch.
#C     RECORD     I   Data record.
#C     STATE      O   State (position and velocity).
#C     MAXTRM     P   Maximum number of terms per difference table
#C                    component.
#C
#C$ Detailed_Input
#C
#C     ET          is an epoch at which a state vector is to be
#C                 computed. The epoch is represented as seconds past
#C                 J2000 TDB.
#C
#C     RECORD      is a data record which, when evaluated at epoch ET,
#C                 will give the state (position and velocity) of an
#C                 ephemeris object, relative to its center of motion,
#C                 in an inertial reference frame.
#C
#C                 The contents of RECORD are as follows:
#C
#C                    RECORD(1):         The difference table size per
#C                                       Cartesian component. Call this
#C                                       size MAXDIM; then the difference
#C                                       line (MDA) size DLSIZE is
#C
#C                                         ( 4 * MAXDIM ) + 11
#C                                    
#C                    RECORD(2)
#C                       ...
#C                    RECORD(1+DLSIZE):  An extended difference line.
#C                                       The contents are:
#C
#C                       Dimension  Description
#C                       ---------  ----------------------------------
#C                       1          Reference epoch of difference line
#C                       MAXDIM     Stepsize function vector
#C                       1          Reference position vector,  x
#C                       1          Reference velocity vector,  x
#C                       1          Reference position vector,  y
#C                       1          Reference velocity vector,  y
#C                       1          Reference position vector,  z
#C                       1          Reference velocity vector,  z
#C                       MAXDIM,3   Modified divided difference
#C                                  arrays (MDAs)
#C                       1          Maximum integration order plus 1
#C                       3          Integration order array
#C
#C$ Detailed_Output
#C
#C     STATE       is the state resulting from evaluation of the input
#C                 record at ET. Units are km and km/sec.
#C
#C$ Parameters
#C
#C     MAXTRM      is the maximum number of terms allowed in
#C                 each component of the difference table 
#C                 contained in the input argument RECORD.
#C                 See the INCLUDE file spk21.inc for the value
#C                 of MAXTRM.
#C                  
#C$ Exceptions
#C
#C     1) If the maximum table size of the input record exceeds 
#C        MAXTRM, the error SPICE(DIFFLINETOOLARGE) is signaled.
#C
#C$ Files
#C
#C     None.
#C
#C$ Particulars
#C
#C     The exact format and structure of type 21 (difference lines)
#C     segments are described in the SPK Required Reading file.
#C
#C     SPKE21 is a modified version of SPKE01. The routine has been
#C     generalized to support variable size difference lines.
#C
#C$ Examples
#C
#C     None.
#C
#C$ Restrictions
#C
#C     Unknown.
#C
#C$ Literature_References
#C
#C     NAIF Document 168.0, "S- and P- Kernel (SPK) Specification and
#C     User's Guide"
#C
#C$ Author_and_Institution
#C
#C     N.J. Bachman    (JPL)
#C     F.T. Krogh      (JPL)
#C     W.L. Taber      (JPL)
#C     I.M. Underwood  (JPL)
#C
#C$ Version
#C
#C-    SPICELIB Version 1.0.0, 03-FEB-2014 (NJB) (FTK) (WLT) (IMU)
#C
#C-&
# 
#C$ Index_Entries
#C
#C     evaluate type_21 spk segment
#C
#C-&

#C
#C     The first element of the input record is the dimension
#C     of the difference table MAXDIM. 
#C

# The FORTRAN source code indicates that RECORD[0] contains MAXDIM, but actual 
# data record does not contain it. MAXDIM is contained in each segment.

        MAXDIM = self.current_segment.MAXDIM

        
        if MAXDIM > MAXTRM:
            mes = ('SPKE21 \nThe input record has a maximum table dimension ' +
                'of {0}, while the maximum supported by this routine is {1}. ' +
                'It is possible that this problem is due to your software ' +
                'beeing out of date.').format(MAXDIM, MAXTRM)
            raise RuntimeError(mes)
            return STATE
        
#C
#C     Unpack the contents of the MDA array.
#C
#C        Name     Dimension  Description
#C        ------   ---------  -------------------------------
#C        TL               1  Reference epoch of record
#C        G           MAXDIM  Stepsize function vector
#C        REFPOS           3  Reference position vector
#C        REFVEL           3  Reference velocity vector
#C        DT      MAXDIM,NTE  Modified divided difference arrays
#C        KQMAX1           1  Maximum integration order plus 1
#C        KQ             NTE  Integration order array
#C
#C     For our purposes, NTE is always 3.
#C

# The FORTRAN source code indicates that RECORD[1] contains TL,  but on the 
# actual data RECORD[0] contains it, and all addresses for following data are 
# shifted forward by one.

        self.TL = RECORD[0]
        self.G = RECORD[1:MAXDIM + 1]

#C     
#C     Collect the reference position and velocity.
#C     
        self.REFPOS[0] = RECORD[MAXDIM + 1]
        self.REFVEL[0] = RECORD[MAXDIM + 2]
        
        self.REFPOS[1] = RECORD[MAXDIM + 3]
        self.REFVEL[1] = RECORD[MAXDIM + 4]
        
        self.REFPOS[2] = RECORD[MAXDIM + 5]
        self.REFVEL[2] = RECORD[MAXDIM + 6]
        
#C
#C     Initializing the difference table is one aspect of this routine
#C     that's a bit different from SPKE01. Here the first dimension of
#C     the table in the input record can be smaller than MAXTRM. So, we
#C     must transfer separately the portions of the table corresponding
#C     to each component.
#C
        self.DT = reshape(RECORD[MAXDIM + 7:MAXDIM * 4 + 7], (MAXDIM, 3), 
                          order='F')
        
        self.KQMAX1 = int(RECORD[4 * MAXDIM +  7])
        self.KQ[0] =  int(RECORD[4 * MAXDIM +  8])
        self.KQ[1] =  int(RECORD[4 * MAXDIM +  9])
        self.KQ[2] =  int(RECORD[4 * MAXDIM + 10])
#C     
#C     Next we set up for the computation of the various differences
#C     
        self.DELTA = ET - self.TL
        self.TP = self.DELTA
        self.MQ2 = self.KQMAX1 - 2
        self.KS = self.KQMAX1 - 1

#C
#C     This is clearly collecting some kind of coefficients.  
#C     The problem is that we have no idea what they are...
#C     
#C     The G coefficients are supposed to be some kind of step size 
#C     vector. 
#C     
#C     TP starts out as the delta t between the request time and the
#C     difference line's reference epoch. We then change it from DELTA
#C     by the components of the stepsize vector G.
#C
        for J in range(1, self.MQ2 + 1):
#C
#C        Make sure we're not about to attempt division by zero.
#C
            if self.G[J-1] == 0.0:
                mes = ('SPKE21\nA value of zero was found at index {0} ' + 
                'of the step size vector.').format(J)
                raise RuntimeError(mes)
                return STATE
                
            self.FC[J] = self.TP / self.G[J-1]
            self.WC[J-1] = self.DELTA / self.G[J-1]
            self.TP = self.DELTA + self.G[J-1]

#C
#C     Collect KQMAX1 reciprocals. 
#C   
        for J in range(1, self.KQMAX1 + 1):
            self.W[J-1] = 1.0 / float(J)

#C
#C     Compute the W(K) terms needed for the position interpolation
#C     (Note,  it is assumed throughout this routine that KS, which 
#C     starts out as KQMAX1-1 (the ``maximum integration'') 
#C     is at least 2.
#C
        self.JX = 0
        self.KS1 = self.KS - 1
        
        while self.KS >= 2:
            
            self.JX = self.JX + 1
            
            for J in range(1, self.JX + 1):
                self.W[J+self.KS-1] = self.FC[J] * self.W[J+self.KS1-1] - self.WC[J-1] * self.W[J+self.KS-1]
            
            self.KS = self.KS1
            self.KS1 = self.KS1 - 1

#C
#C     Perform position interpolation: (Note that KS = 1 right now.
#C     We don't know much more than that.)
#C
        for I in range(1, 3 + 1):
            
            self.KQQ = self.KQ[I-1]
            self.SUM = 0.0
            
            for J in range(self.KQQ, 0, -1):
                self.SUM = self.SUM + self.DT[J-1, I-1] * self.W[J+self.KS-1]
            
            STATE[I-1] = self.REFPOS[I-1] + self.DELTA * (self.REFVEL[I-1] + self.DELTA * self.SUM)

#C
#C     Again we need to compute the W(K) coefficients that are 
#C     going to be used in the velocity interpolation. 
#C     (Note, at this point, KS = 1, KS1 = 0.)
#C      
        for J in range(1, self.JX + 1):
            self.W[J+self.KS-1] = self.FC[J] * self.W[J+self.KS1-1] - self.WC[J-1] * self.W[J+self.KS-1]
        
        self.KS = self.KS - 1
        
#C
#C     Perform velocity interpolation:
#C
        for I in range(1, 3 + 1):
            self.KQQ = self.KQ[I-1]
            self.SUM = 0.0
            
            for J in range(self.KQQ, 0, -1):
                self.SUM = self.SUM + self.DT[J-1, I-1] * self.W[J+self.KS-1]
            
            STATE[I+3-1] = self.REFVEL[I-1] + self.DELTA * self.SUM
        
        return STATE
        
        

class Segment(object):
    """A single segment of a SPK file.

    There are several items of information about each segment that are
    loaded from the underlying SPK file, and made available as object
    attributes:

    segment.source - official ephemeris name, like 'DE-0430LE-0430'
    segment.start_second - initial epoch, as seconds from J2000
    segment.end_second - final epoch, as seconds from J2000
    segment.start_jd - start_second, converted to a Julian Date
    segment.end_jd - end_second, converted to a Julian Date
    segment.center - integer center identifier
    segment.target - integer target identifier
    segment.frame - integer frame identifier
    segment.data_type - integer data type identifier
    segment.start_i - index where segment starts
    segment.end_i - index where segment ends
    """
    def __init__(self, daf, source, descriptor):
        self.daf = daf
        self.source = source
        (self.start_second, self.end_second, self.target, self.center,
         self.frame, self.data_type, self.start_i, self.end_i) = descriptor
        self.start_jd = jd(self.start_second)
        self.end_jd = jd(self.end_second)
        
# 'SPK Required Reading' indicates that the penultimate element of the segment 
# is the difference line size (DLSIZE), but actual data contains there a MAXDIM.
        
        self.MAXDIM = int(self.daf.map_array(self.end_i - 1, self.end_i - 1).item())
        self.DLSIZE = 4 * self.MAXDIM + 11

    def __str__(self):
        return self.describe(verbose=False)

    def describe(self, verbose=True):
        """Return a textual description of the segment.
        """
        center = titlecase(target_names.get(self.center, 'Unknown center'))
        target = titlecase(target_names.get(self.target, 'Unknown target'))
        text = ('{0.start_jd:.2f}..{0.end_jd:.2f}  {1} ({0.center})'
                ' -> {2} ({0.target})'
                ' data_type={0.data_type}'.format(self, center, target))
        if verbose:
            text += ('\n  frame={0.frame} data_type={0.data_type} source={1}'
                     .format(self, self.source.decode('ascii')))
        return text
        
    def get_MDA_record(self, time_sec):
        """Return a Modified Difference Array(MDA) record for the time to 
        evaluate with its effective time boundaries (lower and upper).
        Inputs:
            time_sec - epoch for computation, seconds from J2000
        Returns: mda_record, lower_boundary, upper_boundary
            mda_record: A Modified Difference Array record
            lower_boundary: lower boundary of the record, seconds since J2000
            upper_boundary: upper boundary of the record, seconds since J2000
        """

        # Number of records in this segment
        entry_count = int(self.daf.map_array(self.end_i, self.end_i).item())
        
        # Number of entries in epoch directory 
        epoch_dir_count = entry_count // 100
        
        # serch target epoch in epoch directory to narrow serching aria
        if epoch_dir_count >= 1:
            epoch_dir = self.daf.map_array(self.end_i - epoch_dir_count - 1,
                                            self.end_i - 2)
            found = False
            for i in range(1, epoch_dir_count + 1):
                if epoch_dir[i-1] > time_sec:
                    found = True
                    break
            if found:
                serch_last_index = i * 100
                serch_start_index = (i - 1) * 100 + 1
            else:
                serch_last_index = entry_count
                serch_start_index = epoch_dir_count * 100 + 1
        else:
            serch_last_index = entry_count
            serch_start_index = 1

        # epoch_table contains epochs for all records in this segment        
        epoch_table = self.daf.map_array(self.start_i + (entry_count * self.DLSIZE),
                                       self.start_i + (entry_count * self.DLSIZE) + entry_count - 1)

        # serch target epoch in epoch_table
        found = False
        for i in range(serch_start_index, serch_last_index + 1):
            if epoch_table[i-1] > time_sec:
                found = True
                break
        if not found:
            i = serch_last_index
        record_index = i
        upper_boundary = epoch_table[i-1]
        if i != 1:
            lower_boundary = epoch_table[i-2]
        else:
            lower_boundary = self.start_second
        
        mda_record = self.daf.map_array(self.start_i + ((record_index - 1) * self.DLSIZE),
                                        self.start_i + (record_index * self.DLSIZE) - 1)

        # mda_record : one record of MDA
        # lower_boundary : lower boundary of epoch in this MDA record
        # upper_boundary : upper boundary of epoch in this MDA record
        return mda_record, lower_boundary, upper_boundary

def titlecase(name):
    """Title-case target `name` if it looks safe to do so.
    """
    return name if name.startswith(('1', 'C/', 'DSS-')) else name.title()





