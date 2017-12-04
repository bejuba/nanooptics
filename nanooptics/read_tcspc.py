# PicoHarp 300 file access in python
# This script reads a PicoHarp T2 and T3 Mode data file and writes them to disk.
# Works with file format version 2.0 only!
# Based on the Original matlab and C code by
# Peter Kapusta, Michael Wahl, PicoQuant GmbH 2006, updated May 2007
# Converted to python by Bernd Sontheimer, Humboldt University Berlin, August 2017
# Use at your own risk. No warranties.
# Data is kept in memory despite of the huge amount of memory
# this can take in case of large files. Therefore the file can be split using
# the records_per_split parameter.
# Note that marker events have a lower time resolution and may therefore appear
# in the file slightly out of order with respect to regular (photon) event records.
# This is by design. Markers are designed only for relatively coarse
# synchronization requirements such as image scanning.
import sys
import getopt
import numpy as _np
import string
import struct


def read_picoquant_header_legacy(fid):
    """
    read a picoharp data file header.
    :param fid: input file handle
    :return: file header
    """
    ascii_header = dict()
    binary_header = dict()
    board_header = dict()
    t_mode_header = dict()

    # read ascii header
    ascii_header['Identifier'] = fid.read(16).decode('utf-8')
    ascii_header['Format Version'] = fid.read(6).decode('utf-8')
    if ascii_header['Format Version'][:3] != '2.0':
        print('Warning: This program is for version 2.0 files only.')
    ascii_header['Creator Name'] = fid.read(18).decode('utf-8')
    ascii_header['Creator Version'] = fid.read(12).decode('utf-8')
    ascii_header['File Time'] = fid.read(18).decode('utf-8')
    ascii_header['CRLF'] = fid.read(2).decode('utf-8')
    comment = fid.read(256).decode('utf-8')
    ascii_header['Comment'] = string.replace(comment, '\x00', '')

    # read binary header
    binary_header['Curves'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Bits per Record'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Rounting Channels'], =  _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Number of Boards'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Active Curve'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Measurement Mode'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Sub-Mode'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Range No.'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Offset / ns'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Acquisition Time / ms'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Stop at'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Stop on Overflow'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Restart'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Lin/Log'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Time Axis from'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Time Axis to'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Count Axis from'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Count Axis to'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Display Curve map to'] = _np.fromfile(fid, count=8, dtype=_np.int32)
    binary_header['Display Curve show'] = _np.fromfile(fid, count=8, dtype=_np.int32)
    binary_header['Param0'] = {'Start': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'Stop': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'End': _np.fromfile(fid, count=1, dtype=_np.float32)[0]
                               }
    binary_header['Param1'] = {'Start': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'Stop': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'End': _np.fromfile(fid, count=1, dtype=_np.float32)[0]
                               }
    binary_header['Param2'] = {'Start': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'Stop': _np.fromfile(fid, count=1, dtype=_np.float32)[0],
                               'End': _np.fromfile(fid, count=1, dtype=_np.float32)[0]
                               }
    binary_header['Repeat Mode'], =  _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Repeats per Curve'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Repeat Time'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Repeat Wait Time'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    binary_header['Script Name'] = fid.read(20).decode('utf-8')

    # read board specific header
    board_header['Hardware Identifier'] = fid.read(16).decode('utf-8')
    board_header['Hardware Version'] = fid.read(8).decode('utf-8')
    board_header['Hardware Serial'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['Sync Divider'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['CFD ZeroCross (Ch0) / mV'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['CFD Dischr. (Ch0) / mV'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['CFD ZeroCross (Ch1) / mV'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['CFD Dischr. (Ch1) / mV'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['Resolution / ns'], = _np.fromfile(fid, count=1, dtype=_np.float32)
    board_header['Router Model Code'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['Router Enabled'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    board_header['Rounter Ch1'] = {'Input Type': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Edge': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Present': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD ZeroCross / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   }
    board_header['Rounter Ch2'] = {'Input Type': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Edge': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Present': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD ZeroCross / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   }
    board_header['Rounter Ch3'] = {'Input Type': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Edge': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Present': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD ZeroCross / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   }
    board_header['Rounter Ch4'] = {'Input Type': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input Edge': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Present': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD Level / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   'Input CFD ZeroCross / mV': _np.fromfile(fid, count=1, dtype=_np.int32)[0],
                                   }

    # Router settings are meaningful only for an existing router:
    if board_header['Router Model Code'] == 0:
        del board_header['Rounter Ch1'],
        del board_header['Rounter Ch2'],
        del board_header['Rounter Ch3'],
        del board_header['Rounter Ch4']

    # Read T mode specific header
    t_mode_header['External Devices'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Reserved1'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Reserved2'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Count Rate (Ch0) / Hz'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Count Rate (Ch1) / Hz'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Stop after / ms'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Stop Reason'], = _np.fromfile(fid, count=1, dtype=_np.int32)
    t_mode_header['Number of Records'], = _np.fromfile(fid, count=1, dtype=_np.uint32)

    # Read special imaging header
    imaging_header_size, = _np.fromfile(fid, count=1, dtype=_np.uint32)
    img_header = _np.fromfile(fid, count=imaging_header_size, dtype=_np.int32)

    # combine headers
    header = {
        'ASCII': ascii_header,
        'Binary': binary_header,
        'Board': board_header,
        'T_mode': t_mode_header,
        'Image': img_header
    }
    return header


def read_pt2(s, records_per_split=_np.infty, save_as_npz=False):
    """
    read a picoharp t2 mode .pt2 data file.
    :param s: string of the  input file path, max_memory
    :param records_per_split: split the file every after a number of records and save to disk to conserve memory
    :param save as_as_npz: if True, converted data is saved as compressed and numbered npz file in same directory as
     the input file. If records_per_split is set, save_as_npz defaults to True, else the default is False.
    :return: file header dict, channel numpy array, timestamp numpy array
    """
    with open(s, 'rb') as fid:
        header = read_picoquant_header_legacy(fid)
        number_of_records = header['T_mode']['Number of Records']
        if number_of_records <= records_per_split:
            records_to_read = [number_of_records]
        else:
            splits, remainder = divmod(number_of_records, records_per_split)
            records_to_read = [records_per_split]*splits + [remainder]
            save_as_npz = True
            print('Save output and split file into {} parts'.format(splits))

        # Read the T2 mode event records
        resolution = 4e-12  # 4 ps
        wraparound = _np.uint64(210698240)
        for i, records in enumerate(records_to_read):
            t2_records = _np.fromfile(fid, count=records, dtype=_np.uint32)  # each record is composed of 32 bits
            t2_time = _np.bitwise_and(268435455, t2_records)  # last 28 bits are the time stamp
            t2_channel = _np.right_shift(t2_records, 28)  # first 4 bits are the channel
            t2_channel = _np.uint8(t2_channel)  # first 4 bits are the channel
            del t2_records
            # if the channel is 15, then the record is special
            marker_mask = t2_channel == 15
            # For special records, the lowest 4 bits of the time stamp are actually marker bits
            # by multiplying the mask array with the (marker array + 1) we get 0 where there
            # is no marker, ones where there is an overflow marker and values > 1, where there is
            # a true marker. In this array the originally intended value of the true marker therefore
            # is the true marker -1. We only add one because multiplying arrays in numpy is quite fast.
            t2_marker = _np.multiply(marker_mask, _np.bitwise_and(15, t2_time) + 1)
            t2_time = _np.uint64(t2_time)
            t2_time += _np.cumsum([t2_marker == 1], dtype=_np.uint64) * wraparound
            t2_time = t2_time * resolution
            if save_as_npz:
                _np.savez(s[:-3] + str(i), header=header, t2_channel=t2_channel, t2_time=t2_time)
    return header, t2_channel, t2_time


def read_pt3(s, records_per_split=_np.infty, save_as_npz=False):
    """
    read a picoharp t3 mode .pt3 data file.
    :param s: string of the  input file path
    :param records_per_split: split the file every after a number of records and save to disk to conserve memory
    :param save as_as_npz: if True, converted data is saved as compressed and numbered npz file in same directory as
     the input file. If records_per_split is set, save_as_npz defaults to True, else the default is False.
    :return: file header dict, channel numpy array, timestamp numpy array
    """
    with open(s, 'rb') as fid:
        header = read_picoquant_header_legacy(fid)
        number_of_records = header['T_mode']['Number of Records']
        if number_of_records <= records_per_split:
            records_to_read = [number_of_records]
        else:
            splits, remainder = divmod(number_of_records, records_per_split)
            records_to_read = [records_per_split]*splits + [remainder]
            save_as_npz = True
            print('Save output and split file into {} parts'.format(splits))

        # Read the T3 mode event records
        resolution = 4e-12  # 4 ps is the standard max resolution of the picoharp
        wraparound = _np.uint64(65536)
        for i, records in enumerate(records_to_read):
            t3_records = _np.fromfile(fid, count=records, dtype=_np.uint32)  # each record is composed of 32 bits
            t3_sync = _np.bitwise_and(65535, t3_records)  # last 16 bits are the sync counter stamp
            t3_time = _np.bitwise_and(_np.right_shift(t3_records, 16), 4095)
            t3_channel = _np.right_shift(t3_records, 28)  # first 4 bits are the channel
            t3_channel = _np.uint8(t3_channel)
            del t3_records  # free up unused memory
            # if the channel is 15, then the record is special
            marker_mask = t3_channel == 15
            # For special records, the lowest 4 bits of the time stamp are actually marker bits
            # by multiplying the mask array with the (marker array + 1) we get 0 where there
            # is no marker, ones where there is an overflow marker and values > 1, where there is
            # a true marker. In this array the originally intended value of the true marker therefore
            # is the true marker -1. We only add one because multiplying arrays in numpy is quite fast.
            t3_marker = _np.multiply(marker_mask, _np.bitwise_and(15, t3_time) + 1)
            t3_time = _np.uint64(t3_time)
            t3_time += _np.cumsum([t3_marker == 1], dtype=_np.uint64) * wraparound
            t3_time = t3_time * resolution
            if save_as_npz:
                _np.savez(s[:-4] + str(i), header=header, t3_channel=t3_channel, t3_time=t3_time)
    return header, t3_sync, t3_channel, t3_time


# THIS SHOULD BE IMPLEMENTED SOME TIME IN THE FUTURE :)
#
# def read_ptu(file):
#     def hex2dec(s):
#         """return the integer value of a hexadecimal string"""
#         return int(s, 16)
#     # some constants
#     tyEmpty8      = hex2dec('FFFF0008')
#     tyBool8       = hex2dec('00000008')
#     tyInt8        = hex2dec('10000008')
#     tyBitSet64    = hex2dec('11000008')
#     tyColor8      = hex2dec('12000008')
#     tyFloat8      = hex2dec('20000008')
#     tyTDateTime   = hex2dec('21000008')
#     tyFloat8Array = hex2dec('2001FFFF')
#     tyAnsiString  = hex2dec('4001FFFF')
#     tyWideString  = hex2dec('4002FFFF')
#     tyBinaryBlob  = hex2dec('FFFFFFFF')
#     # RecordTypes
#     rtPicoHarpT3     = hex2dec('00010303')# (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $03 (PicoHarp)
#     rtPicoHarpT2     = hex2dec('00010203')# (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $03 (PicoHarp)
#     rtHydraHarpT3    = hex2dec('00010304')# (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $04 (HydraHarp)
#     rtHydraHarpT2    = hex2dec('00010204')# (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $04 (HydraHarp)
#     rtHydraHarp2T3   = hex2dec('01010304')# (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $03 (T3), HW: $04 (HydraHarp)
#     rtHydraHarp2T2   = hex2dec('01010204')# (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $02 (T2), HW: $04 (HydraHarp)
#     rtTimeHarp260NT3 = hex2dec('00010305')# (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $05 (TimeHarp260N)
#     rtTimeHarp260NT2 = hex2dec('00010205')# (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $05 (TimeHarp260N)
#     rtTimeHarp260PT3 = hex2dec('00010306')# (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $06 (TimeHarp260P)
#     rtTimeHarp260PT2 = hex2dec('00010206')# (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $06 (TimeHarp260P)
#
#     # Globals for subroutines
#     TTResultFormat_TTTRRecType = 0
#     TTResult_NumberOfRecords = 0
#     MeasDesc_Resolution = 0
#     MeasDesc_GlobalResolution = 0
#
#     with open(file, 'rb') as fid:
#         magic = fid.read(8).decode('utf-8')
#         if magic is not 'PQTTR':
#             print('Magic invalid, this is not a PTU file.')
#
#         version = fid.read(8).decode('utf-8')
#         print('Tag version: {}'.format(version))
#
#         #read tag head
#         tag_ident = fid.read(8).decode('utf-8')
#         tag_ident = tag_ident[tag_ident is not 0]

def read_ptu_header(s):

    # some constants
    tyEmpty8 = int('FFFF0008', 16)
    tyBool8 = int('00000008', 16)
    tyInt8 = int('10000008', 16)
    tyBitSet64 = int('11000008', 16)
    tyColor8 = int('12000008', 16)
    tyFloat8 = int('20000008', 16)
    tyTDateTime = int('21000008', 16)
    tyFloat8Array = int('2001FFFF', 16)
    tyAnsiString = int('4001FFFF', 16)
    tyWideString = int('4002FFFF', 16)
    tyBinaryBlob = int('FFFFFFFF', 16)

    # RecordTypes
    rtPicoHarpT3 = int('00010303', 16)  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $03 (PicoHarp)
    rtPicoHarpT2 = int('00010203', 16)  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $03 (PicoHarp)
    rtHydraHarpT3 = int('00010304', 16)  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $04 (HydraHarp)
    rtHydraHarpT2 = int('00010204', 16)  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $04 (HydraHarp)
    rtHydraHarp2T3 = int('01010304', 16)  # (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $03 (T3), HW: $04 (HydraHarp)
    rtHydraHarp2T2 = int('01010204', 16)  # (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $02 (T2), HW: $04 (HydraHarp)
    rtTimeHarp260NT3 = int('00010305', 16)  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $05 (TimeHarp260N)
    rtTimeHarp260NT2 = int('00010205', 16)  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $05 (TimeHarp260N)
    rtTimeHarp260PT3 = int('00010306', 16)  # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $06 (TimeHarp260P)
    rtTimeHarp260PT2 = int('00010206', 16) # (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $06 (TimeHarp260P)

    header = dict()

    with open(s, 'rb') as fid:
        magic = str(fid.read(8))
        if magic[0:6] != "PQTTTR":
            print('Magic invalid, this is not PTU file')
            return
        header['version'] = fid.read(8)
        #print('version', version)

    file_type = {}
    while True:
        # read Tag Head
        TagIdent = fid.read(32)  # TagHead.Ident
        TagIdent = string.replace(TagIdent, '\x00', '')
        header['TagIdent'] = TagIdent
        TagIdx = struct.unpack('i', fid.read(4))[0]  # TagHead.Idx
        TagTyp = _np.array(struct.unpack('i', fid.read(4))[0]).astype(_np.uint32)  # TagHead.Typ
        # TagHead.Value will be read in the right type function
        header['TagIdx']=TagIdx
        if TagIdx > -1:
            EvalName = TagIdent + '(' + str(TagIdx + 1) + ')'
        else:
            EvalName = TagIdent

        # print('eval',str(EvalName))

        if TagTyp == tyEmpty8:
            TagInt = struct.unpack('Q', fid.read(8))[0]
            # print('empty')
        elif TagTyp == tyBool8:
            TagInt = struct.unpack('Q', fid.read(8))[0]
            if TagInt == 0:
                # print('False')
                file_type[EvalName] = False
            else:
                # print('True')
                file_type[EvalName] = True
        elif TagTyp == tyInt8:
            TagInt = struct.unpack('Q', fid.read(8))[0]
            file_type[EvalName] = TagInt
            # print('tyInt8',TagInt)
        elif TagTyp == tyBitSet64:
            TagInt = struct.unpack('Q', fid.read(8))[0]
            file_type[EvalName] = TagInt
            # print('tyBitSet64',TagInt)
        elif TagTyp == tyColor8:
            TagInt = struct.unpack('Q', fid.read(8))[0]
            file_type[EvalName] = TagInt
            # print('tyColor8',TagInt)
        elif TagTyp == tyFloat8:
            TagInt = struct.unpack('d', fid.read(8))[0]
            file_type[EvalName] = TagInt
            # print('tyFloat8',TagInt)
        elif TagTyp == tyFloat8Array:
            TagInt = struct.unpack('Q', fid.read(8))[0]
            file_type[EvalName] = TagInt
            # print '<Float array with'+str(TagInt / 8)+'Entries>'
            # print('tyFloat8Array',TagInt)
            fid.seek(TagInt)
        elif TagTyp == tyTDateTime:
            TagFloat = struct.unpack('d', fid.read(8))[0]
            # print('date'+str(TagFloat))
            file_type[EvalName] = TagFloat
        elif TagTyp == tyAnsiString:
            TagInt = int(struct.unpack('Q', fid.read(8))[0])
            TagString = fid.read(TagInt)
            TagString = string.replace(TagString, '\x00', '')

            # print('tyAnsiString',TagString)
            if TagIdx > -1:
                EvalName = TagIdent + '{' + str(TagIdx + 1) + '}'
            file_type[EvalName] = TagString
        elif TagTyp == tyWideString:
            TagInt = struct.unpack('i', fid.read(4))[0].astype(_np.float64)
            TagString = struct.unpack('i', fid.read(4))[0].astype(_np.float64)

            # print('tyWideString',TagString)
            if TagIdx > -1:
                EvalName = TagIdent + '{' + str(TagIdx + 1) + '}'
            file_type[EvalName] = TagString
        elif TagTyp == tyBinaryBlob:
            TagInt = struct.unpack('i', fid.read(4))[0].astype(_np.float64)
            # print('<Binary Blob with '+str(TagInt)+'Bytes>')
            f.seek(TagInt)
        else:
            print('Illegal Type identifier found! Broken file?', TagTyp)

        if TagIdent == "Header_End":
            break

    TTResultFormat_TTTRRecType = file_type['TTResultFormat_TTTRRecType']
    if TTResultFormat_TTTRRecType == rtPicoHarpT3:
        isT2 = False
        print('PicoHarp T3 data')
    elif TTResultFormat_TTTRRecType == rtPicoHarpT2:
        isT2 = True
        print('PicoHarp T2 data')
    elif TTResultFormat_TTTRRecType == rtHydraHarpT3:
        isT2 = False
        print('HydraHarp V1 T3 data')
    elif TTResultFormat_TTTRRecType == rtHydraHarpT2:
        isT2 = True
        print('HydraHarp V1 T2 data')
    elif TTResultFormat_TTTRRecType == rtHydraHarp2T3:
        isT2 = False
        print('HydraHarp V2 T3 data')
    elif TTResultFormat_TTTRRecType == rtHydraHarp2T2:
        isT2 = True
        print('HydraHarp V2 T2 data')
    elif TTResultFormat_TTTRRecType == rtTimeHarp260NT3:
        isT2 = False
        print('TimeHarp260N T3 data')
    elif TTResultFormat_TTTRRecType == rtTimeHarp260NT2:
        isT2 = True
        print('TimeHarp260P T3 data')
    elif TTResultFormat_TTTRRecType == rtTimeHarp260PT3:
        isT2 = False
        print('TimeHarp260P T3 data')
    elif TTResultFormat_TTTRRecType == rtTimeHarp260PT2:
        isT2 = True
        print('TimeHarp260P T2 data')
    else:
        print('Illegal RecordType')

    if TTResultFormat_TTTRRecType == rtPicoHarpT3:
        return ReadPT3(fid, file_type['TTResult_NumberOfRecords'], file_type['MeasDesc_GlobalResolution'])
    elif TTResultFormat_TTTRRecType == rtPicoHarpT2:  # ReadPT2
        print('this file type of file is not supported yet')
        return False
    elif TTResultFormat_TTTRRecType == rtHydraHarpT3:  # ReadHT3(1)
        return ReadHT3(1, fid, file_type['TTResult_NumberOfRecords'], file_type['MeasDesc_GlobalResolution'])
    elif TTResultFormat_TTTRRecType == rtHydraHarpT2:  # ReadHT3(1)
        return False
    elif TTResultFormat_TTTRRecType == rtHydraHarp2T3:
        return ReadHT3(2, fid, file_type['TTResult_NumberOfRecords'], file_type['MeasDesc_GlobalResolution'])
    elif TTResultFormat_TTTRRecType == rtHydraHarp2T2:  # ReadHT2(2);
        print('this file type of file is not supported yet')
        return False
    elif TTResultFormat_TTTRRecType == rtTimeHarp260NT3:  # ReadHT3(2);
        return ReadHT3(2, fid, file_type['TTResult_NumberOfRecords'], file_type['MeasDesc_GlobalResolution'])
    elif TTResultFormat_TTTRRecType == rtTimeHarp260NT2:  # ReadHT2(2);
        print('this file type of file is not supported yet')
        return False
    elif TTResultFormat_TTTRRecType == rtTimeHarp260PT3:  # ReadHT3(2);
        return ReadHT3(2, fid, file_type['TTResult_NumberOfRecords'], file_type['MeasDesc_GlobalResolution'])
    elif TTResultFormat_TTTRRecType == rtTimeHarp260PT2:  # ReadHT2(2);
        print('this file type of file is not supported yet')
        return False
    else:
        print('Illegal RecordType')
        return False

    return header

