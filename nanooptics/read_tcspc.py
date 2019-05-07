# Read common PicoHarp data files and optionally write them to disk.
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

import numpy as _np
import time

def read_picoquant_legacy_header(fid):
    """
    read a picoharp legacy data file header.
    :param fid: input file handle
    :return: file header
    """
    ascii_header = dict()
    binary_header = dict()
    board_header = dict()
    t_mode_header = dict()

    # read ascii header
    ascii_header['Identifier'] = fid.read(16).decode('utf-8').replace('\x00', '')
    ascii_header['Format Version'] = fid.read(6).decode('utf-8').replace('\x00', '')
    if ascii_header['Format Version'] != '2.0':
        print('Warning: This program is for version 2.0 files only.')
    ascii_header['Creator Name'] = fid.read(18).decode('utf-8').replace('\x00', '')
    ascii_header['Creator Version'] = fid.read(12).decode('utf-8').replace('\x00', '')
    ascii_header['File Time'] = fid.read(18).decode('utf-8').replace('\x00', '')
    ascii_header['CRLF'] = fid.read(2).decode('utf-8').replace('\x00', '')
    ascii_header['Comment'] = fid.read(256).decode('utf-8').replace('\x00', '')

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
    binary_header['Script Name'] = fid.read(20).decode('utf-8').replace('\x00', '')

    # read board specific header
    board_header['Hardware Identifier'] = fid.read(16).decode('utf-8').replace('\x00', '')
    board_header['Hardware Version'] = fid.read(8).decode('utf-8').replace('\x00', '')
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


def read_ptu_header(fid):
    """
    read a picoharp .ptu data file header.
    :param fid: input file handle
    :return: file header
    """
    # Tag Types
    tag_types = {
        0xFFFF0008: 'tyEmpty8',
        0x00000008: 'tyBool8',
        0x10000008: 'tyInt8',
        0x11000008: 'tyBitSet64',
        0x12000008: 'tyColor8',
        0x20000008: 'tyFloat8',
        0x21000008: 'tyTDateTime',
        0x2001FFFF: 'tyFloat8Array',
        0x4001FFFF: 'tyAnsiString',
        0x4002FFFF: 'tyWideString',
        0xFFFFFFFF: 'tyBinaryBlob'
    }

    # Record Types
    rec_types = {
        0x00010303: 'rtPicoHarpT3',
        0x00010203: 'rtPicoHarpT2',
        0x00010304: 'rtHydraHarpT3',
        0x00010204: 'rtHydraHarpT2',
        0x01010304: 'rtHydraHarp2T3',
        0x01010204: 'rtHydraHarp2T2',
        0x00010305: 'rtTimeHarp260NT3',
        0x00010205: 'rtTimeHarp260NT2',
        0x00010306: 'rtTimeHarp260PT3',
        0x00010206: 'rtTimeHarp260PT2',
        0x00010307: 'rtMultiHarpNT3',
        0x00010207: 'rtMultiHarpNT2'
    }

    header = dict()
    magic = fid.read(8).decode('utf-8').replace('\x00', '')
    if magic != 'PQTTTR':
        raise IOError('Not a valid PTU file. '
                      'Magic: '.format(magic))
    header['Magic'] = magic
    header['Version'] = fid.read(8).decode('utf-8').replace('\x00', '')

    while True:
        # read tag header
        tag_ident = fid.read(32).decode('utf-8').replace('\x00', '')
        tag_idx = _np.fromfile(fid, count=1, dtype=_np.int32)[0]
        tag_type = tag_types[_np.fromfile(fid, count=1, dtype=_np.uint32)[0]]
        tag_data = None

        # most tag_value types are int64 but some are double
        if (tag_type == 'tyFloat8') | (type == 'tyFloat8Array'):
            dtype = _np.double
        else:
            dtype = _np.int64

        tag_value = _np.fromfile(fid, count=1, dtype=dtype)[0]

        # Some tag types need additional conversion
        if tag_type == 'tyFloat8':
            tag_value = _np.int64(tag_value).view(_np.float64)
        elif tag_type == 'tyBool8':
            tag_value = bool(tag_value)
        elif tag_type == 'tyTDateTime':
            tag_value = convert_ptu_time(_np.uint64(tag_value).view(_np.float64))

        # Some tag types have additional tag_data
        if tag_type == 'tyAnsiString':
            tag_data = fid.read(tag_value).decode('utf-8').replace('\x00', '')
        elif tag_type == 'tyFloat8Array':
            tag_data = fid.read(tag_value).decode('utf-8').replace('\x00', '')
        elif tag_type == 'tyWideString':
            tag_data = fid.read(tag_value).decode('utf-8').replace('\x00', '')
        elif tag_type == 'tyBinaryBlob':
            tag_data = fid.read(tag_value)
        if tag_ident == "Header_End":
            break
        header[tag_ident] = dict(tag_idx=tag_idx, tag_value=tag_value, tag_data=tag_data)

    # convert to readable rec type
    rec_type = rec_types[header['TTResultFormat_TTTRRecType']['tag_value']]
    header['TTResultFormat_TTTRRecType']['tag_value'] = rec_type
    return header


def read_pt2_records(fid, nrecords):
    """
    read a picoharp T2 Mode records.
    :param fid: input file handle
    :param records: number of records to read
    :return: data: record channel and timestamp numpy arrays
    """

    resolution = 4e-12  # 4ps is the standard max resolution of the picoharp
    wraparound = _np.uint64(210698240)
    
    records = _np.fromfile(fid, count=nrecords, dtype=_np.uint32)  # each record is composed of 32 bits
    # last 28 bits are the time stamp
    time = _np.bitwise_and(int('00001111111111111111111111111111', base=2), records)  
    time = _np.uint64(time)
    # first 4 bits are the channel
    channel = _np.right_shift(records, 28)  
    channel = _np.uint8(channel)  
    
    del records
    
    # if the channel is 15, then the record is special
    marker = channel == int('1111', base=2)
    # For special records, the lowest 4 bits of the time stamp are actually marker bits
    # if this is 0, the marker is an overflow marker
    ofl = marker & (_np.bitwise_and(int('1111', base=2), time) == 0)
    
    marker_channel = _np.bitwise_and(int('1111', base=2), time[marker & ~ofl])
    
    time += _np.cumsum(ofl, dtype=_np.uint64) * wraparound
    time = time * resolution
    
    marker_time = time[marker & ~ofl]
    
    # remove all markers from data:
    channel = channel[~marker]
    time = time[~marker]
    
    records = _np.rec.array([channel, time], 
                            dtype=[('channel', _np.uint8),
                                   ('timestamp', _np.float)])
    markers = _np.rec.array([marker_channel, marker_time],
                            dtype=[('marker', _np.uint8), 
                                   ('timestamp', _np.float)])
    return records, markers

def read_pt3_records(fid, nrecords):
    """
    read a picoharp T3 Mode records.
    :param fid: input file handle
    :param records: number of records to read
    :return: data: record channel and timestamp numpy arrays
    """
        
    resolution = 4e-12  # 4 ps is the standard max resolution of the picoharp
    
    channel_bits = int('11110000000000000000000000000000', base=2)
    dtime_bits   = int('00001111111111110000000000000000', base=2)
    nsync_bits   = int('00000000000000001111111111111111', base=2)
    
    wraparound =  _np.uint64(nsync_bits+1)
    
    records = _np.fromfile(fid, count=nrecords, dtype=_np.uint32)
   
    channel = _np.bitwise_and(channel_bits, records)
    channel = _np.uint8(_np.right_shift(channel, 28))
    
    dtime = _np.bitwise_and(dtime_bits, records)
    dtime = _np.uint64(_np.right_shift(time, 16))
    
    nsync = _np.bitwise_and(nsync_bits, records)
    nsync = _np.uint64(nsync)
    
    wraparound =  nsync_bits + 1

    del records

    # if the channel is 15, then the record is special
    marker = channel == int('1111', base=2)
    # For special records, the lowest 4 bits of the time stamp are actually marker bits
    # if this is 0, the marker is an overflow marker
    ofl = marker & (_np.bitwise_and(int('1111', base=2), nsync) == 0)
    # remove overflows from data:
    marker_channel = _np.bitwise_and(int('1111', base=2), dtime[marker & ~ofl])

    nsync += _np.cumsum(ofl, dtype=_np.uint64) * wraparound
    nsync = nsync * resolution
    dtime = dtime * resolution
    
    channel = channel[~ofl]
    nsync = nsync[~ofl]
    nsync = nsync * resolution
    
    records = _np.rec.array([channel, nsync, dtime], 
                            dtype=[('channel', _np.uint8), 
                                   ('synctime', _np.float), 
                                   ('dtime', _np.float)])
    markers = _np.rec.array([marker_channel, marker_time],
                            dtype=[('marker', _np.uint8),
                                   ('timestamp', _np.float)])
    return records, markers


def read_ht2_records(fid, records, version, resolution):
    """
    read a hydraharp T2 Mode records.
    :param fid: input file handle
    :param records: number of records to read
    :return: data: record channel and timestamp numpy arrays
    """
    # Read the T2 mode event records

    if resolution == 0:
        resolution = 1e-12
    if version == 1:
        wraparound =  _np.uint64(33552000)
    elif version == 2:
        wraparound =  _np.uint64(33554432)
    else: 
        print('Error: Record Version {} not implemented!'.format(version))
        return False
    
    special_bits = int('10000000000000000000000000000000', base=2)
    channel_bits = int('01111110000000000000000000000000', base=2)
    time_bits    = int('00000001111111111111111111111111', base=2)
    
    
    records = _np.fromfile(fid, count=records, dtype=_np.uint32)  # each record is composed of 32 bits
    
    special = _np.bitwise_and(special_bits, records)
    special = special == special_bits
    
    channel = _np.bitwise_and(channel_bits, records)
    channel = _np.uint8(_np.right_shift(channel, 25))
    
    time    = _np.bitwise_and(time_bits, records)
    time = _np.uint64(time)
    
    del records
    ofl = special * (channel == int('111111', base=2))
    if (version == 2):
        ofl = ofl * time
    sync = special & (channel == 0) # all syncs will be in channel 0
    marker = special & (channel >= 1) & (channel <= 15)

    time += _np.cumsum(ofl, dtype=_np.uint64) * wraparound
    time = time * resolution
    
    # extract markers
    marker_channel = channel[marker]
    marker_time = time[marker]
   
    # picoquant adds 1 to all non special event channels
    channel += ~special

    # remove overflows from data:
    channel = channel[~((ofl>0) | marker)]
    time = time[~((ofl>0) | marker)]
        
    print('overflows: {}'.format(np.sum(ofl)))
    print('syncs: {}'.format(np.sum(sync)))
    print('markers: {}'.format(np.sum(marker)))
    
    records = _np.rec.array([channel, time], 
                            dtype=[('channel', _np.uint8),
                                   ('timestamp', _np.float)])
    markers = _np.rec.array([marker_channel, marker_time],
                            dtype=[('marker', _np.uint8), 
                                   ('timestamp', _np.float)])
    return records, markers

def read_ht3_records(fid, records, version, resolution):
    """
    read a hydraharp T3 Mode records.
    :param fid: input file handle
    :param records: number of records to read
    :return: data: record channel and timestamp numpy arrays
    """
    
    if resolution == 0:
        resolution = 1e-12
        
    special_bits = int('10000000000000000000000000000000', base=2)
    channel_bits = int('01111110000000000000000000000000', base=2)
    dtime_bits   = int('00000001111111111111110000000000', base=2)
    nsync_bits   = int('00000000000000000000001111111111', base=2)
    
    wraparound =  _np.uint64(nsync_bits+1)
    
    records = _np.fromfile(fid, count=records, dtype=_np.uint32)  # each record is composed of 32 bits

    special = _np.bitwise_and(special_bits, records)
    special = special == special_bits
    
    channel = _np.bitwise_and(channel_bits, records)
    channel = _np.uint8(_np.right_shift(channel, 25))
    
    dtime = _np.bitwise_and(dtime_bits, records)
    dtime = _np.uint64(_np.right_shift(time, 10))
    
    nsync = _np.bitwise_and(nsync_bits, records)
    nsync = _np.uint64(nsync)

    del records
    
    ofl = special * (channel == int('111111', base=2))
    if (version == 2):
        ofl[(nsync==0) & (ofl>0)] += 1 
        ofl = ofl * nsync 
    marker = special & (channel >= 1) & (channel <= 15)
    
    nsync += _np.cumsum(ofl, dtype=_np.uint64) * wraparound
    nsync = nsync * resolution
    dtime = dtime * resolution
    
    marker_channel = dtime[marker]
    marker_time = nsync[marker]
    
    # picoquant adds 1 to all non special event channels
    channel += ~special
    # remove overflows and markers from data:
    channel = channel[~((ofl>0) | marker)]
    nsync = time[~((ofl>0) | marker)]

    print('overflows: {}'.format(np.sum(ofl)))
    print('syncs: {}'.format(np.sum(sync)))
    print('markers: {}'.format(np.sum(marker)))
    
    records = _np.rec.array([channel, nsync, dtime], 
                            dtype=[('channel', _np.uint8), 
                                   ('synctime', _np.float), 
                                   ('dtime', _np.float)])
    markers = _np.rec.array([marker_channel, marker_time],
                            dtype=[('marker', _np.uint8),
                                   ('timestamp', _np.float)])
    
    return records, markers

    
def read_picoquant(s, records_per_split=_np.infty, save_as_npz=False, ignore_markers=True):
    """
    read a picoharp data files. As of now supported are .pt2 .pt3 and some .ptu files.
    :param s: input file string.
    :param records_per_split: big files may be split if they do not fit into memory.
    if set, save_as_npz defaults to True.
    :param save_as_npz: default is False. if set True the data is written to .npz files in the same directory as a the
    input file.
    :return header: dict of the file header information
    :return data: data converted to list of numpy arrays. Entries are dependent on the input file
    """
    with open(s, 'rb') as fid:
        file_ending = s[-4:]
        # read headers and number of records
        if (file_ending == '.pt2') | (file_ending == '.pt3'):
            header = read_picoquant_legacy_header(fid)
            number_of_records = header['T_mode']['Number of Records']

        elif file_ending == '.ptu':
            header = read_ptu_header(fid)
            number_of_records = header['TTResult_NumberOfRecords']['tag_value']
            resolution = header['MeasDesc_GlobalResolution']['tag_value']
        else:
            print('not a valid file. Only .pt2, .pt3 and .ptu files are supported!')
            return False

        if number_of_records <= records_per_split:
            records_to_read = [number_of_records]
        else:
            splits, remainder = divmod(number_of_records, records_per_split)
            records_to_read = [records_per_split]*splits + [remainder]
            save_as_npz = True
            print('Save output and split file into {} parts'.format(splits))

        data = None
        for i, records in enumerate(records_to_read):
            if file_ending == '.pt2':
                [records, markers] = read_pt2_records(fid, records)
            if file_ending == '.pt3':
                [records, markers] = read_pt3_records(fid, records)
            if file_ending == '.ptu':
                rec_type = header['TTResultFormat_TTTRRecType']['tag_value']
                if rec_type == 'rtPicoHarpT3':
                    [records, markers] = read_pt3_records(fid, records)
                elif rec_type == 'rtPicoHarpT2':
                    [records, markers] = read_pt2_records(fid, records)
                    
                elif rec_type in ['rtHydraHarpT3']:
                    [records, markers] = read_ht3_records(fid, records, 1, resolution)                    
                elif rec_type in ['rtHydraHarp2T3',
                                  'rtTimeHarp260NT3',
                                  'rtTimeHarp260PT3',
                                  'rtMultiHarpNT3']:
                    [records, markers] = read_ht3_records(fid, records, 2, resolution)
                
                elif rec_type in ['rtHydraHarpT2']:
                    [records, markers] = read_ht2_records(fid, records, 1, resolution)
                elif rec_type in ['rtHydraHarp2T2',
                                  'rtTimeHarp260NT2',
                                  'rtTimeHarp260PT2',
                                  'rtMultiHarpNT2']:
                    [records, markers] = read_ht2_records(fid, records, 2, resolution)
                else:
                    print('Illegal or not implemented RecordType')
                    return False
            if save_as_npz:
                _np.savez(s[:-4] + str(i), header=header, records=records, markers=markers)
        if ignore_markers:
            if markers.size > 0:
                print('Warning: There are markers in the file, but they are ignored, \
run with ignore_markers=False to get them')
            return header, records
        else:
            return header, records, markers



def convert_ptu_time(tdatetime):
    """Convert the time used in PTU files to nice time string."""
    epoch_diff = 25569  # days between 30/12/1899 and 01/01/1970
    secs_in_day = 86400  # number of seconds in a day
    t = time.gmtime((tdatetime - epoch_diff) * secs_in_day)
    t = time.strftime("%Y/%m/%d %H:%M:%S", t)
    return t
